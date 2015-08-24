/*
 * Constructs reordered svm_nodes while reading 2 xtc files, then calls train.
 * (Sequential read, non-sequential write)
 */
void xtc_trainA(const char *traj_file1, const char *traj_file2) {
	/* Reordered svm-training trajectory data for 2 trajectories: 3Darray[atom #][frame # x 2][position component]
	 * In the frame # dimension, data for the two trajectories alternate such that:
	 * frame[0] = traj1's pos in frame 0; frame[1] = traj2's pos in frame 0; frame[2] = traj1's pos in frame 1; frame[3] = traj2's pos in frame 1; etc
	 */ 
	struct svm_node ***train_vectors;
	double *targets; // Array of classification labels (Trajectory -1 or 1)

	/* Trajectory data */
	t_fileio *traj1 = NULL, *traj2 = NULL;
	int natoms, step;
	real t, prec;
	matrix box;
	rvec *pos1, *pos2;
	gmx_bool b0k;
	int num_frames = 0, est_frames = FRAMESTEP, cur_atom, cur_frame, cur_fr_index, i;

	/* Open trajectory files and get initial data */
	traj1 = gmx_fio_open(traj_file1, "rb");
	traj2 = gmx_fio_open(traj_file2, "rb");
	read_first_xtc(traj1, &natoms, &step, &t, box, &pos1, &prec, &b0k);
	read_first_xtc(traj2, &natoms, &step, &t, box, &pos2, &prec, &b0k);

	/* Allocate memory for training vectors */
	snew(train_vectors, natoms);
	for(cur_atom = 0; cur_atom < natoms; cur_atom++)
	{
		snew(train_vectors[cur_atom], est_frames * 2);
	}

	/* Read and simultaneously reorder trajectory data from two xtc files into svm training format */ 
	gmx_bool renew = FALSE;

	do {
		cur_frame = num_frames;
		cur_fr_index = cur_frame * 2;
		num_frames++;
		if(num_frames > est_frames) {
			est_frames += FRAMESTEP;
			renew = TRUE;
		}
		for(cur_atom = 0; cur_atom < natoms; cur_atom++) {
			if(renew) {
				srenew(train_vectors[cur_atom], est_frames * 2);
			}
			snew(train_vectors[cur_atom][cur_fr_index], NUMNODE); // Will hold cur_atom's position in current frame in traj1
			snew(train_vectors[cur_atom][cur_fr_index + 1], NUMNODE); // Will hold cur_atom's position in current frame in traj2
			for(i = 0; i < 3; i++) { // Insert position from traj1
				train_vectors[cur_atom][cur_fr_index][i].index = i; // Position components are indexed 0:x, 1:y, 2:z
				train_vectors[cur_atom][cur_fr_index][i].value = pos1[cur_atom][i];
			}
			train_vectors[cur_atom][cur_fr_index][i].index = -1; // -1 index marks end of a data vector
			for(i = 0; i < 3; i++) { // Insert position from traj2
				train_vectors[cur_atom][cur_fr_index + 1][i].index = i;
				train_vectors[cur_atom][cur_fr_index + 1][i].value = pos2[cur_atom][i];
			}
			train_vectors[cur_atom][cur_fr_index + 1][i].index = -1;
		}
		renew = FALSE;
	} while(read_next_xtc(traj1, natoms, &step, &t, box, pos1, &prec, &b0k) 
		&& read_next_xtc(traj2, natoms, &step, &t, box, pos2, &prec, &b0k));

	/* Original vector arrays no longer needed */
	sfree(pos1);
	sfree(pos2);

	/* Build targets array with classification labels */
	int num_data = num_frames * 2;
	snew(targets, num_data);
	for(i = 0; i < num_data; i+=2) {
		targets[i] = LABEL1; // trajectory 1
		targets[i + 1] = LABEL2; // trajectory 2
	}

	/* Verify reordered data */
	//print_train_vecs(train_vectors, natoms, num_data, "train.txt");

	/* Train svm */
	train_trajA(train_vectors, targets, natoms, num_data);

	/* Cleanup */
	gmx_fio_close(traj1);
	gmx_fio_close(traj2);
	sfree(targets);
	// Don't free train_vectors memory because svm_train did something to it? Trying to free causes error
}

/*
 * Reads 2 xtc files first, then reorders data directly into svm_problem and then calls train.
 * (Non-sequential read, sequential write)
 */
void xtc_train(const char *traj_file1, const char *traj_file2) {
	struct svm_problem *probs; // Array of svm problems for training
	double *targets; // Array of classification labels (Trajectory -1 or 1)

	/* Trajectory data */
	t_fileio *traj1 = NULL, *traj2 = NULL;
	int natoms, step;
	real t, prec;
	matrix box;
	rvec **pos1, **pos2;
	gmx_bool b0k;
	int num_frames = 0, est_frames = FRAMESTEP, i;

	snew(pos1, est_frames);
	snew(pos2, est_frames);

	/* Open trajectory files and get initial data */
	traj1 = gmx_fio_open(traj_file1, "rb");
	traj2 = gmx_fio_open(traj_file2, "rb");
	read_first_xtc(traj1, &natoms, &step, &t, box, &(pos1[0]), &prec, &b0k);
	read_first_xtc(traj2, &natoms, &step, &t, box, &(pos2[0]), &prec, &b0k);

	do {
		num_frames++;
		if(num_frames >= est_frames) {
			est_frames += FRAMESTEP;
			srenew(pos1, est_frames);
			srenew(pos2, est_frames);
		}
		snew(pos1[num_frames], natoms);
		snew(pos2[num_frames], natoms);
	} while(read_next_xtc(traj1, natoms, &step, &t, box, pos1[num_frames], &prec, &b0k) 
		&& read_next_xtc(traj2, natoms, &step, &t, box, pos2[num_frames], &prec, &b0k));

	/* Build targets array with classification labels */
	int num_data = num_frames * 2;
	snew(targets, num_data);
	for(i = 0; i < num_frames; i++) {
		targets[i] = LABEL1; // trajectory 1
	}
	for(; i < num_data; i++) {
		targets[i] = LABEL2; // trajectory 2
	}

	/* Construct svm problems */
	snew(probs, natoms);
	int cur_atom, cur_frame, cur_data;
	for(cur_atom = 0; cur_atom < natoms; cur_atom++) {
		probs[cur_atom].l = num_data;
		probs[cur_atom].y = targets;
		snew(probs[cur_atom].x, num_data);
		// Insert positions from traj1
		for(cur_frame = 0; cur_frame < num_frames; cur_frame++) {
			snew(probs[cur_atom].x[cur_frame], NUMNODE);
			for(i = 0; i < 3; i++) {
				probs[cur_atom].x[cur_frame][i].index = i; // Position components are indexed 0:x, 1:y, 2:z
				probs[cur_atom].x[cur_frame][i].value = pos1[cur_frame][cur_atom][i];
			}
			probs[cur_atom].x[cur_frame][i].index = -1; // -1 index marks end of a data vector
		}
		// Insert positions from traj2
		for(cur_frame = 0, cur_data = num_frames; cur_frame < num_frames; cur_frame++, cur_data++) {
			snew(probs[cur_atom].x[cur_data], NUMNODE);
			for(i = 0; i < 3; i++) {
				probs[cur_atom].x[cur_data][i].index = i;
				probs[cur_atom].x[cur_data][i].value = pos2[cur_frame][cur_atom][i];
			}
			probs[cur_atom].x[cur_data][i].index = -1;
		}
	}

	/* Original vector arrays no longer needed */
	for(i = 0; i < num_frames; i++) {
		sfree(pos1[i]);
		sfree(pos2[i]);
	}
	sfree(pos1);
	sfree(pos2);

	/* Verify problem data */
	/*struct svm_node ***train_vectors;
	snew(train_vectors, natoms);
	for(i = 0; i < natoms; i++) {
		train_vectors[i] = probs[i].x;
	}
	print_train_vecs(train_vectors, natoms, targets, num_data, "train.txt");
	sfree(train_vectors);*/

	/* Train svm */
	train_traj(probs, natoms);

	/* Cleanup */
	gmx_fio_close(traj1);
	gmx_fio_close(traj2);
	sfree(probs);
	sfree(targets);
}

void train_trajA(struct svm_node ***data, double *targets, int natoms, int num_data) {
	struct svm_problem *probs; // Array of svm problems for training
	struct svm_parameter param; // Parameters used for training
	struct svm_model **models; // Array of svm models produced by training
	int i;

	snew(probs, natoms);
	snew(models, natoms);

	/* Construct svm problems */
	for(i = 0; i < natoms; i++) {
		probs[i].l = num_data;
		probs[i].y = targets;
		probs[i].x = data[i];
	}
	
	/* Set svm parameters */
	param.svm_type = C_SVC;
	param.kernel_type = RBF;
	param.gamma = 3;
	param.cache_size = 100;
	param.eps = 0.001;
	param.C = 1;
	param.nr_weight = 0;
	param.shrinking = 1;
	param.probability = 0;

	/* Train svm */
	for(i = 0; i < natoms; i++) {
		models[i] = svm_train(&(probs[i]), &param);
	}

	/* Verify data */
	// save_print_models(models, natoms, "modeldata.txt");

	sfree(probs);
	sfree(models);
}