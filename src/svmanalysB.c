void svmanalysB(const char *traj_file1, const char *traj_file2) {
	/* Trajectory data */
	t_fileio *traj1 = NULL, *traj2 = NULL;
	int natoms, step;
	real t, prec;
	matrix box;
	rvec **pos1, **pos2;
	gmx_bool b0k;
	int num_frames = 0, est_frames = FRAMESTEP, cur_atom, cur_frame, cur_fr_index, i;

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
	// int num_data = num_frames * 2;
	// snew(targets, num_data);
	// for(i = 0; i < num_data; i+=2) {
	// 	targets[i] = LABEL1; // trajectory 1
	// 	targets[i + 1] = LABEL2; // trajectory 2
	// }

	/* Verify reordered data */
	print_traj(pos1, num_frames, natoms, "traj1.xtc");
	print_traj(pos2, num_frames, natoms, "traj2.xtc");

	/* Train svm */
	//train_trajB(train_vectors, targets, natoms, num_data);

	gmx_fio_close(traj1);
	gmx_fio_close(traj2);
	for(i = 1; i < num_frames; i++) {
		sfree(pos1[i]);
		sfree(pos2[i]);
	}
	sfree(pos1);
	sfree(pos2);
}

void train_trajB(struct svm_node ***data, double *targets, int natoms, int num_data) {
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

	// /* Verify data */
	save_print_models(models, natoms);

	sfree(probs);
	sfree(models);
}