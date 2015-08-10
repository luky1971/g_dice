/*
 *
 *                This source code is KINDA part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.5.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon

 * svmanalys is a program for analyzing trajectory files produced by GROMACS using support vector machines.
 * Written by Ahnaf Siddiqui and Dr. Sameer Varma.
 * Copyright (c) 2015, University of South Florida.
 */

#include <stdlib.h>
#include "svm.h"
#include "svmutils.h"
#include "xtcio.h"

#define NUMFILES 2 // The number of input files
#define FRAMESTEP 1000 // The number of new elements to reallocate by when expanding an array of length # of trajectory frames
#define NUMNODE 4 // The number of svm_nodes per training vector (4 = 3 xyz pos + 1 for -1 index)
#define MODELFNLEN 20 // The maximum length of an output model file's filename

void svmanalys(const char *traj_file);
void train_traj(struct svm_node ***data, int natoms, int num_frames);
void print_traj(rvec **pos, int nframes, int natoms);
void print_train_vecs(struct svm_node ***train_vectors, int dim1, int dim2);
void save_print_models(struct svm_model **models, int n);

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"svmanalys analyzes trajectory files produced by GROMACS,",
		"using support vector machine algorithms."
	};
	const char *fnames[NUMFILES];
	int files[] = {TRAJ1, TRAJ2};

	/* Initial setup */
	init_log(argv[0]);
	srand(time(NULL));
	get_file_args(argc, argv, desc, asize(desc), files, fnames, NUMFILES);
	
	/* Call analysis function on input files */
	svmanalys(fnames[0]);

	close_log();
}

void svmanalys(const char *traj_file) {
	/* Reordered svm-training trajectory data: 3Darray[atom #][frame #][position component] */
	struct svm_node ***train_vectors;

	/* Trajectory data */
	t_fileio *traj = NULL;
	int natoms, step;
	real t, prec;
	matrix box;
	rvec *pos;
	gmx_bool b0k;
	int num_frames = 0, est_frames = FRAMESTEP, cur_atom, cur_frame, i;

	/* Open trajectory file and get initial data */
	traj = gmx_fio_open(traj_file, "rb");
	read_first_xtc(traj, &natoms, &step, &t, box, &pos, &prec, &b0k);

	/* Allocate memory for training vectors */
	snew(train_vectors, natoms);
	for (cur_atom = 0; cur_atom < natoms; cur_atom++)
	{
		snew(train_vectors[cur_atom], est_frames);
	}
		
	/** Read and simultaneously reorder trajectory data into svm training format */ 
	gmx_bool renew = FALSE;
	do {
		cur_frame = num_frames;
		num_frames++;
		if(num_frames + 1 > est_frames) {
			est_frames += FRAMESTEP;
			renew = TRUE;
		}
		for (cur_atom = 0; cur_atom < natoms; cur_atom++) {
			if(renew) {
				srenew(train_vectors[cur_atom], est_frames);
			}
			snew(train_vectors[cur_atom][cur_frame], NUMNODE);
			for (i = 0; i < 3; ++i) {
				train_vectors[cur_atom][cur_frame][i].index = i; // Position components are indexed 0:x, 1:y, 2:z
				train_vectors[cur_atom][cur_frame][i].value = pos[cur_atom][i]; // Value of node is a position component
			}
			train_vectors[cur_atom][cur_frame][NUMNODE - 1].index = -1; // -1 index marks end of a data vector
		}
		renew = FALSE;
	} while(read_next_xtc(traj, natoms, &step, &t, box, pos, &prec, &b0k));

	train_traj(train_vectors, natoms, num_frames);

	/* Cleanup */
	gmx_fio_close(traj);// Don't free train_vectors memory because svm_train did something to it? Trying to free causes error
}

void train_traj(struct svm_node ***data, int natoms, int num_frames) {
	struct svm_problem *probs; // Array of svm problems for training
	struct svm_parameter *param; // Parameters used for training
	struct svm_model **models; // Array of svm models produed by training
	double *targets; // Array of target classification values 
	int i;

	snew(probs, natoms);
	snew(models, natoms);
	snew(targets, num_frames);

	/* Construct random target values (1 or 2) for testing */
	for(i = 0; i < num_frames; i++) {
		targets[i] = rand() % 2 + 1;
	}

	/* Construct svm problems */
	for(i = 0; i < natoms; i++) {
		probs[i].l = num_frames;
		probs[i].y = targets;
		probs[i].x = data[i];
	}

	/* Set svm parameters */
	param->svm_type = C_SVC;
	param->kernel_type = RBF;
	param->gamma = 3;
	param->cache_size = 100;
	param->eps = 0.001;
	param->C = 1;
	param->shrinking = 1;
	param->probability = 0;

	/* Train svm */
	for(i = 0; i < natoms; i++) {
		models[i] = svm_train(&(probs[i]), param);
	}

	/* Verify data */
	save_print_models(models, natoms);

	sfree(probs);
	sfree(models);
	sfree(targets);
}

/********************************************************
 * Test/debug functions
 ********************************************************/

/*
 * Prints the given trajectory data to a text file called "traj.txt"
 */
void print_traj(rvec **pos, int nframes, int natoms) {
	int fr, x;
	FILE *f = fopen("traj.txt", "w");

	for(fr = 0; fr < nframes; fr++) {
		fprintf(f, "\nFrame %d\n", fr);
		for(x = 0; x < natoms; x++) {
			fprintf(f, "%d: %f %f %f\n", x, pos[fr][x][0], pos[fr][x][1], pos[fr][x][2]);
		}
	}

	fclose(f);
}

/*
 * Prints reordered trajectory training data to a text file called "train.txt" 
 */
void print_train_vecs(struct svm_node ***train_vectors, int dim1, int dim2) {
	int i, j, k;
	FILE *f = fopen("train.txt", "w");

	for(i = 0; i < dim1; i++) {
		fprintf(f, "\nAtom %d:\n", i + 1);
		for(j = 0; j < dim2; j++) {
			fprintf(f, "\nFrame %d: ", j);
			for(k = 0; k < 3; k++) {
				fprintf(f, "%d:%f ", train_vectors[i][j][k].index, train_vectors[i][j][k].value);
			}
			fprintf(f, "%d\n", train_vectors[i][j][3].index);
		}
	}

	fclose(f);
}

/*
 * Saves the given svm models in model files and writes their data to a text file called "modeldata.txt"
 */
void save_print_models(struct svm_model **models, int n) {
	char mod_fn[MODELFNLEN];
	int i, j, nr_class, nr_sv, *labels, *sv_indices;
	FILE *f = fopen("modeldata.txt", "w");

	for(i = 0; i < n; i++) {
		fprintf(f, "Model %d:\n", i + 1);
		fprintf(f, "SVM Type: %d\n", svm_get_svm_type(models[i]));
		fprintf(f, "Number of classes: %d\n", nr_class = svm_get_nr_class(models[i]));
		
		snew(labels, nr_class);
		svm_get_labels(models[i], labels);
		fprintf(f, "Labels: ");
		for(j = 0; j < nr_class; j++) {
			fprintf(f, "%d ", labels[j]);
		}
		fprintf(f, "\n");
		sfree(labels);

		fprintf(f, "Number of support vectors: %d\n", nr_sv = svm_get_nr_sv(models[i]));

		snew(sv_indices, nr_sv);
		svm_get_sv_indices(models[i], sv_indices);
		fprintf(f, "SV Indices: ");
		for(j = 0; j < nr_sv; j++) {
			fprintf(f, "%d ", sv_indices[j]);
		}
		fprintf(f, "\n");
		sfree(sv_indices);

		fprintf(f, "\n");

		sprintf(mod_fn, "trajmodel%d", i + 1);
		svm_save_model(mod_fn, models[i]);
	}

	fclose(f);
}