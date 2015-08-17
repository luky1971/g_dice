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
 * Written by Ahnaf Siddiqui, Mohsen Botlani-Esfahani, and Dr. Sameer Varma.
 * Copyright (c) 2015, University of South Florida.
 */

#include <stdlib.h>
#include <time.h>
#include "svm.h"
#include "svio.h"
#include "svutils.h"
#include "xtcio.h"

#define LABEL1 -1 // Classification label for trajectory 1
#define LABEL2 1 // Classification label for trajectory 2
#define SHRINK 1 // Whether or not to use shrinking heuristics in svm_train

void svmanalys(const char *traj_file1, const char *traj_file2, const char *ndx_file1, const char *ndx_file2);
void train_traj(struct svm_problem *probs, int num_probs);

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"svmanalys analyzes trajectory files produced by GROMACS,",
		"using support vector machine algorithms."
	};
#define NUMFILES 4 // The number of input files
	const char *fnames[NUMFILES];
	int files[] = {TRAJ1, TRAJ2, NDX1, NDX2};

	/* Initial setup */
	init_log(argv[0]);
	get_file_args(argc, argv, desc, asize(desc), files, fnames, NUMFILES);
	
	/* Call analysis function on input files */
	clock_t start = clock();
	svmanalys(fnames[0], fnames[1], fnames[2], fnames[3]);
	clock_t end = clock();
	print_log("Execution time: %d\n", end - start);

	close_log();
}

void svmanalys(const char *traj_file1, const char *traj_file2, const char *ndx_file1, const char *ndx_file2) {
	const char *io_error = "Input trajectory files must be .xtc, .trr, or .pdb!\n";

	/* Training data */
	struct svm_problem *probs; // Array of svm problems for training
	double *targets; // Array of classification labels (Trajectory -1 or 1)

	/* Trajectory data */
	rvec **pos1, **pos2; // Trajectory position vectors
	int nframes, nframes2, natoms, natoms2, num_data, i;

	switch(fn2ftp(traj_file1)) {
		case efXTC:
			read_xtc(traj_file1, &pos1, &nframes, &natoms);
			break;
		case efTRR:
			read_trr(traj_file1, &pos1, &nframes, &natoms);
			break;
		case efPDB:
			// read_pdb(traj_file1, &pos1, &natoms);
			break;
		default:
			log_fatal(FARGS, io_error);
	}

	switch(fn2ftp(traj_file2)) {
		case efXTC:
			read_xtc(traj_file2, &pos2, &nframes2, &natoms2);
			break;
		case efTRR:
			read_trr(traj_file2, &pos2, &nframes2, &natoms2);
			break;
		case efPDB:
			// read_pdb(traj_file2, &pos2, &natoms);
			break;
		default:
			log_fatal(FARGS, io_error);
	}

	/* In case input files have different numbers of frames or atoms */
	nframes = nframes2 < nframes ? nframes2 : nframes;
	natoms = natoms2 < natoms ? natoms2 : natoms;

	/* Verify read data */
	// print_traj(pos1, nframes, natoms, "traj1.txt");
	// print_traj(pos2, nframes, natoms, "traj2.txt");

	/* Index data */
#define NUMGROUPS 1
	int *isize, *isize2;
	atom_id **indx1, **indx2; // Atom indices for the two trajectories
	char **grp_names;

	snew(isize, NUMGROUPS);
	snew(indx1, NUMGROUPS);
	snew(grp_names, NUMGROUPS);

	/* If an index file was given, get atom group with indices that will be trained */
	if(ndx_file1 != NULL) {
		rd_index(ndx_file1, NUMGROUPS, isize, indx1, grp_names);
	}
	else { // If no index file, set default indices as 0 to natoms - 1
		snew(indx1[0], natoms);
		isize[0] = natoms;
		for(i = 0; i < natoms; i++) {
			indx1[0][i] = i;
		}
	}
	if(ndx_file2 != NULL) {
		snew(isize2, NUMGROUPS);
		snew(indx2, NUMGROUPS);
		rd_index(ndx_file2, NUMGROUPS, isize2, indx2, grp_names);
		if(isize2[0] != isize[0]) {
			log_fatal(FARGS, "Given index groups have different numbers of atoms!\n");
		}
	}
	else {
		indx2 = indx1;
	}
	sfree(grp_names);

	/* Verify indices */
	for(i = 0; i < isize[0]; i++) {
		printf("%d: %d and %d\n", i, indx1[0][i], indx2[0][i]);
	}

	num_data = nframes * 2;

	/* Build targets array with classification labels */
	snew(targets, num_data);
	for(i = 0; i < nframes; i++) {
		targets[i] = LABEL1; // trajectory 1
	}
	for(; i < num_data; i++) {
		targets[i] = LABEL2; // trajectory 2
	}

	/* Construct svm problems */
	snew(probs, isize[0]);
	int cur_atom, cur_frame, cur_data;
	for(cur_atom = 0; cur_atom < isize[0]; cur_atom++) {
		probs[cur_atom].l = num_data;
		probs[cur_atom].y = targets;
		snew(probs[cur_atom].x, num_data);
		// Insert positions from traj1
		for(cur_frame = 0; cur_frame < nframes; cur_frame++) {
			snew(probs[cur_atom].x[cur_frame], 4); // (4 = 3 xyz pos + 1 for -1 end index)
			for(i = 0; i < 3; i++) {
				probs[cur_atom].x[cur_frame][i].index = i; // Position components are indexed 0:x, 1:y, 2:z
				probs[cur_atom].x[cur_frame][i].value = pos1[cur_frame][indx1[0][cur_atom]][i];
			}
			probs[cur_atom].x[cur_frame][i].index = -1; // -1 index marks end of a data vector
		}
		// Insert positions from traj2
		for(cur_frame = 0, cur_data = nframes; cur_frame < nframes; cur_frame++, cur_data++) {
			snew(probs[cur_atom].x[cur_data], 4);
			for(i = 0; i < 3; i++) {
				probs[cur_atom].x[cur_data][i].index = i;
				probs[cur_atom].x[cur_data][i].value = pos2[cur_frame][indx2[0][cur_atom]][i];
			}
			probs[cur_atom].x[cur_data][i].index = -1;
		}
	}

	/* Verify problem data */
	print_svm_probs(probs, isize[0], "probs.txt");

	/* No longer need original vectors */
	for(i = 0; i < nframes; i++) {
		sfree(pos1[i]);
		sfree(pos2[i]);
	}
	sfree(pos1);
	sfree(pos2);

	/* Train SVM */
	train_traj(probs, isize[0]);

	sfree(isize);
	sfree(indx1[0]);
	sfree(indx1);
	if(ndx_file2 != NULL) {
		sfree(isize2);
		sfree(indx2[0]);
		sfree(indx2);
	}
	sfree(probs); // Don't free the data within probs, will cause error
	sfree(targets);
}

void train_traj(struct svm_problem *probs, int num_probs) {
	struct svm_parameter param; // Parameters used for training
	struct svm_model **models; // Array of svm models produced by training
	int i;

	snew(models, num_probs);
	
	/* Set svm parameters */
	param.svm_type = C_SVC;
	param.kernel_type = RBF;
	param.gamma = 3;
	param.cache_size = 100;
	param.eps = 0.001;
	param.C = 1;
	param.nr_weight = 0;
	param.shrinking = SHRINK;
	param.probability = 0;

	/* Train svm */
	for(i = 0; i < num_probs; i++) {
		models[i] = svm_train(&(probs[i]), &param);
	}

	/* Verify data */
	save_print_models(models, num_probs, "modeldata.txt");

	sfree(models);
}