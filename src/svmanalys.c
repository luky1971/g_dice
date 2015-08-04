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

#include "svm.h"
#include "svmutils.h"
#include "xtcio.h"

#define FRAMESTEP 1000 // The number of new elements to reallocate by when expanding an array of length # of trajectory frames
#define NUMNODE 4 // The number of svm_nodes per training vector (4 = 3 xyz pos + 1 for -1 index)

void traj2svmprob(const char *traj_file);
void print_traj(rvec **pos, int nframes, int natoms);

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"svmanalys analyzes trajectory files produced by GROMACS,",
		"using support vector machine algorithms."
	};
	
	const char *fnames[1];
	int files[] = {TRAJ1};
	init_log(argv[0]);
	
	get_file_args(argc, argv, desc, asize(desc), files, fnames, NUMFILES);
	
	traj2svmprob(fnames[0]);

	close_log();
}

void traj2svmprob(const char *traj_file) {
	/* Reordered svm-training trajectory data: 3Darray[atom #][frame #][position component] */
	svm_node ***train_vectors;

	/* Trajectory data */
	t_fileio *traj = NULL;
	int natoms, step;
	real t, prec;
	matrix box;
	rvec *pos;
	gmx_bool b0k;
	int num_frames = 0, est_frames = FRAMESTEP, cur_atom, cur_frame;

	/* Open trajectory file and get initial data */
	traj = gmx_fio_open(traj_file, "rb");
	read_first_xtc(traj, &natoms, &step, &t, box, &pos, &prec, &b0k);

	/* Allocate memory for training vectors */
	snew(train_vectors, natoms);
	for (cur_atom = 0; cur_atom < natoms; ++cur_atom)
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
		for (cur_atom = 0; cur_atom < natoms; ++cur_atom) {
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

	/* Cleanup */
	gmx_fio_close(traj);
	for (cur_atom = 0; cur_atom < natoms; ++cur_atom)
	{
		for (cur_frame = 0; cur_frame < num_frames + 1; ++cur_frame)
		{
			sfree(train_vectors[cur_atom][cur_frame]); // Free position vector for each frame
		}
		sfree(train_vectors[cur_atom]); // Free all frames for each atom
	}
	sfree(train_vectors); // Free all atoms
}

void print_traj(rvec **pos, int nframes, int natoms) {
	int f, x;
	for(f = 0; f < nframes; ++f) {
		printf("\nFrame %d\n", f);
		for(x = 0; x < natoms; x++) {
			printf("%d: %f %f %f\n", x, pos[f][x][0], pos[f][x][1], pos[f][x][2]);
		}
	}
}