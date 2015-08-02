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

#include "svmutils.h"
#include "xtcio.h"

void order_xtc(const char *traj_file);
void print_traj(rvec **pos, int nframes, int natoms);

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"svmanalys analyzes trajectory files produced by GROMACS,",
		"using support vector machine algorithms."
	};
	
#define NUMFILES 1
	const char *fnames[NUMFILES];
	int files[] = {TRAJ1};
	init_log(argv[0]);
	
	get_file_args(argc, argv, desc, asize(desc), files, fnames, NUMFILES);
	
	order_xtc(fnames[0]);

	close_log();
}

void order_xtc(const char *traj_file) {
	/* Reordered trajectory: row (outer array index) = atom index, column (inner array index) = frame */
	//rvec **ord_pos;
#define FRAMESTEP 1000
	/* Trajectory data */
	t_fileio *traj = NULL;
	int natoms, step, num_frames = 0, est_frames = FRAMESTEP;
	real t, prec;
	matrix box;
	rvec **pos;
	gmx_bool b0k;
	
	traj = gmx_fio_open(traj_file, "rb");
	
	snew(pos, est_frames);
	
	read_first_xtc(traj, &natoms, &step, &t, box, pos, &prec, &b0k);
		
	do {
		num_frames++;
		if(num_frames + 1 > est_frames) {
			est_frames += FRAMESTEP;
			srenew(pos, est_frames);
		}
		snew(pos[num_frames], natoms);
	} while(read_next_xtc(traj, natoms, &step, &t, box, pos[num_frames], &prec, &b0k));

	print_traj(pos, num_frames, natoms);

	int i;
	for (i = 1; i <= num_frames; ++i) {
		sfree(pos[i]);
	}

	sfree(pos);
	gmx_fio_close(traj);
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