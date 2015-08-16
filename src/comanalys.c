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

 * comanalys is a program for calculating the center of mass of particles in a GROMACS trajectory file.
 * Written by Ahnaf Siddiqui and Dr. Sameer Varma.
 * Copyright (c) 2015, University of South Florida.
 */
 
#include "svutils.h"
#include "xtcio.h"

void trajcom(const char *traj_file, const char *index_file);
void center_of_mass(rvec com, rvec *pos, int natoms, atom_id *indx, int isize);
void log_results(rvec *com, int n, rvec com_avg);

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"comanalys calculates the center of mass of a given group of particles over the frames in a trajectory.",
		"It takes as input a trajectory file specified by -f1,",
		"and an index file specified by -n1.",
		"comanalys produces an ASCII data file."
	};
	
#define NUMFILES 2
	const char *fnames[NUMFILES];
	int files[] = {TRAJ1, NDX1};
	
	init_log(argv[0]);
	
	get_file_args(argc, argv, desc, asize(desc), files, fnames, NUMFILES);
	
	trajcom(fnames[0], fnames[1]);
		
	close_log();
}

void trajcom(const char *traj_file, const char *index_file) {
	/* Index data */
	int isize[1];
	atom_id *indx[1];
	char *grp_names[1];
	
	/* Trajectory data */
	t_fileio *traj = NULL;
	int natoms, step;
	real t, prec;
	matrix box;
	rvec *pos;
	gmx_bool b0k;
	
	/* Center of mass data */
	rvec *com;
	rvec com_avg;
	real x_tot = 0, y_tot = 0, z_tot = 0;
	
	/* Counts */
	int num_frames = 0, com_len;
	
	/* Read index */
	rd_index((char*)index_file, 1, isize, indx, grp_names);
	
	traj = gmx_fio_open(traj_file, "rb");
	
	com_len = 10;
	snew(com, com_len);
	
	read_first_xtc(traj, &natoms, &step, &t, box, &pos, &prec, &b0k);
	
	do {
		num_frames++;
		if(num_frames > com_len) {
			com_len += 10;
			srenew(com, com_len);
		}
		center_of_mass(com[num_frames - 1], pos, natoms, indx[0], isize[0]);
		x_tot += com[num_frames - 1][0];
		y_tot += com[num_frames - 1][1];
		z_tot += com[num_frames - 1][2];
	} while(read_next_xtc(traj, natoms, &step, &t, box, pos, &prec, &b0k));
	
	com_avg[0] = x_tot / num_frames;
	com_avg[1] = y_tot / num_frames;
	com_avg[2] = z_tot / num_frames;
	
	log_results(com, num_frames, com_avg);
	
	sfree(com);
	gmx_fio_close(traj);
}

void center_of_mass(rvec com, rvec *pos, int natoms, atom_id *indx, int isize) {
	real x_tot = 0, y_tot = 0, z_tot = 0;
	int i;
	
	for(i = 0; i < isize; i++) {
		x_tot += pos[indx[i]][0];
		y_tot += pos[indx[i]][1];
		z_tot += pos[indx[i]][2];
	}
	
	com[0] = x_tot / isize;
	com[1] = y_tot / isize;
	com[2] = z_tot / isize;
}

void log_results(rvec *com, int n, rvec com_avg) {
	int i;
	print_log("\nCenters of mass:\n");
	for(i = 0; i < n; i++) {
		print_log("Frame %d: %f %f %f\n", i, com[i][0], com[i][1], com[i][2]);
	}
	print_log("Average: %f %f %f\n", com_avg[0], com_avg[1], com_avg[2]);
	
	printf("\nResults saved in svmlog.txt\n");
}