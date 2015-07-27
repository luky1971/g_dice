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
 
#include "svmutils.h"

void trajcom(const char *traj_file, const char *index_file);

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"comanalys calculates the center of mass of a given group of particles over the frames in a trajectory.",
		"It takes as input a trajectory file specified by -f1,",
		"and an index file specified by -n1.",
		"comanalys produces an ASCII data file."
	};
	
	const char *fnames[2];
	int files[] = {TRAJ1, NDX1};
	
	init_log(argv[0]);
	
	get_file_args(argc, argv, desc, asize(desc), files, fnames, 2);
	
	trajcom(fnames[0], fnames[1]);
		
	close_log();
}

void trajcom(const char *traj_file, const char *index_file) {
	//
}
 
 