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

 * Written by Ahnaf Siddiqui, Mohsen Botlani-Esfahani, and Dr. Sameer Varma.
 * Copyright (c) 2015, University of South Florida.
 * The authors would like to acknowledge the use of the services provided by Research Computing at the University of South Florida.
 */

#include "svio.h"

void read_xtc(const char *traj_fname, rvec ***x, int *nframes, int *natoms) {
	print_log("Reading xtc.\n");
	t_fileio *traj = NULL;
	int step;
	real t, prec;
	matrix box;
	gmx_bool b0k;
	int est_frames = FRAMESTEP;
	*nframes = 0;
	
	traj = gmx_fio_open(traj_fname, "rb");

	snew(*x, est_frames);
	read_first_xtc(traj, natoms, &step, &t, box, &((*x)[0]), &prec, &b0k);
	
	do {
		(*nframes)++;
		if(*nframes >= est_frames) {
			est_frames += FRAMESTEP;
			srenew(*x, est_frames);
		}
		snew((*x)[*nframes], *natoms);
	} while(read_next_xtc(traj, *natoms, &step, &t, box, (*x)[*nframes], &prec, &b0k));
	
	sfree((*x)[*nframes]);
	gmx_fio_close(traj);
}

void read_trr(const char *traj_fname, rvec ***x, int *nframes, int *natoms) {
	print_log("Reading trr.\n");
	t_fileio *traj = NULL;
	t_trnheader header;
	gmx_bool b0k;
	*nframes = 0;
	int est_frames = FRAMESTEP;
	
	traj = gmx_fio_open(traj_fname, "rb");

	snew(*x, est_frames);
	while(fread_trnheader(traj, &header, &b0k)) {
		if(*nframes >= est_frames) {
			est_frames += FRAMESTEP;
			srenew(*x, est_frames);
		}
		snew((*x)[*nframes], header.natoms);
		if(fread_htrn(traj, &header, NULL, (*x)[*nframes], NULL, NULL)) {
			(*nframes)++;
		}
		else {
			break;
		}
	}
	*natoms = header.natoms;
	
	gmx_fio_close(traj);
}

void read_pdb(const char *pdb_fname, rvec **x, int *natoms) {
	FILE *pdb_file;
	char *title;
	t_atoms atoms;
	rvec *pos;
	int ePBC;
	matrix box;
	gmx_bool bChange = FALSE;
	
	pdb_file = fopen(pdb_fname, "r");
	get_pdb_coordnum(pdb_file, natoms);
	fclose(pdb_file);
	
	init_t_atoms(&atoms, *natoms, FALSE);
	snew(pos, *natoms);
	
	read_pdb_conf(pdb_fname, title, &atoms, pos, &ePBC, box, bChange, NULL);
	print_log("%s %d\n", pdb_fname, *natoms);
	sfree(pos);
}

void copy_ndx(const char *in_name, int num_groups) {
	int *isize;
	atom_id **indx;
	char **grp_names;
	t_blocka *blk;
	int i;
	
	snew(isize, num_groups);
	snew(indx, num_groups);
	snew(grp_names, num_groups);
	
	rd_index((char*)in_name, num_groups, isize, indx, grp_names);
	
	blk = new_blocka();
	
	for(i = 0; i < num_groups; i++) {
		add_grp(blk, &grp_names, isize[i], indx[i], grp_names[i]);
	}
	
	write_index("outdex.ndx", blk, grp_names);
	
	sfree(isize);
	sfree(indx);
	sfree(grp_names);
}