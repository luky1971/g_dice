/*
 * Written by Ahnaf Siddiqui.
 * Copyright (c) 2015, University of South Florida.
 *
 * Requires the gromacs library libgmx. Link it with -lgmx when compiling.
 * Header xtcio.h usually located in /usr/local/gromacs/include/gromacs/
 */

#include <stdio.h>
#include "/usr/local/gromacs/include/gromacs/xtcio.h"

int main(int argc, char *argv[]) {
	t_fileio *input = NULL;
	t_fileio *output = NULL;
	int natoms, step;
	real t, prec;
	matrix box;
	rvec *x;
	gmx_bool b0k;

	if(argc < 2) {
		printf("Must specify a .xtc file for input!\n");
		return 1;
	}

	input = open_xtc(argv[1], "rb");
	output = open_xtc("xout.xtc", "wb");
	
	read_first_xtc(input, &natoms, &step, &t, box, &x, &prec, &b0k);
	
	do {
		write_xtc(output, natoms, step, t, box, x, prec);
	} while(read_next_xtc(input, natoms, &step, &t, box, x, &prec, &b0k));
}