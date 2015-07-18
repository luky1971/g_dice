/*
 * Written by Ahnaf Siddiqui.
 * Copyright (c) 2015, University of South Florida.
 *
 * Requires the gromacs library libgmx. Link it with -lgmx when compiling.
 * Header xtcio.h usually located in /usr/local/gromacs/include/gromacs/
 * 
 * TODO: Check if files can be opened!!
 */

#include <stdio.h>
#include <string.h>
#include "/usr/local/gromacs/include/gromacs/macros.h"
#include "/usr/local/gromacs/include/gromacs/statutil.h"
#include "/usr/local/gromacs/include/gromacs/xtcio.h"

void x_io(const char *fn);

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"svanalys analyzes trajectory files.",
		"It takes as input two trajectory files specified by -f1 and -f2,",
		"and two index files specified by -n1 and -n2.",
		"svanalys produces a coordinate pdb file specified by -o_atom,",
		"and two ASCII files specified by -eta_atom and -eta_res."
	};
	
	FILE *out_log;
	
	t_filenm fnm[] = {
		{efTRX, "-f1", "traj1", ffREAD},
		{efTRX, "-f2", "traj2", ffREAD},
		{efNDX, "-n1", "index1", ffOPTRD},
		{efNDX, "-n2", "index2", ffOPTRD},
		{efPDB, "-o_atom", "eta_atom", ffWRITE},
		{efDAT, "-eta_atom", "eta_atom", ffWRITE},
		{efDAT, "-eta_res", "eta_res", ffWRITE}
	};
	
	output_env_t oenv;
	
	if( (out_log = fopen("svlog.txt", "w")) == NULL ) {
		out_log = stdout;
	}
	
	parse_common_args(&argc, argv, 0, 
		asize(fnm), fnm, 0, NULL, asize(desc), desc, 0, NULL, &oenv);
	
	x_io(opt2fn("-f1", asize(fnm), fnm));
	
	fclose(out_log);
	return 0;
}

void x_io(const char *fn) {
	t_fileio *input = NULL;
	t_fileio *output = NULL;
	int natoms, step;
	real t, prec;
	matrix box;
	rvec *x;
	gmx_bool b0k;
	
	input = open_xtc(fn, "rb");
	output = open_xtc("xout.xtc", "wb");
	
	read_first_xtc(input, &natoms, &step, &t, box, &x, &prec, &b0k);
	
	do {
		write_xtc(output, natoms, step, t, box, x, prec);
	} while(read_next_xtc(input, natoms, &step, &t, box, x, &prec, &b0k));
	
	close_xtc(input);
	close_xtc(output);
}