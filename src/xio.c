/*
 * Written by Ahnaf Siddiqui.
 * Copyright (c) 2015, University of South Florida.
 *
 * Requires the gromacs library libgmx. Link it with -lgmx when compiling.
 * Header xtcio.h usually located in /usr/local/gromacs/include/gromacs/
 * 
 */

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include "/usr/local/gromacs/include/gromacs/gmxfio.h"
#include "/usr/local/gromacs/include/gromacs/macros.h"
#include "/usr/local/gromacs/include/gromacs/statutil.h"
#include "/usr/local/gromacs/include/gromacs/xtcio.h"

void copy_xtc(t_fileio *input);
void log_print(FILE *f, char const *fmt, ...);

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"svanalys analyzes trajectory files.",
		"It takes as input two trajectory files specified by -f1 and -f2,",
		"and two index files specified by -n1 and -n2.",
		"svanalys produces a coordinate pdb file specified by -o_atom,",
		"and two ASCII files specified by -eta_atom and -eta_res."
	};
	
	FILE *out_log;
	t_fileio *traj1 = NULL, *traj2 = NULL, *ndx1 = NULL, *ndx2 = NULL;
	t_fileio *pdb_coord = NULL, *dat_coord = NULL, *dat_res = NULL;
	
	t_filenm fnm[] = {
		{efTRX, "-f1", "traj1.xtc", ffREAD},
		{efTRX, "-f2", "traj2.xtc", ffREAD},
		{efNDX, "-n1", "index1.ndx", ffOPTRD},
		{efNDX, "-n2", "index2.ndx", ffOPTRD},
		{efPDB, "-o_atom", "eta_atom.pdb", ffWRITE},
		{efDAT, "-eta_atom", "eta_atom.dat", ffWRITE},
		{efDAT, "-eta_res", "eta_res.dat", ffWRITE}
	};
	
	output_env_t oenv;
	
	out_log = fopen("svlog.txt", "w");
	
	parse_common_args(&argc, argv, 0, asize(fnm), fnm, 0, NULL, asize(desc), desc, 0, NULL, &oenv);
	
	if(opt2bSet("-f1", asize(fnm), fnm))
		traj1 = gmx_fio_open(opt2fn("-f1", asize(fnm), fnm), "rb");
	if(opt2bSet("-f2", asize(fnm), fnm))
		traj2 = gmx_fio_open(opt2fn("-f2", asize(fnm), fnm), "rb");
	if(opt2bSet("-n1", asize(fnm), fnm))
		ndx1 = gmx_fio_open(opt2fn("-n1", asize(fnm), fnm), "r");
	if(opt2bSet("-n2", asize(fnm), fnm))
		ndx2 = gmx_fio_open(opt2fn("-n2", asize(fnm), fnm), "r");
	
	pdb_coord = gmx_fio_open(opt2fn("-o_atom", asize(fnm), fnm), "w");
	dat_coord = gmx_fio_open(opt2fn("-eta_atom", asize(fnm), fnm), "w");
	dat_res = gmx_fio_open(opt2fn("-eta_res", asize(fnm), fnm), "w");
	
	if(traj1 != NULL)
		copy_xtc(traj1);
	
	if(traj1 != NULL)		gmx_fio_close(traj1);
	if(traj2 != NULL)		gmx_fio_close(traj2);
	if(ndx1 != NULL)		gmx_fio_close(ndx1);
	if(ndx2 != NULL)		gmx_fio_close(ndx2);
	if(pdb_coord != NULL)	gmx_fio_close(pdb_coord);
	if(dat_coord != NULL)	gmx_fio_close(dat_coord);
	if(dat_res != NULL)		gmx_fio_close(dat_res);
	fclose(out_log);
	
	return 0;
}

void copy_xtc(t_fileio *input) {
	t_fileio *output = NULL;
	int natoms, step;
	real t, prec;
	matrix box;
	rvec *x;
	gmx_bool b0k;
	
	output = gmx_fio_open("xout.xtc", "wb");
	
	read_first_xtc(input, &natoms, &step, &t, box, &x, &prec, &b0k);
	
	do {
		write_xtc(output, natoms, step, t, box, x, prec);
	} while(read_next_xtc(input, natoms, &step, &t, box, x, &prec, &b0k));
	
	gmx_fio_close(output);
}

void log_print(FILE *f, char const *fmt, ...) {
	va_list arg;
	va_start(arg, fmt);
	vprintf(fmt, arg);
	va_end(arg);
	if(f != NULL) {
		va_start(arg, fmt);
		vfprintf(f, fmt, arg);
		va_end(arg);
	}
}