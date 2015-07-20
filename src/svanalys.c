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
 * Copyright (c) 2001-2008, The GROMACS development team,
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

 * svanalys is a program for analyzing trajectory files produced by GROMACS.
 * Written by Ahnaf Siddiqui and Dr. Sameer Varma.
 * Copyright (c) 2015, University of South Florida.
 */

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include "/usr/local/gromacs/include/gromacs/confio.h"
#include "/usr/local/gromacs/include/gromacs/gmx_fatal.h"
#include "/usr/local/gromacs/include/gromacs/gmxfio.h"
#include "/usr/local/gromacs/include/gromacs/macros.h"
#include "/usr/local/gromacs/include/gromacs/smalloc.h"
#include "/usr/local/gromacs/include/gromacs/statutil.h"
#include "/usr/local/gromacs/include/gromacs/trnio.h"
#include "/usr/local/gromacs/include/gromacs/typedefs.h"
#include "/usr/local/gromacs/include/gromacs/xtcio.h"

void svanalys(int argc, char *argv[]);
void analyze(const char *traj1, const char *traj2, const char *ndx1, const char *ndx2, 
			const char *out_pdb, const char *out_coord_dat, const char *out_res_dat);
void print_log(FILE *f, char const *fmt, ...);
void process_traj(const char *in_name, const char *out_pdb);
void copy_xtc(const char *in_name);
void copy_trr(const char *in_name);
void copy_pdb(const char *in_name, const char *out_name);
void copy_ndx(const char *in_name);

static FILE *out_log;
static output_env_t oenv;

int main(int argc, char *argv[]) {
	svanalys(argc, argv);
	return 0;
}

void svanalys(int argc, char *argv[]) {
	const char *desc[] = {
		"svanalys analyzes trajectory files.",
		"It takes as input two trajectory files specified by -f1 and -f2,",
		"and two index files specified by -n1 and -n2.",
		"svanalys produces a coordinate pdb file specified by -o_atom,",
		"and two ASCII files specified by -eta_atom and -eta_res."
	};
	
	const char *traj1, *traj2, *ndx1, *ndx2;
	const char *out_pdb, *out_coord_dat, *out_res_dat;
	
	t_filenm fnm[] = {
		{efTRX, "-f1", "traj1.xtc", ffREAD},
		{efTRX, "-f2", "traj2.xtc", ffREAD},
		{efNDX, "-n1", "index1.ndx", ffOPTRD},
		{efNDX, "-n2", "index2.ndx", ffOPTRD},
		{efPDB, "-o_atom", "eta_atom.pdb", ffWRITE},
		{efDAT, "-eta_atom", "eta_atom.dat", ffWRITE},
		{efDAT, "-eta_res", "eta_res.dat", ffWRITE}
	};
	
	out_log = fopen("svlog.txt", "w");
	
	parse_common_args(&argc, argv, 0, asize(fnm), fnm, 0, NULL, asize(desc), desc, 0, NULL, &oenv);
	
	traj1 = opt2fn("-f1", asize(fnm), fnm);
	traj2 = opt2fn("-f2", asize(fnm), fnm);
	ndx1 = opt2fn("-n1", asize(fnm), fnm);
	ndx2 = opt2fn("-n2", asize(fnm), fnm);
	out_pdb = opt2fn("-o_atom", asize(fnm), fnm);
	out_coord_dat = opt2fn("-eta_atom", asize(fnm), fnm);
	out_res_dat = opt2fn("-eta_res", asize(fnm), fnm);

	/* Call analysis functions here */
	analyze(traj1, traj2, ndx1, ndx2, out_pdb, out_coord_dat, out_res_dat);
	/***/
	
	fclose(out_log);
}

void analyze(const char *traj1, const char *traj2, const char *ndx1, const char *ndx2, const char *out_pdb, const char *out_coord_dat, const char *out_res_dat) {
	/* Sample analysis code for testing */
	
	FILE *coord_dat, *res_dat;
	
	process_traj(traj1, out_pdb);
	copy_ndx(ndx1);
	
	coord_dat = fopen(out_coord_dat, "w");
	res_dat = fopen(out_res_dat, "w");
	
	/* Write to output dat files */
	
	fclose(coord_dat);
	fclose(res_dat);
	
	/***/
}

/* Prints to both stdout and a given logfile */
void print_log(FILE *f, char const *fmt, ...) {
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

/* Logs and also calls gmx_fatal */
void log_fatal(FILE *f, int fatal_errno, const char *file, int line, char const *fmt, ...) {
	va_list arg;
	if(f != NULL) {
		va_start(arg, fmt);
		fprintf(f, "Fatal error in source file %s line %d:\n", file, line);
		vfprintf(f, fmt, arg);
		va_end(arg);
	}
	va_start(arg, fmt);
	gmx_fatal(fatal_errno, file, line, fmt, arg);
	va_end(arg);
}

/********************************************************
 * Test functions
 ********************************************************/

void process_traj(const char *in_name, const char *out_pdb) {
	switch(fn2ftp(in_name)) {
		case efXTC:
			copy_xtc(in_name);
			break;
		case efTRR:
			copy_trr(in_name);
			break;
		case efPDB:
			copy_pdb(in_name, out_pdb);
			break;
		default:
			log_fatal(out_log, FARGS, "Input trajectory files must be .xtc, .trr, or .pdb!\n");
	}
}

void copy_xtc(const char *in_name) {
	t_fileio *input = NULL, *output = NULL;
	int natoms, step;
	real t, prec;
	matrix box;
	rvec *x;
	gmx_bool b0k;
	
	input = gmx_fio_open(in_name, "rb");
	output = gmx_fio_open("xout.xtc", "wb");
	
	read_first_xtc(input, &natoms, &step, &t, box, &x, &prec, &b0k);
	
	do {
		write_xtc(output, natoms, step, t, box, x, prec);
	} while(read_next_xtc(input, natoms, &step, &t, box, x, &prec, &b0k));
	
	gmx_fio_close(input);
	gmx_fio_close(output);
}

void copy_trr(const char *in_name) {
	t_fileio *input = NULL, *output = NULL;
	t_trnheader header;
	rvec *box, *x, *v;
	gmx_bool b0k;
	
	input = open_trn(in_name, "rb");
	output = open_trn("trout.trr", "wb");
	
	fread_trnheader(input, &header, &b0k);
	
	snew(box, header.box_size);
	snew(x, header.natoms);
	snew(v, header.natoms);
	
	fread_htrn(input, &header, box, x, v, NULL);
	
	do {
		fwrite_trn(output, header.step, header.t, header.lambda, box, header.natoms, x, v, NULL);
	} while(fread_trnheader(input, &header, &b0k) && fread_htrn(input, &header, box, x, v, NULL));
	
	sfree(box);
	sfree(x);
	sfree(v);
	close_trn(input);
	close_trn(output);
}

void copy_pdb(const char *in_name, const char *out_name) {
	FILE *in_pdb;
	char *title;
	t_atoms atoms;
	rvec *x;
	int natoms, ePBC;
	matrix box;
	gmx_bool bChange = FALSE;
	
	in_pdb = fopen(in_name, "r");
	get_pdb_coordnum(in_pdb, &natoms);
	fclose(in_pdb);
	
	init_t_atoms(&atoms, natoms, FALSE);
	snew(x, natoms);
	
	read_pdb_conf(in_name, title, &atoms, x, &ePBC, box, bChange, NULL);
	write_sto_conf(out_name, title, &atoms, x, NULL, ePBC, box);
	
	sfree(x);
}

void copy_ndx(const char *in_name) {
	//
}