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

 * Written by Ahnaf Siddiqui and Dr. Sameer Varma.
 * Copyright (c) 2015, University of South Florida.
 */

#include "svmutils.h"

static FILE *out_log = NULL;
static output_env_t oenv = NULL;

/*
 * Gets file names from command line
 * int files[] are the desired files based on files enum defined in svmutils.h
 * Requested file names will be stored in fnames
 */
void get_file_args(int argc, char *argv[], const char *desc[], int desc_size, int files[], const char *fnames[], int num_files) {
	int filetypes[] = {efTRX, efTRX, efNDX, efNDX, efPDB, efDAT, efDAT};
	const char *options[] = {"-f1", "-f2", "-n1", "-n2", "-o_atom", "-eta_atom", "-eta_res"};
	const char *def_names[] = {"traj1.xtc", "traj2.xtc", "index1.ndx", "index2.ndx", "eta_atom.pdb", "eta_atom.dat", "eta_res.dat"};
	unsigned long modes[] = {ffREAD, ffREAD, ffOPTRD, ffOPTRD, ffWRITE, ffWRITE, ffWRITE};
	int i;
	t_filenm *fnm;
	
	snew(fnm, num_files);
	
	for(i = 0; i < num_files; i++) {
		fnm[i].ftp = filetypes[files[i]];
		fnm[i].opt = options[files[i]];
		fnm[i].fn = def_names[files[i]];
		fnm[i].flag = modes[files[i]];
	}

	parse_common_args(&argc, argv, 0, num_files, fnm, 0, NULL, desc_size, desc, 0, NULL, &oenv);
	
	for(i = 0; i < num_files; i++) {
		fnames[i] = opt2fn(options[files[i]], num_files, fnm);
	}

	sfree(fnm);
}

/* Opens the logfile and logs initial time/date */
void init_log(char *program) {
	out_log = fopen("svmlog.txt", "a");
	
	time_t t = time(NULL);
	struct tm *ltime = localtime(&t);
	fprintf(out_log, "\n%s run: %d-%d-%d %d:%d:%d\n", program, ltime->tm_mon + 1, ltime->tm_mday, ltime->tm_year + 1900, ltime->tm_hour, ltime->tm_min, ltime->tm_sec);
}

/* Closes the logfile */
void close_log(void) {
	fclose(out_log);
}

/* Prints to both stdout and the logfile */
void print_log(char const *fmt, ...) {
	va_list arg;
	va_start(arg, fmt);
	vprintf(fmt, arg);
	va_end(arg);
	if(out_log != NULL) {
		va_start(arg, fmt);
		vfprintf(out_log, fmt, arg);
		va_end(arg);
	}
}

/* 
 * Logs fatal error to logfile and also calls gmx_fatal
 * Hint: Use FARGS for the first 3 arguments.
 */
void log_fatal(int fatal_errno, const char *file, int line, char const *fmt, ...) {
	va_list arg;
	if(out_log != NULL) {
		va_start(arg, fmt);
		fprintf(out_log, "Fatal error in source file %s line %d: ", file, line);
		vfprintf(out_log, fmt, arg);
		va_end(arg);
	}
	va_start(arg, fmt);
	gmx_fatal(fatal_errno, file, line, fmt, arg);
	va_end(arg);
}