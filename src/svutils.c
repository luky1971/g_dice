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

#include "svutils.h"

static FILE *out_log = NULL;
static output_env_t oenv = NULL;

/*
 * Gets file names from command line
 * int files[] are the desired files based on files enum defined in svmutils.h
 * Requested file names will be stored in fnames in same order as files[]
 * opt_set contains whether each corresponding file's option was actually found on the command line.
 */
void get_file_args(int argc, char *argv[], const char *desc[], int desc_size, 
	int files[], const char *fnames[], int num_files) {

	int filetypes[] = {efTRX, efTRX, efNDX, efNDX, efPDB, efDAT, efDAT};
	const char *options[] = {"-f1", "-f2", "-n1", "-n2", "-o_atom", "-eta_atom", "-eta_res"};
	const char *def_names[] = {"traj1.xtc", "traj2.xtc", "index1.ndx", "index2.ndx", 
		"eta_atom.pdb", "eta_atom.dat", "eta_res.dat"};
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
		fnames[i] = opt2fn_null(options[files[i]], num_files, fnm);
	}

	sfree(fnm);
}

/* Opens the logfile and logs initial time/date */
void init_log(const char *program) {
	out_log = fopen("svmlog.txt", "a");
	
	time_t t = time(NULL);
	struct tm *ltime = localtime(&t);
	fprintf(out_log, "\n%s run: %d-%d-%d %d:%d:%d\n", 
		program, ltime->tm_mon + 1, ltime->tm_mday, ltime->tm_year + 1900, 
		ltime->tm_hour, ltime->tm_min, ltime->tm_sec);
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
	if(out_log == NULL) {
		init_log(__FILE__);
	}
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
	if(out_log == NULL) {
		init_log(file);
	}
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

/********************************************************
 * Test/debug functions
 ********************************************************/

/*
 * Prints the given 2D array of vectors to a text file with the given name
 */
void print_traj(rvec **x, int nframes, int natoms, const char *fname) {
	int fr, i;
	FILE *f = fopen(fname, "w");

	for(fr = 0; fr < nframes; fr++) {
		fprintf(f, "\nFrame %d\n", fr);
		for(i = 0; i < natoms; i++) {
			fprintf(f, "%d: %f %f %f\n", i, x[fr][i][0], x[fr][i][1], x[fr][i][2]);
		}
	}

	fclose(f);
}

/*
 * Prints the given array of vectors to a text file with the given name
 */
void print_vecs(rvec *x, int natoms, const char *fname) {
	int i;
	FILE *f = fopen(fname, "w");

	for(i = 0; i < natoms; i++) {
		fprintf(f, "%d: %f %f %f\n", i, x[i][0], x[i][1], x[i][2]);
	}

	fclose(f);
}

/*
 * Prints reordered trajectory training nodes to a text file with the given name
 */
void print_train_vecs(struct svm_node ***train_vectors, 
	int natoms, double *targets, int n_data, const char *fname) {

	int i, j, k;
	FILE *f = fopen(fname, "w");

	for(i = 0; i < natoms; i++) {
		fprintf(f, "\nAtom %d:\n", i + 1);
		for(j = 0; j < n_data; j++) {
			fprintf(f, "\nData %d Traj %f: ", j, targets[j]);
			for(k = 0; k < 3; k++) {
				fprintf(f, "%d:%f ", train_vectors[i][j][k].index, train_vectors[i][j][k].value);
			}
			fprintf(f, "%d\n", train_vectors[i][j][k].index);
		}
	}

	fclose(f);
}

/*
 * Prints svm_problem data to a text file with the given name
 */
void print_svm_probs(struct svm_problem *probs, int nprobs, const char *fname) {
	int i, j, k;
	FILE *f = fopen(fname, "w");

	for(i = 0; i < nprobs; i++) {
		fprintf(f, "\nProblem %d:\n", i);
		for(j = 0; j < probs[i].l; j++) {
			fprintf(f, "\nData %d Target %f: ", j, probs[i].y[j]);
			for(k = 0; k < 3; k++) {
				fprintf(f, "%d:%f ", probs[i].x[j][k].index, probs[i].x[j][k].value);
			}
			fprintf(f, "%d\n", probs[i].x[j][k].index);
		}
	}

	fclose(f);
}

/*
 * Saves the given svm models in model files and writes their data to a text file with the given name
 */
void save_print_models(struct svm_model **models, int n, const char *fname) {
	char mod_fn[20];
	int i, j, nr_class, nr_sv, *labels, *sv_indices;
	FILE *f = fopen(fname, "w");

	for(i = 0; i < n; i++) {
		fprintf(f, "Model %d:\n", i + 1);
		fprintf(f, "SVM Type: %d\n", svm_get_svm_type(models[i]));
		fprintf(f, "Number of classes: %d\n", nr_class = svm_get_nr_class(models[i]));
		
		snew(labels, nr_class);
		svm_get_labels(models[i], labels);
		fprintf(f, "Labels: ");
		for(j = 0; j < nr_class; j++) {
			fprintf(f, "%d ", labels[j]);
		}
		fprintf(f, "\n");
		sfree(labels);

		fprintf(f, "Number of support vectors: %d\n", nr_sv = svm_get_nr_sv(models[i]));

		snew(sv_indices, nr_sv);
		svm_get_sv_indices(models[i], sv_indices);
		fprintf(f, "SV Indices: ");
		for(j = 0; j < nr_sv; j++) {
			fprintf(f, "%d ", sv_indices[j]);
		}
		fprintf(f, "\n");
		sfree(sv_indices);

		fprintf(f, "\n");

		sprintf(mod_fn, "trajmodel%d", i + 1);
		svm_save_model(mod_fn, models[i]);
	}

	fclose(f);
}

// clock_t start = clock();
// svmanalys(fnames[0], fnames[1], fnames[2], fnames[3]);
// clock_t end = clock();
// print_log("Execution time: %d\n", end - start);