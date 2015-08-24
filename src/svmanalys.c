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

 * svmanalys analyzes GROMACS trajectory files using support vector machines and calculates eta values.
 * Written by Ahnaf Siddiqui, Mohsen Botlani-Esfahani, and Dr. Sameer Varma.
 * Copyright (c) 2015, University of South Florida.
 * The authors would like to acknowledge the use of the services provided by Research Computing at the University of South Florida.
 */

#include "svmanalys.h"

#define FRAMESTEP 500 // The number of new frames by which to reallocate an array of length # trajectory frames

static void read_traj(const char *traj_fname, rvec ***x, int *nframes, int *natoms);
static void init_log(const char *program);
static void close_log(void);
static void print_log(char const *fmt, ...);
static void log_fatal(int fatal_errno, const char *file, int line, char const *fmt, ...);

static FILE *out_log = NULL;
static output_env_t oenv = NULL;

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"svmanalys analyzes trajectory files produced by GROMACS using support vector machine algorithms.",
		"It calculates eta values for trained atoms from two trajectories."
	};
	const char *fnames[NUMFILES];
	real gamma = GAMMA, c = COST;
	int nfiles, npa;

	init_log(argv[0]);

	t_filenm fnm[] = {
		{efTRX, "-f1", "traj1.xtc", ffREAD},
		{efTRX, "-f2", "traj2.xtc", ffREAD},
		{efNDX, "-n1", "index1.ndx", ffOPTRD},
		{efNDX, "-n2", "index2.ndx", ffOPTRD},
		{efDAT, "-eta_atom", "eta_atom.dat", ffWRITE}
	};

	t_pargs pa[] = {
		{"-g", FALSE, etREAL, {&gamma}, "gamma parameter for svm-train"},
		{"-c", FALSE, etREAL, {&c}, "C (cost) parameter for svm-train"}
	};

	nfiles = asize(fnm);
	npa = asize(pa);

	parse_common_args(&argc, argv, 0, nfiles, fnm, 
		npa, pa, asize(desc), desc, 0, NULL, &oenv);

	fnames[TRAJ1] = opt2fn("-f1", nfiles, fnm);
	fnames[TRAJ2] = opt2fn("-f2", nfiles, fnm);
	fnames[NDX1] = opt2fn_null("-n1", nfiles, fnm);
	fnames[NDX2] = opt2fn_null("-n2", nfiles, fnm);
	fnames[ETA_DAT] = opt2fn("-eta_atom", nfiles, fnm);

	svmanalys(fnames, gamma, c);

	close_log();
}

void svmanalys(const char *fnames[], real gamma, real c) {
	const char *io_error = "Input trajectory files must be .xtc, .trr, or .pdb!\n";
	const char *ndx_error = "Given index groups have different numbers of atoms!\n";

	/* Trajectory data */
	rvec **x1, **x2; // Trajectory position vectors
	int nframes, nframes2, natoms, natoms2, i;

	/* Training data */
	struct svm_problem *probs; // svm problems for training
	struct svm_model **models; // pointers to models produced by training

	/* eta values */
	double *eta;

	/* Read trajectory files */
	switch(fn2ftp(fnames[TRAJ1])) {
		case efXTC:
		case efTRR:
		case efPDB:
			read_traj(fnames[TRAJ1], &x1, &nframes, &natoms);
			break;
		default:
			log_fatal(FARGS, io_error);
	}
	switch(fn2ftp(fnames[TRAJ2])) {
		case efXTC:
		case efTRR:
		case efPDB:
			read_traj(fnames[TRAJ2], &x2, &nframes2, &natoms2);
			break;
		default:
			log_fatal(FARGS, io_error);
	}

	/* In case traj files have different numbers of frames or atoms */
	nframes = nframes2 < nframes ? nframes2 : nframes;
	natoms = natoms2 < natoms ? natoms2 : natoms;

	/* Index data */
#define NUMGROUPS 1
	int *isize, *isize2;
	atom_id **indx1, **indx2; // Atom indices for the two trajectories
	char **grp_names;

	snew(isize, NUMGROUPS);
	snew(indx1, NUMGROUPS);
	snew(grp_names, NUMGROUPS);

	// If an index file was given, get atom group with indices that will be trained
	if(fnames[NDX1] != NULL) {
		rd_index(fnames[NDX1], NUMGROUPS, isize, indx1, grp_names);
		natoms = isize[0];
	}
	else { // If no index file, set default indices as 0 to natoms - 1
		snew(indx1[0], natoms);
		for(i = 0; i < natoms; i++) {
			indx1[0][i] = i;
		}
	}
	if(fnames[NDX2] != NULL) {
		snew(isize2, NUMGROUPS);
		snew(indx2, NUMGROUPS);
		rd_index(fnames[NDX2], NUMGROUPS, isize2, indx2, grp_names);
		if(isize2[0] != natoms) {
			log_fatal(FARGS, ndx_error);
		}
	}
	else {
		indx2 = indx1;
	}
	sfree(grp_names);

	/* Construct svm problems */
	traj2svm_probs(x1, x2, indx1[0], indx2[0], nframes, natoms, &probs);

	/* No longer need original vectors and index data */
	for(i = 0; i < nframes; i++) {
		sfree(x1[i]);
		sfree(x2[i]);
	}
	sfree(x1);
	sfree(x2);

	sfree(isize);
	sfree(indx1[0]);
	sfree(indx1);
	if(fnames[NDX2] != NULL) {
		sfree(isize2);
		sfree(indx2[0]);
		sfree(indx2);
	}

	/* Train SVM */
	snew(models, natoms);
	train_traj(probs, natoms, gamma, c, models);

	/* Calculate eta values */
	snew(eta, natoms);
	calc_eta(models, natoms, nframes, eta);

	print_eta(eta, natoms, fnames[ETA_DAT]);

	/* Clean up */
	sfree(probs); // Don't free the data within probs, will cause error
	sfree(models);
	sfree(eta);
}

void traj2svm_probs(rvec **x1, rvec **x2, atom_id *indx1, atom_id *indx2, int nframes, int natoms, struct svm_problem **probs) {
	int nvecs = nframes * 2;
	int i;
	double *targets; // classification labels (Trajectory -1 or 1)

	/* Build targets array with classification labels */
	snew(targets, nvecs);
	for(i = 0; i < nframes; i++) {
		targets[i] = LABEL1; // trajectory 1
	}
	for(; i < nvecs; i++) {
		targets[i] = LABEL2; // trajectory 2
	}

	/* Construct svm problems */
	snew(*probs, natoms);
	int cur_atom, cur_frame, cur_data;
	for(cur_atom = 0; cur_atom < natoms; cur_atom++) {
		(*probs)[cur_atom].l = nvecs;
		(*probs)[cur_atom].y = targets;
		snew((*probs)[cur_atom].x, nvecs);
		// Insert positions from traj1
		for(cur_frame = 0; cur_frame < nframes; cur_frame++) {
			snew((*probs)[cur_atom].x[cur_frame], 4); // (4 = 3 xyz pos + 1 for -1 end index)
			for(i = 0; i < 3; i++) {
				(*probs)[cur_atom].x[cur_frame][i].index = i; // Position components are indexed 0:x, 1:y, 2:z
				(*probs)[cur_atom].x[cur_frame][i].value = x1[cur_frame][indx1[cur_atom]][i];
			}
			(*probs)[cur_atom].x[cur_frame][i].index = -1; // -1 index marks end of a data vector
		}
		// Insert positions from traj2
		for(cur_frame = 0, cur_data = nframes; cur_frame < nframes; cur_frame++, cur_data++) {
			snew((*probs)[cur_atom].x[cur_data], 4);
			for(i = 0; i < 3; i++) {
				(*probs)[cur_atom].x[cur_data][i].index = i;
				(*probs)[cur_atom].x[cur_data][i].value = x2[cur_frame][indx2[cur_atom]][i];
			}
			(*probs)[cur_atom].x[cur_data][i].index = -1;
		}
	}
}

void train_traj(struct svm_problem *probs, int num_probs, real gamma, real c, struct svm_model **models) {
	struct svm_parameter param; // Parameters used for training

	/* Set svm parameters */
	param.svm_type = C_SVC;
	param.kernel_type = RBF;
	param.gamma = gamma;
	param.cache_size = 100;
	param.eps = 0.001;
	param.C = c;
	param.nr_weight = 0;
	param.shrinking = 1;
	param.probability = 0;

	/* Train svm */
	int i;
	for(i = 0; i < num_probs; i++) {
		models[i] = svm_train(&(probs[i]), &param);
	}
}

void calc_eta(struct svm_model **models, int num_models, int num_frames, double *eta) {
	int i;
	for(i = 0; i < num_models; i++) {
		eta[i] = 1.0 - svm_get_nr_sv(models[i]) / (2.0 * (double)num_frames);
	}
}

void print_eta(double *eta, int n, const char *fname) {
	int i;
	FILE *f = fopen(fname, "w");

	for(i = 0; i < n; i++) {
		fprintf(f, "%f\n", eta[i]);
	}

	fclose(f);
}

static void read_traj(const char *traj_fname, rvec ***x, int *nframes, int *natoms) {
	t_trxstatus *status = NULL;
	real t;
	matrix box;
	int est_frames = FRAMESTEP;
	*nframes = 0;

	snew(*x, est_frames);
	*natoms = read_first_x(oenv, &status, traj_fname, &t, &((*x)[0]), box);

	do {
		(*nframes)++;
		if(*nframes >= est_frames) {
			est_frames += FRAMESTEP;
			srenew(*x, est_frames);
		}
		snew((*x)[*nframes], *natoms);
	} while(read_next_x(oenv, status, &t, *natoms, (*x)[*nframes], box));

	sfree((*x)[*nframes]);
	close_trx(status);
}

/********************************************************
 * Logging functions
 ********************************************************/

/* Opens the logfile and logs initial time/date */
static void init_log(const char *program) {
	out_log = fopen("svmlog.txt", "a");
	
	time_t t = time(NULL);
	struct tm *ltime = localtime(&t);
	fprintf(out_log, "\n%s run: %d-%d-%d %d:%d:%d\n", 
		program, ltime->tm_mon + 1, ltime->tm_mday, ltime->tm_year + 1900, 
		ltime->tm_hour, ltime->tm_min, ltime->tm_sec);
}

/* Closes the logfile */
static void close_log(void) {
	fclose(out_log);
}

/* Prints to both stdout and the logfile */
static void print_log(char const *fmt, ...) {
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
static void log_fatal(int fatal_errno, const char *file, int line, char const *fmt, ...) {
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