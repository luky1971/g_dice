/*
 * This program uses the GROMACS molecular simulation package API.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.

 * g_ensemble_comp quantifies the difference between two conformational ensembles (two trajectory files) 
 * Quantification is in terms of a true metric, eta=1-Overlap 
 * Leighty and Varma, Quantifying Changes in Intrinsic Molecular Motion Using Support Vector Machines, J. Chem. Theory Comput. 2013, 9, 868-875. 
 * Written by Ahnaf Siddiqui, Mohsen Botlani-Esfahani and Sameer Varma.
 */

#include "ensemble_comp.h"

#define FRAMESTEP 500 // The number of new frames by which to reallocate an array of length # trajectory frames

#define NATOMS 20

static FILE *out_log = NULL;

void ensemble_comp(const char *fnames[], real gamma, real c, 
	real **eta, int *natoms, gmx_bool parallel, output_env_t *oenv) {
	const char *io_error = "Input trajectory files must be .xtc, .trr, or .pdb!\n";
	const char *fr_error = "Input trajectories have differing numbers of frames!\n";
	const char *ndx_error = "Given index groups have differing numbers of atoms!\n";
	const char *natom_error = "Input trajectories have differing numbers of atoms!\n";

	/* Trajectory data */
	rvec **x1, **x2; // Trajectory position vectors
	int nframes, nframes2, natoms2, i;

	/* Training data */
	struct svm_problem *probs; // svm problems for training
	struct svm_model **models; // pointers to models produced by training

	/* Read trajectory files */
	switch(fn2ftp(fnames[eTRAJ1])) {
		case efXTC:
		case efTRR:
		case efPDB:
			read_traj(fnames[eTRAJ1], &x1, &nframes, natoms, oenv);
			break;
		default:
			log_fatal(FARGS, io_error);
	}
	switch(fn2ftp(fnames[eTRAJ2])) {
		case efXTC:
		case efTRR:
		case efPDB:
			read_traj(fnames[eTRAJ2], &x2, &nframes2, &natoms2, oenv);
			break;
		default:
			log_fatal(FARGS, io_error);
	}

	/* In case traj files have different numbers of frames or atoms */
	if(nframes != nframes2) {
		log_fatal(FARGS, fr_error);
	}

	/* Index data */
#define NUMGROUPS 1
	int *isize, *isize2;
	atom_id **indx1, **indx2; // Atom indices for the two trajectories
	char **grp_names;

	snew(isize, NUMGROUPS);
	snew(indx1, NUMGROUPS);
	snew(grp_names, NUMGROUPS);

	// If an index file was given, get atom group with indices that will be trained
	if(fnames[eNDX1] != NULL) {
		rd_index(fnames[eNDX1], NUMGROUPS, isize, indx1, grp_names);
		*natoms = isize[0];
	}
	else { // If no index file, set default indices as 0 to natoms - 1
		snew(indx1[0], *natoms);
		for(i = 0; i < *natoms; i++) {
			indx1[0][i] = i;
		}
	}
	if(fnames[eNDX2] != NULL) {
		snew(isize2, NUMGROUPS);
		snew(indx2, NUMGROUPS);
		rd_index(fnames[eNDX2], NUMGROUPS, isize2, indx2, grp_names);
		if(isize2[0] != *natoms) {
			log_fatal(FARGS, ndx_error);
		}
	}
	else {
		if(natoms2 != *natoms) {
			log_fatal(FARGS, natom_error);
		}
		indx2 = indx1;
	}
	sfree(grp_names);

	/* Construct svm problems */
	traj2svm_probs(x1, x2, indx1[0], indx2[0], nframes, *natoms, &probs);

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
	if(fnames[eNDX2] != NULL) {
		sfree(isize2);
		sfree(indx2[0]);
		sfree(indx2);
	}

	/* Train SVM */
	snew(models, *natoms);
	train_traj(probs, *natoms, gamma, c, parallel, models);

	/* Calculate eta values */
	snew(*eta, *natoms);
	calc_eta(models, *natoms, nframes, *eta);

	/* Clean up */
	sfree(probs); // Don't free the data within probs, will cause error
	for(i = 0; i < *natoms; i++) {
		svm_free_model_content(models[i]);
		svm_free_and_destroy_model(&(models[i]));
	}
	sfree(models);
}

void traj2svm_probs(rvec **x1, rvec **x2, atom_id *indx1, atom_id *indx2, int nframes, int natoms, struct svm_problem **probs) {
	int nvecs = nframes * 2;
	int i;
	double *targets; // trajectory classification labels

	print_log("Constructing svm problems...\n");

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
				(*probs)[cur_atom].x[cur_frame][i].index = i + 1; // Position components are indexed 0:x, 1:y, 2:z
				(*probs)[cur_atom].x[cur_frame][i].value = x1[cur_frame][indx1[cur_atom]][i] * 10.0; // Scaling by 10 gives more accurate results
			}
			(*probs)[cur_atom].x[cur_frame][i].index = -1; // -1 index marks end of a data vector
		}
		// Insert positions from traj2
		for(cur_frame = 0, cur_data = nframes; cur_frame < nframes; cur_frame++, cur_data++) {
			snew((*probs)[cur_atom].x[cur_data], 4);
			for(i = 0; i < 3; i++) {
				(*probs)[cur_atom].x[cur_data][i].index = i + 1;
				(*probs)[cur_atom].x[cur_data][i].value = x2[cur_frame][indx2[cur_atom]][i] * 10.0;
			}
			(*probs)[cur_atom].x[cur_data][i].index = -1;
		}
	}
}

void train_traj(struct svm_problem *probs, int num_probs, real gamma, real c, 
	gmx_bool parallel, struct svm_model **models) {
	struct svm_parameter param; // Parameters used for training

	print_log("svm-training trajectory atoms with gamma = %f and C = %f...\n", gamma, c);

	/* Set svm parameters */
	param.svm_type = C_SVC;
	param.kernel_type = RBF;
	param.degree = 3;
	param.gamma = gamma;
	param.coef0 = 0.0;
	param.cache_size = 100.0;
	param.eps = 0.001;
	param.C = c;
	param.nr_weight = 0;
	param.nu = 0.5;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;

	/* Train svm */
	int i;
#pragma omp parallel for if(parallel) schedule(dynamic) private(i) shared(num_probs,models,probs,param)
	for(i = 0; i < num_probs; i++) {
	#ifdef _OPENMP
		print_log("%d threads running svm-train.\n", omp_get_num_threads());
	#endif
		models[i] = svm_train(&(probs[i]), &param);
	}
}

void calc_eta(struct svm_model **models, int num_models, int num_frames, real *eta) {
	int i;

	print_log("Calculating eta values...\n");

	for(i = 0; i < num_models; i++) {
		eta[i] = 1.0 - svm_get_nr_sv(models[i]) / (2.0 * (real)num_frames);
	}
}

void save_eta(real *eta, int num_etas, const char *eta_fname) {
	int i;
	FILE *f = fopen(eta_fname, "w");

	print_log("Saving eta values to %s...\n", eta_fname);
	fprintf(f, "INDEX\tETA\n");
	for(i = 0; i < num_etas; i++) {
		fprintf(f, "%d\t%f\n", i+1, eta[i]);
	}

	fclose(f);
}

void read_traj(const char *traj_fname, rvec ***x, int *nframes, int *natoms, output_env_t *oenv) {
	t_trxstatus *status = NULL;
	real t;
	matrix box;
	int est_frames = FRAMESTEP;
	*nframes = 0;

	print_log("Reading trajectory %s...\n", traj_fname);

	snew(*x, est_frames);
	*natoms = read_first_x(*oenv, &status, traj_fname, &t, &((*x)[0]), box);

	do {
		(*nframes)++;
		if(*nframes >= est_frames) {
			est_frames += FRAMESTEP;
			srenew(*x, est_frames);
		}
		snew((*x)[*nframes], *natoms);
	} while(read_next_x(*oenv, status, &t,
#ifndef GRO_V5 
		*natoms,
#endif
		(*x)[*nframes], box));

	sfree((*x)[*nframes]);
	close_trx(status);
}

void read_pdb(const char *pdb_fname) {
	char title[256];
	t_atoms atoms;
	rvec x[NATOMS];
	int i;

	atoms.nr = NATOMS;
	snew(atoms.atom, NATOMS);
	snew(atoms.atomname, NATOMS);
	snew(atoms.atomtype, NATOMS);
	snew(atoms.atomtypeB, NATOMS);
	atoms.nres = NATOMS;
	snew(atoms.resinfo, NATOMS);
	snew(atoms.pdbinfo, NATOMS);

	read_pdb_conf(pdb_fname, title, &atoms, x, NULL, NULL, FALSE, NULL);

	for(i = 0; i < atoms.nr; i++) {
		print_log("Atom %d: Residue index %d\n", i+1, atoms.atom[i].resind);
	}

	for(i = 0; i < atoms.nres; i++) {
		print_log("Residue %d: %s\n", atoms.resinfo[i].nr, *(atoms.resinfo[i].name));
	}

	sfree(atoms.atom);
	sfree(atoms.atomname);
	sfree(atoms.atomtype);
	sfree(atoms.atomtypeB);
	sfree(atoms.resinfo);
	sfree(atoms.pdbinfo);
}

// Kill this one
// void read_stx(const char *stx_file) {
// 	char title[256];
// 	t_topology top;
// 	rvec *x;
// 	matrix box; // temp

// 	atoms.nr = NATOMS;
// 	snew(atoms.atom, NATOMS);
// 	snew(atoms.atomname, NATOMS);
// 	snew(atoms.atomtype, NATOMS);
// 	snew(atoms.atomtypeB, NATOMS);
// 	atoms.nres = NATOMS;
// 	snew(atoms.resinfo, NATOMS);
// 	snew(atoms.pdbinfo, NATOMS);

// 	read_tps_conf(stx_file, title, &top, NULL, &x, NULL, box, FALSE);

// 	// sfree(x);
// }

void res_tpx(const char *tpr_fname) {
	// t_tpxheader header;
	t_inputrec ir;
	gmx_mtop_t mtop;
	matrix box;
	// int version, generation;
	int natoms, i;

	// read_tpxheader(tpr_fname, &header, TRUE, &version, &generation);

	// print_log("Topology version: %d\n", version);
	// print_log("Topology generation: $d\n", generation);

	// print_log("\nTPX Header:\n");
	// print_log("bIR: %d\n", header.bIr);
	// print_log("bBox: %d\n", header.bBox);
	// print_log("bTop: %d\n", header.bTop);
	// print_log("bX: %d\n", header.bX);
	// print_log("bV: %d\n", header.bV);
	// print_log("bF: %d\n", header.bF);
	// print_log("natoms: %d\n", header.natoms);
	// print_log("ngtc: %d\n", header.ngtc);
	// print_log("lambda: %f\n", header.lambda);

	// if(header.bX)	snew(*x, header.natoms);
	// else			*x = NULL;

	// if(header.bV)	snew(*v, header.natoms);
	// else			*v = NULL;

	// if(header.bF)	snew(*f, header.natoms);
	// else			*f = NULL;

	read_tpx(tpr_fname, &ir, box, &natoms, NULL, NULL, NULL, &mtop);

	// See inputrec.h
	// print_log("\nInput record:\n");
	// print_log("Integration method: %d\n", ir.eI);
	// print_log("nsteps: %d\n", ir.nsteps);
	// print_log("Simulation part: %d\n", ir.simulation_part);
	// print_log("Init step: %d\n", ir.init_step);
	// print_log("nstcalcenergy: %d\n", ir.nstcalcenergy);
	// print_log("ns_type: %d\n", ir.ns_type);
	// print_log("nstlist: %d\n", ir.nstlist);
	// print_log("ndelta: %d\n", ir.ndelta);
	// print_log("nstcomm: %d\n", ir.nstcomm);
	// print_log("comm_mode: %d\n", ir.comm_mode);
	// print_log("nstcheckpoint: %d\n", ir.nstcheckpoint);
	// print_log("nstlog: %d\n", ir.nstlog);
	// print_log("nstxout: %d\n", ir.nstxout);
	// print_log("nstvout: %d\n", ir.nstvout);
	// print_log("nstfout: %d\n", ir.nstfout);
	// print_log("nstenergy: %d\n", ir.nstenergy);
	// print_log("nstxtcout: %d\n", ir.nstxtcout);
	// print_log("init_t: %f\n", ir.init_t);
	// print_log("delta_t: %f\n", ir.delta_t);

	print_log("\nTopology:\n");
	print_log("Name: %s\n", mtop.name[0]);
	print_log("ffparams:\n");
		print_log("\tntypes: %d\n", mtop.ffparams.ntypes);
		print_log("\tatnr: %d\n", mtop.ffparams.atnr);
		print_log("\treppow: %f\n", mtop.ffparams.reppow);
		print_log("\tfudgeQQ: %f\n", mtop.ffparams.fudgeQQ);
	print_log("Molecular types:\n");
	for(i = 0; i < mtop.nmoltype; i++) {
		print_log("\tName: %s\n", mtop.moltype[i].name[0]);
		print_log("\tNumber of atoms: %d\n", mtop.moltype[i].atoms.nr);
		print_log("\tNumber of charge groups: %d\n", mtop.moltype[i].cgs.nr);
	}
	print_log("Molecule blocks:\n");
	for(i = 0; i < mtop.nmolblock; i++) {
		print_log("\tType index: %d\n", mtop.molblock[i].type);
		print_log("\tNumber of molecules: %d\n", mtop.molblock[i].nmol);
		print_log("\tNumber of atoms in one molecule: %d\n", mtop.molblock[i].natoms_mol);
	}
	print_log("Residue numbering parameter: %d\n", mtop.maxres_renum);
	print_log("Maximum residue number in moltype: %d\n", mtop.maxresnr);
	print_log("Number of atom types: %d\n", mtop.atomtypes.nr);
	print_log("Number of molecules: %d\n", mtop.mols.nr);
	print_log("Number of groups: %d\n", mtop.groups.ngrpname);
	print_log("Number of symbol buffers: %d\n", mtop.symtab.nr);
}

/********************************************************
 * Logging functions
 ********************************************************/

void init_log(const char *logfile, const char *program) {
	out_log = fopen(logfile, "a");
	
	time_t t = time(NULL);
	struct tm *ltime = localtime(&t);
	fprintf(out_log, "\n%s run: %d-%d-%d %d:%d:%d\n", 
		program, ltime->tm_mon + 1, ltime->tm_mday, ltime->tm_year + 1900, 
		ltime->tm_hour, ltime->tm_min, ltime->tm_sec);
}

void close_log() {
	fclose(out_log);
}

void print_log(char const *fmt, ...) {
	va_list arg;
	va_start(arg, fmt);
	vprintf(fmt, arg);
	va_end(arg);
	if(out_log == NULL) {
		init_log("etalog.txt", __FILE__);
	}
	if(out_log != NULL) {
		va_start(arg, fmt);
		vfprintf(out_log, fmt, arg);
		va_end(arg);
	}
}

void log_fatal(int fatal_errno, const char *file, int line, char const *fmt, ...) {
	va_list arg;
	if(out_log == NULL) {
		init_log("etalog.txt", file);
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
