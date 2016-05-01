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

static void read_res_pdb(const char *pdb_fname, int natoms, t_atoms *atoms);
static void read_top_gro(const char *gro_fname, t_topology *top);
static void read_top_tpr(const char *tpr_fname, gmx_mtop_t *mtop);

static void eta_res_pdb(const char *pdb_fname, real *eta, int natoms, eta_res_t *eta_res);
static void eta_res_gro(const char *tps_file, real *eta, eta_res_t *eta_res);
static void eta_res_tpr(const char *tpr_fname, real *eta, eta_res_t *eta_res);
static void f_calc_eta_res(real *eta, t_atoms *atoms, eta_res_t *eta_res);

static void ilist2svm_probs(t_ilist ilist[], rvec **x1, rvec **x2, int nframes, struct svm_problem **probs);
static void ilist_to_internal_coords(t_ilist ilist[], rvec **x, int nframes, int natoms, const char *int_coords_fname);
static real calc_angle(rvec a, rvec b, rvec c);
static real calc_dihedral(rvec a, rvec b, rvec c, rvec d);

static void free_atoms(t_atoms *atoms);
static void free_topology(t_topology *top);

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
	const int NUMGROUPS = 1;
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
	sfree(probs); // Don't free the data within probs, will cause error. FIX THIS (create a free function)
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

	print_log("Constructing svm problems for cartesian coordinates...\n");

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
				(*probs)[cur_atom].x[cur_frame][i].index = i + 1; // Position components are indexed 1:x, 2:y, 3:z
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

	print_log("svm-training trajectory with gamma = %f and C = %f...\n", gamma, c);

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
	#if defined _OPENMP && defined DEBUG
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

void calc_eta_res(const char *res_fname, real *eta, int natoms, eta_res_t *eta_res) {
	print_log("Reading residue info from %s...\n", res_fname);
	switch(fn2ftp(res_fname)) {
		case efPDB:
			eta_res_pdb(res_fname, eta, natoms, eta_res);
			break;
		case efGRO: // try using this for tpr as well, or vice versa?
			eta_res_gro(res_fname, eta, eta_res);
			break;
		case efTPR:
			eta_res_tpr(res_fname, eta, eta_res);
			break;
		default:
			print_log("%s is not a supported filetype for residue information. Skipping eta residue calculation.\n", res_fname);
	}
}

gmx_bool calc_eta_dihedrals(const char *fnames[], real gamma, real c, gmx_bool parallel, output_env_t *oenv, eta_dihedral_t *eta_dih) {
	rvec **x1, **x2;
	int nframes, natoms, ndih;

	struct svm_problem *probs;
	struct svm_model **models;

	t_ilist *ilist = NULL;
	t_topology top;
	gmx_mtop_t mtop;

	switch(fn2ftp(fnames[eTOP1])) {
		case efGRO:
			{
				read_top_gro(fnames[eTOP1], &top);

				if(top.idef.ntypes > 0) {
					ndih = top.idef.il[F_PDIHS].nr;
					ilist = top.idef.il;
				}
				else {
					print_log("No interaction information found in %s. Skipping dihedral eta.\n", fnames[eTOP1]);
					free_topology(&top);
				}
			}
			break;
		case efTPR:
			{
				read_top_tpr(fnames[eTOP1], &mtop);

				ndih = mtop.moltype->ilist[F_PDIHS].nr; // TODO: there are multiple moltypes, so convert this to t_topology to get them all!
				ilist = mtop.moltype->ilist;
			}
			break;
		default:
			print_log("%s is not a supported filetype for topology information. Skipping dihedral eta.\n", fnames[eTOP1]);
	}

	if(ilist) {
		read_traj(fnames[eTRAJ1], &x1, &nframes, &natoms, oenv);
		read_traj(fnames[eTRAJ2], &x2, &nframes, &natoms, oenv);

		ilist2svm_probs(ilist, x1, x2, nframes, &probs);

		snew(models, ndih);
		train_traj(probs, ndih, gamma, c, parallel, models);

		// Construct eta_dih struct
		eta_dih->ndih = ndih;

		snew(eta_dih->atoms, ndih * 4);
		for(int i = 0; i < ndih; ++i) {
			eta_dih->atoms[4*i] = top.idef.il[F_PDIHS].iatoms[i*5+1];
			eta_dih->atoms[4*i+1] = top.idef.il[F_PDIHS].iatoms[i*5+2];
			eta_dih->atoms[4*i+2] = top.idef.il[F_PDIHS].iatoms[i*5+3];
			eta_dih->atoms[4*i+3] = top.idef.il[F_PDIHS].iatoms[i*5+4];
		}

		switch(fn2ftp(fnames[eTOP1])) {
			case efGRO:
				free_topology(&top);
				break;
			case efTPR:
				done_mtop(&mtop, TRUE);
				break;
		}

		// Calculate eta
		snew(eta_dih->eta, ndih);
		calc_eta(models, ndih, nframes, eta_dih->eta);

		// Free memory
		for(int i = 0; i < nframes; i++) {
			sfree(x1[i]);
			sfree(x2[i]);
		}
		sfree(x1);
		sfree(x2);

		sfree(probs);
		for(int i = 0; i < ndih; i++) {
			svm_free_model_content(models[i]);
			svm_free_and_destroy_model(&(models[i]));
		}
		sfree(models);

		return TRUE;
	}

	return FALSE;
}

/********************************************************
 * I/O functions
 ********************************************************/

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

	sfree((*x)[*nframes]); // Nothing was read to the last allocated frame
	close_trx(status);
}

// Reads residue and atom information from a pdb file
// Allocates memory for data within t_atoms *atoms. Use free_atoms(t_atoms *atoms) to free.
static void read_res_pdb(const char *pdb_fname, int natoms, t_atoms *atoms) {
	char title[256];
	rvec *x;

	atoms->nr = natoms;
	snew(atoms->atom, natoms);
	snew(atoms->atomname, natoms);
	snew(atoms->atomtype, natoms);
	snew(atoms->atomtypeB, natoms);
	atoms->nres = natoms;
	snew(atoms->resinfo, natoms);
	snew(atoms->pdbinfo, natoms);

	snew(x, natoms);

	read_pdb_conf(pdb_fname, title, atoms, x, NULL, NULL, FALSE, NULL);

	sfree(x);
}

// Reads topology information from a gro file
// Allocates memory for data within t_topology *top. Use free_topology(t_topology *top) to free.
static void read_top_gro(const char *gro_fname, t_topology *top) {
	char title[256];
	rvec *x = NULL;
	matrix box;
	int ePBC;

	init_top(top);

	read_tps_conf(gro_fname, title, top, &ePBC, &x, NULL, box, FALSE);

	sfree(x);
}

// Reads topology information from a tpr file
// Allocates memory for data within gmx_mtop_t *mtop.
// Use done_mtop(gmx_mtop_t *mtop,gmx_bool bDoneSymtab) to free memory.
static void read_top_tpr(const char *tpr_fname, gmx_mtop_t *mtop) {
	t_inputrec ir;
	matrix box;
	int natoms;

	read_tpx(tpr_fname, &ir, box, &natoms, NULL, NULL, NULL, mtop);
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

void save_eta_res(eta_res_t *eta_res, const char *eta_res_fname) {
	int i;
	FILE *f = fopen(eta_res_fname, "w");

	print_log("Saving residue eta values to %s...\n", eta_res_fname);
	fprintf(f, "RESIDUE\tETA\n");
	for(i = 0; i < eta_res->nres; i++) {
		fprintf(f, "%d%s\t%f\n", eta_res->res_nums[i], eta_res->res_names[i], eta_res->avg_etas[i]);
	}

	fclose(f);
}

/********************************************************
 * Eta residue functions
 ********************************************************/

static void eta_res_pdb(const char *pdb_fname, real *eta, int natoms, eta_res_t *eta_res) {
	t_atoms atoms;

	read_res_pdb(pdb_fname, natoms, &atoms);

	f_calc_eta_res(eta, &atoms, eta_res);

	free_atoms(&atoms);
}

// Used for gro files
// Try res_tpx for gro and tpr instead of this.
static void eta_res_gro(const char *tps_file, real *eta, eta_res_t *eta_res) {
	t_topology top;

	read_top_gro(tps_file, &top);

	f_calc_eta_res(eta, &(top.atoms), eta_res);

	free_topology(&top);
}

// Used for tpr files
// TODO: Does this work for gro files generated by grompp etc?
static void eta_res_tpr(const char *tpr_fname, real *eta, eta_res_t *eta_res) {
	gmx_mtop_t mtop;

	read_top_tpr(tpr_fname, &mtop);

	f_calc_eta_res(eta, &(mtop.moltype->atoms), eta_res);

	done_mtop(&mtop, TRUE);
}

static void f_calc_eta_res(real *eta, t_atoms *atoms, eta_res_t *eta_res) {
	real *sums;
	int *n, i;

	print_log("Calculating residue eta values...\n");

	snew(eta_res->res_nums, atoms->nres);
	snew(eta_res->res_names, atoms->nres);
	snew(eta_res->avg_etas, atoms->nres);

	snew(sums, atoms->nres);
	snew(n, atoms->nres);

	// Add up sums
	for(i = 0; i < atoms->nr; i++) {
		sums[atoms->atom[i].resind] += eta[i];
		n[atoms->atom[i].resind]++;
	}

	// Calculate average etas
	for(i = 0; i < atoms->nres; i++) {
		eta_res->avg_etas[i] = sums[i] / n[i];
	}

	// Store residue info in eta_res_t struct
	eta_res->nres = atoms->nres;
	for(i = 0; i < atoms->nres; i++) {
		eta_res->res_nums[i] = atoms->resinfo[i].nr;
		eta_res->res_names[i] = *(atoms->resinfo[i].name);
	}

	sfree(sums);
	sfree(n);
}

/********************************************************
 * Eta dihedral functions
 ********************************************************/

static void ilist2svm_probs(t_ilist ilist[], rvec **x1, rvec **x2, int nframes, struct svm_problem **probs) {
	int ndata = nframes * 2;
	int i;
	double *targets;

	int nr_dih = ilist[F_PDIHS].nr;
	t_iatom *iatoms = ilist[F_PDIHS].iatoms, a, b, c, d;

	print_log("Constructing svm problems for trajectory coordinates...\n");

	/* Build targets array with classification labels */
	snew(targets, ndata);
	for(i = 0; i < nframes; i++) {
		targets[i] = LABEL1; // trajectory 1
	}
	for(; i < ndata; i++) {
		targets[i] = LABEL2;
	}

	/* Calculate dihedrals and construct svm problems */
	snew(*probs, nr_dih);
	int cur_dih, cur_frame, cur_data;
	for(cur_dih = 0; cur_dih < nr_dih; cur_dih++) {
		(*probs)[cur_dih].l = ndata;
		(*probs)[cur_dih].y = targets;
		snew((*probs)[cur_dih].x, ndata);
		// Calculate and insert dihedrals from traj1
		for(cur_frame = 0; cur_frame < nframes; cur_frame++) {
			snew((*probs)[cur_dih].x[cur_frame], 2); // (2 = 1 dihedral angle + 1 for -1 end index)
			(*probs)[cur_dih].x[cur_frame][0].index = 1;

			// Calculate dihedral angle for this dihedral group and insert into the problem node
			a = iatoms[cur_dih * 5 + 1], b = iatoms[cur_dih * 5 + 2], c = iatoms[cur_dih * 5 + 3], d = iatoms[cur_dih * 5 + 4];
			(*probs)[cur_dih].x[cur_frame][0].value = calc_dihedral(x1[cur_frame][a], x1[cur_frame][b], x1[cur_frame][c], x1[cur_frame][d]);

			(*probs)[cur_dih].x[cur_frame][1].index = -1; // -1 index marks end of a data vector
		}
		// Calculate and insert dihedrals from traj2
		for(cur_frame = 0, cur_data = nframes; cur_frame < nframes; cur_frame++, cur_data++) {
			snew((*probs)[cur_dih].x[cur_data], 2);
			(*probs)[cur_dih].x[cur_data][0].index = 1;

			a = iatoms[cur_dih * 5 + 1], b = iatoms[cur_dih * 5 + 2], c = iatoms[cur_dih * 5 + 3], d = iatoms[cur_dih * 5 + 4];

			(*probs)[cur_dih].x[cur_data][0].value = calc_dihedral(x2[cur_frame][a], x2[cur_frame][b], x2[cur_frame][c], x2[cur_frame][d]);

			(*probs)[cur_dih].x[cur_data][1].index = -1;
		}
	}
}

static void ilist_to_internal_coords(t_ilist ilist[], rvec **x, int nframes, int natoms, const char *int_coords_fname) {
	int i, j;

	FILE *f = fopen(int_coords_fname, "w");

	int nr_angle = ilist[F_ANGLES].nr;
	t_iatom *iatoms_angle = ilist[F_ANGLES].iatoms;

	int nr_dih = ilist[F_PDIHS].nr;
	t_iatom *iatoms_dihedral = ilist[F_PDIHS].iatoms;

	t_iatom a, b, c, d;

	// for(i = 0; i < nr_dih; i++) {
	// 	print_log("\tiatom %d: %d\n", i, iatoms_dihedral[i]);
	// }

	for(i = 0; i < nframes; i++) {
		fprintf(f, "\nFrame %d ANGLES:\n", i);
		fprintf(f, "Atom #s:\tAngle(rad)\n");
		for(j = 0; j < nr_angle; j+=4) {
			a = iatoms_angle[j+1], b = iatoms_angle[j+2], c = iatoms_angle[j+3];
			fprintf(f, "A %d %d %d: %f\n", a, b, c, calc_angle(x[i][a], x[i][b], x[i][c]));
		}

		fprintf(f, "\nFrame %d DIHEDRALS (using proper dihedrals PDIHS)\n", i);
		fprintf(f, "Atom #s:\tDihedral(rad)\n");
		for(j = 0; j < nr_dih; j+=5) {
			a = iatoms_dihedral[j+1], b = iatoms_dihedral[j+2], c = iatoms_dihedral[j+3], d = iatoms_dihedral[j+4];
			fprintf(f, "D %d %d %d %d: %f\n", a, b, c, d, calc_dihedral(x[i][a], x[i][b], x[i][c], x[i][d]));
		}
	}

	fclose(f);
}

static real calc_angle(rvec a, rvec b, rvec c) {
	rvec ba, bc;

	rvec_sub(a, b, ba);
	rvec_sub(c, b, bc);

	return gmx_angle(ba, bc);
}

static real calc_dihedral(rvec a, rvec b, rvec c, rvec d) {
	rvec b1, b2, b3;
	rvec n1, n2;

	rvec_sub(b, a, b1);
	rvec_sub(c, b, b2);
	rvec_sub(d, c, b3);

	cprod(b1, b2, n1);
	cprod(b2, b3, n2);

	return gmx_angle(n1, n2);
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

/********************************************************
 * Cleanup functions
 ********************************************************/

void free_eta_res(eta_res_t *eta_res) {
	sfree(eta_res->res_nums);
	sfree(eta_res->res_names);
	sfree(eta_res->avg_etas);
}

void free_eta_dih(eta_dihedral_t *eta_dih) {
	sfree(eta_dih->atoms);
	sfree(eta_dih->eta);
}

// Frees the dynamic memory in a t_atoms struct
static void free_atoms(t_atoms *atoms) {
	sfree(atoms->atomname);
	sfree(atoms->atomtype);
	sfree(atoms->atomtypeB);
	sfree(atoms->pdbinfo);

	sfree(atoms->atom);
	sfree(atoms->resinfo);
}

// Frees the dynamic memory in a t_topology struct
static void free_topology(t_topology *top) {
	// Cannot use done_top(), causes error- pointer being freed was not allocated. See implementation in typedefs.c
	done_atom(&(top->atoms));
	done_symtab(&(top->symtab));
	done_block(&(top->cgs));
	done_block(&(top->mols));
	done_blocka(&(top->excls));
}
