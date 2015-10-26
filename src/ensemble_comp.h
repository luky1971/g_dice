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

#ifndef _ensemble_comp_h
#define _ensemble_comp_h

#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "macros.h"
#include "smalloc.h"
#include "svm.h"
#include "tpxio.h" // testing
#ifdef GRO_V5
#include "fatalerror.h"
#include "pargs.h"
#include "trxio.h"
#else
#include "gmx_fatal.h"
#include "statutil.h"
#endif

#define LABEL1 -1 // classification label for trajectory 1
#define LABEL2 1 // classification label for trajectory 2
#define GAMMA 0.4 // default gamma parameter for svm_train
#define COST 100.0 // default C parameter for svm_train

/* Indices of filenames */
enum {eTRAJ1, eTRAJ2, eNDX1, eNDX2, eRES1, eETA_ATOM, eETA_RES, eNUMFILES};

/* Struct for holding residue eta data */
typedef struct {
	int nres; // number of residues
	int *res_nums; // array of residue numbers
	const char **res_names; // array of names of the residues
	real *avg_etas; // array of the average eta value of each residue
} eta_res_t;


void ensemble_comp(const char *fnames[], real gamma, real c, 
	real **eta, int *natoms, gmx_bool parallel, output_env_t *oenv);
/* Projects the coordinates in fnames[eTRAJ1] and fnames[eTRAJ2] 
 * in the Hilbert space specified by C and gamma  and calculates discriminability.
 * See enum above for fnames[]. They correspond to the command-line file options described in the README.
 * fnames[eNDX1] and/or fnames[eNDX2] can be NULL.
 * fnames[eETA_ATOM] is not used in this function, and can be NULL or unspecified.
 * Memory is allocated for the eta array. sfree() it after you're done.
 * natoms will hold the number of atoms discriminated by the function.
 * The value of parallel controls whether the svm training routine is parallelized. Recommended value is TRUE. 
 * output_env_t *oenv is needed for reading trajectory files.
 * You can initialize an output_env_t using Gromacs' parse_common_args() function.
 */

void traj2svm_probs(rvec **x1, rvec **x2, atom_id *indx1, atom_id *indx2, 
	int nframes, int natoms, struct svm_problem **probs);
/* Constructs svm problems from the given position vectors.
 * Memory is allocated for the probs array, which can be freed (with sfree()) after use, but don't free data within probs.
 * One problem is generated per atom, containing its positions in all the frames of x1 and then x2.
 * indx1 and indx2 indicate the indices of the atoms in x1 and x2, respectively, that should be included in probs.
 */

void train_traj(struct svm_problem *probs, int num_probs, real gamma, real c, 
	gmx_bool parallel, struct svm_model **models);
/* Calls libsvm's svm_train function with default parameters and given gamma and c parameters.
 */

void calc_eta(struct svm_model **models, int num_models, int num_frames, real *eta);
/* Calculates eta values from the number of support vectors in the given svm models
 * and the given number of frames in each trajectory.
 */

void save_eta(real *eta, int num_etas, const char *eta_fname);
/* Saves the given eta values in a data file with the given name.
 */

void read_traj(const char *traj_fname, rvec ***x, int *nframes, int *natoms, output_env_t *oenv);
/* Reads a trajectory file.
 * 2D memory is allocated for x. sfree() it when you're done.
 */

void res_pdb(const char *pdb_fname, real *eta, int natoms, eta_res_t *eta_res);
/*
 * Memory is allocated for the arrays in eta_res. sfree() them when you're done.
 */

void res_tps(const char *tps_file);
/*
 * TEST
 */

void res_tpx(const char *tpr_fname);
/*
 * WIP
 */

void save_res_eta(eta_res_t *eta_res, const char *eta_res_fname);

/********************************************************
 * Logging functions
 ********************************************************/

void init_log(const char *logfile, const char *program);
/* Opens the logfile and logs initial time/date. Remember to close_log() at end of program.
 */

void close_log();
/* Closes the logfile.
 */

void print_log(char const *fmt, ...);
/* Prints to both stdout and the logfile.
 */

void log_fatal(int fatal_errno, const char *file, int line, char const *fmt, ...);
/* Logs fatal error to logfile and also calls gmx_fatal
 * Hint: Use FARGS for first 3 arguments.
 */

#endif
