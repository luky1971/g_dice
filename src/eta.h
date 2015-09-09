/*
 * This file is part of the GROMACS molecular simulation package.
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

 * etanalys analyzes GROMACS trajectory files using support vector machines and calculates eta values.
 * Written by Ahnaf Siddiqui, Mohsen Botlani-Esfahani, and Dr. Sameer Varma.
 * Copyright (c) 2015, University of South Florida.
 * The authors would like to acknowledge the use of the services provided by Research Computing at the University of South Florida.
 */

#ifndef _eta_h
#define _eta_h

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "macros.h"
#include "smalloc.h"
#include "svm.h"
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
#define GAMMA 0.01 // default gamma parameter for svm_train
#define COST 10.0 // default C parameter for svm_train

/* Indices of filenames */
enum {eTRAJ1, eTRAJ2, eNDX1, eNDX2, eETA_ATOM, eNUMFILES};

void etanalys(const char *fnames[], real gamma, real c, real **eta, int *natoms, output_env_t oenv);
/* Trains the given trajectories in fnames[eTRAJ1] and fnames[eTRAJ2] 
 * using support vector machines and calculates their eta values.
 * See enum above for fnames[]. fnames[eNDX1] and/or fnames[eNDX2] can be NULL.
 * fnames[eETA_ATOM] is not used in this function, and can be NULL or unspecified.
 * Memory is allocated for the eta array.
 * output_env_t oenv is needed for reading trajectory files.
 */

void traj2svm_probs(rvec **x1, rvec **x2, atom_id *indx1, atom_id *indx2, 
	int nframes, int natoms, struct svm_problem **probs);
/* Constructs svm problems from the given position vectors.
 * Memory is allocated for the probs array, which can be freed after use, but don't free data within probs.
 * One problem is generated per atom, containing its positions in all the frames of x1 and then x2.
 * indx1 and indx2 indicate the indices of the atoms in x1 and x2, respectively, that should be included in probs.
 * x1 and x2 should each be size [nframes][natoms]. indx1 and indx2 should each be size [natoms].
 */

void train_traj(struct svm_problem *probs, int num_probs, 
	real gamma, real c, struct svm_model **models);
/* Calls libsvm's svm_train function with default parameters and given gamma and c parameters.
 */

void calc_eta(struct svm_model **models, int num_models, int num_frames, real *eta);
/* Calculates eta values from the number of support vectors in the given svm models
 * and the given number of frames in each trajectory.
 */

void save_eta(real *eta, int num_etas, const char *eta_fname);
/* Saves the given eta values in a data file with the given name.
 */

void read_traj(const char *traj_fname, rvec ***x, int *nframes, int *natoms, output_env_t oenv);
/* Reads a trajectory file.
 */

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