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

 * svmanalys analyzes GROMACS trajectory files using support vector machines and calculates eta values.
 * Written by Ahnaf Siddiqui, Mohsen Botlani-Esfahani, and Dr. Sameer Varma.
 * Copyright (c) 2015, University of South Florida.
 * The authors would like to acknowledge the use of the services provided by Research Computing at the University of South Florida.
 */

#ifndef _svmanalys_h
#define _svmanalys_h

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "macros.h"
#ifdef GRO_V5
#include "pargs.h"
#include "trxio.h"
#else
#include "statutil.h"
#endif
#include "gmx_fatal.h"
#include "smalloc.h"
#include "svm.h"

#define LABEL1 -1 // classification label for trajectory 1
#define LABEL2 1 // classification label for trajectory 2
#define GAMMA 0.01 // default gamma parameter for svm_train
#define COST 1e12 // default C parameter for svm_train

/* Indices of filenames */
enum {eTRAJ1, eTRAJ2, eNDX1, eNDX2, eETA_DAT, eNUMFILES};


void svmanalys(const char *fnames[], real gamma, real c);
/* Trains the given trajectories in fnames[TRAJ1] and fnames[TRAJ2] 
 * using support vector machines and calculates their eta values,
 * which are stored in a file with name indicated by fnames[ETA_DAT].
 * See enum above for fnames[]. fnames[NDX1] and/or fnames[NDX2] can be NULL.
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

void calc_eta(struct svm_model **models, int num_models, int num_frames, double *eta);

void print_eta(double *eta, int num_etas, const char *eta_fname);

void svm_prob2file(const struct svm_problem *prob, const char *fname);

#endif
