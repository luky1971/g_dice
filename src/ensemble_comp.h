/*
 * Copyright 2016 Ahnaf Siddiqui, Mohsen Botlani and Sameer Varma
 *
 * This program uses the GROMACS molecular simulation package API.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed at http://www.gromacs.org.

 * g_ensemble_comp quantifies the difference between two conformational ensembles (two trajectory files)
 * Quantification is in terms of a true metric, eta=1-Overlap
 * Leighty and Varma, Quantifying Changes in Intrinsic Molecular Motion Using Support Vector Machines, J. Chem. Theory Comput. 2013, 9, 868-875.
 */

#ifndef ENSEMBLE_COMP_H
#define ENSEMBLE_COMP_H

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "macros.h"
#include "smalloc.h"
#include "svm.h"
#include "tpxio.h"
#ifdef GRO_V5
#include "atoms.h"
#include "fatalerror.h"
#include "pargs.h"
#include "topology.h"
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

/** Struct for holding eta data */
typedef struct {
    // Input parameters
    // See ensemble_comp function for how they are used.
    const char *fnames[eNUMFILES];
    real gamma;
    real c;
    int nthreads;
    output_env_t oenv;

    // eta output for atoms
    int natoms; // number of atoms
    atom_id *atom_nums; // array of atom IDs, taken from fnames[eNDX1]
    real *eta; // eta value of each atom. array size = natoms

    // eta output for residues
    int nres; // number of residues
    int *res_nums; // array of residue IDs
    const char **res_names; // array of names of the residues
    real *avg_etas; // average eta value of each residue. array size = nres
} eta_dat_t;


void init_eta_dat(eta_dat_t *eta_dat);
/* Initializes an eta_dat_t struct, such as setting pointers to NULL and setting default parameters.
 */

void free_eta_dat(eta_dat_t *eta_dat);
/* Frees the dynaimc memory in an eta_dat_t struct.
 */

void ensemble_comp(eta_dat_t *eta_dat);
/* Projects the coordinates in fnames[eTRAJ1] and fnames[eTRAJ2]
 * in the Hilbert space specified by C and gamma and calculates discriminability (eta).
 * See enum above for fnames[]. They correspond to the command-line file options described in the README.
 * fnames[eNDX1] and/or fnames[eNDX2] can be NULL.
 * Memory is allocated for the eta array.
 * natoms will hold the number of atoms discriminated by the function.
 * nthreads is the number of threads to be used if ensemble_comp was built using openmp.
 * nthreads <= 0 will use all available threads.
 * output_env_t *oenv is needed for reading trajectory files.
 * You can initialize one using output_env_init() in Gromacs's oenv.h.
 *
 * If fnames[eRES1] is not NULL, will calculate average discriminability (eta)
 * per residue using the residue info in the given eRES1 file.
 */

void traj2svm_probs(rvec **x1,
                    rvec **x2,
                    atom_id *indx1,
                    atom_id *indx2,
                    int nframes,
                    int natoms,
                    struct svm_problem **probs);
/* Constructs svm problems from the given position vectors.
 * Memory is allocated for the probs array.
 * One problem is generated per atom, containing its positions in all the frames of x1 and then x2.
 * indx1 and indx2 indicate the indices of the atoms in x1 and x2, respectively,
 * that should be included in probs.
 */

void train_traj(struct svm_problem *probs,
                int num_probs,
                real gamma,
                real c,
                int nthreads,
                struct svm_model **models);
/* Calls libsvm's svm_train function with default parameters and given gamma and c parameters.
 * You can use traj2svm_probs to generate svm_problems.
 * Results are stored in models.
 * Memory for models must be pre-allocated as an array of pointers with length = num_probs.
 * nthreads is the number of threads to be used if ensemble_comp was built using openmp.
 * nthreads <= 0 will use all available threads.
 */

void calc_eta(struct svm_model **models,
              int num_models,
              int num_frames,
              real *eta);
/* Calculates discriminability (eta) values from the number of support vectors in the given svm models
 * and the given number of frames in each trajectory.
 */

// void calc_eta_res(const char *res_fname,
//                   real *eta,
//                   int natoms,
//                   eta_res_t *eta_res);
/* Calculates average discriminability (eta) per residue using residue info in given file
 * and the eta values per atom in the given eta array.
 * Supported file formats include pdb and gro (tpr support is experimental).
 * Memory is allocated for the arrays in eta_res.
 */

void save_eta(eta_dat_t *eta_dat);
/* Saves the given discriminability (eta) values in a text file with the given name.
 */

// void save_eta_res(eta_res_t *eta_res,
//                   const char *eta_res_fname);
/* Saves the given discriminability (eta) residue values in a text file with the given name.
 */

void read_traj(const char *traj_fname,
               rvec ***x,
               int *nframes,
               int *natoms,
               output_env_t *oenv);
/* Reads a trajectory file. rvec **x is position coordinates indexed x[frame #][atom #]
 * 2D memory is allocated for x.
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

void log_fatal(int fatal_errno,
               const char *file,
               int line,
               char const *fmt, ...);
/* Logs fatal error to logfile and also calls gmx_fatal
 * Hint: Use FARGS for first 3 arguments.
 */

#endif // ENSEMBLE_COMP_H
