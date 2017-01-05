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

#include "ensemble_comp.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// The number of new frames by which to reallocate an array of length # trajectory frames
#define FRAMESTEP 500

// static void eta_res_pdb(const char *pdb_fname, real *eta, int natoms, eta_res_t *eta_res);
// static void eta_res_tps(const char *tps_file, real *eta, eta_res_t *eta_res);
// static void eta_res_tpx(const char *tpr_fname, real *eta, eta_res_t *eta_res);
// static void f_calc_eta_res(real *eta, t_atoms *atoms, eta_res_t *eta_res);

static FILE *out_log = NULL;

void init_eta_dat(eta_dat_t *eta_dat) {
    eta_dat->gamma = GAMMA;
    eta_dat->c = COST;
    eta_dat->nthreads = -1;
    eta_dat->oenv = NULL;

    eta_dat->natoms = 0;
    eta_dat->atom_nums = NULL;
    eta_dat->eta = NULL;
    eta_dat->nres = 0;
    eta_dat->res_nums = NULL;
    eta_dat->res_names = NULL;
    eta_dat->avg_etas = NULL;
}

void free_eta_dat(eta_dat_t *eta_dat) {
    if (eta_dat->atom_nums) sfree(eta_dat->atom_nums);
    if (eta_dat->eta)       sfree(eta_dat->eta);
    if (eta_dat->res_nums)  sfree(eta_dat->res_nums);
    if (eta_dat->res_names) sfree(eta_dat->res_names);
    if (eta_dat->avg_etas)  sfree(eta_dat->avg_etas);
}


void ensemble_comp(eta_dat_t *eta_dat) {
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
    switch(fn2ftp(eta_dat->fnames[eTRAJ1])) {
        case efXTC:
        case efTRR:
        case efPDB:
            read_traj(eta_dat->fnames[eTRAJ1], &x1, &nframes, &eta_dat->natoms, &eta_dat->oenv);
            break;
        default:
            log_fatal(FARGS, io_error);
    }
    switch(fn2ftp(eta_dat->fnames[eTRAJ2])) {
        case efXTC:
        case efTRR:
        case efPDB:
            read_traj(eta_dat->fnames[eTRAJ2], &x2, &nframes2, &natoms2, &eta_dat->oenv);
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

    /* If an index file was given, get atom group with indices that will be trained */
    if(eta_dat->fnames[eNDX1] != NULL) {
        rd_index(eta_dat->fnames[eNDX1], NUMGROUPS, isize, indx1, grp_names);
        eta_dat->natoms = isize[0];
    }
    else { // If no index file, set default indices as 0 to natoms - 1
        snew(indx1[0], eta_dat->natoms);
        for(i = 0; i < eta_dat->natoms; ++i) {
            indx1[0][i] = i;
        }
    }
    if(eta_dat->fnames[eNDX2] != NULL) {
        snew(isize2, NUMGROUPS);
        snew(indx2, NUMGROUPS);
        rd_index(eta_dat->fnames[eNDX2], NUMGROUPS, isize2, indx2, grp_names);
        if(isize2[0] != eta_dat->natoms) {
            log_fatal(FARGS, ndx_error);
        }
    }
    else {
        if(natoms2 != eta_dat->natoms) {
            log_fatal(FARGS, natom_error);
        }
        indx2 = indx1;
    }
    eta_dat->atom_nums = indx1[0]; // store atom IDs in output

    /* Construct svm problems */
    traj2svm_probs(x1, x2, indx1[0], indx2[0], nframes, eta_dat->natoms, &probs);

    /* No longer need original vectors */
    for(i = 0; i < nframes; ++i) {
        sfree(x1[i]);
        sfree(x2[i]);
    }
    sfree(x1);
    sfree(x2);

    /* No longer need index junk (except for what we stored in atom_nums) */
    sfree(isize);
    sfree(indx1);
    sfree(grp_names);
    if(eta_dat->fnames[eNDX2] != NULL) {
        sfree(isize2);
        sfree(indx2[0]);
        sfree(indx2);
    }

    /* Train SVM */
    snew(models, eta_dat->natoms);
    train_traj(probs, eta_dat->natoms, eta_dat->gamma, eta_dat->c, eta_dat->nthreads, models);

    /* Calculate eta values */
    snew(eta_dat->eta, eta_dat->natoms);
    calc_eta(models, eta_dat->natoms, nframes, eta_dat->eta);

    /* Clean up */
    if(eta_dat->natoms > 0)
        sfree(probs[0].y);
    for(i = 0; i < eta_dat->natoms; ++i) {
        for(int j = 0; j < nframes * 2; ++j) {
            sfree(probs[i].x[j]);
        }
        sfree(probs[i].x);
    }
    sfree(probs);

    for(i = 0; i < eta_dat->natoms; ++i) {
        svm_free_model_content(models[i]);
        svm_free_and_destroy_model(&(models[i]));
    }
    sfree(models);

    /* If residue information given, calculate eta per residue */
    if (eta_dat->fnames[eRES1] != NULL) {
        // TODO
    }
}

void traj2svm_probs(rvec **x1,
                    rvec **x2,
                    atom_id *indx1,
                    atom_id *indx2,
                    int nframes,
                    int natoms,
                    struct svm_problem **probs) {
    int nvecs = nframes * 2;
    int i;
    double *targets; // trajectory classification labels

    print_log("Constructing svm problems...\n");

    /* Build targets array with classification labels */
    snew(targets, nvecs);
    for(i = 0; i < nframes; ++i) {
        targets[i] = LABEL1; // trajectory 1
    }
    for(; i < nvecs; ++i) {
        targets[i] = LABEL2; // trajectory 2
    }

    /* Construct svm problems */
    snew(*probs, natoms);
    int cur_atom, cur_frame, cur_data;
    for(cur_atom = 0; cur_atom < natoms; ++cur_atom) {
        (*probs)[cur_atom].l = nvecs;
        (*probs)[cur_atom].y = targets;
        snew((*probs)[cur_atom].x, nvecs);
        // Insert positions from traj1
        for(cur_frame = 0; cur_frame < nframes; ++cur_frame) {
            snew((*probs)[cur_atom].x[cur_frame], 4); // (4 = 3 xyz pos + 1 for -1 end index)
            for(i = 0; i < 3; ++i) {
                (*probs)[cur_atom].x[cur_frame][i].index = i + 1; // Position components are indexed 1:x, 2:y, 3:z
                (*probs)[cur_atom].x[cur_frame][i].value = x1[cur_frame][indx1[cur_atom]][i] * 10.0; // Scaling by 10 gives more accurate results
            }
            (*probs)[cur_atom].x[cur_frame][i].index = -1; // -1 index marks end of a data vector
        }
        // Insert positions from traj2
        for(cur_frame = 0, cur_data = nframes; cur_frame < nframes; ++cur_frame, ++cur_data) {
            snew((*probs)[cur_atom].x[cur_data], 4);
            for(i = 0; i < 3; ++i) {
                (*probs)[cur_atom].x[cur_data][i].index = i + 1;
                (*probs)[cur_atom].x[cur_data][i].value = x2[cur_frame][indx2[cur_atom]][i] * 10.0;
            }
            (*probs)[cur_atom].x[cur_data][i].index = -1;
        }
    }
}

void train_traj(struct svm_problem *probs,
                int num_probs,
                real gamma,
                real c,
                int nthreads,
                struct svm_model **models) {
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

#ifdef _OPENMP
    if(nthreads > 0)
        omp_set_num_threads(nthreads);
    if(nthreads > 1 || nthreads <= 0)
        print_log("svm training will be parallelized.\n");
#endif

    /* Train svm */
    int i;
#pragma omp parallel for schedule(dynamic) private(i) shared(num_probs,models,probs,param)
    for(i = 0; i < num_probs; ++i) {
    #if defined _OPENMP && defined EC_DEBUG
        print_log("%d threads running svm-train.\n", omp_get_num_threads());
    #endif
        models[i] = svm_train(&(probs[i]), &param);
    }
}

void calc_eta(struct svm_model **models,
              int num_models,
              int num_frames,
              real *eta) {
    int i;

    print_log("Calculating eta values...\n");

    for(i = 0; i < num_models; ++i) {
        eta[i] = 1.0 - svm_get_nr_sv(models[i]) / (2.0 * (real)num_frames);
    }
}

// void calc_eta_res(const char *res_fname,
//                   real *eta,
//                   int natoms,
//                   eta_res_t *eta_res) {
//     print_log("Reading residue info from %s...\n", res_fname);
//     switch(fn2ftp(res_fname)) {
//         case efPDB:
//             eta_res_pdb(res_fname, eta, natoms, eta_res);
//             break;
//         case efGRO: // TODO: try using this for tpr as well, or vice versa?
//             eta_res_tps(res_fname, eta, eta_res);
//             break;
//         case efTPR:
//             eta_res_tpx(res_fname, eta, eta_res);
//             break;
//         default:
//             print_log("%s is not a supported filetype for residue information. Skipping eta residue calculation.\n", res_fname);
//     }
// }

void save_eta(eta_dat_t *eta_dat) {
    int i;
    FILE *f = fopen(eta_dat->fnames[eETA_ATOM], "w");

    print_log("Saving eta values to %s...\n", eta_dat->fnames[eETA_ATOM]);
    fprintf(f, "# ATOMID\tETA\n");
    for(i = 0; i < eta_dat->natoms; ++i) {
        // TODO: it seems that the atom IDs stored in eta_dat->atom_nums
        // are 1 lower than the IDs stored in the index file.
        // Is this normal?
        fprintf(f, "%d\t%f\n", eta_dat->atom_nums[i], eta_dat->eta[i]);
    }

    fclose(f);

    if (eta_dat->avg_etas) {
        // TODO: save residue etas
    }
}

// void save_eta_res(eta_res_t *eta_res,
//                   const char *eta_res_fname) {
//     int i;
//     FILE *f = fopen(eta_res_fname, "w");
//
//     print_log("Saving residue eta values to %s...\n", eta_res_fname);
//     fprintf(f, "# RES\tETA\n");
//     for(i = 0; i < eta_res->nres; ++i) {
//         fprintf(f, "%d%s\t%f\n", eta_res->res_nums[i], eta_res->res_names[i], eta_res->avg_etas[i]);
//     }
//
//     fclose(f);
// }

void read_traj(const char *traj_fname,
               rvec ***x,
               int *nframes,
               int *natoms,
               output_env_t *oenv) {
    t_trxstatus *status = NULL;
    real t;
    matrix box;
    int est_frames = FRAMESTEP;
    *nframes = 0;

    print_log("Reading trajectory %s...\n", traj_fname);

    snew(*x, est_frames);
    *natoms = read_first_x(*oenv, &status, traj_fname, &t, &((*x)[0]), box);

    do {
        ++(*nframes);
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

// static void eta_res_pdb(const char *pdb_fname,
//                         real *eta,
//                         int natoms,
//                         eta_res_t *eta_res) {
//     char title[256];
//     t_atoms atoms;
//     rvec *x;
//
//     atoms.nr = natoms;
//     snew(atoms.atom, natoms);
//     snew(atoms.atomname, natoms);
//     snew(atoms.atomtype, natoms);
//     snew(atoms.atomtypeB, natoms);
//     atoms.nres = natoms;
//     snew(atoms.resinfo, natoms);
//     snew(atoms.pdbinfo, natoms);
//
//     snew(x, natoms);
//
//     read_pdb_conf(pdb_fname, title, &atoms, x, NULL, NULL, FALSE, NULL);
//
//     sfree(x);
//
//     sfree(atoms.atomname);
//     sfree(atoms.atomtype);
//     sfree(atoms.atomtypeB);
//     sfree(atoms.pdbinfo);
//
//     f_calc_eta_res(eta, &atoms, eta_res);
//
//     sfree(atoms.atom);
//     sfree(atoms.resinfo);
// }
//
// // Try res_tpx for gro and tpr instead of this.
// static void eta_res_tps(const char *tps_file,
//                         real *eta,
//                         eta_res_t *eta_res) {
//     char title[256];
//     t_topology top;
//     rvec *x = NULL;
//     matrix box;
//     int ePBC;
//
//     init_top(&top);
//
//     read_tps_conf(tps_file, title, &top, &ePBC, &x, NULL, box, FALSE);
//
//     f_calc_eta_res(eta, &(top.atoms), eta_res);
//
//     // Cannot use done_top(), causes error- pointer being freed was not allocated. See implementation in typedefs.c
//     done_atom(&(top.atoms));
//     done_symtab(&(top.symtab));
//     done_block(&(top.cgs));
//     done_block(&(top.mols));
//     done_blocka(&(top.excls));
//
//     sfree(x);
// }
//
// // TODO: Does this work for gro files generated by grompp etc?
// static void eta_res_tpx(const char *tpr_fname,
//                         real *eta,
//                         eta_res_t *eta_res) {
//     t_inputrec ir;
//     gmx_mtop_t mtop;
//     matrix box;
//     int natoms, i;
//
//     read_tpx(tpr_fname, &ir, box, &natoms, NULL, NULL, NULL, &mtop);
//
//     f_calc_eta_res(eta, &(mtop.moltype->atoms), eta_res);
// }
//
// static void f_calc_eta_res(real *eta,
//                            t_atoms *atoms,
//                            eta_res_t *eta_res) {
//     real *sums;
//     int *n, i;
//
//     print_log("Calculating residue eta values...\n");
//
//     snew(eta_res->res_nums, atoms->nres);
//     snew(eta_res->res_names, atoms->nres);
//     snew(eta_res->avg_etas, atoms->nres);
//
//     snew(sums, atoms->nres);
//     snew(n, atoms->nres);
//
//     // Add up sums
//     for(i = 0; i < atoms->nr; ++i) {
//         sums[atoms->atom[i].resind] += eta[i];
//         ++(n[atoms->atom[i].resind]);
//     }
//
//     // Calculate average etas
//     for(i = 0; i < atoms->nres; ++i) {
//         eta_res->avg_etas[i] = sums[i] / n[i];
//     }
//
//     // Store residue info in eta_res_t struct
//     eta_res->nres = atoms->nres;
//     for(i = 0; i < atoms->nres; ++i) {
//         eta_res->res_nums[i] = atoms->resinfo[i].nr;
//         eta_res->res_names[i] = *(atoms->resinfo[i].name);
//     }
//
//     sfree(sums);
//     sfree(n);
// }

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

void log_fatal(int fatal_errno,
               const char *file,
               int line,
               char const *fmt, ...) {
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
