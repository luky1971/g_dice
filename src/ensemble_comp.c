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

static void free_svm_model(struct svm_model *model);
static void eta_res_pdb(eta_dat_t *eta_dat);
static void eta_res_tps(eta_dat_t *eta_dat);
static void eta_res_tpx(eta_dat_t *eta_dat);
static void f_calc_eta_res(eta_dat_t *eta_dat, t_atoms *atoms);

static FILE *out_log = NULL;

void init_eta_dat(eta_dat_t *eta_dat) {
    eta_dat->gamma = GAMMA;
    eta_dat->c = COST;
    eta_dat->nthreads = -1;
    eta_dat->oenv = NULL;

    eta_dat->natoms = 0;
    eta_dat->atom_IDs = NULL;
    eta_dat->eta = NULL;
    eta_dat->nres = 0;
    eta_dat->res_IDs = NULL;
    eta_dat->res_names = NULL;
    eta_dat->res_natoms = NULL;
    eta_dat->res_eta = NULL;

    eta_dat->natoms_all = 0;
}

void free_eta_dat(eta_dat_t *eta_dat) {
    if (eta_dat->atom_IDs)   sfree(eta_dat->atom_IDs);
    if (eta_dat->eta)        sfree(eta_dat->eta);
    if (eta_dat->res_IDs)    sfree(eta_dat->res_IDs);
    if (eta_dat->res_names)  sfree(eta_dat->res_names);
    if (eta_dat->res_natoms) sfree(eta_dat->res_natoms);
    if (eta_dat->res_eta)    sfree(eta_dat->res_eta);
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

    /* In case traj files have different numbers of frames */
    if (nframes != nframes2) {
        log_fatal(FARGS, fr_error);
    }

    // Save total natoms before it is potentially changed by index data below.
    // Might be needed, for example, by residue reading functions. */
    eta_dat->natoms_all = eta_dat->natoms;

    /* Index data */
    const int NUMGROUPS = 1;
    int *isize, *isize2;
    atom_id **indx1, **indx2; // Atom indices for the two trajectories
    char **grp_names;

    snew(isize, NUMGROUPS);
    snew(indx1, NUMGROUPS);
    snew(grp_names, NUMGROUPS);

    /* If an index file was given, get atom group with indices that will be trained */
    if (eta_dat->fnames[eNDX1] != NULL) {
        rd_index(eta_dat->fnames[eNDX1], NUMGROUPS, isize, indx1, grp_names);
        eta_dat->natoms = isize[0];
    }
    else { // If no index file, set default indices as 0 to natoms - 1
        snew(indx1[0], eta_dat->natoms);
        for (i = 0; i < eta_dat->natoms; ++i) {
            indx1[0][i] = i;
        }
    }
    if (eta_dat->fnames[eNDX2] != NULL) {
        snew(isize2, NUMGROUPS);
        snew(indx2, NUMGROUPS);
        rd_index(eta_dat->fnames[eNDX2], NUMGROUPS, isize2, indx2, grp_names);
        if (isize2[0] != eta_dat->natoms) {
            log_fatal(FARGS, ndx_error);
        }
    }
    else {
        if (natoms2 != eta_dat->natoms) {
            log_fatal(FARGS, natom_error);
        }
        indx2 = indx1;
    }
    eta_dat->atom_IDs = indx1[0]; // store atom IDs in output

    /* Construct svm problems */
    traj2svm_probs(x1, x2, indx1[0], indx2[0], nframes, eta_dat->natoms, &probs);

    /* No longer need original vectors */
    free_traj(x1, nframes);
    free_traj(x2, nframes);

    /* No longer need index junk (except for what we stored in atom_IDs) */
    sfree(isize);
    sfree(indx1);
    sfree(grp_names);
    if (eta_dat->fnames[eNDX2] != NULL) {
        sfree(isize2);
        sfree(indx2[0]);
        sfree(indx2);
    }

    /* Train SVM */
    snew(models, eta_dat->natoms);
    train_svm_probs(probs, eta_dat->natoms, eta_dat->gamma, eta_dat->c, eta_dat->nthreads, models);

    /* Calculate eta values */
    snew(eta_dat->eta, eta_dat->natoms);
    calc_eta(models, eta_dat->natoms, nframes, eta_dat->eta);

    /* Clean up svm stuff */
    free_svm_probs(probs, eta_dat->natoms, nframes * 2);
    free_svm_models(models, eta_dat->natoms);

    /* If residue information given, calculate eta per residue */
    if (eta_dat->fnames[eRES1] != NULL) {
        calc_eta_res(eta_dat);
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
    double *targets = NULL; // trajectory classification labels
    struct svm_node *nodepool = NULL; // allocated memory for storing svm nodes (feature vectors)

    print_log("Constructing svm problems for %d atoms in %d frames...\n",
        natoms, nframes);
    flush_log();

    /* Build targets array with classification labels */
    snew(targets, nvecs);
    for (i = 0; i < nframes; ++i) {
        targets[i] = LABEL1; // trajectory 1
    }
    for (; i < nvecs; ++i) {
        targets[i] = LABEL2; // trajectory 2
    }

    /* Allocate enough space for storing all svm nodes */
    snew(nodepool, 2 * natoms * nframes * 4);
    if (!nodepool)
        log_fatal(FARGS, "Failed to allocate memory for svm training vectors!\n");

    /* Construct svm problems */
    snew(*probs, natoms);
    int cur_atom, cur_frame, cur_data;
    for (cur_atom = 0; cur_atom < natoms; ++cur_atom) {
        printf("Atom %d...\r", cur_atom);
        fflush(stdout);

        (*probs)[cur_atom].l = nvecs;
        (*probs)[cur_atom].y = targets;
        snew((*probs)[cur_atom].x, nvecs);
        // Insert positions from traj1
        for (cur_frame = 0; cur_frame < nframes; ++cur_frame) {
            // snew((*probs)[cur_atom].x[cur_frame], 4); // (4 = 3 xyz pos + 1 for -1 end index)
            (*probs)[cur_atom].x[cur_frame] = nodepool;
            for (i = 0; i < 3; ++i) {
                (*probs)[cur_atom].x[cur_frame][i].index = i + 1; // Position components are indexed 1:x, 2:y, 3:z
                (*probs)[cur_atom].x[cur_frame][i].value = x1[cur_frame][indx1[cur_atom]][i] * 10.0; // Scaling by 10 gives more accurate results
            }
            (*probs)[cur_atom].x[cur_frame][i].index = -1; // -1 index marks end of a data vector
            nodepool += 4; // (3 nodes for xyz pos + 1 for -1 end index)
        }
        // Insert positions from traj2
        for (cur_frame = 0, cur_data = nframes; cur_frame < nframes; ++cur_frame, ++cur_data) {
            // snew((*probs)[cur_atom].x[cur_data], 4);
            (*probs)[cur_atom].x[cur_data] = nodepool;
            for (i = 0; i < 3; ++i) {
                (*probs)[cur_atom].x[cur_data][i].index = i + 1;
                (*probs)[cur_atom].x[cur_data][i].value = x2[cur_frame][indx2[cur_atom]][i] * 10.0;
            }
            (*probs)[cur_atom].x[cur_data][i].index = -1;
            nodepool += 4;
        }
    }
    printf("\n");
    fflush(stdout);
}

void free_svm_probs(struct svm_problem *probs,
                    int nprobs,
                    int nvecs) {
    if (nprobs > 0) {
        sfree(probs[0].y); // Free target array
        if (nvecs > 0)
            sfree(probs[0].x[0]); // The first atom's first frame's x points to the head of the node space
    }
    for (int i = 0; i < nprobs; ++i) {
        sfree(probs[i].x);
    }
    sfree(probs);
}

void train_svm_probs(struct svm_problem *probs,
                     int num_probs,
                     real gamma,
                     real c,
                     int nthreads,
                     struct svm_model **models) {
    struct svm_parameter param; // Parameters used for training

    print_log("svm-training trajectory atoms with gamma = %f and C = %f...\n", gamma, c);
    flush_log();

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
    if (nthreads > 0)
        omp_set_num_threads(nthreads);
    if (nthreads > 1 || nthreads <= 0)
        print_log("svm training will be parallelized.\n");
#endif

    /* Train svm */
    int i;
#pragma omp parallel for schedule(dynamic) private(i) shared(num_probs,models,probs,param)
    for (i = 0; i < num_probs; ++i) {
    #if defined _OPENMP && defined EC_DEBUG
        print_log("%d threads running svm-train.\n", omp_get_num_threads());
    #endif
        models[i] = svm_train(&(probs[i]), &param);
    }
}

static void free_svm_model(struct svm_model *model) {
    svm_free_model_content(model);
    svm_free_and_destroy_model(&model);
}

void free_svm_models(struct svm_model **models, int num_models) {
    for (int i = 0; i < num_models; ++i) {
        free_svm_model(models[i]);
    }
    sfree(models);
}


void calc_eta(struct svm_model **models,
              int num_models,
              int num_frames,
              real *eta) {
    int i;

    print_log("Calculating eta values...\n");
    flush_log();

    for (i = 0; i < num_models; ++i) {
        eta[i] = 1.0 - svm_get_nr_sv(models[i]) / (2.0 * (real)num_frames);
    }
}

void calc_eta_res(eta_dat_t *eta_dat) {
    print_log("Reading residue info from %s...\n", eta_dat->fnames[eRES1]);
    switch(fn2ftp(eta_dat->fnames[eRES1])) {
        case efPDB:
            eta_res_pdb(eta_dat);
            break;
        case efGRO: // TODO: try using this for tpr as well, or vice versa?
            eta_res_tps(eta_dat);
            break;
        case efTPR:
            eta_res_tpx(eta_dat);
            break;
        default:
            print_log("%s is not a supported filetype for residue information. Skipping eta residue calculation.\n",
                eta_dat->fnames[eRES1]);
            flush_log();
    }
}


void save_eta(eta_dat_t *eta_dat) {
    FILE *f = fopen(eta_dat->fnames[eETA_ATOM], "w");

    // atom etas
    if (f) {
        print_log("Saving eta values to %s...\n", eta_dat->fnames[eETA_ATOM]);
        fprintf(f, "# ATOM\tETA\n");
        for (int i = 0; i < eta_dat->natoms; ++i) {
            // Add 1 to the atom ID because Gromacs's stored ID
            // is 1 lower than in given index files
            fprintf(f, "%d\t%f\n", eta_dat->atom_IDs[i] + 1, eta_dat->eta[i]);
        }

        fclose(f);
        f = NULL;
    }
    else {
        print_log("Failed to open file %s for saving eta values.\n",
            eta_dat->fnames[eETA_ATOM]);
    }

    // residue etas
    if (eta_dat->res_eta) {
        f = fopen(eta_dat->fnames[eETA_RES], "w");

        if (f) {
            print_log("Saving residue eta values to %s...\n",
                eta_dat->fnames[eETA_RES]);

            fprintf(f, "# RES\tNATOMS\tETA\n");
            for (int i = 0; i < eta_dat->nres; ++i) {
                fprintf(f, "%d%s\t%d\t%f\n", eta_dat->res_IDs[i],
                                             eta_dat->res_names[i],
                                             eta_dat->res_natoms[i],
                                             eta_dat->res_eta[i]);
            }

            fclose(f);
            f = NULL;
        }
        else {
            print_log("Failed to open file %s for saving residue eta values.\n",
                eta_dat->fnames[eETA_RES]);
        }
    }
    flush_log();
}

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
    if (!*x)
        log_fatal(FARGS, "Failed to allocate memory for trajectory!\n");

    *natoms = read_first_x(*oenv, &status, traj_fname, &t, *x, box);
    if (!status)
        log_fatal(FARGS, "Failed to open trajectory!\n");

    do {
        ++(*nframes);
        if (*nframes >= est_frames) {
            est_frames += FRAMESTEP;
            srenew(*x, est_frames);
        }
        snew((*x)[*nframes], *natoms);

        if (!(*x)[*nframes]) {
            log_fatal(FARGS, "Failed to allocate memory at frame %d!\n", *nframes);
        }
    } while(read_next_x(*oenv, status, &t,
#ifndef GRO_V5
        *natoms,
#endif
        (*x)[*nframes], box));

    sfree((*x)[*nframes]);
    close_trx(status);
}

void free_traj(rvec **x, int nframes) {
    for (int i = 0; i < nframes; ++i) {
        sfree(x[i]);
    }
    sfree(x);
}

static void eta_res_pdb(eta_dat_t *eta_dat) {
    char title[256];
    t_atoms atoms;
    rvec *x;

    atoms.nr = eta_dat->natoms_all;
    snew(atoms.atom, eta_dat->natoms_all);
    snew(atoms.atomname, eta_dat->natoms_all);
    snew(atoms.atomtype, eta_dat->natoms_all);
    snew(atoms.atomtypeB, eta_dat->natoms_all);
    atoms.nres = eta_dat->natoms_all;
    snew(atoms.resinfo, eta_dat->natoms_all);
    snew(atoms.pdbinfo, eta_dat->natoms_all);

    snew(x, eta_dat->natoms_all);

    read_pdb_conf(eta_dat->fnames[eRES1], title, &atoms, x, NULL, NULL, FALSE, NULL);

    sfree(x);

    sfree(atoms.atomname);
    sfree(atoms.atomtype);
    sfree(atoms.atomtypeB);
    sfree(atoms.pdbinfo);

    f_calc_eta_res(eta_dat, &atoms);

    sfree(atoms.atom);
    sfree(atoms.resinfo);
}

// Try res_tpx for gro and tpr instead of this.
static void eta_res_tps(eta_dat_t *eta_dat) {
    char title[256];
    t_topology top;
    rvec *x = NULL;
    matrix box;
    int ePBC;

    init_top(&top);

    read_tps_conf(eta_dat->fnames[eRES1], title, &top, &ePBC, &x, NULL, box, FALSE);

    f_calc_eta_res(eta_dat, &(top.atoms));

    // Cannot use done_top(), causes error- pointer being freed was not allocated. See implementation in typedefs.c
    done_atom(&(top.atoms));
    done_symtab(&(top.symtab));
    done_block(&(top.cgs));
    done_block(&(top.mols));
    done_blocka(&(top.excls));

    sfree(x);
}

// TODO: Does this work for gro files generated by grompp etc?
static void eta_res_tpx(eta_dat_t *eta_dat) {
    t_inputrec ir;
    gmx_mtop_t mtop;
    matrix box;
    int natoms, i;

    read_tpx(eta_dat->fnames[eRES1], &ir, box, &natoms, NULL, NULL, NULL, &mtop);

    f_calc_eta_res(eta_dat, &(mtop.moltype->atoms));
}

static void f_calc_eta_res(eta_dat_t *eta_dat,
                           t_atoms *atoms) {
    real *sums;

    print_log("Calculating residue eta values...\n");

    snew(eta_dat->res_IDs, atoms->nres);
    snew(eta_dat->res_names, atoms->nres);
    snew(eta_dat->res_natoms, atoms->nres);
    snew(eta_dat->res_eta, atoms->nres);
    snew(sums, atoms->nres);

    // Only sum the atoms that are in the indexes in atom_IDs
    int resid;
    for (int i = 0; i < eta_dat->natoms; ++i) {
        resid = atoms->atom[eta_dat->atom_IDs[i]].resind;
        sums[resid] += eta_dat->eta[i];
        eta_dat->res_natoms[resid] += 1;
    }

    // Calculate average etas
    for (int i = 0; i < atoms->nres; ++i) {
        if (eta_dat->res_natoms[i] > 0)
            eta_dat->res_eta[i] = sums[i] / eta_dat->res_natoms[i];
    }

    // Store residue info in eta_res_t struct
    eta_dat->nres = atoms->nres;
    for (int i = 0; i < atoms->nres; ++i) {
        eta_dat->res_IDs[i] = atoms->resinfo[i].nr;
        eta_dat->res_names[i] = *(atoms->resinfo[i].name);
    }

    sfree(sums);
}

/********************************************************
 * Logging functions
 ********************************************************/

void init_log(const char *logfile, int argc, char *argv[]) {
    out_log = fopen(logfile, "a");

    time_t t = time(NULL);
    struct tm *ltime = localtime(&t);

    fprintf(out_log, "\n");
	for(int i = 0; i < argc; ++i) {
		fprintf(out_log, "%s ", argv[i]);
	}
    fprintf(out_log, "\nRun: %d-%d-%d %d:%d:%d\n",
		ltime->tm_mon + 1, ltime->tm_mday, ltime->tm_year + 1900,
		ltime->tm_hour, ltime->tm_min, ltime->tm_sec);
    flush_log();
}

void close_log() {
    fclose(out_log);
}

void print_log(char const *fmt, ...) {
    va_list arg;
    va_start(arg, fmt);
    vprintf(fmt, arg);
    va_end(arg);
    if (out_log != NULL) {
        va_start(arg, fmt);
        vfprintf(out_log, fmt, arg);
        va_end(arg);
    }
}

void flush_log() {
    if (out_log != NULL)
        fflush(out_log);
}

void log_fatal(int fatal_errno,
               const char *file,
               int line,
               char const *fmt, ...) {
    va_list arg;
    if (out_log != NULL) {
        va_start(arg, fmt);
        fprintf(out_log, "Fatal error in source file %s line %d: ", file, line);
        vfprintf(out_log, fmt, arg);
        va_end(arg);
    }
    va_start(arg, fmt);
    gmx_fatal(fatal_errno, file, line, fmt, arg);
    va_end(arg);
}
