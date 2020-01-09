/*
 * Copyright 2016 Ahnaf Siddiqui, Mohsen Botlani and Sameer Varma
 *
 * This program uses the GROMACS molecular simulation package API.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed at http://www.gromacs.org.

 * g_dice quantifies the difference between two conformational ensembles (two trajectory files)
 * Quantification is in terms of a true metric, eta=1-Overlap
 * Leighty and Varma, Quantifying Changes in Intrinsic Molecular Motion Using Support Vector Machines, J. Chem. Theory Comput. 2013, 9, 868-875.
 */

#include "ensemble_comp.h"

int main(int argc, char *argv[]) {
    const char *desc[] = {
        "g_dice evaluates the difference between two conformational ensembles, R and R'.",
        "Quanitification is in terms of a true metric that satisfies the conditions set forth by the zeroth law of thermodynamics.",
        "The quantification metric eta=1-|Overlap|=|R'|-|Overlap|=DeltaR is normalized, that is, 0<=eta<1,",
        "and takes up a value closer to unity as the difference between the ensembles increases.",
        "For two Gaussian distributions with identical standard deviations of 0.5 A,",
        "eta=0.68 represents a geometric center deviation of 1 A.\n",
        "The two ensembles are provided as two trajectory files specified by the -f1 and -f2 options (supported formats=xtc,trr,pdb).",
        "We recommend that frames numbers in the trajectory files are in the range 2500-5000.",
        "While the speed of the algorithm decreases with increase in ensemble size,",
        "the numerical accuracy of the calculation reduces with decrease in ensemble size,",
        "and a small number of frames may not provide a good representation of the ensemble.",
        "Note that if you are not interested in eta to reflect changes in whole molecule translation/rotation,",
        "then these degress of freedom need to be removed prior to ensemble comparison.",
        "For most cases, this can be accomplished by fitting all conformations",
        "in the two ensembles on to one single representative structure, such as the x-ray structure",
        "For this, we recommend the use of trjconv with the -fit rot+trans option.\n"
        "By default, differences (eta) are estimated for all atoms, but comparisons can be done for a smaller specific group of atoms,",
        "which can be selected from index files -n1 and -n2. ",
        "The average eta can also be calculated for each residue specified in a structure file given by -res (.pdb and .gro are supported).\n",
        "Overlaps are estimated by training a support vector machine in a pre-defined Hilbert space specified by",
        "the width of the RDF Kernel (gamma=0.4) and",
        "the maximum value that can be taken up by the Lagrange multiplier (C=100.0).",
        "The values of C and gamma can be changed with -c and -g,",
        "but such changes may increase mean absolute error (MAE=3.26%) of the method.\n",
        "If g_dice was built with OPENMP, you can set the number of threads to use with -nthreads X,\n",
        "where X is the number of threads to use. The default behavior is to use the maximum number of cores available.\n\n",
        "Methodoligical details and example applications can be found in\n",
        "Leighty and Varma, JCTC, 2013, 9: 868-875.\n",
        "Varma, Botlani and Leighty, Proteins, 2014, 82: 3241-3254.\n",
        "Dutta, Botlani and Varma, JPC B, 2014, 118: 14795-14807.\n"
    };

    eta_dat_t eta_dat;
    init_eta_dat(&eta_dat);

    init_log("eta.log", argc, argv);

    t_filenm fnm[] = {
        {efTRX, "-f1", "traj1.xtc", ffREAD},
        {efTRX, "-f2", "traj2.xtc", ffREAD},
        {efNDX, "-n1", "index1.ndx", ffOPTRD},
        {efNDX, "-n2", "index2.ndx", ffOPTRD},
        {efSTX, "-res", "res.pdb", ffOPTRD},
        {efDAT, "-eta_atom", "eta_atom.dat", ffWRITE},
        {efDAT, "-eta_res", "eta_res.dat", ffWRITE}
    };

    t_pargs pa[] = {
        {"-g", FALSE, etREAL, {&eta_dat.gamma}, "RBD Kernel width (default=0.4)"},
        {"-c", FALSE, etREAL, {&eta_dat.c}, "Max value of Lagrange multiplier (default=100)"},
        {"-nthreads", FALSE, etINT, {&eta_dat.nthreads}, "set the number of parallel threads to use (default is max available)"}
    };

    parse_common_args(&argc, argv, 0, eNUMFILES, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &eta_dat.oenv);

    eta_dat.fnames[eTRAJ1] = opt2fn("-f1", eNUMFILES, fnm);
    eta_dat.fnames[eTRAJ2] = opt2fn("-f2", eNUMFILES, fnm);
    eta_dat.fnames[eNDX1] = opt2fn_null("-n1", eNUMFILES, fnm);
    eta_dat.fnames[eNDX2] = opt2fn_null("-n2", eNUMFILES, fnm);
    eta_dat.fnames[eRES1] = opt2fn_null("-res", eNUMFILES, fnm);
    eta_dat.fnames[eETA_ATOM] = opt2fn("-eta_atom", eNUMFILES, fnm);
    eta_dat.fnames[eETA_RES] = opt2fn("-eta_res", eNUMFILES, fnm);

    // Calculate and output eta
    ensemble_comp(&eta_dat);
    save_eta(&eta_dat);
    free_eta_dat(&eta_dat);

    print_log("%s completed successfully.\n", argv[0]);
    close_log();

    return 0;
}
