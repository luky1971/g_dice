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
 * Copyright (c) 2015, University of South Florida.
 */

#include "ensemble_comp.h"

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"g_ensemble_comp evaluates the difference between two conformational ensembles, R and R'.",
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
                "By default, g_ensemble_comp is parallelized if compiled with OpenMP. To disable parallelization, add the option -nopar.\n",
                "Methodoligical details and example applications can be found in\n",
		"Leighty and Varma, JCTC, 2013, 9: 868-875.\n",
		"Varma, Botlani and Leighty, Proteins, 2014, 82: 3241-3254.\n",
		"Dutta, Botlani and Varma, JPC B, 2014, 118: 14795-14807.\n"
	};
	const char *fnames[eNUMFILES];
	output_env_t oenv = NULL;
	real gamma = GAMMA, c = COST;
	gmx_bool nopar = FALSE;

	/* eta values */
	real *eta;
	int natoms;

	init_log("etalog.txt", argv[0]);

	t_filenm fnm[] = {
		{efTRX, "-f1", "traj1.xtc", ffREAD},
		{efTRX, "-f2", "traj2.xtc", ffREAD},
		{efNDX, "-n1", "index1.ndx", ffOPTRD},
		{efNDX, "-n2", "index2.ndx", ffOPTRD},
		{efSTX, "-res", "res.pdb", ffOPTRD},
		{efSTX, "-top", "top.tpr", ffOPTRD},
		{efDAT, "-eta_atom", "eta_atom.dat", ffWRITE},
		{efDAT, "-eta_res", "eta_res.dat", ffWRITE}
	};

	t_pargs pa[] = {
		{"-g", FALSE, etREAL, {&gamma}, "RBD Kernel width (default=0.4)"},
		{"-c", FALSE, etREAL, {&c}, "Max value of Lagrange multiplier (default=100)"},
		{"-nopar", FALSE, etBOOL, {&nopar}, "Set this option to disable parallelization"}
	};

	parse_common_args(&argc, argv, 0, eNUMFILES, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

	fnames[eTRAJ1] = opt2fn("-f1", eNUMFILES, fnm);
	fnames[eTRAJ2] = opt2fn("-f2", eNUMFILES, fnm);
	fnames[eNDX1] = opt2fn_null("-n1", eNUMFILES, fnm);
	fnames[eNDX2] = opt2fn_null("-n2", eNUMFILES, fnm);
	fnames[eRES1] = opt2fn_null("-res", eNUMFILES, fnm);
	fnames[eTOP1] = opt2fn_null("-top", eNUMFILES, fnm);
	fnames[eETA_ATOM] = opt2fn("-eta_atom", eNUMFILES, fnm);
	fnames[eETA_RES] = opt2fn("-eta_res", eNUMFILES, fnm);

	// ensemble_comp(fnames, gamma, c, &eta, &natoms, !nopar, &oenv);

	// TEST

	if(fnames[eTOP1] != NULL) {
		to_internal_coords(fnames, &oenv, "inter_coords.txt");
	}

	// rvec test[4];
	// test[0][XX] = 2, test[0][YY] = 3, test[0][ZZ] = 0;
	// test[1][XX] = 3, test[1][YY] = -1, test[1][ZZ] = 0;
	// test[2][XX] = 6, test[2][YY] = 0, test[2][ZZ] = 0;
	// test[3][XX] = 8, test[3][YY] = 2.9, test[3][ZZ] = -0.001;

	// calc_angles(test);

	// end TEST

	// save_eta(eta, natoms, fnames[eETA_ATOM]);

	// if(fnames[eRES1] != NULL) {
	// 	eta_res_t eta_res;

	// 	calc_eta_res(fnames[eRES1], eta, natoms, &eta_res);
	// 	save_eta_res(&eta_res, fnames[eETA_RES]);

	// 	sfree(eta_res.res_nums);
	// 	sfree(eta_res.res_names);
	// 	sfree(eta_res.avg_etas);
	// }

	// sfree(eta);

	print_log("%s completed successfully.\n", argv[0]);
	close_log();

	return 0;
}
