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

#include "eta.h"

int main(int argc, char *argv[]) {
	const char *desc[] = {
		"etanalys analyzes trajectory files produced by GROMACS using support vector machine algorithms.\n",
		"It trains and calculates eta values for atoms from two trajectories.\n",
		"-f1 and -f2: Specify the two trajectory files (xtc, trr, and pdb files are supported).\n",
		"-n1 and -n2: Specify optional index files to select atom groups.\n",
		"-eta_atom: Specify the name of the output file (default is eta_atom.dat).\n",
		"-g and -c: Specify your own gamma and C parameters for svm-train (default is gamma = 0.01, c = 10.0).\n"
	};
	const char *fnames[eNUMFILES];
	output_env_t oenv = NULL;
	real gamma = GAMMA, c = COST;

	/* eta values */
	real *eta;
	int natoms;

	init_log("etalog.txt", argv[0]);

	t_filenm fnm[] = {
		{efTRX, "-f1", "traj1.xtc", ffREAD},
		{efTRX, "-f2", "traj2.xtc", ffREAD},
		{efNDX, "-n1", "index1.ndx", ffOPTRD},
		{efNDX, "-n2", "index2.ndx", ffOPTRD},
		{efDAT, "-eta_atom", "eta_atom.dat", ffWRITE}
	};

	t_pargs pa[] = {
		{"-g", FALSE, etREAL, {&gamma}, "gamma parameter for svm-train"},
		{"-c", FALSE, etREAL, {&c}, "C (cost) parameter for svm-train"}
	};

	parse_common_args(&argc, argv, 0, eNUMFILES, fnm, 
		asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

	fnames[eTRAJ1] = opt2fn("-f1", eNUMFILES, fnm);
	fnames[eTRAJ2] = opt2fn("-f2", eNUMFILES, fnm);
	fnames[eNDX1] = opt2fn_null("-n1", eNUMFILES, fnm);
	fnames[eNDX2] = opt2fn_null("-n2", eNUMFILES, fnm);
	fnames[eETA_ATOM] = opt2fn("-eta_atom", eNUMFILES, fnm);

	etanalys(fnames, gamma, c, &eta, &natoms, oenv);

	save_eta(eta, natoms, fnames[eETA_ATOM]);
	print_log("Eta values saved in file %s\n", fnames[eETA_ATOM]);

	sfree(eta);
	close_log();
}
