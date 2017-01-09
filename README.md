### g_ensemble_comp

g_ensemble_comp evaluates the difference between two conformational
ensembles, R and R'. Quanitification is in terms of a true metric that
satisfies the conditions set forth by the zeroth law of thermodynamics. The
quantification metric eta=1-|Overlap|=|R'|-|Overlap|=DeltaR is normalized,
that is, 0<=eta<1, and takes up a value closer to unity as the difference
between the ensembles increases. For two Gaussian distributions with
identical standard deviations of 0.5 A, eta=0.68 represents a geometric
center deviation of 1 A.

The two ensembles are provided as two trajectory files specified
by the `-f1` and `-f2` options (supported formats=xtc,trr,pdb).
We recommend that frame numbers in the trajectory files are in the range
2500-5000. While the speed of the algorithm decreases with increase in
ensemble size, the numerical accuracy of the calculation reduces with
decrease in ensemble size, and a small number of frames may not provide a good
representation of the ensemble.  
Note that if you are not interested in eta to
reflect changes in whole molecule translation/rotation, then these degress
of freedom need to be removed prior to ensemble comparison. For most cases,
this can be accomplished by fitting all conformations in the two ensembles
on to one single representative structure, such as the x-ray structure. For
this, we recommend the use of trjconv with the -fit rot+trans option.

By default, differences (eta) are estimated for all atoms, but comparisons can be done
for a smaller specific group of atoms, which can be selected from index
files `-n1` and `-n2`. The average eta can also be calculated for each residue specified
in a structure file given by `-res` (.pdb and .gro are supported).

Overlaps are estimated by training a support vector
machine in a pre-defined Hilbert space specified by the width of the RDF
Kernel (gamma=0.4) and the maximum value that can be taken up by the
Lagrange multiplier (C=100.0). The values of C and gamma can be changed with `-c` and `-g`,
but such changes may increase mean absolute error (MAE=3.26%) of the method.

By default, g_ensemble_comp is parallelized with OpenMP (see installation instructions below).
To set the number of threads at runtime, set `-nthreads X` where X is the number of threads to use. The default behavior is to use the maximum number of cores available.

Methodoligical details and example applications can be found in
Leighty and Varma, JCTC, 2013, 9: 868-875.
Varma, Botlani and Leighty, Proteins, 2014, 82: 3241-3254.
Dutta, Botlani and Varma, JPC B, 2014, 118: 14795-14807.

### INSTALLATION

The following instructions are for unix-based operating systems such as OSX and Linux. Windows is not currently supported.

1. Install Gromacs version 4.5.x or later from http://www.gromacs.org.

2. `git clone` or otherwise obtain and `cd` to the 'g_ensemble_comp' repository.

3. Run `sudo make install` with the necessary arguments for your environment (see below).

If you do not have Gromacs version 5.x installed, you will need to set the makefile's `VGRO` variable to the root Gromacs version number. If Gromacs is installed in a non-default directory (ie not in /usr/local/gromacs) then you will have to set the `GROMACS` variable to the Gromacs installation directory that contains the 'include' and 'lib' folders. For example, if you are running Gromacs 4.5.3 installed in /home/user/tools/gromacs-4.5.3, then you would run the following command:

``` bash
$ sudo make install VGRO=4 GROMACS=/home/user/tools/gromacs-4.5.3
```

If you must run `make install` without sudo privileges, you will need to set the `INSTALL` variable to a path that you can write to.
The default install path is /usr/local/bin. Depending on your system and chosen installation directory, you may have to add g_ensemble_comp to your PATH.

You will also need to make sure that `CC` and `CXX` are set to the C compiler and C++ compiler commands, respectively, that were used to compile your installation of Gromacs. For example, if your environment's default compiler is clang, but you used gcc to compile Gromacs, you would run the following:

``` bash
$ sudo make install CC=gcc CXX=g++
```

If you want to build without OpenMP, set `PARALLEL=0`. You can also add compilation flags by setting `CFLAGS`, and linker flags/libraries by setting `LIBS`. For example, if you set `LIBS=-static` to statically link g_ensemble_comp's dependencies, you can then run the same binary in a different environment without the same C runtime or Gromacs library present.

### USAGE

``` bash
$ g_ensemble_comp -f1 traj1.xtc -f2 traj2.xtc
```

Run `g_ensemble_comp -h` to get detailed usage instructions. The instructions are also provided in the introductory paragraph above as well as the tutor directory.

### Copyright
(c) 2017 Ahnaf Siddiqui, Mohsen Botlani and Sameer Varma  
The code uses SVM libraries: LIBSVM copyright 2000-2014 Chih-Chung Chang and Chih-Jen Lin.

g_ensemble_comp uses the GROMACS molecular simulation package API.  
Copyright (c) 1991-2000, University of Groningen, The Netherlands.  
Copyright (c) 2001-2004, The GROMACS development team.  
Copyright (c) 2013,2014, by the GROMACS development team,  
led by Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,  
and including many others, as listed at http://www.gromacs.org.
