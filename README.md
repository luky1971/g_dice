# etanalys
etanalys analyzes GROMACS trajectory files using support vector machine algorithms and calculates eta values.

Etanalys takes as input two trajectory files representing two conformations of the same groups of atoms in two different states. It svm-trains each atom from both trajectories and then calculates each atom's eta from the number of support vectors. The eta values calculated for each atom are a quantification of the difference between the two conformations. Therefore, the given trajectories or atom groups within trajectories should have the same number of atoms. Etanalys can also take input index files to select matching atom groups from the two trajectory files. If the resulting eta values are not satisfactory, you can set your own C and gamma parameter values for svm-train with -c and -g.

The following instructions are for unix-based operating systems such as OSX and Linux. Windows is not currently supported.

### To install:

1. Install Gromacs version 4.5.x or later from http://www.gromacs.org.

3. `git clone` or otherwise obtain and `cd` to the 'etanalys' repository.

3. Run `sudo make install` with the necessary arguments for your environment (see below).

If you do not have Gromacs version 5.x installed, you will need to set the makefile's `VGRO` variable to the root Gromacs version number. If Gromacs is installed in a non-default directory (ie not in /usr/local/gromacs) then you will have to set the `GROMACS` variable to the Gromacs installation directory that contains the 'include' and 'lib' folders. For example, if you are running Gromacs 4.5.3 installed in /home/user/tools/gromacs-4.5.3, then you would run the following command:

`sudo make install VGRO=4 GROMACS=/home/user/tools/gromacs-4.5.3`

If you must run `make install` without sudo privileges, you will need to set the `INSTALL` variable to a path that you can write to. The default install path is /usr/local/bin. Depending on your system and chosen installation directory, you may have to add etanalys to your PATH. Google is your friend.

If you are not using gcc, you will also need to set `CC` and `CXX` to your C compiler and C++ compiler commands respectively.

### Usage:

After installing, you can run etanalys with a pair of GROMACS trajectory files with a command such as the following:

`etanalys -f1 <trajectory 1 filename> -f2 <trajectory 2 filename> -eta_atom <output filename>`

These are all the options you can set:

-f1 and -f2: Specify the two trajectory files (xtc, trr, and pdb files are supported).  
-n1 and -n2: Specify optional index (.ndx) files to select atom groups.  
-eta_atom: Specify the name of the output file (default is eta_atom.dat).  
-g and -c: Specify your own gamma and C parameters for svm-train.  

The resulting eta values will be found in the output file.

### Copyright 
(c) 2015, University of South Florida
Authors: Ahnaf Siddiqui, Mohsen Botlani-Esfahani, and Dr. Sameer Varma

### Acknowledgments:

Support vector machine training is provided by libsvm, copyright 2000-2014 Chih-Chung Chang and Chih-Jen Lin.

The authors would like to acknowledge the use of the services provided by Research Computing at the University of South Florida.