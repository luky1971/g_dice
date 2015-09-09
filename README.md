# etanalys
etanalys analyzes GROMACS trajectory files using support vector machine algorithms and calculates eta values.
Support vector machine training is provided by libsvm, Copyright (c) 2000-2014 Chih-Chung Chang and Chih-Jen Lin.

The following instructions are for unix-based operating systems such as OSX and Linux. Windows is not currently supported.

To install:

1. Install Gromacs version 4.5.x or later from http://www.gromacs.org.

3. `git clone` or otherwise obtain and `cd` to the 'etanalys' repository.

3. Run `sudo make install` with the necessary arguments for your environment (see below).

If you do not have Gromacs version 5.x installed, you will need to set the makefile's `VGRO` variable to the root Gromacs version number. If Gromacs is installed in a non-default directory (ie not in /usr/local/gromacs) then you will have to set the `GROMACS` variable to the Gromacs installation directory that contains the 'include' and 'lib' folders. For example, if you are running Gromacs 4.5.3 installed in /home/user/tools/gromacs-4.5.3, then you would run the following command:

`sudo make install VGRO=4 GROMACS=/home/user/tools/gromacs-4.5.3`

If you must run `make install` without sudo privileges, you will need to set the `INSTALL` variable to a path that you can write to. The default install path is /usr/local/bin.

If you are not using gcc, you will also need to set `CC` and `CXX` to your C compiler and C++ compiler commands respectively.

After installing, run `etanalys -h` to see its usage. Depending on your system and chosen installation directory, you may have to add etanalys to your PATH. Google is your friend.