# svmanalys
svmanalys analyzes GROMACS trajectory files using support vector machine algorithms and calculates eta values.

The following instructions are for unix-based operating systems such as OSX and Linux. Windows is not currently supported.

To install:

1. Install Gromacs version 4.5.x or later from http://www.gromacs.org.

2. Install libsvm from https://www.csie.ntu.edu.tw/~cjlin/libsvm/. Version 3.20 is officially supported. 

3. `git clone` and `cd` to the 'svmanalys' repository.

3. Run `sudo make install` with the necessary arguments for your environment (see below).

If you do not have Gromacs version 5.x installed, you will need to set the makefile's `VGRO` variable to the root Gromacs version number. You will also need to specify the installation directory of libsvm in the `SVM` variable. If Gromacs is installed in a non-default directory (ie not in /usr/local/gromacs) then you will have to set the `GROMACS` variable to the Gromacs installation directory containing the 'include' and 'lib' folders. For example, if you are running Gromacs 4.5.3 in its default installation and you have libsvm installed in /usr/local/libsvm-3.20, then you would run the following command:

`sudo make install VGRO=4 SVM=/usr/local/libsvm-3.20`

If you must run `make install` without sudo privileges, you will need to set the `INSTALL` variable to a path that you can write to. The default install path is /usr/local/bin.

If you are not using gcc, you will also need to set `CC` and `CXX` to your C compiler and C++ compiler commands respectively.

After installing, run `svmanalys -h` to see its usage. Depending on your system and chosen installation directory, you may have to add svmanalys to your PATH. Google is your friend.