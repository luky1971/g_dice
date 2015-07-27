comanalys: comanalys.o svmutils.o
	gcc -o comanalys comanalys.o svmutils.o -L/usr/local/gromacs/lib -lgmx

comanalys.o: comanalys.c svmutils.h
	gcc -c comanalys.c -I/usr/local/gromacs/include/gromacs

svmutils.o: svmutils.c svmutils.h
	gcc -c svmutils.c -I/usr/local/gromacs/include/gromacs