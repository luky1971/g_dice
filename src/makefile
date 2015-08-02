INCGRO = /usr/local/gromacs/include/gromacs
LINKGRO = /usr/local/gromacs/lib

svmanalys: svmanalys.o svmutils.o
	gcc -o svmanalys svmanalys.o svmutils.o -L$(LINKGRO) -lgmx

gsvmanalys: svmanalys.c svmutils.c
	gcc -g -o gsvmanalys svmanalys.c svmutils.c -I$(INCGRO) -L$(LINKGRO) -lgmx

comanalys: comanalys.o svmutils.o
	gcc -o comanalys comanalys.o svmutils.o -L$(LINKGRO) -lgmx

svmanalys.o: svmanalys.c svmutils.h
	gcc -c svmanalys.c -I$(INCGRO)

comanalys.o: comanalys.c svmutils.h
	gcc -c comanalys.c -I$(INCGRO)

svmutils.o: svmutils.c svmutils.h
	gcc -c svmutils.c -I$(INCGRO)