CC=gcc
CXX=g++
INSTALL=/usr/local/bin
SVM=../libsvm-3.20
VGRO=5

ifeq ($(VGRO),5)
INCGRO=-I/usr/local/gromacs/include/gromacs/utility -I/usr/local/gromacs/include/gromacs/fileio -I/usr/local/gromacs/include/gromacs/commandline -I/usr/local/gromacs/include/gromacs/legacyheaders
LINKGRO=-L/usr/local/gromacs/lib/i386-linux-gnu
LIBGRO=-lgromacs
DEFV5=-D GRO_V5
else
INCGRO=-I/usr/local/gromacs/include/gromacs
LINKGRO=-L/usr/local/gromacs/lib
LIBGRO=-lgmx
DEFV5=
endif

svmanalys: svmanalys.o
	$(CXX) -o svmanalys svmanalys.o $(SVM)/svm.o $(LINKGRO) $(LIBGRO)

install: svmanalys
	install svmanalys $(INSTALL)

.PHONY: install

svmanalys.o: src/svmanalys.c src/svmanalys.h
	$(CC) -c src/svmanalys.c $(INCGRO) -I$(SVM) $(DEFV5)