CC=gcc
CXX=g++
GROMACS=/usr/local/gromacs
VGRO=5
SVM=../libsvm-3.20
INSTALL=/usr/local/bin

ifeq ($(VGRO),5)
INCGRO=-I$(GROMACS)/include/ -I$(GROMACS)/include/gromacs/utility -I$(GROMACS)/include/gromacs/fileio -I$(GROMACS)/include/gromacs/commandline -I$(GROMACS)/include/gromacs/legacyheaders
LINKGRO=-L$(GROMACS)/lib/i386-linux-gnu
LIBGRO=-lgromacs
DEFV5=-D GRO_V5
else
INCGRO=-I$(GROMACS)/include/gromacs
LINKGRO=-L$(GROMACS)/lib
LIBGRO=-lgmx
DEFV5=
endif

.PHONY: all install clean

all: svmanalys eta_anal 

install: svmanalys
	install svmanalys $(INSTALL)

svmanalys: svmanalys.o
	$(CXX) -o svmanalys svmanalys.o $(SVM)/svm.o $(LINKGRO) $(LIBGRO)

eta_anal: eta_anal.o
	$(CXX) -o eta_anal eta_anal.o $(SVM)/svm.o $(LINKGRO) $(LIBGRO)

svmanalys.o: src/svmanalys.c src/svmanalys.h
	$(CC) -c src/svmanalys.c $(INCGRO) -I$(SVM) $(DEFV5)

eta_anal.o: src/eta_anal.c src/svmanalys.h
	$(CC) -c src/eta_anal.c $(INCGRO) -I$(SVM) $(DEFV5)

clean:
	rm -f svmanalys.o eta_anal.o svmanalys eta_anal