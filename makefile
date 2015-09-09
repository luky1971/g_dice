CC=gcc
CXX=g++
GROMACS=/usr/local/gromacs
VGRO=5
SVM=extern/libsvm-3.20
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

.PHONY: install clean

etanalys: etanalys.o eta.o
	make -C $(SVM) && $(CXX) -o etanalys etanalys.o eta.o $(SVM)/svm.o $(LINKGRO) $(LIBGRO)

install: etanalys
	install etanalys $(INSTALL)

etanalys.o: src/etanalys.c src/eta.h
	$(CC) -c src/etanalys.c $(INCGRO) -I$(SVM) $(DEFV5)

eta.o: src/eta.c src/eta.h
	$(CC) -c src/eta.c $(INCGRO) -I$(SVM) $(DEFV5)

clean:
	make clean -C $(SVM) && rm -f etanalys.o eta.o etanalys