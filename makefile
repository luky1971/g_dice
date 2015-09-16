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

g_ensemble_comp: g_ensemble_comp.o ensemble_comp.o
	make -C $(SVM) && $(CXX) -o g_ensemble_comp g_ensemble_comp.o ensemble_comp.o $(SVM)/svm.o $(LINKGRO) $(LIBGRO)

install: g_ensemble_comp
	install g_ensemble_comp $(INSTALL)

g_ensemble_comp.o: src/g_ensemble_comp.c src/ensemble_comp.h
	$(CC) -c src/g_ensemble_comp.c $(INCGRO) -I$(SVM) $(DEFV5)

ensemble_comp.o: src/ensemble_comp.c src/ensemble_comp.h
	$(CC) -c src/ensemble_comp.c $(INCGRO) -I$(SVM) $(DEFV5)

clean:
	make clean -C $(SVM) && rm -f g_ensemble_comp.o ensemble_comp.o g_ensemble_comp
