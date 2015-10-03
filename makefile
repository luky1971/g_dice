CC=gcc
CXX=g++
PARALLEL=1
GROMACS=/usr/local/gromacs
VGRO=5
SVM=extern/libsvm-3.20
SRC=src
BUILD=build
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
endif

ifneq ($(PARALLEL),0)
CFLAGS+= -fopenmp
endif

.PHONY: install clean

$(BUILD)/g_ensemble_comp: $(BUILD)/g_ensemble_comp.o $(BUILD)/ensemble_comp.o
	make -C $(SVM) && $(CXX) $(CFLAGS) -o $(BUILD)/g_ensemble_comp $(BUILD)/g_ensemble_comp.o $(BUILD)/ensemble_comp.o $(SVM)/svm.o $(LINKGRO) $(LIBGRO)

install: $(BUILD)/g_ensemble_comp
	install $(BUILD)/g_ensemble_comp $(INSTALL)

$(BUILD)/g_ensemble_comp.o: $(SRC)/g_ensemble_comp.c $(SRC)/ensemble_comp.h
	$(CC) $(CFLAGS) -o $(BUILD)/g_ensemble_comp.o -c $(SRC)/g_ensemble_comp.c $(DEFV5) $(INCGRO) -I$(SVM)

$(BUILD)/ensemble_comp.o: $(SRC)/ensemble_comp.c $(SRC)/ensemble_comp.h
	$(CC) $(CFLAGS) -o $(BUILD)/ensemble_comp.o -c $(SRC)/ensemble_comp.c $(DEFV5) $(INCGRO) -I$(SVM)

clean:
	make clean -C $(SVM) && rm -f $(BUILD)/g_ensemble_comp.o $(BUILD)/ensemble_comp.o $(BUILD)/g_ensemble_comp
# clean does not rm -r the BUILD directory for safety purposes.