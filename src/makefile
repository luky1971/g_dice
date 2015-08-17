INSTALL = /usr/local/bin
INCGRO = /usr/local/gromacs/include/gromacs
SVM = ../../libsvm-3.20
LINKGRO = /usr/local/gromacs/lib

svmanalys: svmanalys.o svutils.o svio.o
	g++ -o svmanalys svmanalys.o svutils.o svio.o $(SVM)/svm.o -L$(LINKGRO) -lgmx

install: svmanalys
	install svmanalys $(INSTALL)

.PHONY: install

gsvmanalys: svmanalys.c svutils.c svutils.h svio.c svio.h
	gcc -g -o gsvmanalys svmanalys.c svutils.c svio.c -I$(INCGRO) -I$(SVM) -L$(LINKGRO) -lgmx

comanalys: comanalys.o svutils.o
	gcc -o comanalys comanalys.o svutils.o -L$(LINKGRO) -lgmx

svmanalys.o: svmanalys.c svutils.h svio.h
	gcc -c svmanalys.c -I$(INCGRO) -I$(SVM)

comanalys.o: comanalys.c svutils.h
	gcc -c comanalys.c -I$(INCGRO) -I$(SVM)

svutils.o: svutils.c svutils.h
	gcc -c svutils.c -I$(INCGRO) -I$(SVM)

svio.o: svio.c svio.h
	gcc -c svio.c -I$(INCGRO) -I$(SVM)