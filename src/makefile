INSTALL = /usr/local/bin
INCGRO = /usr/local/gromacs/include/gromacs
SVM = ../../libsvm-3.20
LINKGRO = /usr/local/gromacs/lib

svmanalys: svmanalys.o svmutils.o
	g++ -o svmanalys svmanalys.o svmutils.o $(SVM)/svm.o -L$(LINKGRO) -lgmx

install: svmanalys
	install svmanalys $(INSTALL)

.PHONY: install

gsvmanalys: svmanalys.c svmutils.c
	gcc -g -o gsvmanalys svmanalys.c svmutils.c -I$(INCGRO) -L$(LINKGRO) -lgmx

comanalys: comanalys.o svmutils.o
	gcc -o comanalys comanalys.o svmutils.o -L$(LINKGRO) -lgmx

svmanalys.o: svmanalys.c svmutils.h
	gcc -c svmanalys.c -I$(INCGRO) -I$(SVM)

comanalys.o: comanalys.c svmutils.h
	gcc -c comanalys.c -I$(INCGRO)

svmio.o: svmio.c svmio.h
	gcc -c svmio.c -I$(INCGRO)

svmutils.o: svmutils.c svmutils.h
	gcc -c svmutils.c -I$(INCGRO)