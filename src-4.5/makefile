INSTALL=/usr/local/bin
INCGRO=/usr/local/gromacs/include/gromacs
SVM=../../libsvm-3.20
LINKGRO=/usr/local/gromacs/lib

svmanalys: svmanalys.o
	g++ -o svmanalys svmanalys.o $(SVM)/svm.o -L$(LINKGRO) -lgmx

install: svmanalys
	install svmanalys $(INSTALL)

.PHONY: install

svmanalys.o: svmanalys.c svmanalys.h
	gcc -c svmanalys.c -I$(INCGRO) -I$(SVM)