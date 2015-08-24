INSTALL = /usr/local/bin
INCGRO = -I/usr/local/gromacs/include/gromacs
SVM = ../../libsvm-3.20
LINKGRO = -L/usr/local/gromacs/lib

svmanalys: svmanalys.o
	g++ -o svmanalys svmanalys.o $(SVM)/svm.o $(LINKGRO) -lgmx

install: svmanalys
	install svmanalys $(INSTALL)

.PHONY: install

svmanalys.o: svmanalys.c svmanalys.h
	gcc -c svmanalys.c $(INCGRO) -I$(SVM)