INSTALL=/usr/local/bin
INCGRO=-I/usr/local/gromacs/include/gromacs/utility -I/usr/local/gromacs/include/gromacs/fileio -I/usr/local/gromacs/include/gromacs/commandline
SVM=../../libsvm-3.20
LINKGRO=-L/usr/local/gromacs/lib/i386-linux-gnu
DEFV5=-D GRO_V5  

svmanalys: svmanalys.o
	g++ -o svmanalys svmanalys.o $(SVM)/svm.o $(LINKGRO) -lgromacs

install: svmanalys
	install svmanalys $(INSTALL)

.PHONY: install

svmanalys.o: svmanalys.c svmanalys.h
	gcc -c svmanalys.c $(INCGRO) -I$(SVM) $(DEFV5)