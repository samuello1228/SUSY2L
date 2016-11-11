# Use the following on IRIX 6.5 and Linux 2.4
# Last tested with p15.06.01
#
ROOTLIBS = `root-config --libs`
OPTCOMP = `root-config --cflags` 

SPECIALFLAGS= --exceptions
CFLAGS = -I. 
LFLAGS =

RCXX = $(CFLAGS) $(ROOTCFLAGS)
DBG = -g -L -O -m32
CC = g++ ${DBG} $(RCXX) $(OPTCOMP) 
F77 = f77

all: run

run: run.o analysis.o analysis_Dict.o
	$(CC) run.o analysis.o analysis_Dict.o $(ROOTLIBS) -o run
	\rm -fr *.o
	\rm -fr *Dict*

run.o: run.cpp 
	$(CC) -c run.cpp -o run.o

analysis.o: analysis.cpp analysis.h 
	$(CC) -c analysis.cpp -o analysis.o

analysis_Dict.cxx: analysis.h
	rootcint -f analysis_Dict.cxx -c analysis.h

analysis_Dict.o: analysis_Dict.cxx 
	$(CC) -c analysis_Dict.cxx -o analysis_Dict.o

clean:
	\rm -f *.o
	\rm -fr fortran/*.o
	\rm -fr *~
	\rm -fr analysis_Dict.h
	\rm -fr *.cxx
	\rm -fr run
