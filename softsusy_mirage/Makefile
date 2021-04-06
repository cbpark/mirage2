
# Makefile for SOFTSUSY library and test program
# Ben Allanach

.KEEP_STATE:

# Choose your c++ compiler here: if you want to use the link to fortran,
# it's probably going to have to be g++ on linux
CC=	g++
# fortran compiler
FF=     gfortran

.SUFFIXES:	.o .cpp .h .f
.PRECIOUS:	.cpp .h .f libsoft.a

#Optimisation level, eg: -O2 -march=i686, speeds by a factor of 2 on my machine
OPT=
#OR debug level: -g(n=1,2,3)
DEBUG=
CFLAGS= -I. -L. $(DEBUG) $(OPT) $(shell gsl-config --cflags)

.cpp.o:
	$(CC) -c $(CFLAGS) $<
.cpp.a:
	$(CC) -c $(CFLAGS) $<
	ar r $@ $*.o; rm $*.o
.f.o:
	$(FF) $(OPT) -c $<
.f.a:
	$(FF) $(OPT) -c $<
	ar r $@ $*.o; rm $*.o

# Chooses whether to use ranlib - not if your system is a linux one
all: softsusy.x softpoint.x
	case `uname` in \
	   Linux) RANL=;;\
	   OSF1) CFLAGS="$(CFLAGS) -ieee";;\
	   *) RANL="ranlib libnr.a";;\
	   esac

# Example C++ Softsusy calling program
softsusy.x: Makefile main.o libsoft.a
	$(CC) -o softsusy.x $(CFLAGS) main.o $(shell gsl-config --libs) -lsoft;\
	touch softsusy.x

# Command-line program interface
softpoint.x: Makefile softpoint.o libsoft.a
	$(CC) -o softpoint.x $(CFLAGS) softpoint.o $(shell gsl-config --libs) -lsoft;\
	touch softpoint.x

# twoloop higgs not deleted because it takes so long to compile
clean:
	ar x libsoft.a twoloophiggs.o; \
	rm *.a main.o sabine.o softpoint.o *.x *~; \
	ar r libsoft.a twoloophiggs.o; rm twoloophiggs.o

libsoft.a: libsoft.a(rge.o) libsoft.a(utils.o) \
	libsoft.a(linalg.o) libsoft.a(lowe.o) \
	libsoft.a(susy.o) libsoft.a(softsusy.o) \
	libsoft.a(mycomplex.o) libsoft.a(numerics.o) \
	libsoft.a(physpars.o) libsoft.a(softpars.o) libsoft.a(twoloophiggs.o)
	$(RANL)

backup:
	uufiles < .backupinstructions;

softpoint.o: softpoint.cpp linalg.h utils.h softsusy.h def.h
main.o: main.cpp linalg.h utils.h softsusy.h def.h
libsoft.a(numerics.o): numerics.cpp numerics.h linalg.h def.h
libsoft.a(rge.o): rge.cpp rge.h def.h
libsoft.a(softpars.o): softpars.cpp softpars.h susy.h def.h linalg.h
libsoft.a(physpars.o): physpars.cpp physpars.h
libsoft.a(mycomplex.o): mycomplex.cpp mycomplex.h utils.h
libsoft.a(utils.o): utils.cpp utils.h linalg.h
libsoft.a(linalg.o): linalg.cpp linalg.h utils.h
libsoft.a(lowe.o): lowe.cpp lowe.h linalg.h rge.h def.h
libsoft.a(susy.o): susy.cpp susy.h lowe.h utils.h linalg.h rge.h def.h
libsoft.a(utils.o): utils.cpp utils.h linalg.h def.h
libsoft.a(softsusy.o): softsusy.cpp mycomplex.h softsusy.h susy.h lowe.h \
	linalg.h rge.h def.h utils.h numerics.h softpars.h
libsoft.a(twoloopshiggs.o): twoloophiggs.f twoloophiggs.h
