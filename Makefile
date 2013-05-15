ext=$(shell uname | cut -c1-3)

ifeq ($(ext),Lin)
FC=gfortran
FCFLAGS= -O3 -fopenmp -ffast-math -march=amdfam10  -fprefetch-loop-arrays -funroll-loops
LFLAGS= -lblas -llapack
#FC=ifort
#FCFLAGS = -O3 -xHost -ipo -openmp -mkl
#LFLAGS =
endif 

OBJ=rbfprec.o iorbf.o rbfnd.o rbflike.o inoutfile.o


rbfmain.$(ext): $(OBJ) rbfmain.o
	$(FC) $(FCFLAGS) $(OBJ) rbfmain.o -o $@ $(LFLAGS) 

%.o: %.f90
	$(FC) $(FCFLAGS) -c $*.f90

%.o: %.F90
	$(FC) $(FCFLAGS) -c $*.F90

%.o: %.f95
	$(FC) $(FCFLAGS) -c $*.f95

%.o: %.F95
	$(FC) $(FCFLAGS) -c $*.F95

clean:
	-rm -f *.o *.a *.d core *.mod *.$(ext)
