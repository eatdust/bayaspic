<<<<<<< HEAD
FC = mpif90
CC = mpicc
CXX = mpiCC
FFLAGS += -O3 -DMPI
CFLAGS += -O3 -DMPI

LAPACKLIB = -llapack

NESTLIBDIR = ./

export FC CC CXX FFLAGS CFLAGS LAPACKLIB

 
AR = ar r  
LINKLIB = ld -shared
 
NSOBJECTS = utils.o utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^ 

%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $^ 

default: libnest3.a

all: libnest3.a obj_detect eggboxC eggboxC++ gaussian gauss_shell \
rosenbrock himmelblau ackley
 
libnest3.so: $(NSOBJECTS) 
	$(LINKLIB) -o $(LIBS) $@ $^ 
 
libnest3.a: $(NSOBJECTS) 
	$(AR) $@ $^ 
 
obj_detect:
	make -C example_obj_detect
 
gaussian:
	make -C example_gaussian
 
rosenbrock:
	make -C example_rosenbrock
 
ackley:
	make -C example_ackley
 
himmelblau:
	make -C example_himmelblau
 
gauss_shell:
	make -C example_gauss_shell
	
eggboxC:
	make -C example_eggbox_C
	
eggboxC++:
	make -C example_eggbox_C++

clean: 
	-rm $(NESTLIBDIR)/libnest3.*  *.o *.mod
	
cleanall: clean_exec clean clean_obj_detect clean_gaussian clean_gauss_shell \
clean_example_eggbox_C clean_example_eggbox_C++ clean_rosenbrock clean_himmelblau \
clean_ackley

clean_exec:
	-rm obj_detect gaussian rosenbrock ackley himmelblau gauss_shell eggboxC eggboxC++

clean_obj_detect:
	make -C example_obj_detect clean
	
clean_gaussian:
	make -C example_gaussian clean
	
clean_rosenbrock:
	make -C example_rosenbrock clean
	
clean_ackley:
	make -C example_ackley clean
	
clean_himmelblau:
	make -C example_himmelblau clean
	
clean_gauss_shell:
	make -C example_gauss_shell clean
	
clean_example_eggbox_C:
	make -C example_eggbox_C clean
	
clean_example_eggbox_C++:
	make -C example_eggbox_C++ clean
=======
ext=$(shell uname | cut -c1-3)

ifeq ($(ext),Lin)
FC=gfortran
FCFLAGS= -O3 -fopenmp
LFLAGS= -lblas -llapack
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
>>>>>>> 7a167c0b1ba76de790dfcce4e4aa20ad066ce294
