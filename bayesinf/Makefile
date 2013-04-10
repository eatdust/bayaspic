ext=$(shell uname | cut -c1-3)

ifeq ($(ext),Lin)
FC=gfortran
FCFLAGS= -O -fopenmp 
LFLAGS= -lblas -llapack
endif 

OBJ=rbfprec.o iorbf.o rbfnd.o rbflike.o

RBFDIR = ../rbflikes
NESTDIR = ..

OBJSRBF =  rbfprec.o iorbf.o rbfnd.o rbflike.o
OBJSNEST = utils.o utils1.o priors.o kmeans_clstr.o \
	xmeans_clstr.o posterior.o nested.o


OBJFILES = $(OBJSRBF) $(OBJSNEST) nestprec.o nestparams.o nestwrap.o

infbayes.$(ext): $(OBJFILES) bayesinf.o
	$(FC) $(FCFLAGS) $(OBJFILES) bayesinf.o -o ../$@ $(LFLAGS) 

%.o: $(RBFDIR)/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: $(RBFDIR)/%.F90
	$(FC) $(FCFLAGS)  -c$<

%.o: $(NESTDIR)/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: $(NESTDIR)/%.F90
	$(FC) $(FCFLAGS)  -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod ../bayesinf.$(ext)
