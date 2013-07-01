ext=$(shell uname | cut -c1-3)



ifeq ($(ext),Lin)
FC=mpif90
FCFLAGS= -O -fopenmp -DASPIC -DMPISCHED
F08FLAGS = -ffree-line-length-none
INCLUDE=-I/usr/local/include/aspicq
LFLAGS= -lblas -llapack -L/usr/local/lib -laspicq
endif 


LIKEDIR = ../fastlikes
NESTDIR = ..

OBJSLIKE =  rbfprec.o shepprec.o qshepmdata.o iond.o rbfnd.o rbflike.o \
	qshepmd.o shepnd.o sheplike.o

#with -DASPIC
OBJSASPIC = aspicvars.o aspicmodels.o aspicpriors.o wraspic.o 

OBJSMPI = fifo.o iofifo.o scheduler.o

OBJSNEST = utils.o utils1.o priors.o kmeans_clstr.o \
	xmeans_clstr.o posterior.o nested.o


OBJFILES = nestprec.o $(OBJSASPIC) $(OBJSLIKE) $(OBJSNEST) \
	nestparams.o nestwrap.o

FCFLAGS+=$(INCLUDE)


bayaspic.$(ext): $(OBJFILES) $(OBJSMPI) bayaspic.o
	$(FC) $(FCFLAGS) $(OBJFILES) $(OBJSMPI) bayaspic.o -o ../$@ $(LFLAGS) 

bayslowroll.$(ext): $(OBJFILES) bayslowroll.o
	$(FC) $(FCFLAGS) $(OBJFILES) bayslowroll.o -o ../$@ $(LFLAGS) 

aspicmodels.o: aspicmodels.h aspicnames.pp
aspicpriors.o: aspicpriors.h aspicpriors.pp
nestwrap.o: nestparams.o

%.o: $(LIKEDIR)/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: $(LIKEDIR)/%.F90
	$(FC) $(FCFLAGS)  -c$<

%.o: $(NESTDIR)/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: $(NESTDIR)/%.F90
	$(FC) $(FCFLAGS)  -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F08
	$(FC) $(FCFLAGS) $(F08FLAGS) -c $<

clean:
	rm -f *.o *.mod ../bayaspic.$(ext) ../bayslowroll.$(ext)
