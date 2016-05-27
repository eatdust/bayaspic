ext=$(shell uname | cut -c1-3)



ifeq ($(ext),Lin)
FC=mpif90
FCFLAGS= -O2 -fopenmp -DASPIC -DMPI
#-DMPISCHED -DASPICQ
F08FLAGS = -ffree-line-length-none
INCLUDE=-I/usr/local/include/aspic
#-I/usr/local/include/aspicq
LFLAGS= -lblas -llapack -L/usr/local/lib -laspic
#-laspicq

##ifort specific
#FC=ifort
#FCFLAGS= -O -DASPIC
#F08FLAGS="-DNOF08"
#INCLUDE=-I/opt/intel/aspic/include/aspic
#LFLAGS= -mkl -L/opt/intel/aspic/lib -laspic
#CPP=cpp -P -traditional-cpp $(FCFLAGS)
endif 


LIKEDIR = ../fastlikes
NESTDIR = ../multinest
CHORDDIR = ../polychord/src/polychord

OBJSLIKE =  rbfprec.o shepprec.o qshepmdata.o iond.o rbfnd.o rbflike.o \
	qshepmd.o shepnd.o sheplike.o

#with -DASPIC
OBJSASPIC = aspicvars.o aspicmodels.o aspicpriors.o wraspic.o 

OBJSMPI = fifo.o iofifo.o scheduler.o

OBJSNEST = utilsMN.o utils1.o priorsMN.o kmeans_clstr.o \
	xmeans_clstr.o posterior.o nested.o

OBJSCHORD = utilsPC.o abort.o settings.o mpi_utils.o random_utils.o  \
	array_utils.o calculate.o run_time_info.o params.o priorsPC.o \
	read_write.o feedback.o chordal_sampling.o clustering.o \
	generate.o nested_sampling.o ini.o interfaces.o 

#nestwrap.o
OBJFILES = sampl.o $(OBJSASPIC) $(OBJSLIKE) $(OBJSNEST) $(OBJSCHORD) \
	nestparams.o chordparams.o samplaspic.o

OBJSMIN = sampl.o $(OBJSLIKE) $(OBJSNEST) $(OBJSCHORD) \
	nestparams.o chordparams.o samplsr.o

FCFLAGS+=$(INCLUDE)


bayaspic.$(ext): $(OBJFILES) $(OBJSMPI) bayaspic.o
	$(FC) $(FCFLAGS) $(OBJFILES) $(OBJSMPI) bayaspic.o -o ../$@ $(LFLAGS) 

bayslowroll.$(ext): $(OBJSMIN) bayslowroll.o
	$(FC) $(FCFLAGS) $(OBJSMIN) bayslowroll.o -o ../$@ $(LFLAGS) 

aspicmodels.o: aspicmodels.h aspicmodels.pp
aspicpriors.o: aspicpriors.h aspicpriors.pp
nestwrap.o: nestparams.o wraspic.o
chordparams.o: settings.o
samplsr.o: settings.o interfaces.o

%.o: $(LIKEDIR)/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: $(LIKEDIR)/%.F90
	$(FC) $(FCFLAGS)  -c$<

%.o: $(NESTDIR)/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: $(NESTDIR)/%.F90
	$(FC) $(FCFLAGS) $(F08FLAGS) -c $<

%.o: $(CHORDDIR)/%.f90
	$(FC) $(FCFLAGS) $(F08FLAGS) -c $<

%.o: $(CHORDDIR)/%.F90
	$(FC) $(FCFLAGS) $(F08FLAGS) -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

ifeq ($(FC),ifort)
%.ifort: %.F08
	$(CPP) $< -o $@
%.o: %.ifort
	$(FC) $(FCFLAGS) -c -free -Tf$< 
else
%.o: %.F08
	$(FC) $(FCFLAGS) $(F08FLAGS) -c $<
endif

clean:
	rm -f *.ifort *.o *.mod ../bayaspic.$(ext) ../bayslowroll.$(ext)
