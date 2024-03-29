ext=$(shell uname | cut -c1-3)



ifeq ($(ext),Lin)
FC = mpif90
FCFLAGS = -O2  -fopenmp -DASPIC -DLIKENULL -DMPISCHED
#-DASPICQ -DNOFANN -DLIKENULL -DLIKESHEP -DLIKEFANN -DLIKERBF -DRRAD -DRHOW
#-DMPISCHED or -DMPI

#extra flags for F08 cpp generated source files
F08FLAGS = -ffree-line-length-none

#for multinest
FMNFLAGS = -ffree-line-length-none

#for polychord
FPCFLAGS = -ffree-line-length-none

INCLUDE =-I/usr/local/include/aspic
#-I/usr/local/include/aspicq

#if the fann library is present, build the NN likelihood
LDFANN = -ldoublefann

LFLAGS= -lblas -llapack -L/usr/local/lib -laspic $(LDFANN)
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

ifneq ($(LDFANN),)
FANNLIKE = fnnprec.o fann.o fnnlike.o
FCFLAGS += -DDOUBLEFANN
else
FCFLAGS += -DNOFANN
endif

OBJSLIKE =  mlprec.o ioml.o iocmb.o rbfprec.o iorbf.o rbfuncs.o rbfnd.o rbflike.o \
	shepprec.o qshepmdata.o ioshep.o qshepmd.o shepnd.o sheplike.o $(FANNLIKE)


#with -DASPIC
OBJSASPIC = aspicvars.o aspicnonstd.o aspicmodels.o aspicpriors.o wraspic.o 

OBJSMPI = fifo.o iofifo.o scheduler.o

OBJSNEST = utilsMN.o utils1.o priorsMN.o kmeans_clstr.o \
	xmeans_clstr.o posterior.o nested.o

OBJSCHORD = utilsPC.o abort.o settings.o mpi_utils.o random_utils.o  \
	array_utils.o calculate.o run_time_info.o params.o priorsPC.o \
	read_write.o feedback.o chordal_sampling.o clustering.o \
	generate.o nelder_mead.o maximiser.o nested_sampling.o ini.o interfaces.o 

#nestwrap.o
OBJFILES = sampl.o $(OBJSASPIC) $(OBJSLIKE) $(OBJSNEST) $(OBJSCHORD) \
	nestparams.o chordparams.o samplaspic.o

OBJSMIN = sampl.o $(OBJSLIKE) $(OBJSNEST) $(OBJSCHORD) \
	nestparams.o chordparams.o samplsr.o

FCFLAGS += $(INCLUDE)


bayaspic.$(ext): $(OBJFILES) $(OBJSMPI) bayaspic.o
	$(FC) $(FCFLAGS) $(OBJFILES) $(OBJSMPI) bayaspic.o -o ../$@ $(LFLAGS) 

bayslowroll.$(ext): $(OBJSMIN) bayslowroll.o
	$(FC) $(FCFLAGS) $(OBJSMIN) bayslowroll.o -o ../$@ $(LFLAGS) 

aspicmodels.o: aspicmodels.h aspicmodels.pp
aspicpriors.o: aspicpriors.h aspicpriors.pp
#nestwrap.o: nestparams.o wraspic.o
chordparams.o: settings.o
samplsr.o: settings.o interfaces.o

%.o: $(LIKEDIR)/common/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: $(LIKEDIR)/commom/%.F90
	$(FC) $(FCFLAGS)  -c$<

%.o: $(LIKEDIR)/rbf/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: $(LIKEDIR)/rbf/%.F90
	$(FC) $(FCFLAGS)  -c$<

%.o: $(LIKEDIR)/shep/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: $(LIKEDIR)/shep/%.F90
	$(FC) $(FCFLAGS)  -c$<

%.o: $(LIKEDIR)/fnn/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: $(LIKEDIR)/fnn/%.F90
	$(FC) $(FCFLAGS)  -c$<

%.o: $(NESTDIR)/%.f90
	$(FC) $(FCFLAGS) $(FMNFLAGS) -c $<

%.o: $(NESTDIR)/%.F90
	$(FC) $(FCFLAGS) $(FMNFLAGS) -c $<

%.o: $(CHORDDIR)/%.f90
	$(FC) $(FCFLAGS) $(FPCFLAGS) -c $<

%.o: $(CHORDDIR)/%.F90
	$(FC) $(FCFLAGS) $(FPCFLAGS) -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

ifeq ($(FC),ifort)
%.ifort: $(LIKEDIR)/fnn/%.F03
	$(CPP) $< -o $@
%.ifort: $(LIKEDIR)/fnn/%.f03
	$(CPP) $< -o $@
%.ifort: %.F08
	$(CPP) $< -o $@
%.ifort: %.F03
	$(CPP) $< -o $@
%.ifort: %.f03
	$(CPP) $< -o $@
%.o: %.ifort
	$(FC) $(FCFLAGS) -c -free -Tf$< 
else
%.o: $(LIKEDIR)/fnn/%.F03
	$(FC) $(FCFLAGS) -c $<
%.o: $(LIKEDIR)/fnn/%.f03
	$(FC) $(FCFLAGS) -c $<
%.o: %.F08
	$(FC) $(FCFLAGS) $(F08FLAGS) -c $<
%.o: %.F03
	$(FC) $(FCFLAGS) $(F08FLAGS) -c $<
%.o: %.f03
	$(FC) $(FCFLAGS) $(F08FLAGS) -c $<
endif

clean:
	rm -f *.ifort *.o *.mod ../bayaspic.$(ext) ../bayslowroll.$(ext)
