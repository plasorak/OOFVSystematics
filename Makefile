# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile
DEBUGFLAGS      = -g -O0

# Include machine specific flags and locations (inc. files & libs)
#
include $(T2KREWEIGHT)/make/Make.include

#If mode variable is debug, compile in debug mode
#to set the mode, type 'make mode=debug'
ifeq ($(mode),debug)
CXXFLAGS+=$(DEBUGFLAGS)
endif


#all: CalculateOOFVTheory plotOOFVRateUncertainty
all: plotOOFVRateUncertainty
# Please only add new apps to be compiled by make all if you have checked 
# they still compile with following config:
# ./configure --disable-jnubeam --disable-neut --disable-genie --disable-oaanalysis;
# throwParmsFromJnuBeam 

CalculateOOFVTheory: FORCE
	@echo -e "\n\n***** Compiling CalculateOOFVTheory.cxx\n"
	$(CXX) $(CXXFLAGS) -c CalculateOOFVTheory.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) CalculateOOFVTheory.o $(LIBRARIES) -o CalculateOOFVTheory.exe

plotOOFVRateUncertainty: FORCE
	@echo -e "\n\n***** Compiling plotOOFVRateUncertainty.cxx\n"
	$(CXX) $(CXXFLAGS) -c plotOOFVRateUncertainty.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) plotOOFVRateUncertainty.o $(LIBRARIES) -o plotOOFVRateUncertainty.exe


purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) *.exe

FORCE:

# DO NOT DELETE
