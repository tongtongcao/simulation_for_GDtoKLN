ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = g++
CXXFLAGS += -O -g -Wall -fPIC $(ROOTCFLAGS) -I$(ROOTSYS)/include
LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS)

EXEC = gendpi0
OBJS = JGenPhaseSpace.o JGenPhotonEnergy.o JGenFermiMomentum.o gendpi0.o
OBJSL = JGenPhaseSpace.o JGenPhotonEnergy.o JGenFermiMomentum.o
LDLIBS =$(ROOTLIBS)


all:    $(OBJS) 
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJS) $(LDLIBS)
$(OBJS): %.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	-rm $(OBJS) $(EXEC)
