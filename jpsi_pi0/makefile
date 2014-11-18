INCDIR        = $(shell root-config --incdir)
LIBS          = $(shell root-config --libs)
CXXFLAGS      = -std=c++0x -Wall
FILLER_INCDIR = ${VMCWORKDIR}/src/filler
FILLER_LIBS   = -L${VMCWORKDIR}/build/lib -lfiller

gen: macros/fill.cpp
	g++ macros/fill.cpp -o fill $(CXXFLAGS) -I$(INCDIR) $(LIBS) -I$(FILLER_INCDIR) $(FILLER_LIBS) -lEG
