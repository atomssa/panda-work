ROOT_INC        = $(shell root-config --incdir)
ROOT_LIBS       = $(shell root-config --libs)
CXXFLAGS        = -Wall
CXXFLAGS_GPROF  = -Wall -pg

gen: HBeam.h pion_generator.cpp
	g++ pion_generator.cpp -o gen $(CXXFLAGS) -I$(ROOT_INC) $(ROOT_LIBS) -lMinuit

gen_gprof: HBeam.h pion_generator.cpp
	g++ pion_generator.cpp -o gen_gprof $(CXXFLAGS_GPROF) -I$(ROOT_INC) $(ROOT_LIBS) -lMinuit
