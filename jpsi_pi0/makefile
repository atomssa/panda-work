INCDIR   = $(shell root-config --incdir)
LIBS     = $(shell root-config --libs)
CXXFLAGS = -std=c++0x -Wall

gen: fill.cpp
	g++ fill.cpp -o fill $(CXXFLAGS) -I$(INCDIR) $(LIBS) -lEG
