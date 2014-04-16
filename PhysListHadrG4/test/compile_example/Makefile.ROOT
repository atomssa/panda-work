# Generic makefile for ROOT programs.
# Elwin Dijck, December 2009

# Based on the ROOT Makefile by Fons Rademakers.

#------------------------------------------------------------------------------

# Load proper defaults for this architecture.
include $(ROOTSYS)/test/Makefile.arch

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(LibSuf) .$(DllSuf) .$(ExeSuf)

# Some extra warnings.
CXXFLAGS     += -ansi -pedantic -Wall -Wcast-align -Wcast-qual \
		-Wctor-dtor-privacy -Winit-self -Wnon-virtual-dtor \
		-Wno-long-long -Wpointer-arith -Wshadow -Wstrict-aliasing \
		-Wuninitialized -Wwrite-strings

# Too much warnings from PandaRoot itself
# -Wextra -Wold-style-cast -Woverloaded-virtual

# General compilation rule.
.$(SrcSuf).$(ObjSuf):
	@echo "Compiling $<..."
	@$(CXX) $(CXXFLAGS) $(INCLUDE) $< -MM -MF $(basename $@).dep -MT $@
	@$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<
