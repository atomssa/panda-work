# Makefile for the SimComplete program, example of a macro compiled with g++.
# Elwin Dijck, December 2009

# Based on the ROOT Makefile by Fons Rademakers.

#------------------------------------------------------------------------------

include Makefile.Panda

#------------------------------------------------------------------------------

SIMCOMPLETE := SimComplete$(ExeSuf)
SIMCOMPLETEO := SimComplete.$(ObjSuf) SimCompleteDict.$(ObjSuf)
SIMCOMPLETESO := libSimComplete.$(DllSuf)

DIGICOMPLETE := DigiComplete$(ExeSuf)
DIGICOMPLETEO := DigiComplete.$(ObjSuf) DigiCompleteDict.$(ObjSuf)
DIGICOMPLETESO := libDigiComplete.$(DllSuf)

RECOCOMPLETE := RecoComplete$(ExeSuf)
RECOCOMPLETEO := RecoComplete.$(ObjSuf) RecoCompleteDict.$(ObjSuf)
RECOCOMPLETESO := libRecoComplete.$(DllSuf)

all:            $(SIMCOMPLETE) $(DIGICOMPLETE) $(RECOCOMPLETE)

lib:            $(SIMCOMPLETESO) $(DIGICOMPLETESO) $(RECOCOMPLETESO)

clean:
		@rm -f *.$(ObjSuf) *Dict.h *Dict.cxx *.dep

distclean:      clean
		@rm -f $(SIMCOMPLETE) $(DIGICOMPLETE) $(RECOCOMPLETE) *.$(DllSuf)

$(SIMCOMPLETE): SimCompleteMain.$(ObjSuf) $(SIMCOMPLETEO)
		@echo "Linking $@..."
		@$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(DIGICOMPLETE): DigiCompleteMain.$(ObjSuf) $(DIGICOMPLETEO)
		@echo "Linking $@..."
		@$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(RECOCOMPLETE): RecoCompleteMain.$(ObjSuf) $(RECOCOMPLETEO)
		@echo "Linking $@..."
		@$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(SIMCOMPLETESO): $(SIMCOMPLETEO)
		@echo "Linking $@..."
		@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(LIBS)
		@echo "$@ done"

$(DIGICOMPLETESO): $(DIGICOMPLETEO)
		@echo "Linking $@..."
		@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(LIBS)
		@echo "$@ done"

$(RECOCOMPLETESO): $(RECOCOMPLETEO)
		@echo "Linking $@..."
		@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(LIBS)
		@echo "$@ done"

#------------------------------------------------------------------------------

SimCompleteDict.$(SrcSuf): SimComplete.h SimCompleteLinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $(INCLUDE) $^

DigiCompleteDict.$(SrcSuf): DigiComplete.h DigiCompleteLinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $(INCLUDE) $^

RecoCompleteDict.$(SrcSuf): RecoComplete.h RecoCompleteLinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $(INCLUDE) $^

include $(wildcard *.dep)
