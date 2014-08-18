// $Id: example02LinkDef.h 537 2011-04-20 09:42:14Z ivana $

//------------------------------------------------
// The Virtual Monte Carlo examples
// Copyright (C) 2007, 2008 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: vmc@pcroot.cern.ch
//-------------------------------------------------

/// \file example02LinkDef.h
/// \brief The CINT link definitions for example E02 classes

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
 
#pragma link C++ class  Ex02MCApplication+;
#pragma link C++ class  Ex02MCStack+;
#pragma link C++ class  Ex02Particle+;
#pragma link C++ class  Ex02ChamberParameterisation+;
#pragma link C++ class  Ex02DetectorConstruction+;
#pragma link C++ class  Ex02DetectorConstructionOld+;
#pragma link C++ class  Ex02MagField+;
#pragma link C++ class  Ex02TrackerHit+;
#pragma link C++ class  Ex02TrackerSD+;
#pragma link C++ class  Ex02RootManager+;
#pragma link C++ class  std::stack<Ex02Particle*,deque<Ex02Particle*> >+;

#pragma link C++ enum   FileMode;

#endif





