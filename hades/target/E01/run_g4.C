// $Id: run_g4.C 341 2008-05-26 11:04:57Z ivana $

//------------------------------------------------
// The Virtual Monte Carlo examples
// Copyright (C) 2007, 2008 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: vmc@pcroot.cern.ch
//-------------------------------------------------

/// \ingroup E01
/// \file E01/run_g4.C
/// \brief Macro for running Example01 with Geant4. 

void run_g4(const TString& configMacro = "g4Config.C") 
{
/// Macro function for running Example01 with Geant4 from
/// Root interactive session
/// \param configMacro configuration macro name, default \ref E01/g4Config.C 

  // Load basic libraries
  gROOT->LoadMacro("macro/basiclibs.C");
  basiclibs();
  
  // Load Geant4 libraries
  gROOT->LoadMacro("macro/g4libs.C");
  g4libs();
  
  // Load this example library
  gSystem->Load("lib/tgt_linuxx8664gcc/libexample01");

  // MC application
  Ex01MCApplication* appl 
    = new Ex01MCApplication("Example01", "The example01 MC application");

  // Initialize MC
  appl->InitMC(configMacro);

  // Run MC
  appl->RunMC(2);
  
  delete appl;
}  
