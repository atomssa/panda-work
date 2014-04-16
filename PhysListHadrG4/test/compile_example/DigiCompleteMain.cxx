// Example of a compiled program to perform a full simulation run: a compiled
// version of the digi_complete_tpc.C macro.
//
// Elwin Dijck, December 2009
// JGM, January 2010

#include "DigiComplete.h"

#include <TRint.h>
#include <TROOT.h>

#include <iostream>
#include <sstream>


using namespace std;


template <typename Type>
Type StrTo(char const *str);


int main(int argc, char *argv[])
{
    // Start application.
    gROOT->SetBatch();
    TApplication app("DigiComplete", &argc, argv, 0, -1);

    // Read the command line options or provide default values.
    
    TString inFile   = argc > 1 ? argv[1] : "sim_complete.root";
    TString parFile  = argc > 2 ? argv[2] : "simparams.root";
    TString digiFile = argc > 3 ? argv[3] : "all.par";
    TString outFile  = argc > 4 ? argv[4] : "digi_complete.root";
    
    cout << boolalpha;

    cout << endl << "Starting digi simulation with TPC:" << endl
        << "    input file             : " << inFile << endl
        << "    output file            : " << outFile << endl
        << "    params file            : " << parFile << endl
        << "    input digi params file : " << digiFile << endl << endl;

    // Start the simulation.
    DigiComplete(inFile, parFile, digiFile, outFile);
    
    return 0;
}


// Helper function to convert an argument to any type.
template <typename Type>
Type StrTo(char const *str)
{
    istringstream iss(str);

    iss >> boolalpha;  // Use true/false instead of 1/0 for booleans.

    Type result;
    iss >> result;

    return result;
}
