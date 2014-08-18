#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
 
int NumericalRootFinder()
{
 
   // Create the function and wrap it
  
  TF1 f("ThirdOrderPol", "[0]+[1]*x+[2]*x*x+[3]*x*x", TMath::PiOver2(), TMath::TwoPi() );
   for (int i=0; i<4; ++i) f.SetParameter(
   ROOT::Math::WrappedTF1 wf1(f);
 
   // Create the Integrator
   ROOT::Math::BrentRootFinder brf;
 
   // Set parameters of the method
   brf.SetFunction( wf1, TMath::PiOver2(), TMath::TwoPi() );
   brf.Solve();
 
   cout << brf.Root() << endl;
 
   return 0;
}
