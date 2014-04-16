#ifndef ThreePions_H
#define ThreePions_H 1

#include "FairTask.h"
#include "TLorentzVector.h"

#include <vector>

class TTree;
class TClonesArray;
class TObjectArray;
class TProfile;
class TProfile2D;
class TH2;
class TH1;

class ThreePions : public FairTask
{

public:
	
  // ** Default constructor   
  ThreePions();
	
  // ** Destructor 
  ~ThreePions();	
	
  // ** Virtual method Init 
  virtual InitStatus Init();
	
  // ** Virtual method Exec 
  virtual void Exec(Option_t* opt);
	
  virtual void Finish();

protected:

  int nevt;
  
  ClassDef(ThreePions,1);
  
};

#endif /* ThreePions_H */
