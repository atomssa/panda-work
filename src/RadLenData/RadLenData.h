#ifndef RadLenData_H
#define RadLenData_H 1

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

class RadLenData : public FairTask
{

public:
	
  // ** Default constructor   
  RadLenData();
	
  // ** Destructor 
  ~RadLenData();	
	
  // ** Virtual method Init 
  virtual InitStatus Init();
	
  // ** Virtual method Exec 
  virtual void Exec(Option_t* opt);
	
  virtual void Finish();

protected:

  int nevt;
  static const int ndet =20;
  
  TClonesArray *rlp_array;
  TClonesArray *track_array;

  std::vector<TProfile*> vProf;
  std::vector<TH1*> vHistCumul;

  std::vector<TProfile2D*> vProf2d;
  std::vector<TH2*> vHistCumul2d;  
  
  ClassDef(RadLenData,1);
  
};

#endif /* RadLenData_H */
