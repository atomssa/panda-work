#ifndef AnaTdav2_H
#define AnaTdav2_H

#include "FairTask.h"
#include "PndAnalysis.h"

class RhoCandList;

class AnaTdav2 : public FairTask{

 public:

  AnaTdav2();
  ~AnaTdav2();

  virtual InitStatus Init();

  virtual void Exec(Option_t* opt);

  virtual void Finish();

  void set_verbosity(int v) {verb = v;}

 private:

  int verb;
  int nevt;

  PndAnalysis *fAna;

  ClassDef(AnaTdav2,1);

};

#endif /*AnaTdav2_H*/
