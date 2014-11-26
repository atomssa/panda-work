#include "AnaTdav2.h"

#include "FairTask.h"

using namespace std;

AnaTdav2::AnaTdav2():
  verb(0),
  nevt(0),
  fAna()
{

}

AnaTdav2::~AnaTdav2() {

}

InitStatus AnaTdav2::Init() {

  cout << "AnaTdav2::Init" << endl;
  fAna = new PndAnalysis();

}

void AnaTdav2::Exec(Option_t* opt) {
  if (nevt++%1000==0)
    cout << "======== AnaTdav2::Exec evt " << nevt << " ======== " << endl;

  fAna->GetEvent();

}

void AnaTdav2::Finish() {
  cout << "AnaTdav2::Exec" << endl;
  fAna->Reset();
}
