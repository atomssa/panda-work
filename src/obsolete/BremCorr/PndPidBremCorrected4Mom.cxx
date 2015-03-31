//////////////////////////////////////////////////////////////////////////
//                                                                      //
// PndPidBremCorrected4Mom	                                        //
//                                                                      //
// Container for Bremstrahlung radiaton corrected momentum              //
//                                                                      //
// Author: Klaus Goetzen, GSI, 12.06.08		                        //
// Copyright (C) 2008, GSI Darmstadt.		         		//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "PndPidBremCorrected4Mom.h"

PndPidBremCorrected4Mom::PndPidBremCorrected4Mom(): 
  fXmomentum(0.),
  fYmomentum(0.),
  fZmomentum(0.),
  fEnergy(0.)
{  
}

PndPidBremCorrected4Mom::PndPidBremCorrected4Mom(TLorentzVector &p4) :
  fXmomentum(p4.X()),
  fYmomentum(p4.Y()),  
  fZmomentum(p4.Z()),
  fEnergy(p4.E())
{
}

PndPidBremCorrected4Mom::~PndPidBremCorrected4Mom() {

}

// Do we need to propagate the covariance matrix?

//TMatrixD& PndPidBremCorrected4Mom::P4Cov() const
//{
//  static TMatrixD covP4(4,4);
//	
//  covP4(0,0) = fErrP7[18];  covP4(1,0) = fErrP7[19]; covP4(1,1) = fErrP7[20];
//  covP4(2,0) = fErrP7[21];  covP4(2,1) = fErrP7[22]; covP4(2,2) = fErrP7[23];
//  covP4(3,0) = fErrP7[24];  covP4(3,1) = fErrP7[25]; covP4(3,2) = fErrP7[26]; 
//  covP4(3,3) = fErrP7[27];
//	
//  for (int i=0; i<3; i++)
//    for (int j=i+1; j<4; j++)
//      covP4(i,j)=covP4(j,i);
//			
//  return covP4;
//}

//void PndPidBremCorrected4Mom::SetCov7(const TMatrixD &cov7 )
//{
//  // position error
//    
//  fErrP7[0] = cov7(0,0);  fErrP7[1] = cov7(1,0); fErrP7[2] = cov7(1,1);  
//  fErrP7[3] = cov7(2,0);  fErrP7[4] = cov7(2,1); fErrP7[5] = cov7(2,2);
//    
//  // position-momentum covariance
//    
//  fErrP7[6] = cov7(3,0);   fErrP7[7] = cov7(3,1);  fErrP7[8] = cov7(3,2);
//  fErrP7[9] = cov7(4,0);   fErrP7[10] = cov7(4,1); fErrP7[11] = cov7(4,2);
//  fErrP7[12] = cov7(5,0);  fErrP7[13] = cov7(5,1); fErrP7[14] = cov7(5,2);
//  fErrP7[15] = cov7(6,0);  fErrP7[16] = cov7(6,1); fErrP7[17] = cov7(6,2);
//    
//  // momentum error
//  fErrP7[18] = cov7(3,3);  fErrP7[19] = cov7(4,3); fErrP7[20] = cov7(4,4);
//  fErrP7[21] = cov7(5,3);  fErrP7[22] = cov7(5,4); fErrP7[23] = cov7(5,5);
//  fErrP7[24] = cov7(6,3);  fErrP7[25] = cov7(6,4); fErrP7[26] = cov7(6,5); 
//  fErrP7[27] = cov7(6,6);
//}

//void PndPidBremCorrected4Mom::SetP4Cov(const TMatrixD &covP4 )
//{
//  // position error
//    
//  fErrP7[0] = 0;  fErrP7[1] = 0; fErrP7[2] = 0;  
//  fErrP7[3] = 0;  fErrP7[4] = 0; fErrP7[5] = 0;
//    
//  // position-momentum covariance
//    
//  fErrP7[6] = 0;   fErrP7[7] = 0;  fErrP7[8] = 0;
//  fErrP7[9] = 0;   fErrP7[10] = 0; fErrP7[11] = 0;
//  fErrP7[12] = 0;  fErrP7[13] = 0; fErrP7[14] = 0;
//  fErrP7[15] = 0;  fErrP7[16] = 0; fErrP7[17] = 0;
//    
//  // momentum error
//  fErrP7[18] = covP4(0,0);  fErrP7[19] = covP4(1,0); fErrP7[20] = covP4(1,1);
//  fErrP7[21] = covP4(2,0);  fErrP7[22] = covP4(2,1); fErrP7[23] = covP4(2,2);
//  fErrP7[24] = covP4(3,0);  fErrP7[25] = covP4(3,1); fErrP7[26] = covP4(3,2); 
//  fErrP7[27] = covP4(3,3);
//}


ClassImp(PndPidBremCorrected4Mom)
