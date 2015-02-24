//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RhoCandList                                                          //
//                                                                      //
// Container class for RhoCandidates                                    //
//                                                                      //
// Author List:                                                         //
// Marcel Kunze,  RUB, Feb. 99                                          //
// Copyright (C) 1999-2001, Ruhr-University Bochum.                     //
// Ralf Kliemt, HIM/GSI Feb.2013 (Cleanup & Restructuring)              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include "RhoBase/RhoCandList.h"
#include "RhoBase/RhoFactory.h"
#include "RhoBase/RhoCandidate.h"
#include "RhoBase/RhoParticleSelectorBase.h"
#include "RhoBase/RhoVertexSelectorBase.h"

ClassImp ( RhoCandList )

#include <iostream>
using namespace std;

RhoCandList::RhoCandList ( const char* name, UInt_t capacity ) :
  TNamed ( name,name ), fFast ( kFALSE )
{
  fOwnList = new TObjArray ( capacity );
}

// Perform a deep copy
RhoCandList::RhoCandList (RhoCandList& l)
{
  fFast = l.fFast;
  fOwnList = new TObjArray( l.GetLength() );

  Cleanup();
  const Int_t n = l.GetLength();
  for ( int i=0; i<n; i++ ) {
    Put ( l.Get ( i ) );
  }
}

// Perform a deep copy - const argument version to play well with stl containers
RhoCandList::RhoCandList (const RhoCandList& l)
{
  fFast = l.fFast;
  fOwnList = new TObjArray( l.GetLength() );

  Cleanup();
  const Int_t n = l.GetLength();
  for ( int i=0; i<n; i++ ) {
    Put ( l.GetConst ( i ) );
  }
}

//--------------
// Destructor --
//--------------
// Deletes the list and the owned objects
RhoCandList::~RhoCandList( )
{
  Cleanup();
  delete fOwnList;
}

//----------------------
//-- public Functions --
//----------------------

void RhoCandList::Cleanup( )
{
    fOwnList->Clear(); //Delete() would destruct the candidates which will be cleaned by the TFactory!
}

void RhoCandList::SetNumberOfTracks ( Int_t n )
{
  cerr << "RhoCandList::SetNumberOfTracks is deprecated" << endl;
}

Int_t RhoCandList::GetNumberOfTracks() const
{
  return fOwnList->GetLast() +1;
}

void RhoCandList::Put ( const RhoCandidate* c, Int_t i )
{
  RhoCandidate* newCand = RhoFactory::Instance()->NewCandidate ( c );
  if ( i<0 ) {
    fOwnList->Add ( newCand );
  } else {
    ( *fOwnList ) [i] = newCand;
  }
}

void RhoCandList::InsertAt ( Int_t i, const RhoCandidate* c )
{
  fOwnList->AddAtAndExpand ( ( TObject* ) c,i );
  Put ( c,i );
}


RhoCandidate* RhoCandList::Get ( Int_t i )
{
  if ( i >= GetNumberOfTracks() ) { return 0; }
  return ( RhoCandidate* ) ( fOwnList->UncheckedAt ( i ) );
}

RhoCandidate* RhoCandList::GetConst ( Int_t i ) const
{
  if ( i >= GetNumberOfTracks() ) { return 0; }
  return ( RhoCandidate* ) ( fOwnList->UncheckedAt ( i ) );
}

RhoCandList* RhoCandList::GetFittedList()
{
  // Create and return a new list of the fitted candidates

  // new list with extended name from the unfitted and intelligent starting size
  RhoCandList* fittedlist = new RhoCandList( Form("%s Fitted",fName.Data()) , GetNumberOfTracks() );
  FillFittedList(*fittedlist);
  return fittedlist;
}

void RhoCandList::FillFittedList(RhoCandList &fittedlist)
{
  RhoCandidate* aFit=0;
  fittedlist.Cleanup();
  for(int j=0;j<GetNumberOfTracks();j++) {
    aFit = GetConst(j)->GetFit();
    if(!aFit) Error("RhoCandList::FillFittedList","Fit pointer does not exist!");
    fittedlist.Put(aFit);
  }
  return;
}


// Compare the marker and remove corresponding entry (MK,12/99)
// This allows to remove objects in several lists
Int_t RhoCandList::Remove ( RhoCandidate* c )
{
  Int_t nRemoved = 0;
  Int_t n = GetNumberOfTracks();
  for ( Int_t i=0; i<n; i++ ) {
    RhoCandidate* b = Get ( i );
    if ( b->Equals ( c ) ) {
      fOwnList->RemoveAt ( i );
      nRemoved++;
    }
  }
  fOwnList->Compress();
  return nRemoved;
}

Int_t RhoCandList::RemoveFamily ( RhoCandidate* c )
{
  Int_t nRemoved = 0;
  Int_t n = GetNumberOfTracks();
  for ( Int_t i=0; i<n; i++ ) {
    RhoCandidate* b = Get ( i );
    if ( b->Overlaps ( c ) ) {
      fOwnList->RemoveAt ( i );
      nRemoved++;
    }
  }
  fOwnList->Compress();
  return nRemoved;
}

Int_t RhoCandList::RemoveClones()
{
  Int_t nRemoved = 0;
  Int_t n = GetNumberOfTracks();
  for ( Int_t i=0; i<n-1; i++ ) {
    RhoCandidate* b = Get ( i );
    if ( b==0 ) { continue; }
    for ( Int_t j=i+1; j<n; ++j ) {
      RhoCandidate* c = Get ( j );
      if ( c==0 ) { continue; }
      if ( b->Equals ( c ) ) {
        fOwnList->RemoveAt ( j );
        nRemoved++;
      }
    }
  }
  fOwnList->Compress();
  return nRemoved;
}

Int_t RhoCandList::OccurrencesOf ( RhoCandidate* c )
{
  Int_t nCand = 0;
  const Int_t n = GetNumberOfTracks();
  for ( Int_t i=0; i<n; i++ ) {
    RhoCandidate* b = Get ( i );
    if ( b->Equals ( c ) ) { nCand++; }
  }
  return nCand;
}

Double_t RhoCandList::GetTotalEnergy ( Double_t emin )
{
  Double_t e = 0.0;
  const Int_t n = GetNumberOfTracks();
  for ( Int_t i=0; i<n; i++ ) {
    RhoCandidate* c = Get ( i );
    Double_t energy = c->Energy();
    if ( energy > emin ) { e += energy; }
  }
  return e;
}

TVector3 RhoCandList::GetTotalMomentum ( Double_t pmin )
{
  TVector3 p ( 0.0,0.0,0.0 );
  const Int_t n = GetNumberOfTracks();
  for ( Int_t i=0; i<n; i++ ) {
    RhoCandidate* c = Get ( i );
    TVector3 p3 = c->P3();
    if ( p3.Mag() > pmin ) { p = p + p3; }
  }
  return p;
}

void RhoCandList::Boost ( const TVector3& p )
{
  const Int_t n = GetNumberOfTracks();
  for ( Int_t i=0; i<n; i++ ) {
    Get ( i )->Boost ( p );
  }
}

void RhoCandList::PrintOn ( std::ostream& o ) const
{
  const Int_t n = GetNumberOfTracks();
  o << "content: " << n << " entries." << endl;
  if ( n!=0 ) {
    for ( Int_t i=0; i<n; i++ ) {
      const RhoCandidate* c = GetConst ( i );
      o << "[" << i << "] " << *c;
      o << endl;
    }
  }
}

void RhoCandList::Remainder ( RhoCandList& l )
{
  const Int_t n = l.GetNumberOfTracks();
  for ( int i=0; i<n; i++ ) {
    const RhoCandidate* c = l.Get ( i );
    if ( c->GetMarker() !=0 ) { fOwnList->RemoveAt ( i ); }
  }
  fOwnList->Compress();
}

// Makes now a deep copy rather than just cpoying pointers (MK,12/99)
void RhoCandList::operator = ( const RhoCandList& l )
{
  fFast = l.fFast;

  Cleanup();
  const Int_t n = l.GetNumberOfTracks();
  for ( int i=0; i<n; i++ ) {
    Put ( l.GetConst( i ) );
  }
}

RhoCandidate* RhoCandList::operator[] ( Int_t i )
{
  return Get(i);
}

void RhoCandList::SetType ( const TParticlePDG* pdt )
{
  for ( Int_t i=0; i<GetNumberOfTracks(); i++ ) {
    Get(i)->SetType ( pdt );
  }
}

void RhoCandList::SetType ( const char* name )
{
  for ( Int_t i=0; i<GetNumberOfTracks(); i++ ) {
    Get(i)->SetType ( name );
  }
}

void RhoCandList::SetType ( Int_t pdgcode )
{
  for ( Int_t i=0; i<GetNumberOfTracks(); i++ ) {
    Get(i)->SetType ( pdgcode );
  }
}


void RhoCandList::Combine ( RhoCandList& l1, RhoCandList& l2,RhoVertexSelectorBase* s )
{
  Cleanup();
  CombineAndAppend ( l1, l2, s );
}

void RhoCandList::Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3,RhoVertexSelectorBase* s )
{
  Cleanup();
  CombineAndAppend ( l1, l2, l3, s );
}

void RhoCandList::Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4,RhoVertexSelectorBase* s )
{
  Cleanup();
  CombineAndAppend ( l1, l2, l3, l4, s );
}

void RhoCandList::Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, RhoCandList& l5,RhoVertexSelectorBase* s )
{
  Cleanup();
  CombineAndAppend ( l1, l2, l3, l4, l5, s );
}

void RhoCandList::CombineAndAppend ( RhoCandList& l1, RhoCandList& l2,RhoVertexSelectorBase* selector )
{
  //printf("RhoCandList::CombineAndAppend()\n");
  TLorentzVector vl;
  Double_t charge;
  Bool_t nearby = kTRUE;

  const int len1=l1.GetLength();
  const int len2=l2.GetLength();

  int i1,i2;
  int st2;

  for ( i1=0; i1<len1; i1++ ) {
    st2=0;
    if ( &l1==&l2 ) { st2=i1+1; }

    for ( i2=st2; i2<len2; i2++ ) {
      if ( l1[i1]->Overlaps ( l2[i2] ) ) { continue; }

      vl=l1[i1]->P4() +l2[i2]->P4();
      charge=l1[i1]->Charge() +l2[i2]->Charge();
      if ( selector ) { nearby = selector->Accept ( l1[i1],l2[i2] ); }
      if ( !nearby ) { continue; }

      RhoCandidate c ( vl,charge );
      c.SetCovP4 ( l1[i1]->P4Cov() +l2[i2]->P4Cov() );

      c.SetMarker ( l1[i1]->GetMarker ( 0 ) |l2[i2]->GetMarker ( 0 ),0 );
      c.SetMarker ( l1[i1]->GetMarker ( 1 ) |l2[i2]->GetMarker ( 1 ),1 );
      c.SetMarker ( l1[i1]->GetMarker ( 2 ) |l2[i2]->GetMarker ( 2 ),2 );
      c.SetMarker ( l1[i1]->GetMarker ( 3 ) |l2[i2]->GetMarker ( 3 ),3 );


      if ( selector!=0 ) {
        c.SetPosition ( selector->GetVertex() );
        c.SetVect ( selector->GetMomentum() );
        c.SetEnergy ( c.E() );
      }

      Put ( &c );
      // after putting (does a copy and drops daughter links)
      RhoCandidate* cInList = Get(GetLength()-1);
      cInList->AddDaughterLinkSimple(l1[i1]) ;
      cInList->AddDaughterLinkSimple(l2[i2]) ;
    }
  }

}

void RhoCandList::CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3,RhoVertexSelectorBase* selector )
{
  TLorentzVector vl;
  Double_t charge;
  Bool_t nearby = kTRUE;

  const int len1=l1.GetLength();
  const int len2=l2.GetLength();
  const int len3=l3.GetLength();

  int i1,i2,i3;
  int st2,st3;

  for ( i1=0; i1<len1; i1++ ) {
    st2=0;
    if ( &l2==&l1 ) { st2=i1+1; }

    for ( i2=st2; i2<len2; i2++ ) {
      if ( l1[i1]->Overlaps ( l2[i2] ) ) { continue; }

      st3=0;
      if ( &l3==&l2 ) { st3=i2+1; }
      else if ( &l3==&l1 ) { st3=i1+1; }

      for ( i3=st3; i3<len3; i3++ ) {
        if ( l3[i3]->Overlaps ( l2[i2] ) || l3[i3]->Overlaps ( l1[i1] ) ) { continue; }

        vl=l1[i1]->P4() +l2[i2]->P4() +l3[i3]->P4();
        charge=l1[i1]->Charge() +l2[i2]->Charge() +l3[i3]->Charge();
        if ( selector ) { nearby = selector->Accept ( l1[i1],l2[i2],l3[i3] ); }
        if ( !nearby ) { continue; }

        RhoCandidate c ( vl,charge );
        c.SetCovP4 ( l1[i1]->P4Cov() +l2[i2]->P4Cov() +l3[i3]->P4Cov() );

        c.SetMarker ( l1[i1]->GetMarker ( 0 ) |l2[i2]->GetMarker ( 0 ) |l3[i3]->GetMarker ( 0 ),0 );
        c.SetMarker ( l1[i1]->GetMarker ( 1 ) |l2[i2]->GetMarker ( 1 ) |l3[i3]->GetMarker ( 1 ),1 );
        c.SetMarker ( l1[i1]->GetMarker ( 2 ) |l2[i2]->GetMarker ( 2 ) |l3[i3]->GetMarker ( 2 ),2 );
        c.SetMarker ( l1[i1]->GetMarker ( 3 ) |l2[i2]->GetMarker ( 3 ) |l3[i3]->GetMarker ( 3 ),3 );

        Put ( &c );
        // after putting (does a copy and drops daughter links)
        RhoCandidate* cInList = Get(GetLength()-1);
        cInList->AddDaughterLinkSimple(l1[i1]) ;
        cInList->AddDaughterLinkSimple(l2[i2]) ;
        cInList->AddDaughterLinkSimple(l3[i3]) ;
      }
    }
  }
}

void RhoCandList::CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4,RhoVertexSelectorBase* selector )
{
  TLorentzVector vl;
  Double_t charge;
  Bool_t nearby = kTRUE;

  const int len1=l1.GetLength();
  const int len2=l2.GetLength();
  const int len3=l3.GetLength();
  const int len4=l4.GetLength();

  int i1,i2,i3,i4;
  int st2,st3,st4;

  for ( i1=0; i1<len1; i1++ ) {
    st2=0;
    if ( &l2==&l1 ) { st2=i1+1; }

    for ( i2=st2; i2<len2; i2++ ) {
      if ( l1[i1]->Overlaps ( l2[i2] ) ) { continue; }

      st3=0;
      if ( &l3==&l2 ) { st3=i2+1; }
      else if ( &l3==&l1 ) { st3=i1+1; }

      for ( i3=st3; i3<len3; i3++ ) {
        if ( l3[i3]->Overlaps ( l2[i2] ) || l3[i3]->Overlaps ( l1[i1] ) ) { continue; }

        st4=0;
        if ( &l4==&l3 ) { st4=i3+1; }
        else if ( &l4==&l2 ) { st4=i2+1; }
        else if ( &l4==&l1 ) { st4=i1+1; }

        for ( i4=st4; i4<len4; i4++ ) {
          if ( l4[i4]->Overlaps ( l3[i3] ) || l4[i4]->Overlaps ( l2[i2] ) || l4[i4]->Overlaps ( l1[i1] ) ) { continue; }

          vl=l1[i1]->P4() +l2[i2]->P4() +l3[i3]->P4() +l4[i4]->P4();
          charge=l1[i1]->Charge() +l2[i2]->Charge() +l3[i3]->Charge() +l4[i4]->Charge();
          if ( selector ) { nearby = selector->Accept ( l1[i1],l2[i2],l3[i3],l4[i4] ); }
          if ( !nearby ) { continue; }

          RhoCandidate c ( vl,charge );
          c.SetCovP4 ( l1[i1]->P4Cov() +l2[i2]->P4Cov() +l3[i3]->P4Cov() +l4[i4]->P4Cov() );

          c.SetMarker ( l1[i1]->GetMarker ( 0 ) |l2[i2]->GetMarker ( 0 ) |l3[i3]->GetMarker ( 0 ) |l4[i4]->GetMarker ( 0 ),0 );
          c.SetMarker ( l1[i1]->GetMarker ( 1 ) |l2[i2]->GetMarker ( 1 ) |l3[i3]->GetMarker ( 1 ) |l4[i4]->GetMarker ( 1 ),1 );
          c.SetMarker ( l1[i1]->GetMarker ( 2 ) |l2[i2]->GetMarker ( 2 ) |l3[i3]->GetMarker ( 2 ) |l4[i4]->GetMarker ( 2 ),2 );
          c.SetMarker ( l1[i1]->GetMarker ( 3 ) |l2[i2]->GetMarker ( 3 ) |l3[i3]->GetMarker ( 3 ) |l4[i4]->GetMarker ( 3 ),3 );

          Put ( &c );
          // after putting (does a copy and drops daughter links)
          RhoCandidate* cInList = Get(GetLength()-1);
          cInList->AddDaughterLinkSimple(l1[i1]) ;
          cInList->AddDaughterLinkSimple(l2[i2]) ;
          cInList->AddDaughterLinkSimple(l3[i3]) ;
          cInList->AddDaughterLinkSimple(l4[i4]) ;
        }
      }
    }
  }
}

void RhoCandList::CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, RhoCandList& l5,RhoVertexSelectorBase* selector )
{
  TLorentzVector vl;
  Double_t charge;
  Bool_t nearby = kTRUE;

  const int len1=l1.GetLength();
  const int len2=l2.GetLength();
  const int len3=l3.GetLength();
  const int len4=l4.GetLength();
  const int len5=l5.GetLength();

  int i1,i2,i3,i4,i5;
  int st2,st3,st4,st5;

  for ( i1=0; i1<len1; i1++ ) {
    st2=0;
    if ( &l2==&l1 ) { st2=i1+1; }

    for ( i2=st2; i2<len2; i2++ ) {
      if ( l1[i1]->Overlaps ( l2[i2] ) ) { continue; }

      st3=0;
      if ( &l3==&l2 ) { st3=i2+1; }
      else if ( &l3==&l1 ) { st3=i1+1; }

      for ( i3=st3; i3<len3; i3++ ) {
        if ( l3[i3]->Overlaps ( l2[i2] ) || l3[i3]->Overlaps ( l1[i1] ) ) { continue; }

        st4=0;
        if ( &l4==&l3 ) { st4=i3+1; }
        else if ( &l4==&l2 ) { st4=i2+1; }
        else if ( &l4==&l1 ) { st4=i1+1; }

        for ( i4=st4; i4<len4; i4++ ) {
          if ( l4[i4]->Overlaps ( l3[i3] ) || l4[i4]->Overlaps ( l2[i2] ) || l4[i4]->Overlaps ( l1[i1] ) ) { continue; }

          st5=0;
          if ( &l5==&l4 ) { st5=i4+1; }
          else if ( &l5==&l3 ) { st5=i3+1; }
          else if ( &l5==&l2 ) { st5=i2+1; }
          else if ( &l5==&l1 ) { st5=i1+1; }

          for ( i5=st5; i5<len5; i5++ ) {
            if ( l5[i5]->Overlaps ( l4[i4] ) || l5[i5]->Overlaps ( l3[i3] )
                 || l5[i5]->Overlaps ( l2[i2] ) || l5[i5]->Overlaps ( l1[i1] ) ) { continue; }

            vl=l1[i1]->P4() +l2[i2]->P4() +l3[i3]->P4() +l4[i4]->P4() +l5[i5]->P4();
            charge=l1[i1]->Charge() +l2[i2]->Charge() +l3[i3]->Charge() +l4[i4]->Charge() +l5[i5]->Charge();
            if ( selector ) { nearby = selector->Accept ( l1[i1],l2[i2],l3[i3],l4[i4],l5[i5] ); }
            if ( !nearby ) { continue; }

            RhoCandidate c ( vl,charge );
            c.SetCovP4 ( l1[i1]->P4Cov() +l2[i2]->P4Cov() +l3[i3]->P4Cov() +l4[i4]->P4Cov() +l5[i5]->P4Cov() );

            c.SetMarker ( l1[i1]->GetMarker ( 0 ) |l2[i2]->GetMarker ( 0 )
                          |l3[i3]->GetMarker ( 0 ) |l4[i4]->GetMarker ( 0 ) |l5[i5]->GetMarker ( 0 ),0 );
            c.SetMarker ( l1[i1]->GetMarker ( 1 ) |l2[i2]->GetMarker ( 1 )
                          |l3[i3]->GetMarker ( 1 ) |l4[i4]->GetMarker ( 1 ) |l5[i5]->GetMarker ( 1 ),1 );
            c.SetMarker ( l1[i1]->GetMarker ( 2 ) |l2[i2]->GetMarker ( 2 )
                          |l3[i3]->GetMarker ( 2 ) |l4[i4]->GetMarker ( 2 ) |l5[i5]->GetMarker ( 2 ),2 );
            c.SetMarker ( l1[i1]->GetMarker ( 3 ) |l2[i2]->GetMarker ( 3 )
                          |l3[i3]->GetMarker ( 3 ) |l4[i4]->GetMarker ( 3 ) |l5[i5]->GetMarker ( 3 ),3 );

            Put ( &c );
            // after putting (does a copy and drops daughter links)
            RhoCandidate* cInList = Get(GetLength()-1);
            cInList->AddDaughterLinkSimple(l1[i1]) ;
            cInList->AddDaughterLinkSimple(l2[i2]) ;
            cInList->AddDaughterLinkSimple(l3[i3]) ;
            cInList->AddDaughterLinkSimple(l4[i4]) ;
            cInList->AddDaughterLinkSimple(l5[i5]) ;
          }
        }
      }
    }
  }
}

std::ostream&  operator << ( std::ostream& o, const RhoCandList& a )
{
  a.PrintOn ( o );
  return o;
}

//extern "C" void qsort(void *, size_t, size_t, int (const void *,const void *));

typedef int compare_function ( const void*, const void* );

void RhoCandList::Sort ( int ( *compfunc ) ( const RhoCandidate**, const RhoCandidate** ) )
{
  qsort ( fOwnList, GetNumberOfTracks(), sizeof ( void* ), ( compare_function* ) compfunc );
}


void RhoCandList::Select ( RhoCandList& l, Bool_t ( *selfunc ) ( RhoCandidate* ) )
{
  Cleanup();
  const Int_t n = l.GetNumberOfTracks();
  for ( Int_t i=0; i<n; i++ ) {
    RhoCandidate* c = Get ( i );
    if ( selfunc ( c ) ) {
      Put ( c );
    }
  }
}

// Destructive selection

void RhoCandList::Select ( RhoParticleSelectorBase* pidmgr )
{
  const Int_t n = GetNumberOfTracks();
  for ( Int_t i=0; i<n; i++ ) {
    RhoCandidate* c = Get ( i );
    if ( !pidmgr->Accept ( c ) ) { fOwnList->RemoveAt ( i ); }
  }
  fOwnList->Compress();
}

// Non-destructive selection

void RhoCandList::Select ( RhoCandList& l, RhoParticleSelectorBase* pidmgr )
{
  Cleanup();
  const Int_t n = l.GetNumberOfTracks();
  for ( Int_t i=0; i<n; i++ ) {
    RhoCandidate* c = l.Get ( i );
    if ( pidmgr->Accept ( c ) ) {
      Put ( c );
    }
  }
}

void RhoCandList::Append ( RhoCandList& l, RhoParticleSelectorBase* pidmgr )
{
  const Int_t n = l.GetNumberOfTracks();
  for ( Int_t i=0; i<n; i++ ) {
    RhoCandidate* c = l.Get ( i );
    if ( 0==pidmgr || pidmgr->Accept ( c ) ) {
      Put ( c );
    }
  }
}
