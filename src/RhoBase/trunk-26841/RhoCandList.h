#ifndef RHOCANDLIST_H
#define RHOCANDLIST_H
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
#include <iostream>
#include "TNamed.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TVector.h"
#include "TVector3.h"
#include "TParticlePDG.h"

class RhoCandidate;
class RhoParticleSelectorBase;
class RhoVertexSelectorBase;

class RhoCandList : public TNamed
{

  public:
    //Constructor
    RhoCandList ( const char* name="RhoCandList",UInt_t capacity=512 );
    RhoCandList ( RhoCandList& b );
    RhoCandList ( const RhoCandList&);
    //Destructor
    virtual ~RhoCandList();
    void Cleanup();

    void SetNumberOfTracks ( Int_t );
    void SetNumberOfCandidates ( Int_t n ) {
      SetNumberOfTracks ( n );
    }
    Int_t GetNumberOfTracks() const;
    Int_t GetNumberOfCandidates() const {
      return GetNumberOfTracks();
    }
    Int_t GetLength() const {
      return GetNumberOfTracks();
    }
    void Add ( const RhoCandidate* c ) {
      Put ( c );
    }
    void Append ( const RhoCandidate* c ) {
      Put ( c );
    }
    void Append ( RhoCandList& l, RhoParticleSelectorBase* pidmgr=0 );
    void Put ( const RhoCandidate*, Int_t i = -1 );
    void InsertAt ( Int_t i, const RhoCandidate* c );

	// without type
    void Combine ( RhoCandList& l1, RhoCandList& l2 );
    void Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3 );   //added 06/08 K.Goetzen
    void Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4 );//added 06/08 K.Goetzen
    void Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, RhoCandList& l5 );//added 06/08 K.Goetzen

    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2 );
    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3 );//added 06/08 K.Goetzen
    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4 );//added 06/08 K.Goetzen
    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, RhoCandList& l5 );//added 06/08 K.Goetzen

	// with type as TParticlePDG
	void Combine ( RhoCandList& l1, RhoCandList& l2, const TParticlePDG* pdt ) { Combine(l1, l2); SetType(pdt); }
    void Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, const TParticlePDG* pdt ) { Combine(l1, l2, l3); SetType(pdt); }
    void Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, const TParticlePDG* pdt ) { Combine(l1, l2, l3, l4); SetType(pdt); }
    void Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, RhoCandList& l5, const TParticlePDG* pdt ) { Combine(l1, l2, l3, l4, l5); SetType(pdt); }

    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, const TParticlePDG* pdt ) { int s = GetLength(); CombineAndAppend(l1, l2); SetType(pdt, s); }
    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, const TParticlePDG* pdt ) { int s = GetLength();CombineAndAppend(l1, l2, l3); SetType(pdt, s); }
    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, const TParticlePDG* pdt ) { int s = GetLength(); CombineAndAppend(l1, l2, l3, l4); SetType(pdt, s); }
    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, RhoCandList& l5, const TParticlePDG* pdt )
		{ int s = GetLength(); CombineAndAppend(l1, l2, l3, l4, l5); SetType(pdt, s); }

	// with type as char* (name)
	void Combine ( RhoCandList& l1, RhoCandList& l2, const char* name ) { Combine(l1, l2); SetType(name); }
    void Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, const char* name ) { Combine(l1, l2, l3); SetType(name); }
    void Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, const char* name ) { Combine(l1, l2, l3, l4); SetType(name); }
    void Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, RhoCandList& l5, const char* name ) { Combine(l1, l2, l3, l4, l5); SetType(name); }

    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, const char* name ) { int s = GetLength(); CombineAndAppend(l1, l2); SetType(name, s); }
    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, const char* name ) { int s = GetLength();CombineAndAppend(l1, l2, l3); SetType(name, s); }
    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, const char* name ) { int s = GetLength(); CombineAndAppend(l1, l2, l3, l4); SetType(name, s); }
    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, RhoCandList& l5, const char* name )
		{ int s = GetLength(); CombineAndAppend(l1, l2, l3, l4, l5); SetType(name, s); }

	// with type as int (pdgtype)
	void Combine ( RhoCandList& l1, RhoCandList& l2, Int_t pdgcode ) { Combine(l1, l2); SetType(pdgcode); }
    void Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, Int_t pdgcode ) { Combine(l1, l2, l3); SetType(pdgcode); }
    void Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, Int_t pdgcode ) { Combine(l1, l2, l3, l4); SetType(pdgcode); }
    void Combine ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, RhoCandList& l5, Int_t pdgcode ) { Combine(l1, l2, l3, l4, l5); SetType(pdgcode); }

    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, Int_t pdgcode ) { int s = GetLength(); CombineAndAppend(l1, l2); SetType(pdgcode, s); }
    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, Int_t pdgcode ) { int s = GetLength();CombineAndAppend(l1, l2, l3); SetType(pdgcode, s); }
    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, Int_t pdgcode ) { int s = GetLength(); CombineAndAppend(l1, l2, l3, l4); SetType(pdgcode, s); }
    void CombineAndAppend ( RhoCandList& l1, RhoCandList& l2, RhoCandList& l3, RhoCandList& l4, RhoCandList& l5, Int_t pdgcode )
		{ int s = GetLength(); CombineAndAppend(l1, l2, l3, l4, l5); SetType(pdgcode, s); }

	// set type for complete list or starting from index 'start'
    void SetType ( const TParticlePDG* pdt, int start=0 );
    void SetType ( const char* name, int start=0 );
    void SetType ( Int_t pdgcode, int start=0 );

   void Select ( RhoParticleSelectorBase* pidmgr );
    void Select ( RhoCandList& l, RhoParticleSelectorBase* pidmgr );
    void Select ( RhoCandList& l, Bool_t ( *selfunc ) ( RhoCandidate* ) );
    Int_t OccurrencesOf ( RhoCandidate* );
    Int_t Remove ( RhoCandidate* ); // Returns #removed cands
    Int_t RemoveFamily ( RhoCandidate* ); // Returns #removed cands
    Int_t RemoveClones(); // Returns #removed cands
    void Reset() {
      SetNumberOfTracks ( 0 );
    }
    RhoCandidate* Get ( Int_t );
    RhoCandidate* GetConst ( Int_t ) const;

    RhoCandList* GetFittedList();
    void FillFittedList(RhoCandList &fittedlist);

    // Analysis functions

    void Boost ( const TVector3& );
    Double_t GetTotalEnergy ( Double_t emin=0.0 );
    TVector3 GetTotalMomentum ( Double_t pmin=0.0 );

    //operations which work on existing lists

    void PrintOn ( std::ostream& o=std::cout ) const;
    void Sort ( int ( *compfunc ) ( const RhoCandidate**, const RhoCandidate** ) );

    void Remainder ( RhoCandList& );

    void operator = ( const RhoCandList& );
    RhoCandidate* operator[] ( Int_t );

    void SetFast ( Bool_t yesno=kTRUE ) {
      fFast = yesno;
    }

  private:
    TObjArray*   fOwnList;  // This holds the candidates
    Bool_t   fFast;     // Fast mode = no error matrices

    ClassDef ( RhoCandList,1 ) //Container for RhoCandidates

};

// standalone print
std::ostream&  operator << ( std::ostream& o, const RhoCandList& );

#endif
