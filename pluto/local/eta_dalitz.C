void eta_dalitz() {

  gSystem->Load("libHist.so"); // for TH1
  gSystem->Load("libPhysics.so"); // for TLorentzVector

  gSystem->Load("libPluto.so");

  //Create a file for the NTuple:
  TFile *f = new TFile("ntuple.root", "RECREATE");

  //Create an NTuple with several variables:
  TNtuple *ntuple = new TNtuple("ntuple", "data from Pluto events", "eta_px:eta_py:eta_pz:opang");

  //Create a control histo:
  TH1F * histo1 = new TH1F ("histo1","dilepton mass with opening angle < 9deg",100,0.0,0.7);

  //Define the reaction: pp -> pp eta at 3.5 kinetic beam energy:
  PReaction my_reaction("3.5","p","p", "p p eta [g dilepton [e+ e-]]","eta_dalitz");

  //Inline calculating of a variable (here: Theta of e+ in degree):
  my_reaction.Do("theta_ep = ([e+]->Theta() * 180.)/TMath::Pi()");
  //...same for e-:
  my_reaction.Do("theta_em = ([e-]->Theta() * 180.)/TMath::Pi()");

  //Example of a reaction filter: Put a # in front of the variable name:
  my_reaction.Do("#filter = 1; if theta_ep<18 || theta_ep>85 || theta_em<18 || theta_em>85; #filter = 0");

  //Calculate the variables for the NTuple:
  my_reaction.Do("eta_px = [eta]->Px() ; eta_py = [eta]->Py() ; eta_pz = [eta]->Pz();");

  //Opening angle in lab frame:
  my_reaction.Do("opang = [e+]->Angle([e-])");

  //The NTuple should be written here:
  my_reaction.Output(ntuple);

  //The (conditional) histogram. Prepare the axis (_x)
  my_reaction.Do(histo1,"if opang > (9./180.)*TMath::Pi(); _x = ([e+] + [e-])->M()");

  //Run event loop:
  my_reaction.Loop(1000);

}
