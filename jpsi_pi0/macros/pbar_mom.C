void pbar_mom() {

  double mp=0.937;
  vector<double> pbarp;
  pbarp.push_back(5.519);
  pbarp.push_back(5.513);
  pbarp.push_back(8);
  pbarp.push_back(12);
  for (int i=0; i<pbarp.size(); ++i) {
    double s = 2*mp*(1.0 + sqrt(pow(mp,2)+pow(pbarp[i],2)));
    double sqrts = sqrt(s);
    double epbar = TMath::Hypot(mp,pbarp[i]);
    double beta = pbarp[i]/(epbar+mp);
    double epbar_recalc = (s-2*mp)/2.0/mp;
    double pbar_mom_recalc = sqrt( pow(epbar_recalc,2) - pow(mp,2) );
    cout << "pbar_mom = " << setw(6) << pbarp[i] << " s= " << setw(6) << s << " #sqrt{s}= " << setw(6) << sqrts << " #beta= " << setw(6) << beta
	 << "pbar_mom_rec = " << setw(6) << pbar_mom_recalc << endl;
  }

}
