TPad* padSetup3x2div(TCanvas *tc, int ii) {
  if (tc->GetPad(1) == 0) {
    tc->Divide(2,3);
  }
  return (TPad* ) tc->GetPad(ii+1);
}

TPad* padSetup2x3vert(TCanvas *tc, int ii) {

  int ix= (ii%2);
  int iy= (ii/2);

  double a= 0.15;
  double b= 0.05;
  double _a = a/(1-a);
  double _b = b/(1-b);
  double w = 1/(3+_a+_b);
  double xl = ix==1?0.49:0.0;
  double xh = ix==1?1.0:0.49;
  //double yl[3] = {0.0, (_a+1)*w, (_a+2)*w-0.003};
  //double yh[3] = {(_a+1)*w, (_a+2)*w, (_a+_b+3)*w};
  double yl[3] = {(_a+2)*w-0.003, (_a+1)*w, 0.0 };
  double yh[3] = {(_a+_b+3)*w, (_a+2)*w, (_a+1)*w};

  cout << "p " << ii << "= (" <<  ix << "," << iy << ") = (" << xl << ", " << yl[iy] << ", " << xh << ", " << yh[iy] << ")" << " h= " << yh[iy] - yl[iy] << endl;
  TPad *pad = new TPad(Form("pads%d",ii),Form("pads%d",ii), xl, yl[iy],  xh, yh[iy]);
  pad->SetBorderSize(0);

  tc->cd(0);
  pad->Draw();

  double epsilon=1e-9;
  if (ix==0) {pad->SetRightMargin(epsilon); pad->SetLeftMargin(0.15);}
  if (ix==1) {pad->SetRightMargin(0.1); pad->SetLeftMargin(0.15);}
  //.if (ix==1) {pad->SetRightMargin(epsilon); pad->SetLeftMargin(0.2);}
  //if ((ii%3)==0) {pad->SetRightMargin(epsilon); pad->SetLeftMargin(0.2);}
  //if ((ii%3)==1) {pad->SetLeftMargin(epsilon); pad->SetRightMargin(epsilon); }
  //if ((ii%3)==2) {pad->SetLeftMargin(epsilon); pad->SetRightMargin(0.1);}

  if (iy==0) {
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(epsilon);
  } else if (iy==1) {
    pad->SetTopMargin(epsilon);
    pad->SetBottomMargin(epsilon);
  } else {
    pad->SetTopMargin(epsilon);
    pad->SetBottomMargin(0.15);
  }

  pad->SetTicks(0,1);

  return pad;
}

TPad* padSetup3x2(TCanvas *tc, int ii) {

  double a= 0.15;
  double b= 0.05;
  double _a = a/(1-a);
  double _b = b/(1-b);
  double w = 1/(3+_a+_b);
  double yl = ii<3?0.53:0.0;
  double yh = ii<3?1.0:0.53;
  double xl[3] = {0.0, (_a+1)*w, (_a+2)*w-0.003};
  double xh[3] = {(_a+1)*w, (_a+2)*w, (_a+_b+3)*w};

  cout << "p " << ii << "(" << xl[ii%3] << ", " << yl << ", " << xh[ii%3] << ", " << yh << ")" << " w= " << xh[ii%3] - xl[ii%3] << endl;
  TPad *pad = new TPad(Form("pads%d",ii),Form("pads%d",ii), xl[ii%3], yl,  xh[ii%3], yh);
  pad->SetBorderSize(0);

  tc->cd(0);
  pad->Draw();

  double epsilon=1e-9;
  if ((ii%3)==0) {pad->SetRightMargin(epsilon); pad->SetLeftMargin(0.2);}
  if ((ii%3)==1) {pad->SetLeftMargin(epsilon); pad->SetRightMargin(epsilon); }
  if ((ii%3)==2) {pad->SetLeftMargin(epsilon); pad->SetRightMargin(0.1);}

  if (ii<3) {
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(epsilon);
  } else {
    pad->SetTopMargin(epsilon);
    pad->SetBottomMargin(0.15);
  }

  pad->SetTicks(0,1);

  return pad;
}

TPad* padSetup3x3(TCanvas *tc, int ii) {

  int ix= (ii%3);
  int iy= (ii/3);

  double a= 0.15;
  double b= 0.05;
  double _a = a/(1-a);
  double _b = b/(1-b);
  double w = 1/(3+_a+_b);

  double xl[3] = {0.0, (_a+1)*w, (_a+2)*w-0.003};
  double xh[3] = {(_a+1)*w, (_a+2)*w, (_a+_b+3)*w};

  double yl[3] = {(_a+2)*w-0.003, (_a+1)*w, 0.0};
  double yh[3] = {(_a+_b+3)*w, (_a+2)*w, (_a+1)*w};

  cout << "p " << ii << "(" << xl[ix] << ", " << yl[iy] << ", " << xh[ix] << ", " << yh[iy] << ")" << " w= " << xh[ix] - xl[ix] << endl;
  TPad *pad = new TPad(Form("pads%d",ii),Form("pads%d",ii), xl[ix], yl[iy],  xh[ix], yh[iy]);
  pad->SetBorderSize(0);

  tc->cd(0);
  pad->Draw();

  double epsilon=1e-9;
  if (ix==0) {pad->SetRightMargin(epsilon); pad->SetLeftMargin(0.2);}
  if (ix==1) {pad->SetLeftMargin(epsilon); pad->SetRightMargin(epsilon); }
  if (ix==2) {pad->SetLeftMargin(epsilon); pad->SetRightMargin(0.1);}

  if (iy==0) {pad->SetTopMargin(0.05); pad->SetBottomMargin(epsilon);}
  if (iy==1) {pad->SetTopMargin(epsilon); pad->SetBottomMargin(epsilon);}
  if (iy==2) {pad->SetTopMargin(epsilon); pad->SetBottomMargin(0.2);}
  pad->SetTicks(0,1);

  return pad;
}
