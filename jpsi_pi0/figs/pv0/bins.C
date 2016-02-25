{

  //double tmin[3] = {-0.45, -2.0, -5.0};
  //double tmax[3] = {0.616486, 0.457248, 0.31538};
  //double delta[3] = {0.1, 0.2, 0.5};

  //double tmin[3] = {-0.092, -1.3, -2.85};
  //double ftmin[3] = {-0.1, -1.3, -2.8};
  //double tmax[3] = {0.59, 0.43, 0.3};
  //double delta[3] = {0.08, 0.2, 0.4};
  //for (int ip=0; ip < 3; ++ip) {
  //  cout << "double range[iplab][]= {";
  //  for (int i=0; i<15; ++i) {
  //    double xx = i==0 ? tmin[ip]: (ftmin[ip] + delta[ip]*i);
  //    if (xx>tmax[ip]) {
  //	cout << tmax[ip] << " }; // n=" << i << endl;;
  //	break;
  //    } else {
  //	cout << xx << ",";
  //    }
  //  }
  //}


  double tmin[3] = {-0.092, -1.0, -1.0};
  double tmax[3] = {0.59, 0.43, 0.3};
  double ftmin[3] = {-0.45, -2.0, -5.0};
  double delta[3] = {0.1, 0.2, 0.2};
  //double ftmin[3] = {-0.1, -1.0, -2.0};
  //double delta[3] = {0.08, 0.2, 0.4};
  for (int ip=0; ip < 3; ++ip) {
    cout << "double range[iplab][]= {";
    for (int i=0; i<50; ++i) {
      double xx = ftmin[ip] + delta[ip]*i;
      //double xx = i==0 ? tmin[ip]: (ftmin[ip] + delta[ip]*i);
      if (xx>tmax[ip]) {
	cout << tmax[ip] << " }; // n=" << i << endl;;
	break;
      } else {
	cout << xx << ",";
      }
    }
  }



}
