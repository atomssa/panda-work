{

  double start[3] = {-0.45, -2.0, -5.0};
  double delta[3] = {0.1, 0.2, 0.5};
  double tmax[3] = {0.616486, 0.457248, 0.31538};

  for (int ip=0; ip < 3; ++ip) {
    cout << "double range[iplab][]= {";
    for (int i=0; i<15; ++i) {
      double xx = start[ip] + delta[ip]*i;
      if (xx>tmax[ip]) {
	cout << tmax[ip] << " }; // n=" << i << endl;;
	break;
      } else {
	cout << xx << ",";
      }
    }
  }

}
