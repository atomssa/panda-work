void tmp3() {

  vector<int> v;

  for (int i=0;i<5; ++i) v.push_back(i);
  cout << "v = ";
  for (int i=0;i<5; ++i) cout << v[i] << " ";
  cout << endl;
  
  //int *x = &v[2];
  //*x = 5;

  int &x = v[2];
  x = 5;
    
  cout << "v2 = ";
  for (int i=0;i<5; ++i) cout << v[i] << " ";
  cout << endl;
  
  
}
