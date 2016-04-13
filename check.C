{
  int IntMode;
  double XSecWeight, FluxWeight;
  neut_nominal->SetBranchAddress("IntMode",    &IntMode);
  neut_nominal->SetBranchAddress("XSecWeight", &XSecWeight);
  neut_nominal->SetBranchAddress("FluxWeight", &FluxWeight);

  for(int i = 0;  i < neut_nominal->GetEntries(); i++){
    neut_nominal->GetEntry(i);
    if(XSecWeight > 2 || XSecWeight < 0 // || FluxWeight == 0
       )
      neut_nominal->Show(i);
    
  }
}
