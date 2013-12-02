
{
 
  string cat  = "cat6";
  string type = "DL";

  TH1F* h = new TH1F("h",TString(cat.c_str())+";cut;r",7,0,7);

  vector<string> categories;
  categories.push_back(type+"_"+cat+"_New_v3_reg_0");
  categories.push_back(type+"_"+cat+"_New_v3_reg_1");
  categories.push_back(type+"_"+cat+"_New_v3_reg_2");
  categories.push_back(type+"_"+cat+"_New_v3_reg_3");
  categories.push_back(type+"_"+cat+"_New_v3_reg_4");
  categories.push_back(type+"_"+cat+"_New_v3_reg_5");
  categories.push_back(type+"_"+cat+"_New_v3_reg_6");
  
  for( int b = 0; b < categories.size(); b++){ 
    
    TFile* f = TFile::Open(("higgsCombine"+categories[b]+".Asymptotic.mH120.root").c_str());
    if( f==0 ) continue;
    
    double exp=999;
    
    Double_t r;
    TTree* limit = (TTree*)f->Get("limit");
    limit->SetBranchAddress("limit",&r);
    
    for(int k = 0 ; k< limit->GetEntries() ; k++){
      limit->GetEntry(k);
      if(k==2) exp = r;
    }
    
    cout << categories[b] << ": r<" <<  exp << " @ 95% CL" << endl;

    h->SetBinContent(b+1, exp);
  }
  
  h->GetYaxis()->SetRangeUser(0,15);
  h->Draw();

}
