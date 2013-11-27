

{


  vector<TString> samples;
  vector<TString> systs;

  //samples.push_back("TTH125");
  //samples.push_back("TTJets");
  samples.push_back("TTV");
  samples.push_back("DiBoson");
  samples.push_back("EWK");
  samples.push_back("SingleT");

  systs.push_back("nominal");
  //systs.push_back("csvUp");
  //systs.push_back("csvDown");
  //systs.push_back("JECUp");
  //systs.push_back("JECDown");

  //TString sample = "TTV";
  //TString syst   = "nominal";

  for(int i = 0 ; i < samples.size() ; i++){

    for(int j = 0 ; j < systs.size() ; j++){
      

      float tot1 = 0.;
      float tot2 = 0.;

      TString sample = samples[i];
      TString syst   = systs[j];
      
      cout << string(sample.Data())   << " --- " << string(syst.Data()) << endl;
      TFile* f1 = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/Trees/MEM/Nov21_2013/MEAnalysisNew_DL_"+syst+"_reg_"+sample+".root");
      if( f1!=0 ){
	TTree* tree = (TTree*)f1->Get("tree");
	 tot1 = tree->GetEntries("btag_LR>0.992");
	 //cout << "Old = " << tot << endl;
      }
      TFile* f1 = TFile::Open("./MEAnalysisNew_DL_"+syst+"_v2_reg_"+sample+".root");
      if( f1!=0 ){
	TTree* tree = (TTree*)f1->Get("tree");
	tot2 = tree->GetEntries("btag_LR>0.992");
	//cout << "New = " << tot << endl;
      }

      if( abs(tot1-tot2)>0 ) cout << tot1 << " -- " << tot2 << endl;

    }
  }
  
}
