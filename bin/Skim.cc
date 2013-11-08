#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "TList.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TObjArray.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TKey.h"
#include "TMultiGraph.h"
#include "Bianchi/TTHStudies/interface/Samples.h"
#include "Bianchi/TTHStudies/interface/Test.h"
#include "Math/GenVector/LorentzVector.h"

#include <algorithm>

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"


using namespace std;


int main(int argc, const char* argv[])
{

  std::cout << "Skim" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  PythonProcessDesc builder(argv[1]);
 
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  std::string pathToFile( in.getParameter<std::string>("pathToFile" ) );
  std::string outPath   ( in.getParameter<std::string>("outPath" ) );
  std::string ordering  ( in.getParameter<std::string>("ordering" ) );
  std::string newDir    ( in.getParameter<std::string>("newDir" ) );
  double lumi           ( in.getParameter<double>     ("lumi" ) );
  bool verbose          ( in.getParameter<bool>       ("verbose" ) );


  std::map<string,string> subsets;
  edm::VParameterSet skims;
  if( in.exists("skims") && (argc==2 || argc==4) ) skims = in.getParameter<edm::VParameterSet>("skims") ;
  else{
    if(argc==6){
      string fName = string((argv[4])); 
      string fCut  = string((argv[5]));
      edm::ParameterSet pset;
      pset.addParameter("name",fName);
      pset.addParameter("cut",fCut );
      skims.push_back( pset );
      if( verbose ){
	TCut mycut(fCut.c_str());
	cout << "Cut ==> " << string(mycut.GetTitle()) << endl;
      }
    }
    else{
      cout << "You need to specify name and cut for your sample" << endl;
      return 0;
    }
  }

  BOOST_FOREACH(edm::ParameterSet p, skims) {
    string name = p.getParameter<string>("name");
    string cut  = p.getParameter<string>("cut");
    subsets[name] = cut;
  }

  //const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;
  edm::VParameterSet samples;
  if( in.exists("samples") && argc==2 ) samples = in.getParameter<edm::VParameterSet>("samples") ;
  else{
    if(argc>=4){
      string fName = string((argv[2])); 
      string fNick = string((argv[3]));
      string cut   = "DUMMY";
      edm::ParameterSet pset;
      pset.addParameter("skip",    false);
      pset.addParameter("name",    fName);
      pset.addParameter("nickName",fNick);
      pset.addParameter("color",       1);
      pset.addParameter("xSec",      -1.);
      pset.addParameter("update",  false);
      pset.addParameter("cut",       cut);
      samples.push_back( pset );
    }
    else{
      cout << "You need to specify name and nickname for your sample" << endl;
      return 0;
    }
  }


  bool openAllFiles  = false;
  Samples* mySamples = new Samples(openAllFiles, pathToFile, ordering, samples, lumi, verbose);
  map<string, TH1F*> mapHist;
  vector<string> mySampleFiles;

  if(mySamples->IsOk()){

    cout << "Ok!" << endl;
    mySampleFiles = mySamples->Files();

    for( unsigned int i = 0 ; i < mySampleFiles.size(); i++){
      string sampleName       = mySampleFiles[i];

      if(verbose){
	cout << mySampleFiles[i] << " ==> " << mySamples->GetXSec(sampleName) 
	     << " pb,"
	     << " ==> weight = "            << mySamples->GetWeight(sampleName) << endl;
      }
    }
  }
  else{
    cout << "Problems... leaving" << endl;
    return 0;
  }


  for(unsigned int i = 0 ; i < mySampleFiles.size(); i++){

    string currentName       = mySampleFiles[i];

    cout << "Opening file " << currentName << endl;
    mySamples->OpenFile( currentName );

    string fileName          = mySamples->GetFileName(currentName);
    
    for(std::map<string, string>::iterator it = subsets.begin(); it!=subsets.end(); it++){
      
      string skim_it = it->first;
      TCut cut((it->second).c_str());
   
      string outputName = outPath+"_"+skim_it+"/"+newDir+"/"+fileName+"_"+skim_it+".root"; 
      TFile* fs      = TFile::Open(outputName.c_str(), "RECREATE");
      
      TIter nextkey((mySamples->GetFile( currentName ))->GetListOfKeys());
      TH1F *key;
      while ( (key = (TH1F*)nextkey()) ) {
	string name(key->GetName());
	if( name.find("tree")!=string::npos) continue;
	TH1F* h = (TH1F*)(mySamples->GetFile( currentName ))->FindObjectAny( key->GetName());
	fs->cd();
	h->Write(key->GetName());
      }
      cout << "Start copying skim " << skim_it << endl;
      cout << " ---> " << outputName << endl;
      cout << "Cut: " << it->second << endl;
      TTree* outTree = (mySamples->GetTree( currentName, "tree"))->CopyTree( cut );
      if(!outTree){
	cout << "Null tree... continue" << endl;
	continue;
      }
      cout << "Done..." << endl;
      
      //////////////////////////////////NEW VARIABLES///////////////////////////////
      //////////////////////////////////////////////////////////////////////////////
      
      
      fs->cd();
      outTree->Write("",TObject::kOverwrite );
      fs->Close();
      
      mySamples->GetFile( currentName )->Close();
            
    }

    cout << endl;  
  }

  
  return 0;
  
  
}
