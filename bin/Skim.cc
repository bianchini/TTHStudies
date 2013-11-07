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
  std::string outPath(    in.getParameter<std::string>("outPath" ) );
  std::string ordering(   in.getParameter<std::string>("ordering" ) );
  std::string newDir(     in.getParameter<std::string>("newDir" ) );
  double lumi(            in.getParameter<double>("lumi" ) );
  bool verbose(           in.getParameter<bool>("verbose" ) );


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
      pset.addParameter("skip",false);
      pset.addParameter("name",    fName);
      pset.addParameter("nickName",fNick);
      pset.addParameter("color",1);
      pset.addParameter("xSec",-1.);
      pset.addParameter("update",false);
      pset.addParameter("cut",cut);
      samples.push_back( pset );
    }
    else{
      cout << "You need to specify name and nickname for your sample" << endl;
      return 0;
    }
  }


  bool openAllFiles = false;
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

      //string cleanSE  = "srmrm srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_"+skim_it+"/"+newDir+"/"+fileName+"_"+skim_it+".root"; 
      //cout << cleanSE << endl;
      //gSystem->Exec(cleanSE.c_str());
     
      //string outputName0 = pathToFile+"/"+fileName+"_"+skim_it+".root"; 
      //TFile* fs0      = new TFile(outputName0.c_str(), "RECREATE");
      //string copyToSE0 = " srmcp -2 file:///"+outputName0+
      //" srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_"+skim_it+"/"+fileName+"_"+skim_it+".root"; 
      //cout << copyToSE0 << endl;
      //gSystem->Exec(copyToSE0.c_str());
      //fs0->Close();
      //delete fs0;

      string outputName = outPath+"/AllHDiJetPt_"+skim_it+"/"+newDir+"/"+fileName+"_"+skim_it+".root"; 
      TFile* fs      = TFile::Open(outputName.c_str(), "RECREATE");
     
      //fwlite::TFileService fs = fwlite::TFileService( ("TThNTuples_"+currentName+".root").c_str());
      //TTree* outTree = fs.make<TTree>("outTree","TTH Tree");
      //TTree* outTree = new TTree("tree","TTH tree");
      
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
      
      //string moveToScratch = "mv ./"+outputName+" "+pathToFile;
      //cout << moveToScratch << endl;
      //gSystem->Exec(moveToScratch.c_str());

      //string cleanSE  = "srmrm srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_"+skim_it+"/"+fileName+"_"+skim_it+".root"; 
      //cout << cleanSE << endl;
      //gSystem->Exec(cleanSE.c_str());
      //string copyToSE = " srmcp -2 file:///"+outputName+
      //" srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_"+skim_it+"/"+fileName+"_"+skim_it+".root"; 
      //cout << copyToSE << endl;
      //gSystem->Exec(copyToSE.c_str());
      
    }


    cout << endl;  

    //string copyToSE = " srmcp -2 file:///"+string(mySamples->GetFile( currentName )->GetName())+
    //" srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_Step2/"+
    //mySamples->GetFileName(currentName);

    //cout << copyToSE << endl;
    //gSystem->Exec(copyToSE.c_str());

    //string removeFromScratch = "rm "+string(mySamples->GetFile( currentName )->GetName());
    //cout << removeFromScratch << endl;
    //gSystem->Exec(removeFromScratch.c_str());
  }

  
  return 0;
  
  
}
