#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TPRegexp.h"  
#include "TSystem.h"
#include "TROOT.h"


#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

#include "HelperTMVA.h"
#include "HelperFunctions.h"


struct Label
{
    std::string xlabel;
    std::string unit;
    char type;
};

Label MakeLabel(TString str)
{
    if (str.CountChar(';') != 2)
        std::cout << "ERROR: Must contain 2 ';': " << str << std::endl;
    str.ReplaceAll(";;", "; ;");
    TStringToken token(str, ";");
    TString labeltype;
    Label label;

    if(token.NextToken())  label.xlabel = token;
    if(token.NextToken())  label.unit = token;
    if(token.NextToken())  labeltype = token;
    
    if (label.unit == " ")
        label.unit = "";
    if (labeltype != "F" && labeltype != "I")
        std::cout << "ERROR: Must be either 'F' or 'I': " << labeltype << std::endl;
    else
        label.type = (labeltype == "F") ? 'F' : 'I';
    
    return label;
}


////////////////////////////////////////////////////////////////////////////////
/// Main                                                                     ///
////////////////////////////////////////////////////////////////////////////////

void TrainRegression(TString myMethodList="BDTG", TString outfileName="TMVAReg.root" , TString target = "pt_part")
{
    gROOT->SetBatch(1);
    gROOT->LoadMacro("HelperFunctions.h" );  // make functions visible to TTreeFormula

    if (!TString(gROOT->GetVersion()).Contains("5.34")) {
        std::cout << "INCORRECT ROOT VERSION! Please use 5.34:" << std::endl;
        std::cout << "source /uscmst1/prod/sw/cms/slc5_amd64_gcc462/lcg/root/5.34.02-cms/bin/thisroot.csh" << std::endl;
        std::cout << "Return without doing anything." << std::endl;
        return;
    }
    
 
    TMVA::Tools::Instance();


    // Default MVA methods to be trained + tested
    std::map<std::string, int> Use;

    Use["BDTG"]            = 1;
   
    std::cout << std::endl;
    std::cout << "==> Start TMVARegression" << std::endl;

    if (myMethodList != "") {
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

        std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
        for (UInt_t i=0; i<mlist.size(); i++) {
            std::string regMethod(mlist[i]);

            if (Use.find(regMethod) == Use.end()) {
                std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
                for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
                std::cout << std::endl;
                return;
            }
            Use[regMethod] = 1;
        }
    }

    //--------------------------------------------------------------------------
    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    TMVA::Factory *factory = new TMVA::Factory( "TMVARegression_target-"+target, outputFile, 
                                                "!V:!Silent:!Color:!DrawProgressBar:Transformations=I:AnalysisType=Regression" );
    
    const std::vector<std::string> & inputExpressions      = GetInputExpressionsReg();
    const std::vector<std::string> & inputExpressionLabels = GetInputExpressionLabelsReg();
    assert(inputExpressions.size() == inputExpressionLabels.size());

    
    for (UInt_t iexpr=0; iexpr!=inputExpressions.size(); iexpr++){
        Label label = MakeLabel(inputExpressionLabels.at(iexpr));
        TString expr = inputExpressions.at(iexpr);
        factory->AddVariable(expr, label.xlabel, label.unit, label.type);
    }

    factory->AddTarget( target );

    //factory->AddSpectator( "pt_part", "parton pt",  "GeV", 'F');
    //factory->AddSpectator( "pt_gen",  "genjet pt",  "GeV", 'F');
    //factory->AddSpectator( "e_part", "parton e",  "GeV", 'F');
    //factory->AddSpectator( "e_gen",  "genjet e",  "GeV", 'F');
    //factory->AddSpectator( "eta_part", "parton eta",  "", 'F');
    //factory->AddSpectator( "eta_rec",  "parton eta",  "", 'F');

    TFile *input(0);
    TString dirname = "/shome/bianchi/CMSSW_5_3_3_patch2_New/src/Bianchi/TTHStudies/bin/root/";
      //  "gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V2/";
    TString prefix  = "";
    TString suffix  = ".root";

    TTree *regTrainTree(0), *regTestTree(0);
    
    std::vector<std::string> processes;
    //processes.push_back("TTH_HToBB_M-125_8TeV-pythia6");
    processes.push_back("treeProducerTRAIN");

    std::vector<TFile *> files;

    for (UInt_t i=0; i<processes.size(); i++){

      std::string process = processes.at(i);
      input = (TFile*) TFile::Open(dirname + prefix + process + suffix, "READ");
      if (!input) {
	std::cout << "ERROR: Could not open input file." << std::endl;
	exit(1);
      }
      std::cout << "--- TMVARegression           : Using input file: " << input->GetName() << std::endl;
      files.push_back(input);
        
      // --- Register the regression tree
      regTrainTree = (TTree*) input->Get("genJetHeavyTree"); //tree
      regTestTree  = (TTree*) input->Get("genJetHeavyTree"); //tree
      
      // Global event weights per tree (see below for setting event-wise weights)
      Double_t regWeight = 1.0;
      
      // You can add an arbitrary number of regression trees
      factory->AddRegressionTree(regTrainTree, regWeight, TMVA::Types::kTraining);
      factory->AddRegressionTree(regTestTree , regWeight, TMVA::Types::kTesting );
    }


    //factory->SetWeightExpression( "var1", "Regression" );

    //TCut mycut = "hJet_genPt[0] > 0. && hJet_genPt[1] > 0.2 && hJet_csv_nominal[0] > 0.2 && hJet_csv_nominal[1] > 0.2 && hJet_pt[0] > 30. && hJet_pt[1] > 30. && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && abs(hJet_flavour[0])==5 && abs(hJet_flavour[1])==5";
    TCut mycut = "hJet_genPt>20. && csv_rec>0.2 && hJet_pt>30. && TMath::Abs(hJet_eta)<2.5";

    factory->PrepareTrainingAndTestTree( mycut, "V:nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents" );

    if (Use["BDTG"])
        factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                             "!H:V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.7:nCuts=200:MaxDepth=3:NNodesMax=15" );
    //--------------------------------------------------------------------------
    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // --- Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // --- Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    //--------------------------------------------------------------------------
    // Save the output
    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVARegression is done!" << std::endl;

    for (UInt_t i=0; i<files.size(); i++)
        files.at(i)->Close();

    delete outputFile;
    delete factory;

}
