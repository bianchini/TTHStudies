#!/usr/bin/env python
"""
Calculate jet transfer function for use with the Matrix Element.
Based on TransferFunctions.cc 

IMPORTANT:
  At the moment only "doOnlyJetTF" is reproduced.

TESTS:
  Made sure on March 20 that the produced
  ControlPlotsTEST_std_part.root file agrees with:
    TransferFunctions ../python/transferFunctions_SL.py
  (2ee5731c35e21a5387ffed651169a3bad04ec985)

TODO:
  - also implement CSV TFs

  - when calling from command line argv[2]=reg, argv[3]=gen

  - document caller argument of calcTF

  - throw exceptions instead of returns

  - pythonize - right now it's a straight port

"""

########################################
# Imports
########################################

import sys
import math

import ROOT
ROOT.gROOT.SetBatch(True)    


########################################
# calcTF
########################################

def calcTF(
        outFileName   =  "transferFunctions",
        pathToFile    = "./root/treeProducerTESTpt_gen.root",
        doOnlyJetTF   = 1,
        relWeightBin0 = 0.65,
        relWeightBin1 = 0.65,
        relShiftBin0  = 0.92,
        reg           = "_std",
        gen           = "_part"):
    """ Main function work calculation of transfer functions 

    Arguments:
      outFileName     [string]   Name given to the output file (if doOnlyJetTF=False?)
      pathToFile      [string]   Full path and filename for the input tree root-file
      doOnlyJetTF     [bool]     If False, also do CSV probability
                                   density. Not implemented yet
      relWeightBin0   [float]    Relative weight given to first (of two) Gaussians in 
                                   eta bin0 for b-jet response                          
      relWeightBin1   [float]    Relative weight given to first (of two) Gaussians in 
                                   eta bin0 for b-jet response                          
      relShiftBin0    [float]    Shift of b-jet response in eta bin0 
                                   (understand what exactly this is used for)
      reg             [string]   Use for fit of b-jet energy regression
      gen             [string]   Suffix of the truth tree: should 
                                   be _part (partons) or _gen (genJets)

    Return value: None
    """
    
    extraname =  "TEST"+reg+gen;

    # Open the output file for control plots
    fout = ROOT.TFile("./root/ControlPlots"+extraname+".root","RECREATE");
    fout.cd()

    if "std" in reg:
        reg = ""

        
    # Open and check input file
    input_file = ROOT.TFile.Open( pathToFile,"READ");    
    if ( not input_file or input_file.IsZombie()):
        print "No input file..."
        return

    # Get the trees from the input file
    tree = input_file.Get("genTree")
    if not tree:
        print "The tree called 'genTree' is not there..."
        return

    treeEvent = input_file.Get("genEventTree")
    if not treeEvent:
        print "The tree called 'genEventTree' is not there..."
        return

    treeJetsLight = input_file.Get("genJetLightTree")
    if not treeJetsLight:
        print "The tree called 'genJetLightTree' is not there..."
        return
    
    treeJetsHeavy = input_file.Get("genJetHeavyTree")
    if not treeJetsHeavy:
        print "The tree called 'genJetHeavyTree' is not there..."
        return 

    
    # Define eta-binning and number of pt-bins
    etaBinning = [0.0, 1.0, 2.5]
    nBinsY     = 58 

    for jetType in ["Light", "Heavy"]:

        # RESOLUTION (GAUSSIAN)
        resolBin0  = ROOT.TH1F("resol"+jetType+"Bin0","",nBinsY,30,300)
        respBin0   = ROOT.TH1F("resp" +jetType+"Bin0","",nBinsY,30,300)
        resolBin1  = ROOT.TH1F("resol"+jetType+"Bin1","",nBinsY,30,300)
        respBin1   = ROOT.TH1F("resp" +jetType+"Bin1","",nBinsY,30,300)

        # RESOLUTION (DOUBLE GAUSSIAN)
        respG1Bin0     = ROOT.TH1F("respG1" +jetType+"0","",  nBinsY,30,300)
        respG2Bin0     = ROOT.TH1F("respG2" +jetType+"0","",  nBinsY,30,300)  
        resolG1Bin0    = ROOT.TH1F("resolG1"+jetType+"0","",  nBinsY,30,300)
        resolG2Bin0    = ROOT.TH1F("resolG2"+jetType+"0","",  nBinsY,30,300)
        respG1Bin1     = ROOT.TH1F("respG1" +jetType+"1","",  nBinsY,30,300)
        respG2Bin1     = ROOT.TH1F("respG2" +jetType+"1","",  nBinsY,30,300)  
        resolG1Bin1    = ROOT.TH1F("resolG1"+jetType+"1","",  nBinsY,30,300)
        resolG2Bin1    = ROOT.TH1F("resolG2"+jetType+"1","",  nBinsY,30,300)
         
        # Define the fitting functions
        # SG = Single Gaussian
        # DG = Double Gaussian        
        SG_formula_Bin0 = ("[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])");
        SG_formula_Bin1 = ("[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])");

        DG_formula_Bin0 = "[0]*( {0}*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-{1})*exp(-0.5*(x-[1]*{2})*(x-[1]*{3})/[3]/[3]))".format( relWeightBin0, 
                                                                                                                                 relWeightBin0, 
                                                                                                                                 relShiftBin0, 
                                                                                                                                 relShiftBin0 )

        DG_formula_Bin1 = "[0]*( {0}*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-{1})*exp(-0.5*(x-[4])*(x-[4])/[3]/[3]))".format(relWeightBin1, 
                                                                                                                                relWeightBin1)

        # Loop over eta bins
        for k in range(2):

            
            if jetType == "Light":
                cut_Bin = "(pt_rec>20)  *(TMath::Abs(eta{0})>{1} && TMath::Abs(eta{2})<{3})".format( gen, 
                                                                                                     etaBinning[k],  
                                                                                                     gen, 
                                                                                                     etaBinning[k+1] )
            elif jetType == "Heavy":
                cut_Bin = "(pt_rec{0}>20)*(TMath::Abs(eta{1})>{2} && TMath::Abs(eta{3})<{4})".format( reg, 
                                                                                                      gen, 
                                                                                                      etaBinning[k],  
                                                                                                      gen, 
                                                                                                      etaBinning[k+1])
                
            
            hCorr = ROOT.TH2F("hCorr"+jetType,"", nBinsY, 30, 300, 100, 0, 500)
            if jetType == "Light":
                treeJetsLight.Draw("e_rec:e"+gen+">>hCorr"+jetType, cut_Bin )
            elif jetType == "Heavy":
                treeJetsHeavy.Draw("e_rec"+reg+":e"+gen+">>hCorr"+jetType, cut_Bin )


            # Loop over bins in hCorr        
            for j in range(1, hCorr.GetNbinsX()+1):


                hSlice = hCorr.ProjectionY("_py{0}_{1}".format(j,k), j, j)

                print hSlice.GetName()

                
                if jetType == "Light":
                    if (hSlice.GetEntries()<10):
                        continue
                elif jetType == "Heavy":
                    if (hSlice.GetEntries()<8):
                        continue

                print hSlice.GetMean(), hSlice.GetRMS()

                resol = ROOT.TF1("resol","gaus", 0,500)

                
                hSlice.Fit(resol)

                fout.mkdir( "Bin{0}{1}_{2}".format( jetType, j, k))
                fout.cd(    "Bin{0}{1}_{2}".format( jetType, j, k))
                hSlice.Write("",ROOT.TObject.kOverwrite)


                
                if jetType == "Light":
                    if (k==0):
                        resolBin0.SetBinContent(j, resol.GetParameter(2))
                        resolBin0.SetBinError  (j, resol.GetParError (2))
                        respBin0.SetBinContent (j, resol.GetParameter(1))
                        respBin0.SetBinError   (j, resol.GetParError (1))

                    if(k==1):
                        resolBin1.SetBinContent(j, resol.GetParameter(2))
                        resolBin1.SetBinError  (j, resol.GetParError (2))
                        respBin1.SetBinContent (j, resol.GetParameter(1))
                        respBin1.SetBinError   (j, resol.GetParError (1))  

                
                elif jetType == "Heavy":
                    if k==0:	
                        resolDG = ROOT.TF1("resolDG", DG_formula_Bin0, 0,500)
                        resolDG.SetParameter(1, hSlice.GetMean())
                        resolDG.SetParameter(2, hSlice.GetRMS())
                        resolDG.SetParameter(3, hSlice.GetRMS())

                    if k==1:
                        resolDG = ROOT.TF1("resolDG", DG_formula_Bin1, 0,500)
                        resolDG.SetParameter(1, hSlice.GetMean())
                        resolDG.SetParameter(2, hSlice.GetRMS())
                        resolDG.SetParameter(3, hSlice.GetRMS())
                        resolDG.SetParLimits(4, hSlice.GetMean()-2*hSlice.GetRMS() ,hSlice.GetMean())

                    hSlice.Fit(resolDG)
                    hSlice.Fit(resolDG)

                    hSlice.Write("",ROOT.TObject.kOverwrite)

                    if k==0:
                        resolBin0.SetBinContent(j, resol.GetParameter(2))
                        resolBin0.SetBinError  (j, resol.GetParError (2))
                        respBin0.SetBinContent (j, resol.GetParameter(1))
                        respBin0.SetBinError   (j, resol.GetParError (1))
                        respG1Bin0.SetBinContent   (j,  resolDG.GetParameter(1) )
                        respG1Bin0.SetBinError     (j,  resolDG.GetParError (1) )
                        respG2Bin0.SetBinContent   (j,  relShiftBin0 * resolDG.GetParameter(1) )
                        respG2Bin0.SetBinError     (j,  relShiftBin0 * resolDG.GetParError (1) )
                        #respG2Bin0.SetBinContent   (j,  resolDG.GetParameter(4) )
                        #respG2Bin0.SetBinError     (j,  resolDG.GetParError (4) )
                        resolG1Bin0.SetBinContent  (j,  resolDG.GetParameter(2) )
                        resolG1Bin0.SetBinError    (j,  resolDG.GetParError (2) )
                        resolG2Bin0.SetBinContent  (j,  resolDG.GetParameter(3) ) 
                        resolG2Bin0.SetBinError    (j,  resolDG.GetParError (3) )

                    if k==1:
                        resolBin1.SetBinContent(j, resol.GetParameter(2))
                        resolBin1.SetBinError  (j, resol.GetParError (2))
                        respBin1.SetBinContent (j, resol.GetParameter(1))
                        respBin1.SetBinError   (j, resol.GetParError (1))
                        respG1Bin1.SetBinContent   (j,  resolDG.GetParameter(1) )
                        respG1Bin1.SetBinError     (j,  resolDG.GetParError (1) )
                        respG2Bin1.SetBinContent   (j,  resolDG.GetParameter(4) )
                        respG2Bin1.SetBinError     (j,  resolDG.GetParError (4) )
                        resolG1Bin1.SetBinContent  (j,  resolDG.GetParameter(2) )
                        resolG1Bin1.SetBinError    (j,  resolDG.GetParError (2) )
                        resolG2Bin1.SetBinContent  (j,  resolDG.GetParameter(3) )
                        resolG2Bin1.SetBinError    (j,  resolDG.GetParError (3) )

                # End treatment of Heavy Jets                
            # End of loop over bins in hCorr
        # End of loop over eta bins



        sigmaBin0 = ROOT.TF1("sigma"+jetType+"Bin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)
        sigmaBin1 = ROOT.TF1("sigma"+jetType+"Bin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)

        meanBin0  = ROOT.TF1("mean"+jetType+"Bin0", "x*[0]",0,400)
        meanBin1  = ROOT.TF1("mean"+jetType+"Bin1", "x*[0]",0,400)

        if jetType == "Heavy":
            meanG1Bin0  = ROOT.TF1("meanG1"+jetType+"Bin0", "x*[0]+[1]",0,400)
            meanG2Bin0  = ROOT.TF1("meanG2"+jetType+"Bin0", "x*[0]+[1]",0,400)
            sigmaG1Bin0 = ROOT.TF1("sigmaG1"+jetType+"Bin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)
            sigmaG2Bin0 = ROOT.TF1("sigmaG2"+jetType+"Bin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)

            meanG1Bin1  = ROOT.TF1("meanG1"+jetType+"Bin1", "x*[0]+[1]",0,400)
            meanG2Bin1  = ROOT.TF1("meanG2"+jetType+"Bin1", "x*[0]+[1]",0,400)
            sigmaG1Bin1 = ROOT.TF1("sigmaG1"+jetType+"Bin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)
            sigmaG2Bin1 = ROOT.TF1("sigmaG2"+jetType+"Bin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)
        # End of jetType == "Heavy"

        # Single-Gauss Resolutions
        resolBin0.Fit(sigmaBin0,"","",60,160)
        resolBin0.Fit(sigmaBin0,"","",40,200)
        # extra fit for heavy jets
        if jetType == "Heavy":
            resolBin0.Fit(sigmaBin0,"","",60,300)

        resolBin1.Fit(sigmaBin1,"","",60,160)
        resolBin1.Fit(sigmaBin1,"","",40,200)
        # extra fit for heavy jets        
        if jetType == "Heavy":
            resolBin1.Fit(sigmaBin1,"","",40,300)
                    
        # Single-Gauss Responses
        respBin0.Fit (meanBin0, "","",60,160)
        respBin0.Fit (meanBin0, "","",40,200)
        respBin1.Fit (meanBin1, "","",60,160)
        respBin1.Fit (meanBin1, "","",40,200)

        if jetType == "Heavy":
            # Double Gauss Responses
            respG1Bin0.Fit  (meanG1Bin0, "","",60,160)
            respG1Bin0.Fit  (meanG1Bin0, "","",40,200)

            respG2Bin0.Fit  (meanG2Bin0, "","",60,160)
            respG2Bin0.Fit  (meanG2Bin0, "","",40,200)

            respG1Bin1.Fit  (meanG1Bin1, "","",80,200)
            respG1Bin1.Fit  (meanG1Bin1, "","",60,300)
            respG1Bin1.Fit  (meanG1Bin1, "","",40,300)

            respG2Bin1.Fit  (meanG2Bin1, "","",80,200)
            respG2Bin1.Fit  (meanG2Bin1, "","",60,300)
            respG2Bin1.Fit  (meanG2Bin1, "","",40,300)

            # Double Gauss Resolutions
            if "part" in gen:
                resolG1Bin0.Fit (sigmaG1Bin0,"","",100,300)
            else:
                resolG1Bin0.Fit (sigmaG1Bin0,"","",150,300)

            resolG2Bin0.Fit (sigmaG2Bin0,"","",100,300)

            if "part" in gen:
                resolG1Bin1.Fit (sigmaG1Bin1,"","",30,  300)
                resolG1Bin1.Fit (sigmaG1Bin1,"","",60,  300)
                resolG1Bin1.Fit (sigmaG1Bin1,"","",80,  300)
            else:
                resolG1Bin1.Fit (sigmaG1Bin1,"","",30,  300)

            resolG2Bin1.Fit (sigmaG2Bin1,"","",130,300)
            resolG2Bin1.Fit (sigmaG2Bin1,"","",130,300)
        # End of jetType == "Heavy"


        

        # Loop over eta bins 
        for k in range(2):

            hCorr = ROOT.TH2F("hCorr"+jetType,"",  nBinsY, 30,300 , 100, 0,500)

            # loop over bins in hCorr
            for j in range(1,hCorr.GetNbinsX()+1):

                fout.cd()  
                en = hCorr.GetBinCenter(j)

                hCorr_py = fout.Get("Bin{0}{1}_{2}/_py{1}_{2}".format(jetType, j, k))

                if not hCorr_py:
                    print "No histo", jetType, j, k
                    sys.exit()
                    continue
                
                if jetType == "Light":                
                    resolDG_old =  hCorr_py.GetFunction("resol")
                elif jetType == "Heavy":                
                    resolDG_old =  hCorr_py.GetFunction("resolDG")
                if (not resolDG_old):
                    print "No func"
                    continue
                    
                resolDG = 0
                resolSG = 0

                if k==0:
                    
                    if jetType == "Light":
                        resolSG = ROOT.TF1("resolSG_fit", SG_formula_Bin0 , -800, 800)
                        resolSG.SetNpx(1000)
                        resolSG.SetLineColor(ROOT.kMagenta)
                        norm = resolDG_old.Integral(-800,800)	
                        m = meanBin0.Eval(en)
                        s = sigmaBin0.Eval(en)
                        resolSG.SetParameter(0, norm/math.sqrt(2*ROOT.TMath.Pi())/s )
                        resolSG.SetParameter(1, m)
                        resolSG.SetParameter(2, s)

                        fout.cd( "BinLight{0}_{1}".format(j, k) )
                        resolDG_old.Write("resol{0}_{1}".format( j,k),ROOT.TObject.kOverwrite)
                        resolSG    .Write("resolSG_fit{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)

                    elif jetType == "Heavy":

                        resolDG = ROOT.TF1("resolDG_fit", DG_formula_Bin0 , -800,800)
                        resolDG.SetNpx(1000)
                        resolDG.SetLineColor(ROOT.kBlue)

                        norm = resolDG_old.Integral(-800,800)
                        m1   = meanG1Bin0.Eval(en)
                        m2   = meanG2Bin0.Eval(en)
                        s1   = sigmaG1Bin0.Eval(en)
                        s2   = sigmaG2Bin0.Eval(en)

                        resolDG.SetParameter(0, norm/(relWeightBin0*math.sqrt(2*ROOT.TMath.Pi())*s1 + (1-relWeightBin0)*math.sqrt(2*ROOT.TMath.Pi())*s2 )  )
                        resolDG.SetParameter(1, m1)
                        resolDG.SetParameter(2, s1)
                        resolDG.SetParameter(3, s2)
                        resolDG.SetParameter(4, m2)

                        resolSG = ROOT.TF1("resolSG_fit", SG_formula_Bin0 , -800, 800)
                        resolSG.SetNpx(1000)
                        resolSG.SetLineColor(ROOT.kMagenta)	
                        m = meanBin0.Eval(en)
                        s = sigmaBin0.Eval(en)
                        resolSG.SetParameter(0, norm/math.sqrt(2*ROOT.TMath.Pi())/s )
                        resolSG.SetParameter(1, m)
                        resolSG.SetParameter(2, s)

                        fout.cd( "BinHeavy{0}_{1}".format(j, k) )
                        resolDG_old.Write("resolDG{0}_{1}".format(j,k),    ROOT.TObject.kOverwrite)
                        resolDG    .Write("resolDG_fit{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)
                        resolSG    .Write("resolSG_fit{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)

                
                if k==1:
                    if jetType == "Light":
                        resolSG = ROOT.TF1("resolSG_fit",SG_formula_Bin1, -800, 800)
                        resolSG.SetNpx(1000)
                        resolSG.SetLineColor(ROOT.kMagenta)	
                        norm = resolDG_old.Integral(-800,800)
                        m = meanBin1.Eval(en)
                        s = sigmaBin1.Eval(en)
                        resolSG.SetParameter(0, norm/math.sqrt(2*ROOT.TMath.Pi())/s )
                        resolSG.SetParameter(1, m)
                        resolSG.SetParameter(2, s)

                        fout.cd( "BinLight{0}_{1}".format(j, k) )
                        resolDG_old.Write( "resol{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)
                        resolSG    .Write( "resolSG_fit{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)

                    elif jetType == "Heavy":
                        
                        resolDG = ROOT.TF1("resolDG_fit", DG_formula_Bin1 , -800,800)
                        resolDG.SetNpx(1000)
                        resolDG.SetLineColor(ROOT.kBlue)

                        norm = resolDG_old.Integral(-800,800)
                        m1 = meanG1Bin1.Eval(en)
                        m2 = meanG2Bin1.Eval(en)
                        s1 = sigmaG1Bin1.Eval(en)
                        s2 = sigmaG2Bin1.Eval(en)
                        resolDG.SetParameter(0, norm/(relWeightBin1*math.sqrt(2*ROOT.TMath.Pi())*s1 + (1-relWeightBin1)*math.sqrt(2*ROOT.TMath.Pi())*s2 )  )
                        resolDG.SetParameter(1, m1)
                        resolDG.SetParameter(2, s1)
                        resolDG.SetParameter(3, s2)
                        resolDG.SetParameter(4, m2)

                        resolSG = ROOT.TF1("resolSG_fit", SG_formula_Bin1 , -800, 800)
                        resolSG.SetNpx(1000)
                        resolSG.SetLineColor(ROOT.kMagenta)	
                        m = meanBin1.Eval(en)
                        s = sigmaBin1.Eval(en)
                        resolSG.SetParameter(0, norm/math.sqrt(2*ROOT.TMath.Pi())/s )
                        resolSG.SetParameter(1, m)
                        resolSG.SetParameter(2, s)

                        fout.cd( "BinHeavy{0}_{1}".format(j, k) )
                        resolDG_old.Write( "resolDG{0}_{1}".format(    j,k),ROOT.TObject.kOverwrite)
                        resolDG    .Write( "resolDG_fit{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)
                        resolSG    .Write( "resolSG_fit{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)
                    # End of jetType == Heavy
                # End of eta_bin == 1
            # end of loop over bins in hCorr
        # end of loop over eta bins

        if doOnlyJetTF:
            fout.cd()  
            resolBin0.Write("",ROOT.TObject.kOverwrite)
            respBin0. Write("",ROOT.TObject.kOverwrite)
            resolBin1.Write("",ROOT.TObject.kOverwrite)
            respBin1. Write("",ROOT.TObject.kOverwrite)

            if jetType == "Heavy":
                respG1Bin0. Write("",ROOT.TObject.kOverwrite)
                respG2Bin0. Write("",ROOT.TObject.kOverwrite)
                resolG1Bin0.Write("",ROOT.TObject.kOverwrite)
                resolG2Bin0.Write("",ROOT.TObject.kOverwrite)
                respG1Bin1. Write("",ROOT.TObject.kOverwrite)
                respG2Bin1. Write("",ROOT.TObject.kOverwrite)
                resolG1Bin1.Write("",ROOT.TObject.kOverwrite)
                resolG2Bin1.Write("",ROOT.TObject.kOverwrite)
            # End of jetType == Heavy
        # End of doOnlyJetTF
    # End of loop over jetType

    hCorr_py = fout.Get("Bin{0}{1}_{2}/_py{1}_{2}".format("Light", 1, 0))
    resolDG_old =  hCorr_py.GetFunction("resol")                    


  
    if doOnlyJetTF:
        fout.Close()
        return
# end of calcTF


#########################################
# Allow use from command line and as import
#########################################

# If callled from the command-line: execute calcTF
if __name__ == "__main__":
    calcTF()

# end of handling command line
