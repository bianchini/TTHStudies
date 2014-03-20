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
    """ Main function work calculation of transfer functions """
    
    extraname =  "TEST"+reg+gen;

    # Open the output file
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


    # RESOLUTION (GAUSSIAN)
    resolLightBin0  = ROOT.TH1F("resolLightBin0","",nBinsY,30,300)
    respLightBin0   = ROOT.TH1F("respLightBin0","", nBinsY,30,300)
    resolLightBin1  = ROOT.TH1F("resolLightBin1","",nBinsY,30,300)
    respLightBin1   = ROOT.TH1F("respLightBin1","", nBinsY,30,300)

    resolHeavyBin0  = ROOT.TH1F("resolHeavyBin0","",nBinsY,30,300)
    respHeavyBin0   = ROOT.TH1F("respHeavyBin0","", nBinsY,30,300)
    resolHeavyBin1  = ROOT.TH1F("resolHeavyBin1","",nBinsY,30,300)
    respHeavyBin1   = ROOT.TH1F("respHeavyBin1","", nBinsY,30,300)
    

    # RESOLUTION (DOUBLE GAUSSIAN)
    respG1HeavyBin0     = ROOT.TH1F("respG1HeavyBin0","",   nBinsY,30,300)
    respG2HeavyBin0     = ROOT.TH1F("respG2HeavyBin0","",   nBinsY,30,300)  
    resolG1HeavyBin0    = ROOT.TH1F("resolG1HeavyBin0","",  nBinsY,30,300)
    resolG2HeavyBin0    = ROOT.TH1F("resolG2HeavyBin0","",  nBinsY,30,300)
    respG1HeavyBin1     = ROOT.TH1F("respG1HeavyBin1","",   nBinsY,30,300)
    respG2HeavyBin1     = ROOT.TH1F("respG2HeavyBin1","",   nBinsY,30,300)  
    resolG1HeavyBin1    = ROOT.TH1F("resolG1HeavyBin1","",  nBinsY,30,300)
    resolG2HeavyBin1    = ROOT.TH1F("resolG2HeavyBin1","",  nBinsY,30,300)

    
    # Define the fitting functions
    # SG = Single Gaussian
    # DG = Double Gaussian
    SG_formula = ("[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])");
    DG_formula_Bin0 = "[0]*( {0}*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-{1})*exp(-0.5*(x-[1]*{2})*(x-[1]*{3})/[3]/[3]))".format( relWeightBin0, 
                                                                                                                                 relWeightBin0, 
                                                                                                                                 relShiftBin0, 
                                                                                                                                 relShiftBin0 )

    DG_formula_Bin1 = "[0]*( {0}*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-{1})*exp(-0.5*(x-[4])*(x-[4])/[3]/[3]))".format(relWeightBin1, 
                                                                                                                        relWeightBin1)

    # Loop over eta bins
    for k in range(2):

        
        cut_Bin_Light = "(pt_rec>20)  *(TMath::Abs(eta{0})>{1} && TMath::Abs(eta{2})<{3})".format( gen, 
                                                                                                   etaBinning[k],  
                                                                                                   gen, 
                                                                                                   etaBinning[k+1] )

        cut_Bin_Heavy = "(pt_rec{0}>20)*(TMath::Abs(eta{1})>{2} && TMath::Abs(eta{3})<{4})".format( reg, 
                                                                                                    gen, 
                                                                                                    etaBinning[k],  
                                                                                                    gen, 
                                                                                                    etaBinning[k+1])


        hCorrLight = ROOT.TH2F("hCorrLight","", nBinsY, 30,300 , 100, 0,500)
        treeJetsLight.Draw("e_rec:e"+gen+">>hCorrLight", cut_Bin_Light )
        
        # Loop over bins in hCorrLight        
        for j in range(1, hCorrLight.GetNbinsX()+1):


            hSlice = hCorrLight.ProjectionY("_py", j, j)
            if (hSlice.GetEntries()<10):
                continue
                
            print hSlice.GetMean(), hSlice.GetRMS()

            resol = ROOT.TF1("resol","gaus", 0,500)
            hSlice.Fit(resol)

            fout.mkdir( "BinLight{0}_{1}".format( j, k))
            fout.cd(    "BinLight{0}_{1}".format( j, k))
            hSlice.Write("",ROOT.TObject.kOverwrite)
     
            if (k==0):
                resolLightBin0.SetBinContent(j, resol.GetParameter(2))
                resolLightBin0.SetBinError  (j, resol.GetParError (2))
                respLightBin0.SetBinContent (j, resol.GetParameter(1))
                respLightBin0.SetBinError   (j, resol.GetParError (1))
            
            if(k==1):
                resolLightBin1.SetBinContent(j, resol.GetParameter(2))
                resolLightBin1.SetBinError  (j, resol.GetParError (2))
                respLightBin1.SetBinContent (j, resol.GetParameter(1))
                respLightBin1.SetBinError   (j, resol.GetParError (1))  
        
        # End of loop over bins in hCorrLight


        
        hCorrHeavy = ROOT.TH2F("hCorrHeavy","",  nBinsY, 30,300 , 100, 0,500)
        treeJetsHeavy.Draw("e_rec"+reg+":e"+gen+">>hCorrHeavy", cut_Bin_Heavy )

        # Loop over bins in hCorrHeavy
        for j in range(1, hCorrHeavy.GetNbinsX()+1):


            hSlice = hCorrHeavy .ProjectionY("_py",j,j)
            if(hSlice.GetEntries()<8):
                continue

            print hSlice.GetMean(), hSlice.GetRMS()
            resol = ROOT.TF1("resol","gaus", 0,500)
            hSlice.Fit(resol)

            fout.mkdir( "BinHeavy{0}_{1}".format(j,k))
            fout.cd(  "BinHeavy{0}_{1}".format(j,k))
            hSlice.Write("", ROOT.TObject.kOverwrite)
  
            
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
                resolHeavyBin0.SetBinContent(j, resol.GetParameter(2))
                resolHeavyBin0.SetBinError  (j, resol.GetParError (2))
                respHeavyBin0.SetBinContent (j, resol.GetParameter(1))
                respHeavyBin0.SetBinError   (j, resol.GetParError (1))
                respG1HeavyBin0.SetBinContent   (j,  resolDG.GetParameter(1) )
                respG1HeavyBin0.SetBinError     (j,  resolDG.GetParError (1) )
                respG2HeavyBin0.SetBinContent   (j,  relShiftBin0 * resolDG.GetParameter(1) )
                respG2HeavyBin0.SetBinError     (j,  relShiftBin0 * resolDG.GetParError (1) )
                #respG2HeavyBin0.SetBinContent   (j,  resolDG.GetParameter(4) )
                #respG2HeavyBin0.SetBinError     (j,  resolDG.GetParError (4) )
                resolG1HeavyBin0.SetBinContent  (j,  resolDG.GetParameter(2) )
                resolG1HeavyBin0.SetBinError    (j,  resolDG.GetParError (2) )
                resolG2HeavyBin0.SetBinContent  (j,  resolDG.GetParameter(3) ) 
                resolG2HeavyBin0.SetBinError    (j,  resolDG.GetParError (3) )

            
            if k==1:
                resolHeavyBin1.SetBinContent(j, resol.GetParameter(2))
                resolHeavyBin1.SetBinError  (j, resol.GetParError (2))
                respHeavyBin1.SetBinContent (j, resol.GetParameter(1))
                respHeavyBin1.SetBinError   (j, resol.GetParError (1))
                respG1HeavyBin1.SetBinContent   (j,  resolDG.GetParameter(1) )
                respG1HeavyBin1.SetBinError     (j,  resolDG.GetParError (1) )
                respG2HeavyBin1.SetBinContent   (j,  resolDG.GetParameter(4) )
                respG2HeavyBin1.SetBinError     (j,  resolDG.GetParError (4) )
                resolG1HeavyBin1.SetBinContent  (j,  resolDG.GetParameter(2) )
                resolG1HeavyBin1.SetBinError    (j,  resolDG.GetParError (2) )
                resolG2HeavyBin1.SetBinContent  (j,  resolDG.GetParameter(3) )
                resolG2HeavyBin1.SetBinError    (j,  resolDG.GetParError (3) )
                
      # Enf of Loop over bins in hCorrHeavy

    # End of loop over eta bins

    sigmaLightBin0 = ROOT.TF1("sigmaLightBin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)
    sigmaLightBin1 = ROOT.TF1("sigmaLightBin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)
    sigmaHeavyBin0 = ROOT.TF1("sigmaHeavyBin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)
    sigmaHeavyBin1 = ROOT.TF1("sigmaHeavyBin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)

    meanLightBin0  = ROOT.TF1("meanLightBin0", "x*[0]",0,400)
    meanLightBin1  = ROOT.TF1("meanLightBin1", "x*[0]",0,400)
    meanHeavyBin0  = ROOT.TF1("meanHeavyBin0", "x*[0]",0,400)
    meanHeavyBin1  = ROOT.TF1("meanHeavyBin1", "x*[0]",0,400)

    meanG1HeavyBin0  = ROOT.TF1("meanG1HeavyBin0", "x*[0]+[1]",0,400)
    meanG2HeavyBin0  = ROOT.TF1("meanG2HeavyBin0", "x*[0]+[1]",0,400)
    sigmaG1HeavyBin0 = ROOT.TF1("sigmaG1HeavyBin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)
    sigmaG2HeavyBin0 = ROOT.TF1("sigmaG2HeavyBin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)

    meanG1HeavyBin1  = ROOT.TF1("meanG1HeavyBin1", "x*[0]+[1]",0,400)
    meanG2HeavyBin1  = ROOT.TF1("meanG2HeavyBin1", "x*[0]+[1]",0,400)
    sigmaG1HeavyBin1 = ROOT.TF1("sigmaG1HeavyBin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)
    sigmaG2HeavyBin1 = ROOT.TF1("sigmaG2HeavyBin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)


    resolLightBin0.Fit(sigmaLightBin0,"","",60,160)
    resolLightBin1.Fit(sigmaLightBin1,"","",60,160)
    resolHeavyBin0.Fit(sigmaHeavyBin0,"","",60,160)
    resolHeavyBin1.Fit(sigmaHeavyBin1,"","",60,160)
    respLightBin0.Fit (meanLightBin0, "","",60,160)
    respLightBin1.Fit (meanLightBin1, "","",60,160)
    respHeavyBin0.Fit (meanHeavyBin0, "","",60,160)
    respHeavyBin1.Fit (meanHeavyBin1, "","",60,160)

    resolLightBin0.Fit(sigmaLightBin0,"","",40,200)
    resolLightBin1.Fit(sigmaLightBin1,"","",40,200)

    resolHeavyBin0.Fit(sigmaHeavyBin0,"","",40,200)
    resolHeavyBin0.Fit(sigmaHeavyBin0,"","",60,300)

    resolHeavyBin1.Fit(sigmaHeavyBin1,"","",40,200)
    resolHeavyBin1.Fit(sigmaHeavyBin1,"","",40,300)

    respLightBin0.Fit (meanLightBin0, "","",40,200)
    respLightBin1.Fit (meanLightBin1, "","",40,200)
    respHeavyBin0.Fit (meanHeavyBin0, "","",40,200)
    respHeavyBin1.Fit (meanHeavyBin1, "","",40,200)

    # responses
    respG1HeavyBin0.Fit  (meanG1HeavyBin0, "","",60,160)
    respG1HeavyBin0.Fit  (meanG1HeavyBin0, "","",40,200)

    respG2HeavyBin0.Fit  (meanG2HeavyBin0, "","",60,160)
    respG2HeavyBin0.Fit  (meanG2HeavyBin0, "","",40,200)

    respG1HeavyBin1.Fit  (meanG1HeavyBin1, "","",80,200)
    respG1HeavyBin1.Fit  (meanG1HeavyBin1, "","",60,300)
    respG1HeavyBin1.Fit  (meanG1HeavyBin1, "","",40,300)

    respG2HeavyBin1.Fit  (meanG2HeavyBin1, "","",80,200)
    respG2HeavyBin1.Fit  (meanG2HeavyBin1, "","",60,300)
    respG2HeavyBin1.Fit  (meanG2HeavyBin1, "","",40,300)
  
    # resolutions
    if "part" in gen:
        resolG1HeavyBin0.Fit (sigmaG1HeavyBin0,"","",100,300)
    else:
        resolG1HeavyBin0.Fit (sigmaG1HeavyBin0,"","",150,300)


    resolG2HeavyBin0.Fit (sigmaG2HeavyBin0,"","",100,300)
 
    if "part" in gen:
        resolG1HeavyBin1.Fit (sigmaG1HeavyBin1,"","",30,  300)
        resolG1HeavyBin1.Fit (sigmaG1HeavyBin1,"","",60,  300)
        resolG1HeavyBin1.Fit (sigmaG1HeavyBin1,"","",80,  300)
    else:
        resolG1HeavyBin1.Fit (sigmaG1HeavyBin1,"","",30,  300)

    resolG2HeavyBin1.Fit (sigmaG2HeavyBin1,"","",130,300)
    resolG2HeavyBin1.Fit (sigmaG2HeavyBin1,"","",130,300)
    #resolG2HeavyBin1.Fit (sigmaG2HeavyBin1,"","",50, 180)
    #resolG2HeavyBin1.Fit (sigmaG2HeavyBin1,"","",50, 250)
    #resolG2HeavyBin1.Fit (sigmaG2HeavyBin1,"","",100,300)
    #resolG2HeavyBin1.Fit (sigmaG2HeavyBin1,"","",120,300)
    #resolG2HeavyBin1.Fit (sigmaG2HeavyBin1,"","",80, 300)
    #resolG2HeavyBin1.Fit (sigmaG2HeavyBin1,"","",40, 250)


    # Loop over eta bins 
    for k in range(2):

        hCorrHeavy = ROOT.TH2F("hCorrHeavy","",  nBinsY, 30,300 , 100, 0,500)

        # loop over bins in hCorrHeavy
        for j in range(1,hCorrHeavy.GetNbinsX()+1):

            en = hCorrHeavy.GetBinCenter(j)

            hCorrHeavy_py = fout.Get( "BinHeavy{0}_{1}/hCorrHeavy_py".format(j, k))

            if not hCorrHeavy_py:
                print "No histo"
                continue

            resolDG_old =  hCorrHeavy_py.GetFunction("resolDG")
            if not resolDG_old:
                print "No func"
                continue

            resolDG = 0
            resolSG = 0
            
            if k==0:
                resolDG = ROOT.TF1("resolDG_fit", DG_formula_Bin0 , -800,800)
                resolDG.SetNpx(1000)
                resolDG.SetLineColor(ROOT.kBlue)

                norm = resolDG_old.Integral(-800,800)
                m1   = meanG1HeavyBin0.Eval(en)
                m2   = meanG2HeavyBin0.Eval(en)
                s1   = sigmaG1HeavyBin0.Eval(en)
                s2   = sigmaG2HeavyBin0.Eval(en)

                resolDG.SetParameter(0, norm/(relWeightBin0*math.sqrt(2*ROOT.TMath.Pi())*s1 + (1-relWeightBin0)*math.sqrt(2*ROOT.TMath.Pi())*s2 )  )
                resolDG.SetParameter(1, m1)
                resolDG.SetParameter(2, s1)
                resolDG.SetParameter(3, s2)
                resolDG.SetParameter(4, m2)

                resolSG = ROOT.TF1("resolSG_fit", SG_formula , -800, 800)
                resolSG.SetNpx(1000)
                resolSG.SetLineColor(ROOT.kMagenta)	
                m = meanHeavyBin0.Eval(en)
                s = sigmaHeavyBin0.Eval(en)
                resolSG.SetParameter(0, norm/math.sqrt(2*ROOT.TMath.Pi())/s )
                resolSG.SetParameter(1, m)
                resolSG.SetParameter(2, s)

                fout.cd( "BinHeavy{0}_{1}".format(j, k) )
                resolDG_old.Write("resolDG{0}_{1}".format(j,k),    ROOT.TObject.kOverwrite)
                resolDG    .Write("resolDG_fit{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)
                resolSG    .Write("resolSG_fit{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)

            if k==1:
                resolDG = ROOT.TF1("resolDG_fit", DG_formula_Bin1 , -800,800)
                resolDG.SetNpx(1000)
                resolDG.SetLineColor(ROOT.kBlue)

                norm = resolDG_old.Integral(-800,800)
                m1 = meanG1HeavyBin1.Eval(en)
                m2 = meanG2HeavyBin1.Eval(en)
                s1 = sigmaG1HeavyBin1.Eval(en)
                s2 = sigmaG2HeavyBin1.Eval(en)
                resolDG.SetParameter(0, norm/(relWeightBin1*math.sqrt(2*ROOT.TMath.Pi())*s1 + (1-relWeightBin1)*math.sqrt(2*ROOT.TMath.Pi())*s2 )  )
                resolDG.SetParameter(1, m1)
                resolDG.SetParameter(2, s1)
                resolDG.SetParameter(3, s2)
                resolDG.SetParameter(4, m2)

                resolSG = ROOT.TF1("resolSG_fit", SG_formula , -800, 800)
                resolSG.SetNpx(1000)
                resolSG.SetLineColor(ROOT.kMagenta)	
                m = meanHeavyBin1.Eval(en)
                s = sigmaHeavyBin1.Eval(en)
                resolSG.SetParameter(0, norm/math.sqrt(2*ROOT.TMath.Pi())/s )
                resolSG.SetParameter(1, m)
                resolSG.SetParameter(2, s)
	
                fout.cd( "BinHeavy{0}_{1}".format(j, k) )
                resolDG_old.Write( "resolDG{0}_{1}".format(    j,k),ROOT.TObject.kOverwrite)
                resolDG    .Write( "resolDG_fit{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)
                resolSG    .Write( "resolSG_fit{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)
                 
        # end of loop over bins in hCorrHeavy

                
        hCorrLight = ROOT.TH2F("hCorrLight","",  nBinsY, 30,300 , 100, 0,500)

        # loop over bins in hCorrLight
        for j in range(1,hCorrLight.GetNbinsX()+1):

            en = hCorrLight.GetBinCenter(j)

            hCorrLight_py = fout.Get("BinLight{0}_{1}/hCorrLight_py".format(j, k))

            if (not hCorrLight_py):
                print "No histo"
                continue


            resolDG_old =  hCorrLight_py.GetFunction("resol")
            if (not resolDG_old):
                print "No func"
                continue

            
            if (k==0):

                resolSG = ROOT.TF1("resolSG_fit", SG_formula , -800, 800)
                resolSG.SetNpx(1000)
                resolSG.SetLineColor(ROOT.kMagenta)
                norm = resolDG_old.Integral(-800,800)	
                m = meanLightBin0.Eval(en)
                s = sigmaLightBin0.Eval(en)
                resolSG.SetParameter(0, norm/math.sqrt(2*ROOT.TMath.Pi())/s )
                resolSG.SetParameter(1, m)
                resolSG.SetParameter(2, s)

                fout.cd( "BinLight{0}_{1}".format(j, k) )
                resolDG_old.Write("resol{0}_{1}".format( j,k),ROOT.TObject.kOverwrite)
                resolSG    .Write("resolSG_fit{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)
                
            if (k==1):
                resolSG = ROOT.TF1("resolSG_fit",SG_formula, -800, 800)
                resolSG.SetNpx(1000)
                resolSG.SetLineColor(ROOT.kMagenta)	
                norm = resolDG_old.Integral(-800,800)
                m = meanLightBin1.Eval(en)
                s = sigmaLightBin1.Eval(en)
                resolSG.SetParameter(0, norm/math.sqrt(2*ROOT.TMath.Pi())/s )
                resolSG.SetParameter(1, m)
                resolSG.SetParameter(2, s)

                fout.cd( "BinLight{0}_{1}".format(j, k) )
                resolDG_old.Write( "resol{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)
                resolSG    .Write( "resolSG_fit{0}_{1}".format(j,k),ROOT.TObject.kOverwrite)


        # loop over bins in hCorrLight
                
    # end of loop over eta bins

  
    if doOnlyJetTF:
        fout.cd()  
        resolLightBin0.Write("",ROOT.TObject.kOverwrite)
        respLightBin0.Write("",ROOT.TObject.kOverwrite)
        resolLightBin1.Write("",ROOT.TObject.kOverwrite)
        respLightBin1.Write("",ROOT.TObject.kOverwrite)
        resolHeavyBin0.Write("",ROOT.TObject.kOverwrite)
        respHeavyBin0.Write("",ROOT.TObject.kOverwrite)
        resolHeavyBin1.Write("",ROOT.TObject.kOverwrite)
        respHeavyBin1.Write("",ROOT.TObject.kOverwrite)
        respG1HeavyBin0.Write("",ROOT.TObject.kOverwrite)
        respG2HeavyBin0.Write("",ROOT.TObject.kOverwrite)
        resolG1HeavyBin0.Write("",ROOT.TObject.kOverwrite)
        resolG2HeavyBin0.Write("",ROOT.TObject.kOverwrite)
        respG1HeavyBin1.Write("",ROOT.TObject.kOverwrite)
        respG2HeavyBin1.Write("",ROOT.TObject.kOverwrite)
        resolG1HeavyBin1.Write("",ROOT.TObject.kOverwrite)
        resolG2HeavyBin1.Write("",ROOT.TObject.kOverwrite)
        fout.Close()
        return
    # end of doOnlyJetTF
# end of calcTF


#########################################
# Allow use from command line and as import
#########################################

# If callled from the command-line: execute calcTF
if __name__ == "__main__":
    calcTF()

# end of handling command line
