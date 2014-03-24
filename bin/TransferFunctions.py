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
        gen           = "_part",
        li_jet_types  = ["Light", "Heavy"],                
):
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
      li_jet_types    [list of strings] Jet types for which to calculate the TF.
                                          Available: Light, Heavy, BDRSSub

    Return value: None
    """

    ########################################
    # Validate arguments
    ########################################

    li_available_jet_types = ["Light", "Heavy", "BDRSSub"]

    for jt in li_jet_types:
        if jt not in li_available_jet_types:
            print "Not able to process: ", jt
            print "Available are: ", " ".join(li_jet_types)
            print "Exiting..."
            return


    ########################################
    # Prepare Output
    ########################################
    
    extraname =  "TEST"+reg+gen;

    # Open the output file for control plots
    fout = ROOT.TFile("./root/ControlPlots"+extraname+".root","RECREATE");
    fout.cd()

    if "std" in reg:
        reg = ""


    ########################################
    # Get Input Trees
    ########################################

    # Open and check input file
    input_file = ROOT.TFile.Open( pathToFile,"READ");    
    if ( not input_file or input_file.IsZombie()):
        print "No input file..."
        return

    if "Light" in li_jet_types:
        treeJetsLight = input_file.Get("genJetLightTree")
        if not treeJetsLight:
            print "Could not get treee for Light jets"
            print "Exiting..."

    if "Heavy" in li_jet_types:
        treeJetsHeavy = input_file.Get("genJetHeavyTree")
        if not treeJetsHeavy:
            print "Could not get treee for Heavy jets"
            print "Exiting..."

    if "BDRSSub" in li_jet_types:
        treeJetsBDRSSub = input_file.Get("pfBDRS15xar03")
        if not treeJetsBDRSSub:
            print "Could not get treee for BDRSSub jets"
            print "Exiting..."
                
    
    # Define eta-binning and number of pt-bins
    etaBinning = [0.0, 1.0, 2.5]
    nBinsY     = 58 


    ########################################
    # Loop over jet types
    ########################################

    for jetType in li_jet_types:
        
        ########################################
        # Loop eta bins
        ########################################

        for i_eta in range( len(etaBinning)-1 ):

            # RESOLUTION (GAUSSIAN)
            resolBin  = ROOT.TH1F("resol"+jetType+"Bin"+str(i_eta),"", nBinsY,30,300)
            respBin   = ROOT.TH1F("resp" +jetType+"Bin"+str(i_eta),"", nBinsY,30,300)

            # RESOLUTION (DOUBLE GAUSSIAN)
            respG1Bin     = ROOT.TH1F("respG1" +jetType+str(i_eta),"", nBinsY,30,300)
            respG2Bin     = ROOT.TH1F("respG2" +jetType+str(i_eta),"", nBinsY,30,300)  
            resolG1Bin    = ROOT.TH1F("resolG1"+jetType+str(i_eta),"", nBinsY,30,300)
            resolG2Bin    = ROOT.TH1F("resolG2"+jetType+str(i_eta),"", nBinsY,30,300)

            # Define the fitting functions
            # SG = Single Gaussian
            # DG = Double Gaussian        
            SG_formula_Bin = ("[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])")

            if i_eta==0:
                DG_formula_Bin = "[0]*( {0}*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-{1})*exp(-0.5*(x-[1]*{2})*(x-[1]*{3})/[3]/[3]))".format( relWeightBin0, 
                                                                                                                                            relWeightBin0, 
                                                                                                                                            relShiftBin0, 
                                                                                                                                            relShiftBin0 )
            elif i_eta==1:
                DG_formula_Bin = "[0]*( {0}*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-{1})*exp(-0.5*(x-[4])*(x-[4])/[3]/[3]))".format(relWeightBin1, 
                                                                                                                                   relWeightBin1)


            ########################################
            # Prepare 2d histogram hCorr
            ########################################

            hCorr = ROOT.TH2F("hCorr"+jetType,"", nBinsY, 30, 300, 100, 0, 500)

            # Light
            if jetType == "Light":
                minEntriesPerSlice = 10
                cut_Bin = "(pt_rec>20)  *(TMath::Abs(eta{0})>{1} && TMath::Abs(eta{2})<{3})".format( gen, 
                                                                                                     etaBinning[i_eta],  
                                                                                                     gen, 
                                                                                                     etaBinning[i_eta+1] )
                treeJetsLight.Draw("e_rec:e"+gen+">>hCorr"+jetType, cut_Bin )


            # Heavy
            elif jetType == "Heavy":
                minEntriesPerSlice = 8
                cut_Bin = "(pt_rec{0}>20)*(TMath::Abs(eta{1})>{2} && TMath::Abs(eta{3})<{4})".format( reg, 
                                                                                                      gen, 
                                                                                                      etaBinning[i_eta],  
                                                                                                      gen, 
                                                                                                      etaBinning[i_eta+1])
                treeJetsHeavy.Draw("e_rec"+reg+":e"+gen+">>hCorr"+jetType, cut_Bin )


            # BDRSSub
            elif jetType == "BDRSSub":                

                minEntriesPerSlice = 10
                
                # Loop over the two subjets and draw them both into
                # the same 2D histogram
                first_sj = True                
                for sj in ["sj1", "sj2"]:
                    cut_Bin  = "((mass>0)" 
                    cut_Bin +=   "&& ({0}_drtrue<0.3)"
                    cut_Bin +=   "&& ({0}_pt>20)"
                    cut_Bin +=   "&& (TMath::Abs( {0}_trueeta ) > {1})"
                    cut_Bin +=   "&& (TMath::Abs( {0}_trueeta ) < {2}))"
                    cut_Bin = cut_Bin.format( sj, etaBinning[i_eta], etaBinning[i_eta+1] )

                    # Just draw the first
                    if first_sj:
                        target_string = ">>"
                        first_sj = False
                    # Add all later subjets
                    else:
                        target_string = ">>+"

                    draw_string = "{0}_pt:{0}_truept{1}hCorr{2}".format( sj, target_string, jetType)
                    treeJetsBDRSSub.Draw( draw_string, cut_Bin)                

                # end of loop over subjets
            # End of treating BDRSSub jets



            # Loop over bins in hCorr        
            for j in range(1, hCorr.GetNbinsX()+1):


                hSlice = hCorr.ProjectionY("_py{0}_{1}".format(j,i_eta), j, j)
                
                if (hSlice.GetEntries()<minEntriesPerSlice):
                    continue

                print hSlice.GetMean(), hSlice.GetRMS()

                resol = ROOT.TF1("resol","gaus", 0,500)
                
                hSlice.Fit(resol)

                fout.mkdir( "Bin{0}{1}_{2}".format( jetType, j, i_eta))
                fout.cd(    "Bin{0}{1}_{2}".format( jetType, j, i_eta))
                hSlice.Write("",ROOT.TObject.kOverwrite)
                
                if jetType == "Light":
                    resolBin.SetBinContent(j, resol.GetParameter(2))
                    resolBin.SetBinError  (j, resol.GetParError (2))
                    respBin.SetBinContent (j, resol.GetParameter(1))
                    respBin.SetBinError   (j, resol.GetParError (1))
                
                elif jetType == "Heavy":

                    resolDG = ROOT.TF1("resolDG", DG_formula_Bin, 0,500)
                    resolDG.SetParameter(1, hSlice.GetMean())
                    resolDG.SetParameter(2, hSlice.GetRMS())
                    resolDG.SetParameter(3, hSlice.GetRMS())
                    if i_eta==1:
                        resolDG.SetParLimits(4, hSlice.GetMean()-2*hSlice.GetRMS() ,hSlice.GetMean())

                    hSlice.Fit(resolDG)
                    hSlice.Fit(resolDG)

                    hSlice.Write("",ROOT.TObject.kOverwrite)

                    resolBin.SetBinContent(j, resol.GetParameter(2))
                    resolBin.SetBinError  (j, resol.GetParError (2))
                    respBin.SetBinContent (j, resol.GetParameter(1))
                    respBin.SetBinError   (j, resol.GetParError (1))

                    respG1Bin.SetBinContent   (j,  resolDG.GetParameter(1) )
                    respG1Bin.SetBinError     (j,  resolDG.GetParError (1) )

                    resolG1Bin.SetBinContent  (j,  resolDG.GetParameter(2) )
                    resolG1Bin.SetBinError    (j,  resolDG.GetParError (2) )
                    resolG2Bin.SetBinContent  (j,  resolDG.GetParameter(3) ) 
                    resolG2Bin.SetBinError    (j,  resolDG.GetParError (3) )

                    if i_eta==0:
                        respG2Bin.SetBinContent   (j,  relShiftBin0 * resolDG.GetParameter(1) )
                        respG2Bin.SetBinError     (j,  relShiftBin0 * resolDG.GetParError (1) )
                    elif i_eta==1:
                        respG2Bin.SetBinContent   (j,  resolDG.GetParameter(4) )
                        respG2Bin.SetBinError     (j,  resolDG.GetParError (4) )


                elif jetType == "BDRSSub":
                    resolBin.SetBinContent(j, resol.GetParameter(2))
                    resolBin.SetBinError  (j, resol.GetParError (2))
                    respBin.SetBinContent (j, resol.GetParameter(1))
                    respBin.SetBinError   (j, resol.GetParError (1))


                # End treatment of BDRSSub jets
            # End of loop over bins in hCorr

            sigmaBin = ROOT.TF1("sigma"+jetType+"Bin"+str(i_eta),"sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)
            meanBin  = ROOT.TF1("mean"+jetType+"Bin"+str(i_eta), "x*[0]",0,400)

            if jetType == "Heavy":
                meanG1Bin  = ROOT.TF1("meanG1"+jetType+"Bin"+str(i_eta), "x*[0]+[1]",0,400)
                meanG2Bin  = ROOT.TF1("meanG2"+jetType+"Bin"+str(i_eta), "x*[0]+[1]",0,400)
                sigmaG1Bin = ROOT.TF1("sigmaG1"+jetType+"Bin"+str(i_eta),"sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)
                sigmaG2Bin = ROOT.TF1("sigmaG2"+jetType+"Bin"+str(i_eta),"sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400)
            # End of jetType == "Heavy"

            # Single-Gauss Resolutions
            resolBin.Fit(sigmaBin,"","",60,160)
            resolBin.Fit(sigmaBin,"","",40,200)

            # extra fit for heavy jets
            if jetType == "Heavy":
                if i_eta==0:
                    resolBin.Fit(sigmaBin,"","",60,300)
                elif i_eta==1:
                    resolBin.Fit(sigmaBin,"","",40,300)

            # Single-Gauss Responses
            respBin.Fit (meanBin, "","",60,160)
            respBin.Fit (meanBin, "","",40,200)

            if jetType == "Heavy":
                # Double Gauss Responses
                if i_eta==0:
                    respG1Bin.Fit  (meanG1Bin, "","",60,160)
                    respG1Bin.Fit  (meanG1Bin, "","",40,200)
                    respG2Bin.Fit  (meanG2Bin, "","",60,160)
                    respG2Bin.Fit  (meanG2Bin, "","",40,200)

                    if "part" in gen:
                        resolG1Bin.Fit (sigmaG1Bin,"","",100,300)
                    else:
                        resolG1Bin.Fit (sigmaG1Bin,"","",150,300)
                    resolG2Bin.Fit (sigmaG2Bin,"","",100,300)

                elif i_eta==1:
                    respG1Bin.Fit  (meanG1Bin, "","",80,200)
                    respG1Bin.Fit  (meanG1Bin, "","",60,300)
                    respG1Bin.Fit  (meanG1Bin, "","",40,300)
                    respG2Bin.Fit  (meanG2Bin, "","",80,200)
                    respG2Bin.Fit  (meanG2Bin, "","",60,300)
                    respG2Bin.Fit  (meanG2Bin, "","",40,300)

                    if "part" in gen:
                        resolG1Bin.Fit (sigmaG1Bin,"","",30,  300)
                        resolG1Bin.Fit (sigmaG1Bin,"","",60,  300)
                        resolG1Bin.Fit (sigmaG1Bin,"","",80,  300)
                    else:
                        resolG1Bin.Fit (sigmaG1Bin,"","",30,  300)

                    resolG2Bin.Fit (sigmaG2Bin,"","",130,300)
                    resolG2Bin.Fit (sigmaG2Bin,"","",130,300)
                # End of i_eta==1
            # End of jetType == "Heavy"

            # loop over bins in hCorr
            for j in range(1,hCorr.GetNbinsX()+1):

                fout.cd()  
                en = hCorr.GetBinCenter(j)

                hCorr_py = fout.Get("Bin{0}{1}_{2}/_py{1}_{2}".format(jetType, j, i_eta))

                if not hCorr_py:
                    print "No histo", jetType, j, i_eta
                    continue
                
                if jetType == "Light":                
                    resolDG_old =  hCorr_py.GetFunction("resol")
                elif jetType == "Heavy":                
                    resolDG_old =  hCorr_py.GetFunction("resolDG")
                elif jetType == "BDRSSub":                
                    resolDG_old =  hCorr_py.GetFunction("resol")


                if (not resolDG_old):
                    print "No func"
                    continue
                    
                resolDG = 0
                resolSG = 0

                if jetType == "Light":
                    resolSG = ROOT.TF1("resolSG_fit", SG_formula_Bin , -800, 800)
                    resolSG.SetNpx(1000)
                    resolSG.SetLineColor(ROOT.kMagenta)
                    norm = resolDG_old.Integral(-800,800)	
                    m = meanBin.Eval(en)
                    s = sigmaBin.Eval(en)
                    resolSG.SetParameter(0, norm/math.sqrt(2*ROOT.TMath.Pi())/s )
                    resolSG.SetParameter(1, m)
                    resolSG.SetParameter(2, s)

                    fout.cd( "BinLight{0}_{1}".format(j, i_eta) )
                    resolDG_old.Write("resol{0}_{1}".format( j,i_eta),ROOT.TObject.kOverwrite)
                    resolSG    .Write("resolSG_fit{0}_{1}".format(j,i_eta),ROOT.TObject.kOverwrite)


                elif jetType == "Heavy":
                    resolDG = ROOT.TF1("resolDG_fit", DG_formula_Bin , -800,800)
                    resolDG.SetNpx(1000)
                    resolDG.SetLineColor(ROOT.kBlue)


                    norm = resolDG_old.Integral(-800,800)
                    m1   = meanG1Bin.Eval(en)
                    m2   = meanG2Bin.Eval(en)
                    s1   = sigmaG1Bin.Eval(en)
                    s2   = sigmaG2Bin.Eval(en)

                    if i_eta==0:
                        resolDG.SetParameter(0, norm/(relWeightBin0*math.sqrt(2*ROOT.TMath.Pi())*s1 + (1-relWeightBin0)*math.sqrt(2*ROOT.TMath.Pi())*s2 )  )
                    elif i_eta==1:
                        resolDG.SetParameter(0, norm/(relWeightBin1*math.sqrt(2*ROOT.TMath.Pi())*s1 + (1-relWeightBin1)*math.sqrt(2*ROOT.TMath.Pi())*s2 )  )

                    resolDG.SetParameter(1, m1)
                    resolDG.SetParameter(2, s1)
                    resolDG.SetParameter(3, s2)
                    resolDG.SetParameter(4, m2)

                    resolSG = ROOT.TF1("resolSG_fit", SG_formula_Bin , -800, 800)
                    resolSG.SetNpx(1000)
                    resolSG.SetLineColor(ROOT.kMagenta)	
                    m = meanBin.Eval(en)
                    s = sigmaBin.Eval(en)
                    resolSG.SetParameter(0, norm/math.sqrt(2*ROOT.TMath.Pi())/s )
                    resolSG.SetParameter(1, m)
                    resolSG.SetParameter(2, s)

                    fout.cd( "BinHeavy{0}_{1}".format(j, i_eta) )
                    resolDG_old.Write("resolDG{0}_{1}".format(j,i_eta),    ROOT.TObject.kOverwrite)
                    resolDG    .Write("resolDG_fit{0}_{1}".format(j,i_eta),ROOT.TObject.kOverwrite)
                    resolSG    .Write("resolSG_fit{0}_{1}".format(j,i_eta),ROOT.TObject.kOverwrite)


                elif jetType == "BDRSSub":
                    resolSG = ROOT.TF1("resolSG_fit", SG_formula_Bin , -800, 800)
                    resolSG.SetNpx(1000)
                    resolSG.SetLineColor(ROOT.kMagenta)
                    norm = resolDG_old.Integral(-800,800)	
                    m = meanBin.Eval(en)
                    s = sigmaBin.Eval(en)
                    resolSG.SetParameter(0, norm/math.sqrt(2*ROOT.TMath.Pi())/s )
                    resolSG.SetParameter(1, m)
                    resolSG.SetParameter(2, s)

                    fout.cd( "BinBDRSSub{0}_{1}".format(j, i_eta) )
                    resolDG_old.Write("resol{0}_{1}".format( j,i_eta),ROOT.TObject.kOverwrite)
                    resolSG    .Write("resolSG_fit{0}_{1}".format(j,i_eta),ROOT.TObject.kOverwrite)

                # End of jetType==BDRSSub
            # End of loop over bins in hCorr

            if doOnlyJetTF:
                fout.cd()  
                resolBin.Write("",ROOT.TObject.kOverwrite)
                respBin. Write("",ROOT.TObject.kOverwrite)

                if jetType == "Heavy":
                    respG1Bin. Write("",ROOT.TObject.kOverwrite)
                    respG2Bin. Write("",ROOT.TObject.kOverwrite)
                    resolG1Bin.Write("",ROOT.TObject.kOverwrite)
                    resolG2Bin.Write("",ROOT.TObject.kOverwrite)
                # End of jetType == Heavy
            # End of doOnlyJetTF
        # End of loop over eta bins
    # End of loop over jetType
  
    if doOnlyJetTF:
        fout.Close()
        return
# end of calcTF


#########################################
# Allow use from command line and as import
#########################################

# If callled from the command-line: execute calcTF
if __name__ == "__main__":
    
    # Default
    # calcTF()
    
    # Testing of BDRSSub
    calcTF(
        outFileName   =  "transferFunctions",
        pathToFile    = "/shome/gregor/data/tth_htobb_m125_8tev_v93.root",
        li_jet_types  = ["BDRSSub"]
    )



# end of handling command line
