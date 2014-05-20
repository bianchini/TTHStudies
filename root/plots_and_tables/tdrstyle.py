#!/usr/bin/env python

from ROOT import gROOT, gStyle

def tdrstyle():
    gROOT.SetStyle("Plain")
    
    gStyle.SetAxisColor(1, "XYZ")
    
    gStyle.SetCanvasColor(0)
    #gStyle.SetCanvasBorderSize(10)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetCanvasDefH(700)
    gStyle.SetCanvasDefW(700)
    gStyle.SetCanvasDefX(0)
    gStyle.SetCanvasDefY(0)
    
    gStyle.SetFitFormat("5.4g")
    gStyle.SetFuncColor(2)
    gStyle.SetFuncStyle(1)
    gStyle.SetFuncWidth(1)
    
    gStyle.SetFrameBorderMode(0)
    gStyle.SetFrameBorderSize(1)
    gStyle.SetFrameFillStyle(0)
    gStyle.SetFrameFillColor(0)
    gStyle.SetFrameLineColor(1)
    gStyle.SetFrameLineStyle(1)  # 0?
    gStyle.SetFrameLineWidth(1)  # 1?
    
    gStyle.SetGridColor(0)
    gStyle.SetGridStyle(3)
    gStyle.SetGridWidth(1)
    
    #gStyle.SetHistFillColor(1)
    #gStyle.SetHistFillStyle(0)
    gStyle.SetHistLineColor(1)
    gStyle.SetHistLineStyle(0)
    gStyle.SetHistLineWidth(1)
    
    gStyle.SetLabelColor(1, "XYZ")
    gStyle.SetLabelFont(42,"XYZ")
    gStyle.SetLabelOffset(0.007,"XYZ")  # 0.010?
    gStyle.SetLabelSize(0.05,"XYZ")  # 0.04?
    
    gStyle.SetLegendBorderSize(0)
    gStyle.SetLegendFillColor(0)
    gStyle.SetLegendFont(42)
    
    gStyle.SetMarkerSize(1.0)
    gStyle.SetMarkerStyle(20)
    
    gStyle.SetLineColor(1)
    gStyle.SetLineWidth(2)
    #gStyle.SetLineScalePS(2)
    
    gStyle.SetOptDate(0)
    gStyle.SetOptFile(0)
    gStyle.SetOptFit(1)
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    #gStyle.SetOptLogx(0)
    #gStyle.SetOptLogy(0)
    #gStyle.SetOptLogz(0)
    
    gStyle.SetPadColor(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetPadBorderSize(10)
    gStyle.SetPadTopMargin(0.05)  # 0.08?
    gStyle.SetPadBottomMargin(0.13)
    gStyle.SetPadLeftMargin(0.16)
    gStyle.SetPadRightMargin(0.03)  # 0.05?
    gStyle.SetPadGridX(0)
    gStyle.SetPadGridY(0)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    
    gStyle.SetStatColor(0)
    gStyle.SetStatFont(42)
    gStyle.SetStatFontSize(0.025)
    gStyle.SetStatTextColor(1)
    gStyle.SetStatFormat("6.4g")
    gStyle.SetStatBorderSize(1)
    gStyle.SetStatH(0.1)
    gStyle.SetStatW(0.15)
    #gStyle.SetStatX(0)
    #gStyle.SetStatY(0)
    
    #gStyle.SetTextSize(0.055)
    gStyle.SetTextFont(42)
    
    gStyle.SetTitleBorderSize(0)
    gStyle.SetTitleColor(1)
    gStyle.SetTitleFont(42)
    gStyle.SetTitleColor(1,"XYZ")
    gStyle.SetTitleFont(42,"XYZ")
    gStyle.SetTitleSize(0.06,"XYZ")  # 0.05?
    #gStyle.SetTitleOffset(1.4,"XYZ")
    gStyle.SetTitleOffset(0.9,"X")
    gStyle.SetTitleOffset(1.20,"Y")
    gStyle.SetTitleFillColor(10)
    gStyle.SetTitleFontSize(0.05)
    gStyle.SetTitleTextColor(1)
    #gStyle.SetTitleH(0)
    #gStyle.SetTitleW(0)
    #gStyle.SetTitleX(0)
    #gStyle.SetTitleY(0.985)
    #gStyle.SetTitleStyle(1001)
    
    gStyle.SetPalette(1)
    #gStyle.SetNdivisions(510, "XYZ")  # 505?
    gStyle.SetNdivisions(505, "XYZ")
    gStyle.SetEndErrorSize(0)  # 2?
    #gStyle.SetErrorMarker(20)
    #gStyle.SetErrorX(0.)
    #gStyle.SetPaperSize(20.,20.)
    gStyle.SetStripDecimals(1)
    gStyle.SetTickLength(0.03, "XYZ")
    return 1

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i+lv/3], 16) for i in range(0, lv, lv/3))

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb
