import os
from math import sqrt
from collections import OrderedDict
from tdrstyle_JERC import *
import tdrstyle_JERC as TDR
rt.gROOT.SetBatch(rt.kTRUE)
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptFit(0)


def ominus(a,b):
    return sqrt(max(a*a-b*b,0.0))

def oplus(a,b):
    return sqrt(a*a+b*b)

def HistToGraph(hist):
    graph = rt.TGraphErrors()
    for bin in range(1,hist.GetNbinsX()+1):
        x = hist.GetBinCenter(bin)
        y = hist.GetBinContent(bin)
        ex = hist.GetBinWidth(bin)/2
        ey = hist.GetBinError(bin)
        graph.SetPoint(bin-1, x, y)
        graph.SetPointError(bin-1, ex, ey)
    return graph
    


def MakeRatioHistograms(h1, h2, name):
    hratio = h1.Clone(name)
    for bin in range(1,h1.GetNbinsX()+1):
        r, dr = (0,0)
        num = h1.GetBinContent(bin)
        den = h2.GetBinContent(bin) 
        if den!=0 and num!=0:
            r = num/den
            dr = r*oplus(h1.GetBinError(bin)/num,h2.GetBinError(bin)/den)
        hratio.SetBinContent(bin, r)
        hratio.SetBinError(bin, dr)
    return hratio


def MakeRatioGraphs(g1, g2, name):
    gratio = g1.Clone(name)
    for bin in range(0,g1.GetN()):
        r, dr = (0,0)
        x = list(gratio.GetX())[bin]
        num = list(g1.GetY())[bin]
        den = list(g2.GetY())[bin]
        if den!=0 and num!=0:
            r = num/den
            # dr = r*oplus(list(g1.GetEY())[bin]/num,list(g2.GetEY())[bin]/den)
        gratio.SetPoint(bin, x, r)
        # gratio.SetPointError(bin, x, dr)
    return gratio