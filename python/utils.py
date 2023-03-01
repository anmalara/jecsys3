import os
from math import sqrt
from collections import OrderedDict
from tdrstyle_JERC import *
import tdrstyle_JERC as TDR
rt.gROOT.SetBatch(rt.kTRUE)
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptFit(0)

from graph_utils import *


def ominus(a,b):
    return sqrt(max(a*a-b*b,0.0))

def oplus(a,b):
    return sqrt(a*a+b*b)

def FixXAsisPartition(canv, shift=None, textsize=0.05, bins=[30,100,300,1000, 3000]):
        canv.SetLogx(True)
        GettdrCanvasHist(canv).GetXaxis().SetNoExponent(True)
        latex = rt.TLatex()
        latex.SetTextFont(42)
        latex.SetTextSize(textsize)
        latex.SetTextAlign(23)
        if shift==None:
            YMin, YMax = (GettdrCanvasHist(canv).GetYaxis().GetXmin(), GettdrCanvasHist(canv).GetYaxis().GetXmax())
            shift = YMin-0.018*(YMax-YMin)
        for xbin in bins:
            latex.DrawLatex(xbin,shift,str(xbin))


if __name__ == "__main__":
    # do nothing
    pass