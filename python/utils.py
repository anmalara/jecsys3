import os, ctypes
from collections import OrderedDict
from tdrstyle_JERC import *
import tdrstyle_JERC as TDR
rt.gROOT.SetBatch(rt.kTRUE)
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptFit(0)

from graph_utils import *
from printing_utils import *

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

def GetEmptyMarker(marker):
    if marker == rt.kFullTriangleUp:   return rt.kOpenTriangleUp
    if marker == rt.kFullCircle:       return rt.kOpenCircle
    if marker == rt.kFullSquare:       return rt.kOpenSquare
    if marker == rt.kFullTriangleDown: return rt.kOpenTriangleDown
    if marker == rt.kFullCross:        return rt.kOpenCross
    if marker == rt.kFullDiamond:      return rt.kOpenDiamond

def PrintFuncParameters(func, info=None, color=yellow):
    if info is None:
        info = f'  --> Fit results for {func.GetExpFormula()}:'
    print(color(info))
    for i in range(func.GetNpar()):
        text = f'    --> {i}: {func.GetParameter(i):.3f} +- {func.GetParError(i):.3f}'
        p_min,p_max = ctypes.c_double(0), ctypes.c_double(0)
        func.GetParLimits(i, p_min, p_max)
        if p_min.value!=p_max.value:
            text += f' in [{p_min.value:.3f},{p_max.value:.3f}]'
        elif p_min.value!=0:
            text += ' fixed'
        print(color(text))

if __name__ == "__main__":
    # do nothing
    pass