import ROOT as rt
from array import array
from math import sqrt

def ominus(a,b):
    return sqrt(max(a*a-b*b,0.0))

def oplus(a,b):
    return sqrt(a*a+b*b)

def RemoveFunc(graph):
    graph.GetListOfFunctions().Remove(list(graph.GetListOfFunctions())[0])

def HistToGraph(hist):
    graph = rt.TGraphErrors()
    for bin in range(1,hist.GetNbinsX()+1):
        x = hist.GetBinCenter(bin)
        y = hist.GetBinContent(bin)
        ex = hist.GetBinWidth(bin)/2
        ey = hist.GetBinError(bin)
        if (y==0 and ey==0):
            continue
        graph.SetPoint(bin-1, x, y)
        graph.SetPointError(bin-1, ex, ey)
    return graph
    
def Hist2DToGraphAsymm(h2d, xmin, xmax, yerr=None):
    y_ax = h2d.GetYaxis()
    x_vals, y_vals, x_errs_lo, x_errs_hi, y_errs = [], [], [], [], []
    for x_bin in range(1, h2d.GetNbinsX()+2):
        if h2d.GetXaxis().GetBinLowEdge(x_bin)!= float(xmin): continue
        if h2d.GetXaxis().GetBinLowEdge(x_bin+1)!= float(xmax): continue
        for y_bin in range(1, h2d.GetNbinsY()+2):
            y_ = h2d.GetBinContent(x_bin,y_bin)
            if y_==0: continue
            x_vals.append(y_ax.GetBinCenter(y_bin))
            y_vals.append(y_)
            x_errs_lo.append(y_ax.GetBinCenter(y_bin)-y_ax.GetBinLowEdge(y_bin))
            x_errs_hi.append(y_ax.GetBinLowEdge(y_bin+1)- y_ax.GetBinCenter(y_bin))
            y_errs.append(0.2 if yerr else h2d.GetBinError(x_bin,y_bin))
    graph = rt.TGraphAsymmErrors(len(x_vals), array('d',x_vals), array('d', y_vals),  array('d', x_errs_lo),  array('d', x_errs_hi),  array('d', y_errs), array('d', y_errs))
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
    center_x = list(g1.GetX())
    vals_num = list(g1.GetY())
    vals_den = list(g2.GetY())
    for bin in range(0,g1.GetN()):
        r, dr = (0,0)
        x = center_x[bin]
        j = g2.GetXaxis().FindBin(x)
        if j < 1 or j > g2.GetN():
            if x == list(g2.GetX())[bin]:
                j = bin
            else:
                continue
        num = vals_num[bin]
        den = vals_den[bin]
        if den!=0 and num!=0:
            r = num/den
            dr = r*oplus(g1.GetErrorY(bin)/num,g2.GetErrorY(bin)/den)
        gratio.SetPoint(bin, x, r)
        if type(gratio) is rt.TGraphAsymmErrors:
            gratio.SetPointError(bin, g1.GetErrorXlow(bin), g1.GetErrorXhigh(bin), dr, dr)
        else:
            gratio.SetPointError(bin, g1.GetErrorX(bin), dr)
    return gratio

def MakeRatioGraphFunc(graph, func, isAsymmError=True, scale = lambda x: 1):
    ratio = rt.TGraphAsymmErrors() if isAsymmError else rt.TGraphErrors()
    for i in range(graph.GetN()):
        x = list(graph.GetX())[i]
        y = list(graph.GetY())[i]
        ex, ey = graph.GetErrorX(i), graph.GetErrorY(i)
        r, dr = (0,0)
        if func.Eval(x) != 0:
            r= y/(func.Eval(x)*scale(x))
            dr= ey/(func.Eval(x)*scale(x))
        ratio.SetPoint(i, x, r)
        if isAsymmError:
            ratio.SetPointError(i, graph.GetErrorXlow(i), graph.GetErrorXhigh(i), dr, dr)
        else:
            ratio.SetPointError(i, ex, dr)
    return ratio

def TruncateGraph(graph, min, max):
    for i in reversed(range((graph.GetN()))):
        x = list(graph.GetX())[i]
        if (min>0 and x<min): graph.RemovePoint(i)
        if (max>0 and x>max): graph.RemovePoint(i)
    return graph

def MergeGraphs(graph1, graph2, scale_err=None):
    points = []
    for i in range(graph1.GetN()):
        x = list(graph1.GetX())[i]
        y = list(graph1.GetY())[i]
        ex, ey = graph1.GetErrorX(i), graph1.GetErrorY(i)
        if scale_err:
            ey *= scale_err
        points.append((x, y, ex, ey*1.10+(0.04 if y>0.9 else 0.002)))
    for i in range(graph2.GetN()):
        x = list(graph2.GetX())[i]
        y = list(graph2.GetY())[i]
        ex, ey = graph2.GetErrorX(i), graph2.GetErrorY(i)
        if scale_err:
            ey *= scale_err
        points.append((x, y, ex, ey*1.10+(0.04 if y>0.9 else 0.002)))
    
    points.sort(key=lambda p: p[0])
    x = array('d',  [p[0] for p in points])
    y = array('d',  [p[1] for p in points])
    ex = array('d', [p[2] for p in points])
    ey = array('d', [p[3] for p in points])
    combined_graph = rt.TGraphErrors(len(points), x,y,ex,ey)
    return combined_graph

def ShiftHist(hist, shift):
    new_hist = hist.Clone(hist.GetName()+'_shifted')
    for bin in range(1,hist.GetNbinsX()+1):
        new_hist.SetBinContent(bin, hist.GetBinContent(bin) + shift)
    return new_hist


if __name__ == "__main__":
    # do nothing
    pass