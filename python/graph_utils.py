import ROOT as rt
from array import array

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
        x = list(g1.GetX())[bin]
        j = g2.GetXaxis().FindBin(x)
        if j < 1 or j > g2.GetN():
            continue
        num = list(g1.GetY())[bin]
        den = list(g2.GetY())[bin]
        if den!=0 and num!=0:
            r = num/den
            # dr = r*oplus(list(g1.GetEY())[bin]/num,list(g2.GetEY())[bin]/den)
        gratio.SetPoint(bin, x, r)
        # gratio.SetPointError(bin, x, dr)
    return gratio

def MakeRatioGraphFunc(graph, func, scale = lambda x: 1):
    ratio = rt.TGraphErrors()
    for i in range(graph.GetN()):
        x = list(graph.GetX())[i]
        y = list(graph.GetY())[i]
        ex, ey = graph.GetErrorX(i), graph.GetErrorY(i)
        r, dr = (0,0)
        if func.Eval(x) != 0:
            r= y/(func.Eval(x)*scale(x))
            dr= ey/(func.Eval(x)*scale(x))
        ratio.SetPoint(i, x, r)
        ratio.SetPointError(i, ex, dr)
    return ratio

def TruncateGraph(graph, min, max):
    for i in reversed(range((graph.GetN()))):
        x = list(graph.GetX())[i]
        if (min>0 and x<min): graph.RemovePoint(i)
        if (max>0 and x>max): graph.RemovePoint(i)
    return graph

def MergeGraphs(graph1, graph2):
    points = []
    for i in range(graph1.GetN()):
        x = list(graph1.GetX())[i]
        y = list(graph1.GetY())[i]
        ex, ey = graph1.GetErrorX(i), graph1.GetErrorY(i)
        points.append((x, y, ex, ey*1.10+(0.04 if y>0.9 else 0.002)))
    for i in range(graph2.GetN()):
        x = list(graph2.GetX())[i]
        y = list(graph2.GetY())[i]
        ex, ey = graph2.GetErrorX(i), graph2.GetErrorY(i)
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