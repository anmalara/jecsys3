import ROOT as rt
import numpy as np
from array import array
from math import sqrt


def ominus(a, b):
    return sqrt(max(a * a - b * b, 0.0))


def oplus(a, b):
    return sqrt(a * a + b * b)


def RemoveFunc(graph):
    graph.GetListOfFunctions().Remove(list(graph.GetListOfFunctions())[0])


def FuncToGraph(func, bins):
    graph = rt.TGraph()
    for bin, x in enumerate(bins):
        y = func.Eval(x)
        graph.SetPoint(bin, x, y)
    return graph


def HistToGraph(hist):
    graph = rt.TGraphErrors()
    for bin in range(1, hist.GetNbinsX() + 1):
        x = hist.GetBinCenter(bin)
        y = hist.GetBinContent(bin)
        ex = hist.GetBinWidth(bin) / 2
        ey = hist.GetBinError(bin)
        if y == 0 and ey == 0:
            continue
        graph.SetPoint(bin - 1, x, y)
        graph.SetPointError(bin - 1, ex, ey)
    return graph


def Hist2DToGraphAsymm(h2d, xmin, xmax, yerr=None):
    y_ax = h2d.GetYaxis()
    x_vals, y_vals, x_errs_lo, x_errs_hi, y_errs = [], [], [], [], []
    for x_bin in range(1, h2d.GetNbinsX() + 2):
        if h2d.GetXaxis().GetBinLowEdge(x_bin) != float(xmin):
            continue
        if h2d.GetXaxis().GetBinLowEdge(x_bin + 1) != float(xmax):
            continue
        for y_bin in range(1, h2d.GetNbinsY() + 2):
            y_ = h2d.GetBinContent(x_bin, y_bin)
            if y_ == 0:
                continue
            x_vals.append(y_ax.GetBinCenter(y_bin))
            y_vals.append(y_)
            x_errs_lo.append(y_ax.GetBinCenter(y_bin) - y_ax.GetBinLowEdge(y_bin))
            x_errs_hi.append(y_ax.GetBinLowEdge(y_bin + 1) - y_ax.GetBinCenter(y_bin))
            y_errs.append(0.2 if yerr else h2d.GetBinError(x_bin, y_bin))
    graph = rt.TGraphAsymmErrors(
        len(x_vals), array("d", x_vals), array("d", y_vals), array("d", x_errs_lo), array("d", x_errs_hi), array("d", y_errs), array("d", y_errs)
    )
    return graph


def Hist2DToGraphAsymmAverage(h2d, xmin, xmax, yerr=None):
    y_ax = h2d.GetYaxis()
    x_vals, y_vals, x_errs_lo, x_errs_hi, y_errs = [], [], [], [], []
    x_bin_min, x_bin_max = None, None
    import pdb

    # pdb.set_trace()
    for x_bin in range(1, h2d.GetNbinsX() + 2):
        if h2d.GetXaxis().GetBinLowEdge(x_bin) == float(xmin):
            x_bin_min = x_bin
        if h2d.GetXaxis().GetBinLowEdge(x_bin + 1) == float(xmax):
            x_bin_max = x_bin
    # pdb.set_trace()
    for y_bin in range(1, h2d.GetNbinsY() + 2):
        y_, y_err = [], []
        for x_bin in range(x_bin_min, x_bin_max + 1, 1):
            bin_content = h2d.GetBinContent(x_bin, y_bin)
            bin_error = h2d.GetBinError(x_bin, y_bin)
            if bin_content == 0 or bin_error == 0:
                continue
            y_.append(bin_content)
            y_err.append(bin_error)
        if len(y_) == 0:
            continue
        # pdb.set_trace()
        y_ = np.array(y_)
        y_err = np.array(y_err)
        y_err = np.mean(y_err) + np.std(y_)
        y_err = np.std(y_)
        y_ = np.mean(y_)
        x_vals.append(y_ax.GetBinCenter(y_bin))
        y_vals.append(y_)
        x_errs_lo.append(y_ax.GetBinCenter(y_bin) - y_ax.GetBinLowEdge(y_bin))
        x_errs_hi.append(y_ax.GetBinLowEdge(y_bin + 1) - y_ax.GetBinCenter(y_bin))
        y_errs.append(0.2 if yerr else h2d.GetBinError(x_bin, y_bin))
    graph = rt.TGraphAsymmErrors(
        len(x_vals), array("d", x_vals), array("d", y_vals), array("d", x_errs_lo), array("d", x_errs_hi), array("d", y_errs), array("d", y_errs)
    )
    return graph


def MakeRatioHistograms(h1, h2, name):
    hratio = h1.Clone(name)
    for bin in range(1, h1.GetNbinsX() + 1):
        r, dr = (0, 0)
        num = h1.GetBinContent(bin)
        den = h2.GetBinContent(bin)
        if den != 0 and num != 0:
            r = num / den
            dr = r * oplus(h1.GetBinError(bin) / num, h2.GetBinError(bin) / den)
        hratio.SetBinContent(bin, r)
        hratio.SetBinError(bin, dr)
    return hratio


def MakeRatioGraphs(g1, g2, name):
    gratio = g1.Clone(name)
    center_x1 = list(g1.GetX())
    center_x2 = list(g2.GetX())
    common_indices_1 = [i for i, x in enumerate(center_x1) if x in center_x2]
    common_indices_2 = [i for i, x in enumerate(center_x2) if x in center_x1]
    if len(common_indices_1) != len(common_indices_2):
        print("center_x1", center_x1)
        print("center_x2", center_x2)
        print("common_indices_1", common_indices_1)
        print("common_indices_2", common_indices_2)
        raise RuntimeError("Impossible to make ratio. Please fix this")
    vals_num = list(g1.GetY())
    vals_den = list(g2.GetY())
    for index in range(len(common_indices_1)):
        r, dr = (0, 0)
        x = center_x1[index]
        num = vals_num[index]
        den = vals_den[index]
        index1 = common_indices_1[index]
        index2 = common_indices_2[index]
        if den != 0 and num != 0:
            r = num / den
            dr = r * oplus(g1.GetErrorY(index1) / num, g2.GetErrorY(index2) / den)
        gratio.SetPoint(index, x, r)
        if type(gratio) is rt.TGraphAsymmErrors:
            gratio.SetPointError(index, g1.GetErrorXlow(index1), g1.GetErrorXhigh(index2), dr, dr)
        else:
            gratio.SetPointError(index, g1.GetErrorX(index1), dr)
    return gratio


def MakeRatioGraphFunc(graph, func, isAsymmError=True, scale=lambda x: 1):
    ratio = rt.TGraphAsymmErrors() if isAsymmError else rt.TGraphErrors()
    for i in range(graph.GetN()):
        x = list(graph.GetX())[i]
        y = list(graph.GetY())[i]
        ex, ey = graph.GetErrorX(i), graph.GetErrorY(i)
        r, dr = (0, 0)
        if func.Eval(x) != 0:
            r = y / (func.Eval(x) * scale(x))
            dr = ey / (func.Eval(x) * scale(x))
        ratio.SetPoint(i, x, r)
        if isAsymmError:
            ratio.SetPointError(i, graph.GetErrorXlow(i), graph.GetErrorXhigh(i), dr, dr)
        else:
            ratio.SetPointError(i, ex, dr)
    return ratio


def TruncateGraph(graph, min, max):
    for i in reversed(range((graph.GetN()))):
        x = list(graph.GetX())[i]
        if min > 0 and x < min:
            graph.RemovePoint(i)
        if max > 0 and x > max:
            graph.RemovePoint(i)
    return graph


def RemovePointFromGraph(graph, min, max):
    for i in reversed(range((graph.GetN()))):
        x = list(graph.GetX())[i]
        if min <= x and x <= max:
            graph.RemovePoint(i)
    return graph


def MergeGraphs(graph1, graph2, scale_err=None):
    points = []
    for i in range(graph1.GetN()):
        x = list(graph1.GetX())[i]
        y = list(graph1.GetY())[i]
        ex, ey = graph1.GetErrorX(i), graph1.GetErrorY(i)
        if scale_err:
            ey *= scale_err
        points.append((x, y, ex, ey * 1.10 + (0.04 if y > 0.9 else 0.002)))
    for i in range(graph2.GetN()):
        x = list(graph2.GetX())[i]
        y = list(graph2.GetY())[i]
        ex, ey = graph2.GetErrorX(i), graph2.GetErrorY(i)
        if scale_err:
            ey *= scale_err
        points.append((x, y, ex, ey * 1.10 + (0.04 if y > 0.9 else 0.002)))

    points.sort(key=lambda p: p[0])
    x = array("d", [p[0] for p in points])
    y = array("d", [p[1] for p in points])
    ex = array("d", [p[2] for p in points])
    ey = array("d", [p[3] for p in points])
    combined_graph = rt.TGraphErrors(len(points), x, y, ex, ey)
    return combined_graph


def ShiftHist(hist, shift):
    new_hist = hist.Clone(hist.GetName() + "_shifted")
    for bin in range(1, hist.GetNbinsX() + 1):
        new_hist.SetBinContent(bin, hist.GetBinContent(bin) + shift)
    return new_hist


if __name__ == "__main__":
    # do nothing
    pass
