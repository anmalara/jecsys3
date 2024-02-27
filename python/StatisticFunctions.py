import ROOT
from array import array

def lnN(nominal, variation, minAllowed=0.001):
    if variation==0: raise Exception("var==0")
    # return 1+abs(nominal-variation)/nominal
    var = (nominal-variation)/nominal
    if abs(var)<minAllowed:
        var = minAllowed*var/abs(var)
    return 1+var

def CountBinsMinMax(hist, min, max):
    Nbins = 0
    xmin = 1e6
    xmax = 0
    for i in range(hist.GetNbinsX()+1):
        if hist.GetXaxis().GetBinLowEdge(i)>=min:
            if (xmin>1e5): xmin = hist.GetXaxis().GetBinLowEdge(i)
            if hist.GetXaxis().GetBinUpEdge(i)<=max:
                xmax = hist.GetXaxis().GetBinUpEdge(i)
                Nbins += 1
    return (Nbins,xmin,xmax)

def GetBinEdges(hist, fmin, fmax):
    edges = []
    for i in range(1, hist.GetNbinsX()+1):
        edge = hist.GetXaxis().GetBinLowEdge(i)
        if edge>=fmin and edge<=fmax:
            edges.append(edge)
    return edges

def GetBinEdgesGraph(graph, fmin, fmax):
    edges = []
    x_values = list(graph.GetX())
    for i in range(graph.GetN()):
        x = x_values[i]
        ex = graph.GetErrorX(i)
        edge = abs(x-ex)
        if edge>=fmin and edge<=fmax:
            edges.append(edge)
    edges.append(x_values[-1]+ex)
    return edges


def GetConfidenceIntervals(func, fitRes, Nbins, xcenters, ci, cl = 0.68, func2=None, fitRes2=None):
    def get_func_info(func, fitRes):
        npar = func.GetNumberFreeParameters()
        npar_real = func.GetNpar()
        if npar_real != npar:
            npar = npar_real
        chi2 = func.GetChisquare()
        ndf = func.GetNDF()
        covmatr = fitRes.GetCovarianceMatrix()
        return npar, npar_real, chi2, ndf, covmatr
    
    isfit2 = func2 is not None
    fixed = 0
    npar, npar_real, chi2, ndf, covmatr = get_func_info(func, fitRes)
    if isfit2:
        npar2, npar_real2, chi2_2, ndf2, covmatr2 = get_func_info(func2, fitRes2)
        ndf += ndf2
        chi2 += chi2_2

    tStudent = ROOT.TMath.StudentQuantile(0.5 + cl/2, ndf)
    chindf = ROOT.TMath.Sqrt(chi2/ndf)

    def setup_vectors(npar_real, npar):
        grad = array('d',[0.]*npar_real)
        sum_vector = array('d',[0.]*npar)
        return grad, sum_vector

    grad, sum_vector = setup_vectors(npar_real, npar)
    if isfit2:
        grad2, sum_vector2 = setup_vectors(npar_real2, npar2)
    
    def update_vectors(ipoint, npar, func, covmatr, grad, sum_vector):
        func.GradientPar(array('d',[xcenters[ipoint]]), grad)
        # multiply the covariance matrix by gradient
        for irow in range(npar):
            sum_vector[irow]=0
            for icol in range(npar):
                igrad=0
                ifree=0
                if fixed:
                    print("ERROR. Check other function.")
                else:
                    igrad = icol
                sum_vector[irow] += covmatr[irow][icol]*grad[igrad]

    def calculate_ci(npar, grad, sum_vector):
        ci_=0
        igrad = 0
        for i_ in range(npar):
            igrad=0
            ifree=0
            if fixed:
                print("ERROR. Check other function")
            else:
                igrad = i_
            ci_ += grad[igrad]*sum_vector[i_]
        ci_ = ROOT.TMath.Sqrt(ci_)
        return ci_

    for ipoint in range(Nbins):
        update_vectors(ipoint=ipoint, npar=npar, func=func, covmatr=covmatr, grad=grad, sum_vector=sum_vector)
        if isfit2:
            update_vectors(ipoint=ipoint, npar=npar2, func=func2, covmatr=covmatr2, grad=grad2, sum_vector=sum_vector2)
        ci_ = calculate_ci(npar, grad, sum_vector)
        ci_ /= func.Eval(xcenters[ipoint])
        ci_ = ci_**2
        if isfit2:
            ci2_ = calculate_ci(npar2, grad2, sum_vector2)
            ci2_ /= func2.Eval(xcenters[ipoint])
            ci2_ = ci2_**2
            ci_ += ci2_
        ref = func.Eval(xcenters[ipoint])
        if isfit2:
            ref /= func2.Eval(xcenters[ipoint])
        ci[ipoint] = ROOT.TMath.Sqrt(ci_)*tStudent*chindf * ref

def ComputeHistWithCL(name, func, fitRes, hist, cl=0.68, do_poisson=False, func2=None, fitRes2=None):
    # create a histogram for the fit region only
    # Nbins,xmin,xmax = CountBinsMinMax(hist, func.GetXmin(), func.GetXmax())
    # band = ROOT.TH1F("band"+name+str(cl),name+str(cl), Nbins, xmin, xmax)
    # band_pull  = ROOT.TH1F("band_pull"+name+str(cl),name+str(cl), Nbins, xmin, xmax)
    # pull = ROOT.TH1F("pull"+name+str(cl),name+str(cl), Nbins, xmin, xmax)

    isGraph = isinstance(hist, ROOT.TGraph)
    if isGraph:
        graph = hist
        edges = GetBinEdgesGraph(graph, func.GetXmin(), func.GetXmax())
        Nbins = len(edges)
    else:
        edges = GetBinEdges(hist, func.GetXmin(), func.GetXmax())
        Nbins = len(edges)
    edges = array('d',edges)
    ci = array('d',[0]*Nbins)

    band = ROOT.TH1F("band"+name+str(cl),name+str(cl), Nbins-1, edges)
    band_pull  = ROOT.TH1F("band_pull"+name+str(cl),name+str(cl), Nbins-1, edges)
    if do_poisson:
        pull = ROOT.TGraphAsymmErrors()
    else:
        pull = ROOT.TH1F("pull"+name+str(cl),name+str(cl), Nbins-1, edges)

    # now get an array with the bin centers and compute the CLs
    GetConfidenceIntervals(func, fitRes, Nbins, edges, ci, cl, func2=func2, fitRes2=fitRes2)
    
    if isGraph:
        x_values = list(graph.GetX())
        y_values = list(graph.GetY())
    bin_loop = range(1,Nbins)
    for i in bin_loop:
        if isGraph:
            x_ = x_values[i-1]
            y_hist = y_values[i-1]
            y_err = graph.GetErrorY(i-1)
        else:
            x_ = band.GetXaxis().GetBinCenter(i)
            ibin = hist.GetXaxis().FindBin(x_)
            y_hist = hist.GetBinContent(ibin)
            y_err =  hist.GetBinError(ibin)
        y_ = func.Eval(x_)
        if func2 is not None:
            y_ /= func2.Eval(x_)
        if y_==0: y_ = 1e-07
        band.SetBinContent(i, y_)
        band.SetBinError(i, ci[i-1])
        band_pull.SetBinContent(i, 1)
        band_pull.SetBinError(i, ci[i-1]/y_)
        if do_poisson:
            alpha = 1 - 0.6827
            width = hist.GetBinWidth(ibin)
            N = y_hist*width
            L =  0  if N==0 else ROOT.Math.gamma_quantile(alpha/2,N,1.)
            U =  ROOT.Math.gamma_quantile_c(alpha/2,N+1,1)
            pull.SetPoint(i-1, x_, N/width/y_)
            pull.SetPointEYlow(i-1, (N-L)/width/y_)
            pull.SetPointEYhigh(i-1, (U-N)/width/y_)
        else:
            pull.SetBinContent(i, y_hist/y_)
            pull.SetBinError(i, y_err/y_)
    # for i in range(1, band.GetNbinsX()+1):
    #     print(name, i, band.GetBinCenter(i), band.GetBinContent(i), band.GetBinError(i))
    return band, band_pull, pull



def CrystalBall_Fit(x, par):
    mean = par[0]
    sigma = par[1]
    alpha = par[2]
    n = par[3]
    norm = par[4]
    std =(x[0]-mean)/sigma
    if (alpha < 0): std *= -1
    alpha = ROOT.TMath.Abs(alpha)
    result = 0
    if std >= -alpha:
        result = ROOT.TMath.Exp(-0.5*std*std)
    else:
        A = ROOT.TMath.Power(n/alpha, n)*ROOT.TMath.Exp(-0.5*alpha*alpha)
        B = n/alpha-alpha
        result = A/ROOT.TMath.Power(B-std, n)
    return norm*result


def ExpGaussExp_Fit(x, par):
    mean = par[0]
    sigma = par[1]
    kl = par[2]
    kh = par[3]
    norm = par[4]
    std =(x[0]-mean)/sigma
    result = 0
    if std < -kl:
        result = ROOT.TMath.Exp(kl*kl*0.5+kl*std)
    elif -kl <= std < kh:
        result = ROOT.TMath.Exp(-0.5*std*std)
    else:
        result = ROOT.TMath.Exp(kh*kh*0.5-kh*std)
    return norm*result