from utils import *


def ProfileRMS(h2d, hname, scale=sqrt(2)*0.9):
    hproj = h2d.ProjectionX(hname)
    for bin in range(1,hproj.GetNbinsX()+1):
        htmp = h2d.ProjectionY('temp',bin,bin)
        rms = htmp.GetRMS()
        hproj.SetBinContent(bin, scale*rms)
        # hproj.SetBinError(bin, scale*htmp.GetRMSError()+0.0005+(0.1*htmp.GetRMSError()/rms if rms!=0 else 0))
        hproj.SetBinError(bin, scale*htmp.GetRMSError())
    return hproj

# def CalculateMPFX(h1, h2, name, RC):
#     hnew = h1.Clone(name)
#     for bin in range(1,h1.GetNbinsX()+1):
#         hnew.SetBinContent(bin, ominus(h1.GetBinContent(bin),h2.GetBinContent(bin)))
#         hnew.SetBinContent(bin, oplus(hnew.GetBinContent(bin),RC/hnew.GetBinCenter(bin)))
#         hnew.SetBinError(bin, oplus(h1.GetBinError(bin),h2.GetBinError(bin)))
#     return hnew

def CalculateMPFX(g1, g2, name, RC):
    gnew = g1.Clone(name)
    for bin in range(0,g1.GetN()):
        x = list(gnew.GetX())[bin]
        y1 = list(g1.GetY())[bin]
        y2 = list(g2.GetY())[bin]
        ex = list(g1.GetEX())[bin]
        ey = oplus(list(g1.GetEY())[bin], list(g2.GetEY())[bin])
        if (x==0): continue
        gnew.SetPoint(bin, x, oplus(ominus(y1,y2), RC/x) )
        gnew.SetPointError(bin, ex, ey)
    return gnew


def CalculateMPFX2(g1, g2, name, RC):
    gnew = g1.Clone(name)
    for bin in range(0,g1.GetN()):
        x = list(gnew.GetX())[bin]
        y1 = list(g1.GetY())[bin]
        y2 = list(g2.GetY())[bin]
        ex = list(g1.GetEX())[bin]
        ey = oplus(list(g1.GetEY())[bin], list(g2.GetEY())[bin])
        if (x==0): continue
        gnew.SetPoint(bin, x, oplus(oplus(y1,y2), RC/x) )
        gnew.SetPointError(bin, ex, ey)
    return gnew




class CombineJER():
    def __init__(self, run='Run2', year='UL18', algo='AK4 CHS', etamin = '0.0', etamax = '0.5'):
        self.run = run
        self.year = year
        self.algo = algo
        self.etamin = etamin
        self.etamax = etamax
        self.inputPath = './rootfiles/CombineJER/'
        self.outputPath = './pdfs/CombineJER/'
        if not os.path.exists(self.inputPath):
            self.inputPath = '../rootfiles/CombineJER/'
            self.outputPath = '../pdfs/CombineJER/'
        os.system('mkdir -p '+self.outputPath)
        TDR.extraText  = 'Preliminary'
        TDR.cms_lumi = TDR.commonScheme['legend'][self.year]+', '+TDR.commonScheme['lumi'][self.year]+' fb^{-1}'
        TDR.cms_energy = TDR.commonScheme['energy'][self.year]
        TDR.extraText3 = []
        TDR.extraText3.append(self.algo)
        TDR.extraText3.append(self.etamin+' < |#eta| < '+self.etamax)

        self.mu_bins = ['Mu0to10','Mu10to20','Mu20to30','Mu30to40','Mu40to50','Mu50to60']
        self.mu_bins_short = ['Mu0to10', 'Mu30to40', 'Mu50to60']
        self.mpfx_samples = ['low-PU', 'high-PU']
        self.balance_samples = ['dijet SM', 'dijet FE', 'zjet']
        self.all_types = ['Data', 'MC', 'Ratio']
        self.pdfextraname = ''

    def LoadMPFX(self):
        self.files = OrderedDict()
        self.graphs = OrderedDict()
        self.files['high-PU MC'] = rt.TFile(self.inputPath+'output-P8CP55to7000-2b-UL18V5V2_ABCD-noSF.root')
        self.files['high-PU Data'] = rt.TFile(self.inputPath+'output-DATA-2b-UL18V5V2_ABCD-noSF.root')
        self.files['low-PU MC'] = rt.TFile(self.inputPath+'Sevgi-output-MC-mpf_mpfx_2Dand3D.root')
        self.files['low-PU Data'] = rt.TFile(self.inputPath+'Sevgi-output-DATA-mpf_mpfx_2Dand3D.root')
        folder = 'Standard/Eta_'+self.etamin+'-'+self.etamax
        for fname, f_ in self.files.items():
            folder_ = folder
            for mode in ['MPF','MPFX']:
                name = fname+' '+mode
                h2d = f_.Get(folder_+'/h2'+mode.lower())
                h2d.RebinX(2)
                # self.graphs[name] = rt.TGraphErrors(ProfileRMS(h2d,name))
                self.graphs[name] = HistToGraph(ProfileRMS(h2d,name))
                # self.graphs[name].SetDirectory(0)
                if 'Data' in name:
                    self.graphs[name.replace('Data', 'Ratio')] = MakeRatioGraphs(self.graphs[name], self.graphs[name.replace('Data', 'MC')], name.replace('Data', 'Ratio'))
                    # self.graphs[name.replace('Data', 'Ratio')].SetDirectory(0)
            name = fname+' JER'
            RC = self.data_RC if 'Data' in name else self.mc_RC
            if 'low-PU' in fname: RC=0
            self.graphs[name] = CalculateMPFX(self.graphs[name.replace('JER', 'MPF')], self.graphs[name.replace('JER', 'MPFX')], name, RC)
            # self.graphs[name].SetDirectory(0)
            if 'Data' in name:
                    self.graphs[name.replace('Data', 'Ratio')] = MakeRatioGraphs(self.graphs[name], self.graphs[name.replace('Data', 'MC')], name.replace('Data', 'Ratio'))
                    # self.graphs[name.replace('Data', 'Ratio')].SetDirectory(0)
        

        for type in ['MC','Data']:
            name1 = ' '.join(['high-PU',type,'JER'])
            name2 = ' '.join(['low-PU',type,'JER'])
            name  = ' '.join(['MPFX',type,'Ratio'])
            self.graphs[name] = MakeRatioGraphs(self.graphs[name1], self.graphs[name2], name)
    
    def LoadRC(self):
        self.data_RC = 0
        self.mc_RC = 0
        f_RC = rt.TFile(self.inputPath+'RC.root')
        graphs_RC= {}
        for type in ['Data', 'MC']:
            graphs_RC[type] = f_RC.Get(type+'/RMS')
        count = 0
        for bin in range(1,graphs_RC['MC'].GetN()):
            eta = list(graphs_RC['Data'].GetX())[bin]
            if eta< float(self.etamin) or eta> float(self.etamax): continue
            dtrms = list(graphs_RC['Data'].GetY())[bin]
            mcrms = list(graphs_RC['MC'].GetY())[bin]
            self.data_RC = oplus(self.data_RC, dtrms)
            self.mc_RC = oplus(self.mc_RC, mcrms)
            count += 1
        self.data_RC /=count
        self.mc_RC /=count
        # self.data_RC *=1.5
        # self.mc_RC = 0
        # self.data_RC = 0
        f_RC.Close()
    
    def LoadDijet(self):
        self.files['dijet'] = rt.TFile(self.inputPath+'dijet_balance_UL18.root')
        for mode in ['SM', 'FE']:
            for type in ['Data', 'MC']:
                name = ' '.join(['dijet',mode,type,'JER'])
                self.graphs[name] = self.files['dijet'].Get('dijet_balance_jer_'+type+'_0p00000_0p261_'+mode+'_nominal')
                # self.graphs[name] = self.files['dijet'].Get('dijet_balance_jer_'+type+'_2p65_2p2_'+mode+'_nominal')
                if type=='MC':
                    self.graphs[name.replace('MC', 'Ratio')] = MakeRatioGraphs(self.graphs[name.replace('MC', 'Data')], self.graphs[name], name.replace('MC', 'Ratio'))
        self.files['dijet CHS'] = rt.TFile(self.inputPath+'dijet_balance_UL18_Summer19UL18_V5_AK4CHS.root')
        for mode in ['SM', 'FE']:
            for type in ['Data', 'MC']:
                name = ' '.join(['dijet CHS',mode,type,'JER'])
                self.graphs[name] = self.files['dijet CHS'].Get('dijet_balance_jer_'+type+'_0p261_0p522_'+mode+'_nominal')
                if type=='MC':
                    self.graphs[name.replace('MC', 'Ratio')] = MakeRatioGraphs(self.graphs[name.replace('MC', 'Data')], self.graphs[name], name.replace('MC', 'Ratio'))
        
        self.files['yannick'] = rt.TFile(self.inputPath+'MPF+MPFx.root')
        for type in ['MC']:
            name = ' '.join(['yannick',type,'JER'])
            self.graphs[name.replace('JER', 'MPF')] = HistToGraph(self.files['yannick'].Get('MPF'))
            self.graphs[name] = HistToGraph(self.files['yannick'].Get('MPFx'))
            self.graphs[name.replace('JER', 'MPFX')] = CalculateMPFX(self.graphs[name.replace('JER', 'MPF')], self.graphs[name], name.replace('JER', 'MPFX'), 0)
    
    def LoadZjet(self):
        self.files['zjet'] = rt.TFile(self.inputPath+'zjet_balance_UL2018_jetpt_nominal_small.root')
        for type in ['Data', 'MC','Ratio']:
            name = ' '.join(['zjet',type,'JER'])
            self.graphs[name] = self.files['zjet'].Get(type.lower())
    
    def LoadMCTruth(self):
        self.funcs={}
        name = 'MC Truth avg-mu JER'
        self.files[name] = rt.TFile(self.inputPath+'JER_MCtruth_avg_mu_UL18.root')
        self.graphs[name] = self.files[name].Get('ak4pfchsl1l2l3/RelResVsJetPt_JetEta0.261to0.522_Mu0to60')
        self.funcs[name] = list(self.graphs[name].GetListOfFunctions())[0]
        self.graphs[name].GetListOfFunctions().Remove(self.funcs[name])
        name = 'MC Truth bin-mu JER'
        self.files[name] = rt.TFile(self.inputPath+'JER_MCtruth_UL18.root')
        for mu in self.mu_bins:
            self.graphs[name+mu] = self.files[name].Get('ak4pfchsl1l2l3/RelResVsJetPt_JetEta0.261to0.522_'+mu)
            self.funcs[name+mu] = list(self.graphs[name+mu].GetListOfFunctions())[0]
            self.graphs[name+mu].GetListOfFunctions().Remove(self.funcs[name+mu])

    def LoadInputs(self):
        self.LoadRC()
        self.LoadMPFX()
        self.LoadDijet()
        self.LoadZjet()
        self.LoadMCTruth()
        for mu in ['Mu0to10','Mu10to20','Mu20to30','Mu30to40','Mu40to50','Mu50to60']:
            for pu in ['high','low']:
                name  = pu+'-PU MC JER'
                # scale = (lambda x: 1) if pu =='low' else (lambda x: (1.05 if x<100 or x>1000 else 1.))
                scale = (lambda x: 1)
                self.graphs[(name+mu).replace('JER', 'Ratio')] = MakeRatioGraphFunc(self.graphs[name], self.funcs['MC Truth bin-mu JER'+mu], scale=scale)
            # scale =lambda x: (1.05 if x<100 or x>1000 else 1.)
            name  = 'dijet SM MC JER'
            self.graphs[(name+mu).replace('JER', 'Ratio')] = MakeRatioGraphFunc(self.graphs[name], self.funcs['MC Truth bin-mu JER'+mu], scale=scale)
            name  = 'dijet FE MC JER'
            self.graphs[(name+mu).replace('JER', 'Ratio')] = MakeRatioGraphFunc(self.graphs[name], self.funcs['MC Truth bin-mu JER'+mu], scale=scale)
            name  = 'dijet CHS SM MC JER'
            self.graphs[(name+mu).replace('JER', 'Ratio')] = MakeRatioGraphFunc(self.graphs[name], self.funcs['MC Truth bin-mu JER'+mu], scale=scale)
            name  = 'dijet CHS FE MC JER'
            self.graphs[(name+mu).replace('JER', 'Ratio')] = MakeRatioGraphFunc(self.graphs[name], self.funcs['MC Truth bin-mu JER'+mu], scale=scale)
            name  = 'zjet MC JER'
            self.graphs[(name+mu).replace('JER', 'Ratio')] = MakeRatioGraphFunc(self.graphs[name], self.funcs['MC Truth bin-mu JER'+mu], scale=scale)

    def CreateCanvas(self, canvName='', zoom=False, nEntries=8,nFunc=2):
        if 'dicanv' in self.__dict__: self.dicanv.Close()
        if 'dicanv_comb' in self.__dict__: self.dicanv_comb.Close()
        if 'canv' in self.__dict__: self.canv.Close()
        XMin, XMax = (15, 4500)
        YMin, YMax = (0.85,1.2) if zoom else (0.85,1.6)
        xName, yName = ('p_{T,jet} [GeV]', 'Jet Energy Resolution')
        RName = 'MC/Truth' if zoom else 'Data/MC'
        self.dicanv = tdrDiCanvas('dicanvas_JER'+self.year+canvName, XMin, XMax, 0.0001,0.45, YMin, YMax, xName, yName, RName)
        self.dicanv.cd(1).SetLogx(True)
        self.dicanv.cd(2).SetLogx(True)
        self.dicanv_comb = tdrDiCanvas('dicanvas_JER_Comb'+self.year+canvName, XMin, XMax, 0.0001,0.45, YMin,YMax, xName, yName, RName)
        self.dicanv_comb.cd(1).SetLogx(True)
        self.dicanv_comb.cd(2).SetLogx(True)

        self.canv = tdrCanvas('JER'+self.year+canvName, XMin, XMax, 0, 1.2, xName, yName, square=kSquare, isExtraSpace=True)
        self.canv.SetLogx(True)
        self.leg = tdrLeg(0.60,0.90-(nEntries+1)*0.045,0.90,0.90)
        self.leg_func = tdrLeg(0.45,0.90-(nFunc+1)*0.045,0.65,0.90)
        FixXAsisPartition(self.canv)
        FixXAsisPartition(self.dicanv.cd(2), shift=0.825 , textsize=0.11, bins=[30,300, 3000])
        FixXAsisPartition(self.dicanv_comb.cd(2), shift=0.825 , textsize=0.11, bins=[30,300, 3000])
        self.lines= {}
        for y in [1,1.1]:
            self.lines[y] = rt.TLine(XMin, y, XMax, y)
            self.dicanv.cd(2)
            tdrDrawLine(self.lines[y], lcolor=rt.kBlack, lstyle=rt.kDashed, lwidth=1)
            self.dicanv_comb.cd(2)
            tdrDrawLine(self.lines[y], lcolor=rt.kBlack, lstyle=rt.kDashed, lwidth=1)
        

    def PlotMPFX(self):
        styles = OrderedDict([
            ('low-PU Data MPF',   {'mcolor':rt.kRed,    'marker':rt.kFullCircle, 'msize':0.8}),
            ('low-PU MC MPF',     {'mcolor':rt.kRed,    'marker':rt.kOpenCircle, 'msize':0.8}),
            ('low-PU Data MPFX',  {'mcolor':rt.kBlue,   'marker':rt.kFullCircle, 'msize':0.8}),
            ('low-PU MC MPFX',    {'mcolor':rt.kBlue,   'marker':rt.kOpenCircle, 'msize':0.8}),
            ('high-PU Data MPF',  {'mcolor':rt.kRed+2,  'marker':rt.kFullSquare, 'msize':0.6}),
            ('high-PU MC MPF',    {'mcolor':rt.kRed+2,  'marker':rt.kOpenSquare, 'msize':0.6}),
            ('high-PU Data MPFX', {'mcolor':rt.kBlue+2, 'marker':rt.kFullSquare, 'msize':0.6}),
            ('high-PU MC MPFX',   {'mcolor':rt.kBlue+2, 'marker':rt.kOpenSquare, 'msize':0.6}),

            # ('low-PU Data JER',   {'mcolor':rt.kOrange+1, 'marker':rt.kFullCircle, 'msize':0.4}),
            # ('low-PU MC JER',     {'mcolor':rt.kOrange+1, 'marker':rt.kOpenCircle, 'msize':0.4}),
            # ('high-PU Data JER',  {'mcolor':rt.kGreen+2,  'marker':rt.kFullSquare, 'msize':0.4}),
            # ('high-PU MC JER',    {'mcolor':rt.kGreen+2,  'marker':rt.kOpenSquare, 'msize':0.4}),
            # ('yannick MC JER',         {'mcolor':rt.kBlack,  'marker':rt.kFullCross,       'msize':0.8}),
            ('yannick MC MPF',     {'mcolor':rt.kGray,   'marker':rt.kFullCross,       'msize':0.8}),
            ('yannick MC MPFX',    {'mcolor':rt.kGray+1, 'marker':rt.kFullCross,       'msize':0.8}),
        ])
        self.CreateCanvas()
        self.canv.cd()
        for name,style in styles.items():
            hist = self.graphs[name]
            tdrDraw(hist, 'Pz', **style)
            self.leg.AddEntry(hist, name, 'PLE')
        fixOverlay()
        self.canv.SaveAs(os.path.join(self.outputPath, 'MPFX.pdf'))
        self.canv.Close()
    
    def CleanGraphs(self):
        self.pdfextraname = '_cleaned'+self.pdfextraname
        ranges = {
            'low-PU':       (30, 2000),
            'high-PU':      (30, 2000),
            'zjet':         (30, 2000),
            'dijet SM':     (30, 2000),
            'dijet FE':     (30, 2000),
            'dijet CHS SM': (30, 2000),
            'dijet CHS FE': (30, 2000),
            }
        for type in self.all_types:
            combined_name = ' '.join(['combined',type,'JER'])
            combined_mpfx_name = ' '.join(['combined MPFX',type,'JER'])
            self.graphs[combined_name]= rt.TGraphErrors()
            self.graphs[combined_mpfx_name]= rt.TGraphErrors()
            for sample, (min_,max_) in ranges.items():
                name = ' '.join([sample,type, 'JER'])
                self.graphs[name] = TruncateGraph(self.graphs[name],min_,max_)
                if any([x in sample for x in ['low-PU', 'CHS']]): continue
                if any([x in sample for x in ['high-PU']]):
                    self.graphs[combined_mpfx_name] = MergeGraphs(self.graphs[combined_mpfx_name], self.graphs[name])
                else:
                    self.graphs[combined_name] = MergeGraphs(self.graphs[combined_name], self.graphs[name])
    
    def FitJER(self):
        fit_k = False
        fix_N = False
        fit_min, fit_max = 20, 3000
        initial_pars = [2,0.5,0.05,-1,1,1,1]
        functional_form_MC = "TMath::Sqrt([0]*[0]/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])"
        functional_form_Data = "TMath::Sqrt([0]*[0]*[4]*[4]/(x*x)+[1]*[1]*[5]*[5]*pow(x,[3])+[2]*[2]*[6]*[6])"
        functional_form_Ratio, npars_ratio = functional_form_Data+"/"+functional_form_MC, 7
        if not fit_k:
            functional_form_Ratio, npars_ratio = "TMath::Sqrt([4]*[4]/(x*x)+[5]*[5]*pow(x,[7])+[6]*[6])/TMath::Sqrt([0]*[0]/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])", 8
        
        for mode in ['combined', 'combined MPFX']:
            for type in ['MC','Data']:
                func_name = ' '.join([mode, 'func', type,'JER'])
                if fit_k and type=='Data':
                    functional_form, npars = functional_form_Data, 7
                else:
                    functional_form, npars = functional_form_MC, 4
                self.funcs[func_name] = rt.TF1(func_name, functional_form,fit_min, fit_max, npars)
                func = self.funcs[func_name]
                for i in range(func.GetNpar()):
                    func.SetParameter(i,initial_pars[i])
                if fix_N and 'MC' in type:
                    func.SetParLimits(0,5,100)

                graph_name = ' '.join([mode, type ,'JER'])
                if 'Data' in type:
                    func_MC = self.funcs[func_name.replace('Data','MC')]
                    for i in range(func_MC.GetNpar()):
                        if not fit_k: continue
                        func.FixParameter(i,func_MC.GetParameter(i))
                self.graphs[graph_name].Fit(func, 'Q+', '', 50, 2000)
                self.graphs[graph_name].GetListOfFunctions().Remove(list(self.graphs[graph_name].GetListOfFunctions())[0])
        
            func_name = ' '.join([mode, 'func', 'Ratio','JER'])
            self.funcs[func_name] = rt.TF1(func_name, functional_form_Ratio,fit_min, fit_max,npars_ratio)
            func = self.funcs[func_name]
            if fit_k:
                func_Data = self.funcs[func_name.replace('Ratio','Data')]
                for i in range(npars_ratio):
                    # print(mode, i,func_Data.GetParameter(i))
                    func.FixParameter(i,func_Data.GetParameter(i))
            else:
                func_MC = self.funcs[func_name.replace('Ratio','MC')]
                func_Data = self.funcs[func_name.replace('Ratio','Data')]
                for i in range(func_MC.GetNpar()):
                    # print(mode, i,func_MC.GetParameter(i))
                    func.FixParameter(i,func_MC.GetParameter(i))
                for i in range(func_Data.GetNpar()):
                    # print(mode, i+func_MC.GetNpar(),func_Data.GetParameter(i))
                    func.FixParameter(i+func_MC.GetNpar(),func_Data.GetParameter(i))
        
    def Plot(self, to_plot, pdfname, nEntries=8,nFunc=2):
        self.CreateCanvas(zoom=any(['RatioMu' in x for x in to_plot]), nEntries=nEntries,nFunc=nFunc)
        func_styles = OrderedDict([
            ('combined func Data JER',       {'label':'Data',  'lcolor':rt.kAzure+2, 'lstyle':rt.kSolid}),
            ('combined func MC JER',         {'label':'MC',    'lcolor':rt.kRed+1,   'lstyle':rt.kSolid}),
            ('combined func Ratio JER',      {'label':'Ratio', 'lcolor':rt.kBlack,   'lstyle':rt.kSolid}),
            ('combined MPFX func Data JER',  {'label':'Data',  'lcolor':rt.kAzure+2, 'lstyle':rt.kSolid}),
            ('combined MPFX func MC JER',    {'label':'MC',    'lcolor':rt.kRed+1,   'lstyle':rt.kSolid}),
            ('combined MPFX func Ratio JER', {'label':'Ratio', 'lcolor':rt.kBlack,   'lstyle':rt.kSolid}),
        ])
        styles = OrderedDict([
            ('low-PU Data JER',               {'mcolor':rt.kRed+1,     'marker':rt.kFullCircle, 'msize':0.8}),
            ('low-PU MC JER',                 {'mcolor':rt.kRed+1,     'marker':rt.kOpenCircle, 'msize':0.8}),
            ('low-PU Ratio JER',              {'mcolor':rt.kRed+1,     'marker':rt.kFullCircle, 'msize':0.8}),
            ('low-PU MC RatioMu0to10',        {'mcolor':rt.kRed+1,     'marker':rt.kFullCircle, 'msize':0.9}),
            ('low-PU MC RatioMu10to20',       {'mcolor':rt.kRed+1,     'marker':rt.kFullCircle, 'msize':0.9}),

            ('high-PU Data JER',              {'mcolor':rt.kMagenta+2, 'marker':rt.kFullSquare, 'msize':0.8}),
            ('high-PU MC JER',                {'mcolor':rt.kMagenta+2, 'marker':rt.kOpenSquare, 'msize':0.8}),
            ('high-PU Ratio JER',             {'mcolor':rt.kMagenta+2, 'marker':rt.kFullSquare, 'msize':0.8}),
            ('high-PU MC RatioMu50to60',      {'mcolor':rt.kMagenta+2, 'marker':rt.kFullSquare, 'msize':0.9}),
            ('high-PU MC RatioMu30to40',      {'mcolor':rt.kMagenta+2, 'marker':rt.kFullSquare, 'msize':0.9}),

            ('dijet SM Data JER',             {'mcolor':rt.kOrange+1,  'marker':rt.kFullTriangleDown, 'msize':0.8}),
            ('dijet SM MC JER',               {'mcolor':rt.kOrange+1,  'marker':rt.kOpenTriangleDown, 'msize':0.8}),
            ('dijet SM Ratio JER',            {'mcolor':rt.kOrange+1,  'marker':rt.kFullTriangleDown, 'msize':0.8}),
            ('dijet SM MC RatioMu50to60',     {'mcolor':rt.kOrange+1,  'marker':rt.kFullTriangleDown, 'msize':0.9}),
            ('dijet SM MC RatioMu30to40',     {'mcolor':rt.kOrange+1,  'marker':rt.kFullTriangleDown, 'msize':0.9}),

            ('dijet FE Data JER',             {'mcolor':rt.kGreen+2,   'marker':rt.kFullTriangleUp,  'msize':0.8}),
            ('dijet FE MC JER',               {'mcolor':rt.kGreen+2,   'marker':rt.kOpenTriangleUp,  'msize':0.8}),
            ('dijet FE Ratio JER',            {'mcolor':rt.kGreen+2,   'marker':rt.kFullTriangleUp,  'msize':0.8}),
            ('dijet FE MC RatioMu50to60',     {'mcolor':rt.kGreen+2,   'marker':rt.kFullTriangleUp,  'msize':0.9}),
            ('dijet FE MC RatioMu30to40',     {'mcolor':rt.kGreen+2,   'marker':rt.kFullTriangleUp,  'msize':0.9}),

            ('dijet CHS SM Data JER',         {'mcolor':rt.kOrange-1,  'marker':rt.kFullCross,       'msize':0.8}),
            ('dijet CHS SM MC JER',           {'mcolor':rt.kOrange-1,  'marker':rt.kOpenCross,       'msize':0.8}),
            ('dijet CHS SM Ratio JER',        {'mcolor':rt.kOrange-1,  'marker':rt.kFullCross,       'msize':0.8}),
            ('dijet CHS SM MC RatioMu50to60', {'mcolor':rt.kOrange-1,  'marker':rt.kOpenCross,       'msize':0.8}),
            ('dijet CHS SM MC RatioMu30to40', {'mcolor':rt.kOrange-1,  'marker':rt.kOpenCross,       'msize':0.8}),
            
            ('dijet CHS FE Data JER',         {'mcolor':rt.kGreen+4,   'marker':rt.kFullDiamond,     'msize':0.8}),
            ('dijet CHS FE MC JER',           {'mcolor':rt.kGreen+4,   'marker':rt.kOpenDiamond,     'msize':0.8}),
            ('dijet CHS FE Ratio JER',        {'mcolor':rt.kGreen+4,   'marker':rt.kFullDiamond,     'msize':0.8}),
            ('dijet CHS FE MC RatioMu50to60', {'mcolor':rt.kGreen+4,   'marker':rt.kOpenDiamond,     'msize':0.8}),
            ('dijet CHS FE MC RatioMu30to40', {'mcolor':rt.kGreen+4,   'marker':rt.kOpenDiamond,     'msize':0.8}),

            ('zjet Data JER',                 {'mcolor':rt.kAzure+2,   'marker':rt.kFullTriangleUp,  'msize':0.8}),
            ('zjet MC JER',                   {'mcolor':rt.kAzure+2,   'marker':rt.kOpenTriangleUp,  'msize':0.8}),
            ('zjet Ratio JER',                {'mcolor':rt.kAzure+2,   'marker':rt.kFullTriangleUp,  'msize':0.8}),
            ('zjet MC RatioMu50to60',         {'mcolor':rt.kAzure+2,   'marker':rt.kFullTriangleUp,  'msize':0.9}),
            ('zjet MC RatioMu30to40',         {'mcolor':rt.kAzure+2,   'marker':rt.kFullTriangleUp,  'msize':0.9}),

            ('combined Data JER',             {'mcolor':rt.kAzure+2,   'marker':rt.kFullCircle,      'msize':0.8, 'label':'Data'  }),
            ('combined MC JER',               {'mcolor':rt.kRed+1,     'marker':rt.kOpenCircle,      'msize':0.8, 'label':'MC'    }),
            ('combined Ratio JER',            {'mcolor':rt.kBlack,     'marker':rt.kFullCircle,      'msize':0.8, 'label':'Ratio' }),

            ('combined MPFX Data JER',        {'mcolor':rt.kAzure+2,   'marker':rt.kFullCircle,      'msize':0.8, 'label':'Data'  }),
            ('combined MPFX MC JER',          {'mcolor':rt.kRed+1,     'marker':rt.kOpenCircle,      'msize':0.8, 'label':'MC'    }),
            ('combined MPFX Ratio JER',       {'mcolor':rt.kBlack,     'marker':rt.kFullCircle,      'msize':0.8, 'label':'Ratio' }),

            ('yannick MC JER',         {'mcolor':rt.kBlack,  'marker':rt.kFullCross,       'msize':0.8}),
            ('yannick MC MPF',         {'mcolor':rt.kGray,   'marker':rt.kFullCross,       'msize':0.8}),
            ('yannick MC MPFX',        {'mcolor':rt.kGray+1, 'marker':rt.kFullCross,       'msize':0.8}),
        ])

        # to_plot.append('yannick MC JER')
        # to_plot.append('yannick MC MPF')
        # to_plot.append('yannick MC MPFX')
        self.leg.Clear()
        self.leg_func.Clear()
        self.dicanv.cd(1)
        for mu in self.mu_bins:
            name = 'MC Truth bin-mu JER'+mu
            bin = mu.replace('Mu','[').replace('to',',')+']'
            lstyle = rt.kSolid
            if mu== 'Mu30to40': lstyle= rt.kDotted
            if mu== 'Mu50to60': lstyle= rt.kDashed
            if not name in to_plot: continue
            tdrDrawLine(self.funcs[name], lcolor=rt.kViolet+7, lstyle=lstyle)
            self.leg.AddEntry(self.funcs[name], 'MC Truth  #mu='+bin, 'l')
        self.dicanv.cd(1)
        self.leg.Draw('same')
        for name,style in func_styles.items():
            if not name in to_plot: continue
            if not name in self.funcs: continue
            if 'Ratio' in name: self.dicanv.cd(2)
            else: self.dicanv.cd(1)
            legName = style.get('label', name)
            style.pop('label', None)
            tdrDrawLine(self.funcs[name], **style)
            if not 'Ratio' in name:
                self.leg_func.AddEntry(self.funcs[name], legName, 'l')
        self.dicanv.cd(1)
        self.leg_func.Draw('same')
        for name,style in styles.items():
            if not name in to_plot: continue
            if not name in self.graphs: continue
            if 'Ratio' in name: self.dicanv.cd(2)
            else: self.dicanv.cd(1)
            legName = style.get('label', name)
            style.pop('label', None)
            graph = self.graphs[name]
            tdrDraw(graph, 'Pz', **style)
            if not 'Ratio' in name:
                self.leg.AddEntry(graph, legName, 'PLE')
        self.dicanv.SaveAs(os.path.join(self.outputPath, pdfname+self.pdfextraname+'.pdf'))

    def PlotCombination(self):
        to_plot = []
        for type in self.all_types:
            to_plot.append(' '.join(['combined',type,'JER']))
            to_plot.append(' '.join(['combined','func', type,'JER']))
        self.Plot(to_plot=to_plot, pdfname='JER_Combination', nEntries=2,nFunc=2)

        to_plot = []
        for type in self.all_types:
            to_plot.append(' '.join(['combined MPFX',type,'JER']))
            to_plot.append(' '.join(['combined MPFX','func', type,'JER']))
        self.Plot(to_plot=to_plot, pdfname='JER_Combination_MPFX', nEntries=2,nFunc=2)
    
    def PlotAll(self):
        self.PlotMPFX()
        to_plot = []
        for mu in self.mu_bins_short:
            to_plot.append('MC Truth bin-mu JER'+mu)
        self.Plot(to_plot=to_plot, pdfname='JER_MC_Truth')

        for sample in self.mpfx_samples+self.balance_samples:
            for type in ['MC']:
                to_plot.append(' '.join([sample,type,'JER']))
        self.Plot(to_plot=to_plot, pdfname='JER_MC')

        for sample in self.mpfx_samples+self.balance_samples:
            mu = 'RatioMu50to60' if not 'low-PU' in sample else 'RatioMu0to10'
            to_plot.append(' '.join([sample,'MC', mu]))
        self.Plot(to_plot=to_plot, pdfname='JER_MC_ratio')

        to_plot = list(filter(lambda x: not 'RatioMu' in x, to_plot))
        for sample in self.mpfx_samples+self.balance_samples:
            mu = 'RatioMu30to40' if not 'low-PU' in sample else 'RatioMu10to20'
            to_plot.append(' '.join([sample,'MC', mu]))
        self.Plot(to_plot=to_plot, pdfname='JER_MC_ratio2')

        to_plot = []
        for sample in self.balance_samples:
            for type in self.all_types:
                to_plot.append(' '.join([sample,type,'JER']))
        self.Plot(to_plot=to_plot, pdfname='JER_Data_noMPFX')

        to_plot = []
        for sample in self.mpfx_samples:
            for type in self.all_types:
                to_plot.append(' '.join([sample,type,'JER']))
        self.Plot(to_plot=to_plot, pdfname='JER_Data_MPFX')

        to_plot = []
        for sample in self.mpfx_samples+self.balance_samples:
            for type in self.all_types:
                to_plot.append(' '.join([sample,type,'JER']))
        self.Plot(to_plot=to_plot, pdfname='JER_Data')

    def Store(self):
        self.files['output'] = rt.TFile(self.outputPath+'output_CombineJER.root', 'RECREATE')
        self.files['output'].cd()
        for pu in self.mpfx_samples:
            for type in self.all_types:
                for mode in ['JER']:
                    name = ' '.join([pu,type,mode])
                    self.graphs[name].Write() 

    def Close(self):
        for f_ in self.files.values():
            f_.Close()


def main():
    CJ = CombineJER()
    CJ.LoadInputs()
    CJ.PlotAll()
    CJ.CleanGraphs()
    CJ.PlotAll()
    CJ.FitJER()
    CJ.PlotCombination()
    CJ.Store()
    CJ.Close()

if __name__ == '__main__':
    main()
