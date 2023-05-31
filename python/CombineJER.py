from utils import *

def AddRC(graph, name, RC, RC_err):
    gnew = graph.Clone(name)
    xpos = list(gnew.GetX())
    vals = list(graph.GetY())
    for bin in range(0,graph.GetN()):
        x  = xpos[bin]
        if (x==0): continue
        err = oplus(graph.GetErrorY(bin), RC_err/x)
        gnew.SetPoint(bin, x, oplus(vals[bin], RC/x))
        gnew.SetPointError(bin, graph.GetErrorXlow(bin), graph.GetErrorXhigh(bin), err, err)
    return gnew


class CombineJER():
    def __init__(self, year='UL18', algo='AK4 CHS', etamin = '0.0', etamax = '0.261'):
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

        self.mpfx_samples = ['MPFX RC', 'MPFX']
        self.mpfx_samples = ['MPFX']
        # self.mpfx_samples = []
        # self.balance_samples = ['dijet FE', 'dijet SM', 'dijet v19 FE', 'dijet v19 SM', 'zjet']
        # self.balance_samples = ['dijet FE', 'dijet SM', 'dijet v19 FE', 'dijet v19 SM']
        # self.balance_samples = ['dijet SM', 'dijet v19 SM']
        self.balance_samples = ['dijet FE', 'dijet SM']
        # self.balance_samples = ['zjet']
        # self.mc_truth_modes = ['CI', 'Gauss']
        self.mc_truth_modes = ['CI']
        self.combinations = {'mpfx': self.mpfx_samples, 'balance': self.balance_samples, 'all': self.balance_samples+self.mpfx_samples+['MC Truth CI']}

        self.mc_truth_fits = ['', ' N fix', ' P fix', ' NP fix']
        # self.mc_truth_fits = ['']
        self.all_types = ['Data', 'MC', 'Ratio']
        self.pdfextraname = ''
        self.functional_form_MC = "TMath::Sqrt([0]*[0]/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])"
        self.func_min, self.func_max = 10, 4500
        self.fit_min, self.fit_max = 50, 1000
        self.fit_k = True # Fix same exponent for S term for Data and MC
        self.FixD = False
        self.FixN = True # Fix N term to RC values for Data and MC
        self.FixMC = False
        self.setks_max  = True # Fix max k factor for the S-Term
        self.setkc_max  = False # Fix max k factor for the C-Term
        self.SetStyle()

    def LoadMPFX(self):
        modes = {'MPFX RC': 'h2jer', 'MPFX': 'h2mpfx'}
        for mode, hname in modes.items():
            if not mode in self.mpfx_samples: continue
            for type in ['Data', 'MC']:
                name = f'{mode} {type}'
                self.files[name] = rt.TFile(os.path.join(self.inputPath,f'jmenano_{type.lower()}_cmb_Run2_v26.root'))
                h2d = self.files[name].Get('Dijet2/'+hname)
                self.graphs[name] = Hist2DToGraphAsymm(h2d, xmin=self.etamin, xmax=self.etamax)
                if not 'RC' in mode:
                    graph = self.graphs[f'Noise Term {type}']
                    RC = self.funcs[f'Noise Term {type}'].Eval(30)
                    RC_err = graph.GetErrorY(int(graph.GetN()/2))
                    self.graphs[name] = AddRC(self.graphs[name], name, RC, RC_err)
            name = f'{mode} Ratio'
            self.graphs[name] = MakeRatioGraphs(self.graphs[name.replace('Ratio', 'Data')], self.graphs[name.replace('Ratio', 'MC')], name)
    
    def LoadRC(self, mode = 'nominal'):
        #mode = 'sigrc'
        self.files['RC'] = rt.TFile(self.inputPath+'RC_noise_UL18.root')
        for type in ['MC', 'Data']:
            name = f'Noise Term {type}'
            h2d = self.files['RC'].Get(f'rc_noiseterm_vs_npu_jer_{type}_{mode}')
            self.graphs[name] = Hist2DToGraphAsymm(h2d, xmin=self.etamin, xmax=self.etamax, yerr=0.2)
            self.funcs[name] = rt.TF1(name, 'TMath::Sqrt([0]*[0]+[1]*[1]*x)', 0, 45 if 'UL16' in self.year else 65)
            self.graphs[name].Fit(self.funcs[name],'RQMS')
            RemoveFunc(self.graphs[name])
        name = 'Noise Term Ratio'
        self.graphs[name] = MakeRatioGraphs(self.graphs[name.replace('Ratio', 'Data')], self.graphs[name.replace('Ratio', 'MC')], name)
        self.funcs[name+'pol0'] = rt.TF1(name+'pol0', 'pol0', 0, 45 if 'UL16' in self.year else 65)
        self.graphs[name].Fit(self.funcs[name+'pol0'],'RQMS')
        RemoveFunc(self.graphs[name])
        self.funcs[name] = rt.TF1(name, 'TMath::Sqrt([0]*[0]+[1]*[1]*x)/TMath::Sqrt([2]*[2]+[3]*[3]*x)', 0, 45 if 'UL16' in self.year else 65, 4)
        for par in range(0,2):
            self.funcs[name].SetParameter(par,   self.funcs[name.replace('Ratio', 'Data')].GetParameter(par))
            self.funcs[name].SetParameter(par+1, self.funcs[name.replace('Ratio', 'MC')].GetParameter(par))
    
    def LoadDijet(self):
        etamin = self.etamin.replace('.','p').replace('0p0', '0p00000')
        etamax = self.etamax.replace('.','p')
        vers = {'dijet v19': 'dijet_balance_UL18',  'dijet': 'dijet_balance_UL18_Summer20UL18_V2_AK4CHS'}
        for ver, fname in vers.items():
            self.files[ver] = rt.TFile(self.inputPath+fname+'.root')
            for mode in ['SM', 'FE']:
                for type in ['Data', 'MC']:
                    name = f'{ver} {mode} {type}'
                    self.graphs[name] = self.files[ver].Get(f'dijet_balance_jer_{type}_{etamin}_{etamax}_{mode}_nominal')
                    if type=='MC':
                        self.graphs[name.replace('MC', 'Ratio')] = MakeRatioGraphs(self.graphs[name.replace('MC', 'Data')], self.graphs[name], name.replace('MC', 'Ratio'))

    def LoadZjet(self):
        etamin = self.etamin.replace('.','p').replace('0p0', '0p000')
        etamax = self.etamax.replace('.','p')
        modes = {'zjet': 'zjet_balance_UL2018_jetpt_nominal'}
        for mode, fname in modes.items():
            self.files[mode] = rt.TFile(self.inputPath+fname+'.root')
            for type in ['Data', 'MC','Ratio']:
                name = f'{mode} {type}'
                if 'Ratio' in name:
                    self.graphs[name] = MakeRatioGraphs(self.graphs[name.replace('Ratio', 'Data')], self.graphs[name.replace('Ratio', 'MC')], name)
                else:
                    self.graphs[name] = self.files[mode].Get(f'zjet_balance_jer_{type}_{etamin}to{etamax}_nominal')
    
    def LoadMCTruth(self):
        modes = {'CI': 'JER_MCtruth_avg_mu_UL18_CI', 'Gauss': 'JER_MCtruth_avg_mu_UL18_Gauss', 'RMS': 'JER_MCtruth_avg_mu_UL18'}
        for mode, fname in modes.items():
            if not mode in self.mc_truth_modes: continue
            name = f'MC Truth {mode}'
            self.files[name] = rt.TFile(self.inputPath+fname+'.root')
            self.graphs[name] = self.files[name].Get('ak4pfchsl1l2l3/RelResVsJetPt_JetEta0.261to0.522_Mu0to60')
            self.funcs[name+'default'] = list(self.graphs[name].GetListOfFunctions())[0]
            RemoveFunc(self.graphs[name])
            self.funcs[name] = rt.TF1(name, self.functional_form_MC, self.func_min, self.func_max, 4)
            if self.FixN:
                self.funcs[name].FixParameter(0,self.funcs['Noise Term MC'].Eval(30))
            if self.FixD:
                self.funcs[name].FixParameter(3,-1)
            self.graphs[name].Fit(self.funcs[name], 'Q+', '', self.fit_min, self.fit_max)
            RemoveFunc(self.graphs[name])
            PrintFuncParameters(func = self.funcs[name], color=orange)
            self.graphs[name+' Ratio'] = MakeRatioGraphFunc(self.graphs[name], self.funcs[name])

            fname = name+ ' N fix'
            self.funcs[fname] = rt.TF1(fname, self.functional_form_MC, self.func_min, self.func_max, 4)
            self.funcs[fname].FixParameter(0,self.funcs['Noise Term MC'].Eval(30))
            self.graphs[name].Fit(self.funcs[fname], 'Q+', '', 50, 1000)
            self.graphs[fname+' Ratio'] = MakeRatioGraphFunc(self.graphs[name], self.funcs[fname])
            RemoveFunc(self.graphs[name])

            fname = name+ ' P fix'
            self.funcs[fname] = rt.TF1(fname, self.functional_form_MC, self.func_min, self.func_max, 4)
            self.funcs[fname].FixParameter(3,-1)
            self.graphs[name].Fit(self.funcs[fname], 'Q+', '', 50, 1000)
            self.graphs[fname+' Ratio'] = MakeRatioGraphFunc(self.graphs[name], self.funcs[fname])
            RemoveFunc(self.graphs[name])

            fname = name+ ' NP fix'
            self.funcs[fname] = rt.TF1(fname, self.functional_form_MC, self.func_min, self.func_max, 4)
            self.funcs[fname].FixParameter(0,self.funcs['Noise Term MC'].Eval(30))
            self.funcs[fname].FixParameter(3,-1)
            self.graphs[name].Fit(self.funcs[fname], 'Q+', '', 50, 1000)
            self.graphs[fname+' Ratio'] = MakeRatioGraphFunc(self.graphs[name], self.funcs[fname])
            RemoveFunc(self.graphs[name])

    def LoadInputs(self):
        self.files = OrderedDict()
        self.graphs = OrderedDict()
        self.funcs = OrderedDict()
        self.LoadRC()
        self.LoadMCTruth()
        self.LoadMPFX()
        self.LoadDijet()
        self.LoadZjet()
        for name in self.mpfx_samples+self.balance_samples:
            for mode in self.mc_truth_modes:
                num = f'{name} MC'
                den = f'MC Truth {mode}'
                ratio = f'{name} MC Ratio {mode}'
                self.graphs[ratio] = MakeRatioGraphFunc(self.graphs[num], self.funcs[den])

    def CreateCanvas(self, canvName='', zoom=True, nEntries=8,nFunc=2):
        if 'leg' in self.__dict__: self.leg.Clear()
        if 'leg_func' in self.__dict__: self.leg_func.Clear()
        if 'dicanv' in self.__dict__: self.dicanv.Close()
        if 'dicanv_comb' in self.__dict__: self.dicanv_comb.Close()
        if 'canv' in self.__dict__: self.canv.Close()
        XMin, XMax = (15, 4500)
        YMin, YMax, shift = (0.85,1.2, 0.835) if zoom else (0.95,1.3, 0.935) #(0.55,1.6, 0.505)
        xName, yName = ('p_{T,jet} [GeV]', 'Jet Energy Resolution')
        RName = 'MC/Truth' if 'MC' in canvName else 'Data/MC'
        self.dicanv = tdrDiCanvas('dicanvas_JER'+self.year+canvName, XMin, XMax, 0.0001,0.45, YMin, YMax, xName, yName, RName)
        self.dicanv.cd(1).SetLogx(True)
        self.dicanv.cd(2).SetLogx(True)
        self.dicanv_comb = tdrDiCanvas('dicanvas_JER_Comb'+self.year+canvName, XMin, XMax, 0.0001,0.45, YMin,YMax, xName, yName, RName)
        self.dicanv_comb.cd(1).SetLogx(True)
        self.dicanv_comb.cd(2).SetLogx(True)

        self.canv = tdrCanvas('JER'+self.year+canvName, XMin, XMax, 0, 1.2, xName, yName, square=kSquare, isExtraSpace=True)
        self.canv.SetLogx(True)
        self.leg = tdrLeg(0.65,0.90-(nEntries+1)*0.045,0.90,0.90)
        self.leg_func = tdrLeg(0.40,0.90-(nFunc+1)*0.045,0.65,0.90)
        if 'JER_comb' in canvName:
            self._JER_DATA = rt.TLine()
            self._JER_MC = rt.TLine()
            tdrDrawLine(self._JER_DATA, lcolor=rt.kGray, lstyle=rt.kSolid)
            tdrDrawLine(self._JER_MC, lcolor=rt.kGray, lstyle=rt.kDashed)
            self.leg_func.AddEntry(self._JER_DATA, 'Data', 'l')
            self.leg_func.AddEntry(self._JER_MC, 'MC', 'l')
        FixXAsisPartition(self.canv)
        FixXAsisPartition(self.dicanv.cd(2), shift=shift , textsize=0.11, bins=[30, 300, 3000])
        FixXAsisPartition(self.dicanv_comb.cd(2), shift=shift , textsize=0.11, bins=[30,300, 3000])
        self.lines= {}
        for y in [1,1.1]:
            self.lines[y] = rt.TLine(XMin, y, XMax, y)
            self.dicanv.cd(2)
            tdrDrawLine(self.lines[y], lcolor=rt.kBlack, lstyle=rt.kDashed, lwidth=1)
            self.dicanv_comb.cd(2)
            tdrDrawLine(self.lines[y], lcolor=rt.kBlack, lstyle=rt.kDashed, lwidth=1)
        
    def PlotRC(self):
        if 'dicanv' in self.__dict__: self.dicanv.Close()
        XMin, XMax = -3, 63
        YMin, YMax = 0.1, 5
        self.dicanv = tdrDiCanvas('dicanvas_NoiseTerm'+self.year, XMin, XMax, YMin, YMax, 0.90, 1.2, 'N_{PU}', 'Noise Term', 'Data/MC')
        for name, graph in self.graphs.items():
            if not 'Noise Term' in name: continue
            color, marker = rt.kOrange+1, rt.kFullTriangleDown
            if 'MC' in name:
                color, marker = rt.kAzure+2, rt.kFullTriangleUp
            self.dicanv.cd(2 if 'Ratio' in name else 1)
            tdrDraw(graph, 'Pz', mcolor=color, marker=marker)
            if 'Ratio' in name:
                tdrDrawLine(self.funcs[name+'pol0'], lcolor=color+1)
            tdrDrawLine(self.funcs[name], lcolor=color)
        self.dicanv.SaveAs(os.path.join(self.outputPath, 'NoiseTerm.pdf'))
    
    def CleanGraphs(self):
        self.pdfextraname = '_cleaned'+self.pdfextraname
        ranges = {
            'MPFX':         (40, 2000),
            'MPFX RC':      (40, 2000),
            'zjet':         (30, 2000),
            'dijet SM':     (300, 2000),
            'dijet FE':     (300, 2000),
            'dijet v19 SM': (80, 2000),
            'dijet v19 FE': (80, 2000),
            }
        for name in self.graphs.keys():
            for sample, (min_,max_) in ranges.items():
                if sample in name:
                    self.graphs[name] = TruncateGraph(self.graphs[name],min_,max_)
    
    def CombineDatasets(self):
        for type in self.all_types:
            for combination, samples in self.combinations.items():
                combined_name = f'combined {combination} {type}'
                self.graphs[combined_name]= rt.TGraphErrors()
                print(blue(f'Combining samples for {combination} {type}'))
                for sample in samples:
                    name = f'{sample} {type}'
                    if 'Truth' in sample:
                        if 'MC' in type:
                            name = f'{sample}'
                        else: continue
                    if 'RC' in sample: continue
                    if 'SM' in sample and any(['SM' in x for x in self.balance_samples])and any(['FE' in x for x in self.balance_samples]): continue
                    # if any([x in sample for x in ['RC']]): continue
                    print(cyan(f'  --> {sample}'))
                    self.graphs[combined_name] = MergeGraphs(self.graphs[combined_name], self.graphs[name], scale_err=0.1 if 'mpfx' in sample.lower() else None)
    
    def FitJER(self):
        initial_pars = [3.3,0.9,0.03,-1, 1.1,1.1,1.1]
        functional_form_Data = "TMath::Sqrt([0]*[0]*[4]*[4]/(x*x)+[1]*[1]*[5]*[5]*pow(x,[3])+[2]*[2]*[6]*[6])"
        functional_form_Ratio, npars_ratio = functional_form_Data+"/"+self.functional_form_MC, 7
        if not self.fit_k:
            functional_form_Ratio, npars_ratio = "TMath::Sqrt([4]*[4]/(x*x)+[5]*[5]*pow(x,[3]*[7])+[6]*[6])/TMath::Sqrt([0]*[0]/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])", 8
        
        # for mode in ['combined', 'combined MPFX']:
        for mode in [f'combined {x}' for x in self.combinations.keys()]:
            for type in ['MC','Data']:
                func_name = f'{mode} func {type}'
                if self.fit_k and type=='Data':
                    functional_form, npars = functional_form_Data, 7
                else:
                    functional_form, npars = self.functional_form_MC, 4
                self.funcs[func_name] = rt.TF1(func_name, functional_form, self.func_min, self.func_max, npars)
                func = self.funcs[func_name]
                for i in range(func.GetNpar()):
                    func.SetParameter(i,initial_pars[i])
                if self.FixD:
                    func.FixParameter(3,-1)
                    func.FixParameter(7, 1)

                graph_name = f'{mode} {type}'
                if 'Data' in type:
                    func_MC = self.funcs[func_name.replace('Data','MC')]
                    for i in range(func_MC.GetNpar()):
                        if self.fit_k:
                            func.FixParameter(i,func_MC.GetParameter(i))
                            # if i>=4:
                            #     func.SetParLimits(i,1,10)
                        # else:
                        #     func.SetParameter(i,func_MC.GetParameter(i))
                        #     func.SetParLimits(i,func_MC.GetParameter(i),func_MC.GetParameter(i)*10)
                    n = 1 if self.fit_k else func_MC.GetParameter(2)
                    s = 1 if self.fit_k else func_MC.GetParameter(2)
                    c = 1 if self.fit_k else func_MC.GetParameter(2)
                    if self.setks_max:
                        func.SetParLimits(5 if self.fit_k else 1, 1, 100*s)
                    if self.setkc_max:
                        func.SetParameter(6 if self.fit_k else 7,    1.2*c)
                        func.SetParLimits(6 if self.fit_k else 7, 1, 1.5*c)
                if self.FixN:
                    if self.fit_k and 'Data' in type:
                        func.FixParameter(0,self.funcs[f'Noise Term MC'].Eval(30))
                        func.FixParameter(4,self.funcs[f'Noise Term {type}'].Eval(30)/func.GetParameter(0))
                    else:
                        func.FixParameter(0,self.funcs[f'Noise Term {type}'].Eval(30))
                if self.FixMC and 'MC' in type:
                    func.FixParameter(0,self.funcs[f'MC Truth CI'].GetParameter(0))
                    func.FixParameter(1,self.funcs[f'MC Truth CI'].GetParameter(1))
                    func.FixParameter(2,self.funcs[f'MC Truth CI'].GetParameter(2))
                    func.FixParameter(3,self.funcs[f'MC Truth CI'].GetParameter(3))
                print(blue(f'Fitting {type} with {functional_form}'))
                PrintFuncParameters(func=func, info='  --> Pre-fit:', color=orange)
                self.graphs[graph_name].Fit(func, 'Q+', '', self.fit_min, self.fit_max)
                self.graphs[graph_name+' Ratio'] = MakeRatioGraphFunc(self.graphs[graph_name], func)
                PrintFuncParameters(func=func, info='  --> Post-fit:', color=green)
                print(yellow(f'Chi^{2}/ndf = {func.GetChisquare()/func.GetNDF():.3f}'))
                RemoveFunc(self.graphs[graph_name])                
                self.funcs[f'{mode} N-Term {type}'] = rt.TF1(f'N-Term {type}', "TMath::Sqrt([0]*[0]/(x*x))", self.func_min, self.func_max, 1)
                self.funcs[f'{mode} S-Term {type}'] = rt.TF1(f'S-Term {type}', "TMath::Sqrt([0]*[0]*pow(x,[1]))", self.func_min, self.func_max, 2)
                self.funcs[f'{mode} C-Term {type}'] = rt.TF1(f'C-Term {type}', "TMath::Sqrt([0]*[0])", self.func_min, self.func_max, 1)
                if self.fit_k and 'Data' in type:
                    self.funcs[f'{mode} N-Term {type}'].FixParameter(0, func.GetParameter(0)*func.GetParameter(4))
                    self.funcs[f'{mode} S-Term {type}'].FixParameter(0, func.GetParameter(1)*func.GetParameter(5))
                    self.funcs[f'{mode} S-Term {type}'].FixParameter(1, func.GetParameter(3))
                    self.funcs[f'{mode} C-Term {type}'].FixParameter(0, func.GetParameter(2)*func.GetParameter(6))
                else:
                    self.funcs[f'{mode} N-Term {type}'].FixParameter(0, func.GetParameter(0))
                    self.funcs[f'{mode} S-Term {type}'].FixParameter(0, func.GetParameter(1))
                    self.funcs[f'{mode} S-Term {type}'].FixParameter(1, func.GetParameter(3))
                    self.funcs[f'{mode} C-Term {type}'].FixParameter(0, func.GetParameter(2))
        
            func_name = f'{mode} func Ratio'
            print(blue(f'Output ratio for {functional_form_Ratio}'))
            self.funcs[func_name] = rt.TF1(func_name, functional_form_Ratio, self.func_min, self.func_max, npars_ratio)
            func = self.funcs[func_name]
            if self.fit_k:
                func_Data = self.funcs[func_name.replace('Ratio','Data')]
                for i in range(npars_ratio):
                    func.FixParameter(i,func_Data.GetParameter(i))
            else:
                func_MC = self.funcs[func_name.replace('Ratio','MC')]
                func_Data = self.funcs[func_name.replace('Ratio','Data')]
                for i in range(func_MC.GetNpar()):
                    func.FixParameter(i,func_MC.GetParameter(i))
                for i in range(func_Data.GetNpar()):
                    func.FixParameter(i+func_MC.GetNpar(),func_Data.GetParameter(i))
                    if not self.fit_k and i==3:
                        func.FixParameter(i+func_MC.GetNpar(),func_Data.GetParameter(i)/func_MC.GetParameter(i))
            PrintFuncParameters(func=func, info='  --> Fixed at:', color=green)
            self.funcs[f'{mode} N-Term Ratio'] = rt.TF1('N-Term Ratio', "pol0", self.func_min, self.func_max, 1)
            self.funcs[f'{mode} S-Term Ratio'] = rt.TF1('S-Term Ratio', "pol0", self.func_min, self.func_max, 1)
            self.funcs[f'{mode} C-Term Ratio'] = rt.TF1('C-Term Ratio', "pol0", self.func_min, self.func_max, 1)
            if self.fit_k:
                self.funcs[f'{mode} N-Term Ratio'].SetParameter(0, func.GetParameter(4))
                self.funcs[f'{mode} S-Term Ratio'].SetParameter(0, func.GetParameter(5))
                self.funcs[f'{mode} C-Term Ratio'].SetParameter(0, func.GetParameter(6))
                print(blue(f'SF k-parameters:'))
            else:
                self.funcs[f'{mode} N-Term Ratio'].SetParameter(0, func.GetParameter(4)/func.GetParameter(0))
                self.funcs[f'{mode} S-Term Ratio'].SetParameter(0, func.GetParameter(5)/func.GetParameter(1))
                self.funcs[f'{mode} C-Term Ratio'].SetParameter(0, func.GetParameter(6)/func.GetParameter(2))
                print(cyan(f'SF NSC ratio:'))
            print(cyan(f'  --> kN {self.funcs[f"{mode} N-Term Ratio"].GetParameter(0)}'))
            print(cyan(f'  --> kS {self.funcs[f"{mode} S-Term Ratio"].GetParameter(0)}'))
            print(cyan(f'  --> kC {self.funcs[f"{mode} C-Term Ratio"].GetParameter(0)}'))

        
    def Plot(self, to_plot, pdfname, nEntries=8,nFunc=2, zoom=True):
        self.CreateCanvas(canvName=pdfname, nEntries=nEntries, nFunc=nFunc, zoom=zoom)
        self.dicanv.cd(1)
        self.leg.Draw('same')
        for name,style in self.func_style.items():
            if not name in to_plot: continue
            if not name in self.funcs: continue
            self.dicanv.cd(2 if 'Ratio' in name else 1)
            legName = style.get('label', None)
            style.pop('label', None) # label needed for legend but not for tdrDraw
            tdrDrawLine(self.funcs[name], **style)
            if not 'Ratio' in name and legName:
                if 'JER_comb' in pdfname and not 'Term' in legName:
                    if 'Data' in legName:
                        legName = 'JER'
                    else: continue
                self.leg_func.AddEntry(self.funcs[name], legName, 'l')
        
        self.dicanv.cd(1)
        self.leg_func.Draw('same')
        for name,style in self.graph_style.items():
            if not name in to_plot: continue
            if not name in self.graphs: continue
            if 'Ratio' in name: self.dicanv.cd(2)
            else: self.dicanv.cd(1)
            legName = style.get('label', name)
            style.pop('label', None) # label needed for legend but not for tdrDraw
            graph = self.graphs[name]
            tdrDraw(graph, 'Pz', **style)
            if not 'Ratio' in name:
                self.leg.AddEntry(graph, legName, 'PLE')
        self.dicanv.SaveAs(os.path.join(self.outputPath, pdfname+self.pdfextraname+'.pdf'))

    def PlotCombination(self):
        for mode in self.combinations.keys():
            to_plot = []
            to_plot.append(f'MC Truth CI')
            for type in self.all_types:
                to_plot.append(f'combined {mode} {type}')
                to_plot.append(f'combined {mode} {type} Ratio')
                to_plot.append(f'combined {mode} func {type}')
                to_plot.append(f'combined {mode} N-Term {type}')
                to_plot.append(f'combined {mode} S-Term {type}')
                to_plot.append(f'combined {mode} C-Term {type}')
                to_plot.append(f'MC Truth CI')
            self.Plot(to_plot=to_plot, pdfname=f'JER_comb_{mode}', nEntries=3,nFunc=5, zoom=False)
    
    def PlotAll(self):
        self.PlotRC()
        
        to_plot = []
        for mode in self.mc_truth_modes:
            for fit in self.mc_truth_fits:
                to_plot.append(f'MC Truth {mode}{fit}')
                to_plot.append(f'MC Truth {mode}{fit} Ratio')
        self.Plot(to_plot=to_plot, pdfname='JER_MC_truth', nEntries=len(self.mc_truth_modes), nFunc=len(self.mc_truth_modes)*len(self.mc_truth_fits))

        to_plot = []
        type = 'MC'
        for mode in self.mc_truth_modes:
            to_plot.append(f'MC Truth {mode}')
            to_plot.append(f'MC Truth {mode} Ratio')
        for sample in self.mpfx_samples+self.balance_samples:
            to_plot.append(f'{sample} {type}')
            for mode in self.mc_truth_modes:
                to_plot.append(f'{sample} {type} Ratio {mode}')
        self.Plot(to_plot=to_plot, pdfname='JER_MC')

        # for sample in self.mpfx_samples+self.balance_samples:
        #     mu = 'RatioMu50to60' if not 'low-PU' in sample else 'RatioMu0to10'
        #     to_plot.append(' '.join([sample,'MC', mu]))
        # self.Plot(to_plot=to_plot, pdfname='JER_MC_ratio')

        # to_plot = list(filter(lambda x: not 'RatioMu' in x, to_plot))
        # for sample in self.mpfx_samples+self.balance_samples:
        #     mu = 'RatioMu30to40' if not 'low-PU' in sample else 'RatioMu10to20'
        #     to_plot.append(' '.join([sample,'MC', mu]))
        # self.Plot(to_plot=to_plot, pdfname='JER_MC_ratio Gauss')

        # to_plot = []
        # for sample in self.balance_samples:
        #     for type in self.all_types:
        #         to_plot.append(' '.join([sample,type]))
        # self.Plot(to_plot=to_plot, pdfname='JER_Data_noMPFX')

        # to_plot = []
        # for sample in self.mpfx_samples:
        #     for type in self.all_types:
        #         to_plot.append(' '.join([sample,type,'JER']))
        # self.Plot(to_plot=to_plot, pdfname='JER_Data_MPFX')

        to_plot = []
        for sample in self.mpfx_samples+self.balance_samples:
            for type in self.all_types:
                to_plot.append(f'{sample} {type}')
        self.Plot(to_plot=to_plot, pdfname='JER_Data')

    def SetStyle(self):
        self.style =  OrderedDict([
            ('MC Truth Gauss', {'color': rt.kMagenta+2, 'marker': rt.kFullCross,        'msize':1.2, }),
            ('MC Truth Gauss N fix', {'color': rt.kRed+1, 'marker': rt.kFullCross,        'msize':1.2, }),
            ('MC Truth Gauss P fix', {'color': rt.kOrange+1, 'marker': rt.kFullCross,        'msize':1.2, }),
            ('MC Truth Gauss NP fix', {'color': rt.kGreen+1, 'marker': rt.kFullCross,        'msize':1.2, }),
            ('MC Truth CI',    {'color': rt.kViolet+1,  'marker': rt.kFullCross,        'msize':1.2, }),
            ('MC Truth CI N fix',    {'color': rt.kRed+1,  'marker': rt.kFullCross,        'msize':1.2, }),
            ('MC Truth CI P fix',    {'color': rt.kOrange+1,  'marker': rt.kFullCross,        'msize':1.2, }),
            ('MC Truth CI NP fix',    {'color': rt.kGreen+1,  'marker': rt.kFullCross,        'msize':1.2, }),
            ('N-Term',         {'color': rt.kRed+1,     'marker': rt.kFullCircle,       'msize':0.8, }),
            ('S-Term',         {'color': rt.kOrange+1,  'marker': rt.kFullCircle,       'msize':0.8, }),
            ('C-Term',         {'color': rt.kAzure+2,   'marker': rt.kFullCircle,       'msize':0.8, }),
            ('MPFX',           {'color': rt.kOrange+1,  'marker': rt.kFullSquare,       'msize':0.8, }),
            ('MPFX RC',        {'color': rt.kMagenta+7, 'marker': rt.kFullDiamond,      'msize':1.2, }),
            ('dijet v19 SM',   {'color': rt.kAzure+2,   'marker': rt.kFullCross,        'msize':1.2, }),
            ('dijet v19 FE',   {'color': rt.kBlue+2,    'marker': rt.kFullDiamond,      'msize':1.2, }),
            ('dijet SM',       {'color': rt.kGreen+2,   'marker': rt.kFullTriangleDown, 'msize':0.9, }),
            ('dijet FE',       {'color': rt.kAzure+2,   'marker': rt.kFullTriangleUp,   'msize':0.9, }),
            ('zjet',           {'color': rt.kRed+1,     'marker': rt.kFullDiamond,      'msize':1.2, }),
            ('combined',       {'color': rt.kGreen+2,   'marker': rt.kFullCircle,       'msize':0.8, }),
        ])

        self.graph_style = OrderedDict([
            ('MC Truth Gauss',              {'label':'MC Truth Gauss'}),
            ('MC Truth Gauss Ratio',        {'label':''}),
            ('MC Truth CI',                 {'label':'MC Truth CI'}),
            ('MC Truth CI Ratio',           {'label':''}),
            ('combined all Data',           {'mcolor':rt.kAzure+2,   'marker':rt.kFullCircle,       'msize':0.8, 'label':'Data'  }),
            ('combined all MC',             {'mcolor':rt.kRed+1,     'marker':rt.kOpenCircle,       'msize':0.8, 'label':'MC'    }),
            ('combined all Ratio',          {'mcolor':rt.kBlack,     'marker':rt.kFullCircle,       'msize':0.8, 'label':'Ratio' }),

            ('combined all Data Ratio',     {'mcolor':rt.kGray+1,       'marker':rt.kFullCircle,       'msize':0.8, 'label':''}),
            ('combined all MC Ratio',       {'mcolor':rt.kGray+2,       'marker':rt.kOpenCircle,       'msize':0.8, 'label':''}),

            ('combined balance Data',       {'mcolor':rt.kAzure+2,   'marker':rt.kFullCircle,       'msize':0.8, 'label':'Data'  }),
            ('combined balance MC',         {'mcolor':rt.kRed+1,     'marker':rt.kOpenCircle,       'msize':0.8, 'label':'MC'    }),
            ('combined balance Ratio',      {'mcolor':rt.kBlack,     'marker':rt.kFullCircle,       'msize':0.8, 'label':'Ratio' }),

            ('combined balance Data Ratio', {'mcolor':rt.kGray+1,      'marker':rt.kFullCircle,       'msize':0.8, 'label':''  }),
            ('combined balance MC Ratio',   {'mcolor':rt.kGray+2,      'marker':rt.kOpenCircle,       'msize':0.8, 'label':''    }),

            ('combined mpfx Data',          {'mcolor':rt.kAzure+2,   'marker':rt.kFullCircle,       'msize':0.8, 'label':'Data'  }),
            ('combined mpfx MC',            {'mcolor':rt.kRed+1,     'marker':rt.kOpenCircle,       'msize':0.8, 'label':'MC'    }),
            ('combined mpfx Ratio',         {'mcolor':rt.kBlack,     'marker':rt.kFullCircle,       'msize':0.8, 'label':'Ratio' }),

            ('combined mpfx Data Ratio',    {'mcolor':rt.kGray+1,      'marker':rt.kFullCircle,       'msize':0.8, 'label':''  }),
            ('combined mpfx MC Ratio',      {'mcolor':rt.kGray+2,      'marker':rt.kOpenCircle,       'msize':0.8, 'label':''    }),
        ])

        for sample in self.mpfx_samples+self.balance_samples:
            for type in self.all_types+['MC Ratio CI', 'MC Ratio Gauss']:
                self.graph_style[f'{sample} {type}'] = {}
            
        self.func_style = OrderedDict()
        for type in self.all_types:
            for mode in self.combinations.keys():
                for func in ['func', 'N-Term','S-Term','C-Term']:
                    self.func_style[f'combined {mode} {func} {type}'] = {}
        
        for mode in ['CI', 'Gauss']:
            for fit in self.mc_truth_fits:
                self.func_style[f'MC Truth {mode}{fit}'] = {}
                self.graph_style[f'MC Truth {mode}{fit} Ratio'] = {}

        for func in self.func_style.keys():
            for name, info in self.style.items():
                # if name in func:
                if name in func and not ('combined' in name and 'Term' in func):
                    self.func_style[func]['lcolor'] = info['color']
            if 'Data' in func:
                label='Data'
                if 'Term' in func:
                    label = func.replace(' Data', '')
                    for comb in self.combinations:
                        label = label.replace(f'{comb} ','')
            elif 'Term' in func: label= ''
            elif 'Ratio' in func: label= 'Ratio'
            elif 'MC' in func: label= 'MC'
            if 'Truth' in func: label = func
            self.func_style[func]['lstyle'] = rt.kDashed if ('MC' in func and not 'Truth' in func) else rt.kSolid 
            self.func_style[func]['label'] = label
        
        for func in self.graph_style.keys():
            for name, info in self.style.items():
                if name in func:
                    if 'Ratio' in func and 'combined' in name: continue
                    if 'Ratio' in func and 'combined' in name: continue
                    marker = GetEmptyMarker(info['marker']) if 'MC' in func and not 'CI' in func else info['marker']
                    self.graph_style[func]['mcolor'] = info['color']
                    self.graph_style[func]['msize'] = info['msize']
                    self.graph_style[func]['marker'] = marker
            


    def Store(self):
        self.files['output'] = rt.TFile(self.outputPath+'output_CombineJER.root', 'RECREATE')
        self.files['output'].cd()
        for pu in self.mpfx_samples:
            for type in self.all_types:
                name = f'{pu} {type}'
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
    CJ.CombineDatasets()
    CJ.FitJER()
    CJ.PlotCombination()
    CJ.Store()
    CJ.Close()

if __name__ == '__main__':
    main()
