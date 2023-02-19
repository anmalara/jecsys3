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
    print(list(gnew.GetX()))
    for bin in range(0,g1.GetN()):
        x = list(gnew.GetX())[bin]
        y1 = list(g1.GetY())[bin]
        y2 = list(g2.GetY())[bin]
        ex = list(g1.GetEX())[bin]
        ey = oplus(list(g1.GetEY())[bin], list(g2.GetEY())[bin])
        gnew.SetPoint(bin, x, oplus(ominus(y1,y2), RC/x) )
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
        f_RC.Close()
    
    def LoadDijet(self):
        self.files['dijet'] = rt.TFile(self.inputPath+'dijet_balance_UL18.root')
        for type in ['Data', 'MC']:
            for mode in ['SM', 'FE']:
                name = ' '.join(['dijet',type,mode,'JER'])
                self.graphs[name] = self.files['dijet'].Get('dijet_balance_jer_'+type+'_0p00000_0p261_'+mode+'_nominal')
                if type=='MC':
                    self.graphs[name.replace('MC', 'Ratio')] = MakeRatioGraphs(self.graphs[name.replace('MC', 'Data')], self.graphs[name], name.replace('MC', 'Ratio'))
    
    def LoadMCTruth(self):
        name = 'MC Truth avg-mu JER'
        self.files[name] = rt.TFile(self.inputPath+'JER_MCtruth_avg_mu_UL18.root')
        self.graphs[name] = self.files[name].Get('ak4pfchsl1l2l3/RelResVsJetPt_JetEta0.261to0.522_Mu0to60')
        name = 'MC Truth bin-mu JER'
        self.files[name] = rt.TFile(self.inputPath+'JER_MCtruth_UL18.root')
        self.graphs[name+'Mu0to10'] = self.files[name].Get('ak4pfchsl1l2l3/RelResVsJetPt_JetEta0.261to0.522_Mu0to10')
        self.graphs[name+'Mu10to20'] = self.files[name].Get('ak4pfchsl1l2l3/RelResVsJetPt_JetEta0.261to0.522_Mu10to20')
        self.graphs[name+'Mu20to30'] = self.files[name].Get('ak4pfchsl1l2l3/RelResVsJetPt_JetEta0.261to0.522_Mu20to30')
        self.graphs[name+'Mu30to40'] = self.files[name].Get('ak4pfchsl1l2l3/RelResVsJetPt_JetEta0.261to0.522_Mu30to40')
        self.graphs[name+'Mu40to50'] = self.files[name].Get('ak4pfchsl1l2l3/RelResVsJetPt_JetEta0.261to0.522_Mu40to50')
        self.graphs[name+'Mu50to60'] = self.files[name].Get('ak4pfchsl1l2l3/RelResVsJetPt_JetEta0.261to0.522_Mu50to60')


    def LoadInputs(self):
        self.LoadRC()
        self.LoadMPFX()
        self.LoadDijet()
        self.LoadMCTruth()


    def FixXAsisPartition(self, canv, shift=None):
        canv.SetLogx(True)
        GettdrCanvasHist(canv).GetXaxis().SetNoExponent(True)
        latex = rt.TLatex()
        latex.SetTextFont(42)
        latex.SetTextSize(0.05)
        latex.SetTextAlign(23)
        if shift==None:
            YMin, YMax = (GettdrCanvasHist(canv).GetYaxis().GetXmin(), GettdrCanvasHist(canv).GetYaxis().GetXmax())
            shift = YMin-0.018*(YMax-YMin)
        for xbin in [30,100,300,1000, 3000]:
            latex.DrawLatex(xbin,shift,str(xbin))

    def CreateCanvas(self, canvName=''):
        XMin, XMax = (15, 2500)
        YMin, YMax = (0,1)
        xName, yName = ('p_{T,ave} [GeV]', 'RMS')
        self.dicanv = tdrDiCanvas('dicanvas_JER'+self.year+canvName, XMin, XMax, 0,0.35, 0.85,1.35, xName, yName, 'Data/MC')
        self.dicanv.cd(1).SetLogx(True)
        self.dicanv.cd(2).SetLogx(True)
        self.canv = tdrCanvas('JER'+self.year+canvName, XMin, XMax, YMin, YMax, xName, yName, square=kSquare, isExtraSpace=True)
        self.canv.SetLogx(True)
        self.leg = tdrLeg(0.50,0.90-8*0.045,0.70,0.90)
        self.FixXAsisPartition(self.canv)


    def Plot(self):
        self.CreateCanvas()
        styles = OrderedDict([
            ('low-PU Data MPF',   {'mcolor':rt.kRed,    'marker':rt.kFullCircle, 'msize':0.8}),
            ('low-PU MC MPF',     {'mcolor':rt.kRed,    'marker':rt.kOpenCircle, 'msize':0.8}),
            ('low-PU Data MPFX',  {'mcolor':rt.kBlue,   'marker':rt.kFullCircle, 'msize':0.8}),
            ('low-PU MC MPFX',    {'mcolor':rt.kBlue,   'marker':rt.kOpenCircle, 'msize':0.8}),
            ('high-PU Data MPF',  {'mcolor':rt.kRed+2,  'marker':rt.kFullSquare, 'msize':0.5}),
            ('high-PU MC MPF',    {'mcolor':rt.kRed+2,  'marker':rt.kOpenSquare, 'msize':0.5}),
            ('high-PU Data MPFX', {'mcolor':rt.kBlue+2, 'marker':rt.kFullSquare, 'msize':0.5}),
            ('high-PU MC MPFX',   {'mcolor':rt.kBlue+2, 'marker':rt.kOpenSquare, 'msize':0.5}),
        ])
        for name,style in styles.items():
            hist = self.graphs[name]
            tdrDraw(hist, 'Pz', **style)
            self.leg.AddEntry(hist, name, 'PLE')
        fixOverlay()
        self.canv.SaveAs(os.path.join(self.outputPath, 'MPFX.pdf'))
        self.canv.Close()

        styles = OrderedDict([
            ('low-PU Data JER',   {'mcolor':rt.kRed,    'marker':rt.kFullCircle, 'msize':0.8}),
            ('low-PU MC JER',   {'mcolor':rt.kRed,    'marker':rt.kOpenCircle, 'msize':0.8}),
            ('low-PU Ratio JER',   {'mcolor':rt.kRed,    'marker':rt.kFullCircle, 'msize':0.8}),
            ('high-PU Data JER',  {'mcolor':rt.kBlue+2,  'marker':rt.kFullSquare, 'msize':0.5}),
            ('high-PU MC JER',  {'mcolor':rt.kBlue+2,  'marker':rt.kOpenSquare, 'msize':0.5}),
            ('high-PU Ratio JER',  {'mcolor':rt.kBlue+2,  'marker':rt.kFullSquare, 'msize':0.5}),

            ('dijet Data SM JER',  {'mcolor':rt.kOrange+1,  'marker':rt.kFullTriangleDown, 'msize':0.5}),
            ('dijet MC SM JER',  {'mcolor':rt.kOrange+1,  'marker':rt.kOpenTriangleDown, 'msize':0.5}),
            ('dijet Ratio SM JER',  {'mcolor':rt.kOrange+1,  'marker':rt.kFullTriangleDown, 'msize':0.5}),

            ('dijet Data FE JER',  {'mcolor':rt.kGreen+2,  'marker':rt.kFullTriangleUp, 'msize':0.5}),
            ('dijet MC FE JER',  {'mcolor':rt.kGreen+2,  'marker':rt.kOpenTriangleUp, 'msize':0.5}),
            ('dijet Ratio FE JER',  {'mcolor':rt.kGreen+2,  'marker':rt.kFullTriangleUp, 'msize':0.5}),


            # ('MC Truth avg-mu JER',  {'mcolor':rt.kViolet+7,  'marker':rt.kFullStar, 'msize':0.9}),
            # ('MC Truth bin-mu JERMu0to10',  {'mcolor':rt.kViolet+7,  'marker':rt.kOpenStar, 'msize':0.9}),
            ('MC Truth bin-mu JERMu10to20',  {'mcolor':rt.kViolet+7,  'marker':rt.kOpenStar, 'msize':0.9}),
            # ('MC Truth bin-mu JERMu20to30',  {'mcolor':rt.kViolet+7,  'marker':rt.kOpenStar, 'msize':0.9}),
            # ('MC Truth bin-mu JERMu30to40',  {'mcolor':rt.kViolet+7,  'marker':rt.kOpenStar, 'msize':0.9}),
            # ('MC Truth bin-mu JERMu40to50',  {'mcolor':rt.kViolet+7,  'marker':rt.kOpenStar, 'msize':0.9}),
            ('MC Truth bin-mu JERMu50to60',  {'mcolor':rt.kViolet+7,  'marker':rt.kOpenStar, 'msize':0.9}),

        ])
        self.leg.Clear()
        self.dicanv.cd(1)
        self.leg.Draw('same')
        for name,style in styles.items():
            # if 'Data' in name: continue
            if 'Ratio' in name: self.dicanv.cd(2)
            else: self.dicanv.cd(1)
            hist = self.graphs[name]
            tdrDraw(hist, 'Pz', **style)
            if not 'Ratio' in name: self.leg.AddEntry(hist, name, 'PLE')
        self.dicanv.SaveAs(os.path.join(self.outputPath, 'JER.pdf'))
    
    def Store(self):
        self.files['output'] = rt.TFile(self.outputPath+'output_CombineJER.root', 'RECREATE')
        self.files['output'].cd()
        for pu in ['low-PU','high-PU']:
            for type in ['Data', 'MC', 'Ratio']:
                for mode in ['JER']:
                    name = ' '.join([pu,type,mode])
                    self.graphs[name].Write() 

    def Close(self):
        for f_ in self.files.values():
            f_.Close()


def main():
    CJ = CombineJER()
    CJ.LoadInputs()
    CJ.Plot()
    CJ.Store()
    CJ.Close()

if __name__ == '__main__':
    main()
