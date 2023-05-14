import os, ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
import ctypes, math
from collections import OrderedDict
from tdrstyle_JERC import *
import tdrstyle_JERC as TDR
TDR.extraText  = 'Preliminary'

extraname  = '_zghm'
# extraname  = '_z'
# extraname  = '_g'
# extraname  = '_h'

extraname += '_mpfdb'
# extraname += '_mpf'
# extraname += '_db'

extraname += '_pf'
# extraname += '_nopf'

extraname += '_old'
# extraname += '_new'

class PlotGlobalFit():
    def __init__(self, run='Run2', year='Run2', algo='AK4 CHS'):
        self.run = run
        self.year = year
        self.algo = algo
        self.inputPath = './rootfiles/'
        self.outputPath = './pdfs/'
        if not os.path.exists(self.inputPath):
            self.inputPath = '../rootfiles/'
            self.outputPath = '../pdfs/'
        os.system('mkdir -p '+self.outputPath)
        # self.filename = 'output'+run+'.root'
        self.filename = 'output'+run+'_old.root'
        # self.filename = 'output'+run+'_new.root'
        TDR.cms_lumi = TDR.commonScheme['legend'][self.year]+', '+TDR.commonScheme['lumi'][self.year]+' fb^{-1}'
        TDR.cms_energy = TDR.commonScheme['energy'][self.year]
        TDR.extraText3 = []
        TDR.extraText3.append(self.algo)
        TDR.extraText3.append('|#eta| < 1.3')

    def LoadInputs(self):
        self.inputfile = ROOT.TFile(os.path.join(self.inputPath,self.filename))

        self.functionforms = OrderedDict([
            ('track',   {'leg': 'tracks (p0)',    'color': rt.kGreen+2,  'npar':4, 'func': 'ftd'}),
            ('pho',     {'leg': 'photons (p1)',   'color': rt.kCyan+1,   'npar':1, 'func': 'fp'}),
            ('had',     {'leg': 'hadrons (p2)',   'color': rt.kRed+1,    'npar':4, 'func': 'fhx'}),
            ('hcal',    {'leg': 'had. hcal (p3)', 'color': rt.kAzure-2,  'npar':4, 'func': 'fhh'}),
            ('ecal',    {'leg': 'had. ecal (p4)', 'color': rt.kViolet+2, 'npar':4, 'func': 'feh'}),
            ('shower',  {'leg': 'shower (p5)',    'color': rt.kOrange,   'npar':6, 'func': 'fhw'}),
            ('offset',  {'leg': 'offset (p6)',    'color': rt.kTeal+2,   'npar':3, 'func': '(fl1-1)*100'}),
            ('mctrack', {'leg': 'MC tracks (p7)', 'color': rt.kGreen-1,  'npar':4, 'func': 'ftd-ftm'}),
            ('flav',    {'leg': 'flavor (p8)',    'color': rt.kOrange+1, 'npar':3, 'func': 'f1q3-1'}),

            # ('offset',  {'leg': 'offset',                          'color': rt.kTeal+2,   'npar':3, 'func': 'fl1'}),
            # ('shower',  {'leg': 'shower',                          'color': rt.kOrange,   'npar':6, 'func': 'fhw'}),
            # ('flav',    {'leg': 'flavor',                          'color': rt.kOrange+1, 'npar':3, 'func': 'f1q3-1'}),
            # ('hcal',    {'leg': 'HCAL +3%',                        'color': rt.kRed+1,    'npar':5, 'func': 'hadHcalp3'}),
            # ('ZB3',     {'leg': 'ZB -3%',                          'color': rt.kCyan+1,   'npar':5, 'func': 'hadHcalZB097'}),
            # ('ZB0',     {'leg': 'ZB +0%',                          'color': rt.kAzure-2,  'npar':5, 'func': 'hadHcalZB100'}),
            # ('ZB6',     {'leg': 'ZB +6%',                          'color': rt.kViolet+2, 'npar':5, 'func': 'hadHcalZB106'}),
            # ('ecal',    {'leg': 'ECAL +3%',                        'color': rt.kGreen+2,  'npar':5, 'func': 'ecalm3'}),
            # ('gain1',   {'leg': 'gain 1 +3%',                      'color': rt.kTeal+2,   'npar':5, 'func': 'ecalGain1p3'}),
            # ('gain6',   {'leg': 'gain 6 +3%',                      'color': rt.kGreen-1,  'npar':5, 'func': 'ecalGain6p3'}),
            # ('gain12',  {'leg': 'gain 12 +3%',                     'color': rt.kGreen+1,  'npar':5, 'func': 'ecalGain12p3'}),
            # ('TkrEff1', {'leg': '#varepsilon_{trk}=0.999^{N-1}',   'color': rt.kOrange-1, 'npar':5, 'func': 'trkEff0999Nm1'}),
            # ('TkrEff2', {'leg': '#varepsilon_{trk}=0.998^{N-1}',   'color': rt.kOrange+0, 'npar':5, 'func': 'trkEff0998Nm1'}),
            # ('TrkN1',   {'leg': 'N_{trk}=1 +3%',                   'color': rt.kOrange+1, 'npar':5, 'func': 'trkEffNtrk1m3'}),
            # ('TrkN2',   {'leg': 'N_{trk}>1 +3%',                   'color': rt.kOrange+2, 'npar':5, 'func': 'trkEffNtrk2ToInfm3'}),

            ('const',   {'leg': 'const (p0)',     'color': rt.kRed+1,    'npar':1, 'func': 'const'}),
            ])

        self.shapes = OrderedDict()
        for shape, info in self.functionforms.items():
            for mode in ['input','prefit','postfit']:
                for pf in ['','_chf','_nhf','_nef']:
                    self.shapes[shape+mode+pf] = self.inputfile.Get('shape_'+mode+'_'+info['func']+pf)

        self.infos = OrderedDict()
        self.infos['herr']                = {'objName':'herr',                     'legName': 'herr',                                        'legType':'extra', 'legStyle': 'f',   'plotinfo': {'opt':'e3', 'fcolor':ROOT.kCyan+1, 'lcolor':ROOT.kCyan+1, 'msize':0} }
        self.infos['hjesref']             = {'objName':'herr_ref',                 'legName': 'herr_ref',                                    'legType':'extra', 'legStyle': 'f',   'plotinfo': {'opt':'e3', 'fcolor':ROOT.kYellow+1, 'lcolor':ROOT.kYellow+1, 'msize':0} }
        self.infos['jesFit']              = {'objName':'jesFit_Resp',              'legName': '#chi^{2}/n.d.f.',                             'legType':'extra', 'legStyle': 'l',   'plotinfo': { 'lcolor' : ROOT.kBlack, 'lstyle' : ROOT.kSolid,  'lwidth' : 1,} }
        self.infos['jesFit_down']         = {'objName':'jesFit_down_Resp',         'legName': None,                                          'legType':'',      'legStyle': 'l',   'plotinfo': { 'lcolor' : ROOT.kBlack, 'lstyle' : ROOT.kDashed, 'lwidth' : 1,} }
        self.infos['jesFit_up']           = {'objName':'jesFit_up_Resp',           'legName': None,                                          'legType':'',      'legStyle': 'l',   'plotinfo': { 'lcolor' : ROOT.kBlack, 'lstyle' : ROOT.kDashed, 'lwidth' : 1,} }

        for mode in ['mpf','db']:
            self.infos['zjet_'+mode]            = {'objName':'Resp_zjet_'+mode,            'legName': 'Z+jet',                                       'legType':mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'msize':1.2, 'marker': ROOT.kFullStar if mode!='db' else ROOT.kOpenStar,                 'mcolor':ROOT.kRed+1} }
            self.infos['gjet_'+mode]            = {'objName':'Resp_gamjet_'+mode,          'legName': '#gamma+jet',                                  'legType':mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'msize':0.8, 'marker': ROOT.kFullSquare if mode!='db' else ROOT.kOpenSquare,             'mcolor':ROOT.kBlue+1} }
            self.infos['hadw_'+mode]            = {'objName':'Resp_hadw_'+mode,            'legName': 'W#rightarrow qq',                             'legType':mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'msize':0.9, 'marker': ROOT.kFullCircle if mode!='db' else ROOT.kOpenCircle,             'mcolor':ROOT.kGreen+2} }
            self.infos['multijet_'+mode]        = {'objName':'Resp_multijet_'+mode,        'legName': 'Multijet (p_{T}^{'+ScaleText('leading')+'})', 'legType':mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'msize':0.9, 'marker': ROOT.kFullTriangleUp if mode!='db' else ROOT.kOpenTriangleUp,     'mcolor':ROOT.kBlack} }
            self.infos['multijet_recoil_'+mode] = {'objName':'Resp_multijet_recoil_'+mode, 'legName': 'Multijet (p_{T}^{'+ScaleText('recoil')+'})',  'legType':mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'msize':0.9, 'marker': ROOT.kFullTriangleDown if mode!='db' else ROOT.kOpenTriangleDown, 'mcolor':ROOT.kGray+1} }
            self.infos['incljet_'+mode]         = {'objName':'Resp_incljet_'+mode,         'legName': 'Incl. jet',                                   'legType':mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'msize':0.8, 'marker': ROOT.kFullDiamond if mode!='db' else ROOT.kOpenDiamond,           'mcolor':ROOT.kOrange+2} }

        self.pf_style = OrderedDict([('chf',(ROOT.kRed+1,ROOT.kFullTriangleUp)), ('nhf',(ROOT.kGreen+2,ROOT.kFullTriangleDown)), ('nef',(ROOT.kAzure+2,ROOT.kFullSquare))])
        self.infos_PF = OrderedDict()
        mode, color, marker = ('Resp', ROOT.kBlack, ROOT.kFullCircle)
        self.infos_PF[mode+'variation_output'] = {'objName':'jesFit_graph_'+mode, 'legName': mode, 'legStyle': 'lep', 'plotinfo': {'opt':'e3', 'fcolor':color, 'lcolor':color, 'alpha':0.5} }
        self.infos_PF[mode+'prefit'] = {'objName':mode+'_comb_shift', 'legName': mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'mcolor':color, 'marker':marker, 'msize':0.9} }
        # self.infos_PF[mode+'postfit'] = {'objName':mode+'_comb', 'legName': mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'mcolor':color, 'marker':marker} }

        for mode, (color,marker) in self.pf_style.items():
            self.infos_PF[mode+'variation_input'] = {'objName':mode+'_variation_input', 'legName': mode, 'legStyle': 'lep', 'plotinfo': {'opt':'e3', 'fcolor':color, 'lcolor':ROOT.kGray, 'alpha':0.3} }
            self.infos_PF[mode+'variation_output'] = {'objName':mode+'_variation_output', 'legName': mode, 'legStyle': 'lep', 'plotinfo': {'opt':'e3', 'fcolor':color, 'lcolor':color, 'alpha':0.5} }
            self.infos_PF[mode+'prefit'] = {'objName':mode+'_prefit', 'legName': mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'mcolor':color, 'marker':marker} }


        # Load objects
        for name, info in self.infos.items():
            if 'Resp_' in info['objName']:
                for mode in ['prefit', 'postfit', 'raw']:
                    info['obj_'+mode] = self.inputfile.Get(info['objName']+'_'+mode)
            else:
                info['obj'] = self.inputfile.Get(info['objName'])

        for name, info in self.infos_PF.items():
            info['obj'] = self.inputfile.Get(info['objName'])

        # Remove null objects
        for name in list(self.shapes.keys()):
            if not self.shapes[name]:
                del self.shapes[name]

        for name, info in self.infos.copy().items():
            if ('obj_raw' in info and not info['obj_raw']):
                del self.infos[name]
            if ('recoil' in name and not name.replace('_recoil','') in self.infos.keys()):
                del self.infos[name]

        for name, info in self.infos_PF.copy().items():
            if not info['obj']:
                del self.infos_PF[name]


    def CreateCanvasGlobalFit(self, canvName):
        XMin, XMax = (15, 4500)
        YMin, YMax = (0.95,1.07)
        # YMin, YMax = (0.99,1.02)
        yName = 'jet response ratio'
        if canvName =='raw': yName = 'raw '+yName
        if canvName =='postfit': yName = 'post-fit '+yName
        yName = yName.capitalize()
        self.canv = tdrCanvas('PlotGlobalFit'+self.year+canvName, XMin, XMax, YMin, YMax, 'p_{T} [GeV]', yName, square=kSquare, isExtraSpace=True)
        FixXAsisPartition(self.canv)
        l,r,t,b = (self.canv.GetLeftMargin(),self.canv.GetRightMargin(),self.canv.GetTopMargin(),self.canv.GetBottomMargin())
        nSamples_mpf   = sum(info['legType']=='mpf'   for info in self.infos.values())
        nSamples_db    = sum(info['legType']=='db'    for info in self.infos.values())
        nSamples_extra = sum(info['legType']=='extra' for info in self.infos.values())
        xpos = l+0.04*(1-l-r)+0.02
        ypos = b+0.04*(1-t-b)
        self.legs = {}
        self.legs['extra'] = tdrLeg(xpos,ypos,xpos+0.2,ypos+0.04*(nSamples_extra+1), textSize=0.04)
        xpos = 1-r-0.04*(1-l-r)-0.07
        ypos = 1-t-0.04*(1-t-b)
        self.legs['mpf'] = tdrLeg(xpos-0.2,ypos-0.04*(nSamples_mpf+1),xpos,ypos, textSize=0.04)
        self.legs['db'] = tdrLeg(xpos-0.24,ypos-0.04*(nSamples_db+1),xpos-0.2,ypos, textSize=0.04)
        tdrHeader(self.legs['mpf'], 'MPF', textAlign=13)
        tdrHeader(self.legs['db'], 'DB ', textAlign=33)
        self.line = ROOT.TLine(XMin, 1, XMax, 1)
        tdrDrawLine(self.line, lcolor=ROOT.kBlack, lstyle=ROOT.kDashed, lwidth=1)

    def PlotResponse(self, canvName):
        self.CreateCanvasGlobalFit(canvName)
        for name, info in self.infos.items():
            obj = info['obj_'+canvName if 'obj_'+canvName in info else 'obj']
            if not 'jesFit' in info['objName']:
                tdrDraw(obj, **info['plotinfo'])
            elif canvName=='postfit':
                tdrDrawLine(obj, **info['plotinfo'])
            leg_type = 'extra'
            if 'db' in name:
                self.legs['db'].AddEntry(obj, ' ', info['legStyle'])
            elif 'mpf' in name:
                self.legs['mpf'].AddEntry(obj, info['legName'], info['legStyle'])
            else:
                if info['legName']:
                    self.legs['extra'].AddEntry(obj, info['legName'], info['legStyle'])

        # info = self.infos_PF['Resppostfit']
        # obj = info['obj']
        # tdrDraw(obj, **info['plotinfo'])
        self.line.Draw('same')
        fixOverlay()
        self.canv.SaveAs(os.path.join(self.outputPath, self.filename.replace('.root', '_'+canvName+extraname+'.pdf')))
        self.canv.Close()


    def PlotPFVariations(self):
        XMin, XMax = (15, 4500)
        YMin, YMax = (-2+0.001,2.5-0.001)
        canvName = 'PFVariations'+self.year
        yName = 'PF composition change (10^{-2})'
        self.canv = tdrCanvas(canvName, XMin, XMax, YMin, YMax, 'p_{T} [GeV]', yName, square=kSquare, isExtraSpace=True)
        FixXAsisPartition(self.canv)
        line = ROOT.TLine(XMin,0, XMax, 0)
        tdrDrawLine(line, lcolor=ROOT.kBlack, lstyle=ROOT.kDashed, lwidth=1)
        for name, info in self.infos_PF.items():
            if 'Resppostfit' ==name: continue
            obj = info['obj']
            tdrDraw(obj, **info['plotinfo'])
            # if info['legName']:
            #     self.leg.AddEntry(obj, info['legName'], info['legStyle'])
        for pf in ['_chf','_nhf','_nef']:
            prefit, postfit = ('0','0')
            for name, shape in self.shapes.items():
                if 'input' in name or not pf in name: continue
                par_ = self.shapes[name.replace(pf,'')].GetParameter(0)
                par_ = str(par_) if par_<0 else '+'+str(par_)
                func_ = shape.GetTitle().replace('[0]',par_)
                if 'postfit' in name: postfit += func_
                elif 'prefit' in name: prefit += func_
                else: raise RuntimeError('Unexpected function')
            # self.shapes['prefit'+pf] = rt.TF1('prefit'+pf, prefit,15,4500)
            self.shapes['postfit'+pf] = rt.TF1('postfit'+pf, postfit,15,4500)
            # tdrDrawLine(self.shapes['prefit'+pf], lcolor=self.pf_style[pf.replace('_','')][0], lstyle=rt.kDashed, lwidth=2)
            tdrDrawLine(self.shapes['postfit'+pf], lcolor=self.pf_style[pf.replace('_','')][0], lstyle=rt.kSolid, lwidth=2)
        
        fixOverlay()
        self.canv.SaveAs(os.path.join(self.outputPath, canvName+extraname+'.pdf'))
        self.canv.Close()

    def PlotShapes(self, mode):
        XMin, XMax = (15, 4500)
        # YMin, YMax = (-2+0.001,2.5-0.001)
        YMin, YMax = (-2,2)
        canvName = 'GlobalFitShapes_'+mode+self.year
        self.canv = tdrCanvas(canvName, XMin, XMax, YMin, YMax, 'p_{T} [GeV]', 'JES change (%)', square=kSquare, isExtraSpace=True)
        FixXAsisPartition(self.canv)
        leg = tdrLeg(0.4, 0.9- 0.035*(len(self.functionforms)+(3 if mode=='postfit' else 1))/2, 0.92, 0.9, 0.035)
        leg.SetNColumns(2)
        sum='0'
        for shape, info in self.functionforms.items():
            if not shape+mode in self.shapes: continue
            func = self.shapes[shape+mode]
            sum += '+('+func.GetTitle()+')'
            if mode !='input':
                for par in range(func.GetNpar()):
                    par_ = func.GetParameter(par)
                    sum = sum.replace('['+str(par)+']',str(par_) if par_<0 else '+'+str(par_))
            tdrDrawLine(func, lcolor=info['color'], lstyle=rt.kSolid, lwidth=2)
            leg.AddEntry(func, info['leg'], 'l')
        if mode=='prefit': self.sum = sum
        self.shapes['sum'] = rt.TF1('sum', sum,15,4500)
        tdrDrawLine(self.shapes['sum'], lcolor=rt.kBlack, lstyle=rt.kSolid, lwidth=2)
        leg.AddEntry(self.shapes['sum'], 'sum', 'l')
        if mode=='postfit':
            var_sum = sum+'-('+self.sum+')'
            self.shapes['var_sum'] = rt.TF1('var_sum', var_sum,15,4500)
            tdrDrawLine(self.shapes['var_sum'], lcolor=rt.kBlack, lstyle=rt.kDashed, lwidth=2)
            leg.AddEntry(self.shapes['var_sum'], 'post-pre', 'l')
        fixOverlay()
        self.canv.SaveAs(os.path.join(self.outputPath, canvName+extraname+'.pdf'))
        self.canv.Close()

    def PlotCorrelation(self):
        print('Correlation2')

    def PlotAll(self):
        self.LoadInputs()
        self.PlotResponse('raw')
        self.PlotResponse('prefit')
        self.PlotResponse('postfit')
        self.PlotShapes('input')
        self.PlotShapes('prefit')
        self.PlotShapes('postfit')
        self.PlotPFVariations()
        self.PlotCorrelation()
        self.inputfile.Close()

def main():
    PlotGlobalFit().PlotAll()

if __name__ == '__main__':
    main()
