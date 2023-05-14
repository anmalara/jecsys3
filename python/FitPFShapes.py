from utils import *

def FTest(chi2_1, chi2_2, npar_1, npar_2, n_):
    return 1. - rt.TMath.FDistI(((chi2_1-chi2_2)/(npar_2-npar_1))/(chi2_2/(n_-npar_2-1)), npar_2-npar_1, n_-npar_2)

class FitPFShapes():
    def __init__(self, year='UL18'):
        self.year = year
        self.types = ['Rjet', 'chf', 'nhf', 'gammaf']
        self.shape_names = ['hadHcalp3', 'hadHcalZB097', 'hadHcalZB100', 'hadHcalZB106',
            'ecalm3', 'ecalGain1p3', 'ecalGain6p3', 'ecalGain12p3',
            'trkEff0999Nm1', 'trkEff0998Nm1', 'trkEffNtrk1m3', 'trkEffNtrk2ToInfm3'
            ]
        self.inputPath = './rootfiles/PFShapes/'
        self.outputPath = './pdfs/PFShapes/'
        if not os.path.exists(self.inputPath):
            self.inputPath = '.'+self.inputPath
            self.outputPath = '.'+self.outputPath
        os.system('mkdir -p '+self.outputPath)
        TDR.extraText   = 'Simulation'
        TDR.extraText2  = 'Preliminary'
        TDR.cms_lumi = TDR.commonScheme['legend'][self.year]
        TDR.cms_energy = TDR.commonScheme['energy'][self.year]
        TDR.extraText3 = []
        TDR.extraText3.append('AK4 CHS')
        TDR.extraText3.append('|#eta| < 1.3, #mu=0')
        self.pdfextraname = ''

        self.info = {
            'Rjet':               ('1-Resp', rt.kBlack),
            'chf':                ('chf',  rt.kRed+1),
            'nhf':                ('nhf',  rt.kGreen+2),
            'gammaf':             ('nef',  rt.kBlue+2),
            
            'hadHcalp3':          ('HCAL scale +3%', rt.kRed+1),
            'hadHcalZB097':       ('HCAL ZB -3%', rt.kCyan+1),
            'hadHcalZB100':       ('HCAL ZB +0%', rt.kAzure-2),
            'hadHcalZB106':       ('HCAL ZB +6%', rt.kViolet+2),
            
            'ecalm3':             ('ECAL scale +3%', rt.kGreen+2),
            'ecalGain1p3':        ('ECAL gain 1 +3%', rt.kTeal+2),
            'ecalGain6p3':        ('ECAL gain 6 +3%', rt.kGreen-1),
            'ecalGain12p3':       ('ECAL gain 12 +3%', rt.kGreen+1),
            
            'trkEff0999Nm1':      ('#varepsilon_{trk}=0.999^{N-1}', rt.kOrange-1),
            'trkEff0998Nm1':      ('#varepsilon_{trk}=0.998^{N-1}', rt.kOrange),
            'trkEffNtrk1m3':      ('N_{trk}=1 +3%', rt.kOrange+1),
            'trkEffNtrk2ToInfm3': ('N_{trk}>1 +3%', rt.kOrange+2),
            }
    
    def LoadInputs(self):
        self.files = OrderedDict()
        self.hists = OrderedDict()
        self.funcs = OrderedDict()
        self.files[self.year] = rt.TFile(self.inputPath+'pf_shape_variations_'+self.year+'.root')
        self.CreateCanvas() # Needed to get rid of annoying stat box
        for type in self.types:
            for shape in self.shape_names:
                hname = (type,shape)
                self.hists[hname] = self.files[self.year].Get(type+'_'+shape)
                self.hists[hname].SetDirectory(0)
                if type =='Rjet':
                    self.hists[hname] = ShiftHist(self.hists[hname], -1)
                self.hists[hname].Scale(100)
    
    def Fit(self):
        formulas = {
            2: "+[0]+[1]*log(x)",
            3: "+[0]+[1]*log(x)+[2]*pow(log(x),2)",
            4: "+[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)",
            5: "+[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)",
            6: "+[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)+[5]*pow(log(x),5)"
        }
        deg  = 5
        for hname, hist in self.hists.items():
            form = formulas.get(deg)
            self.funcs[hname] = rt.TF1(hname[0]+hname[1],form,15,7000)
            hist.Fit(self.funcs[hname], "QRS+")
            for par in range(deg+1):
                form = form.replace('+['+str(par)+']', "{:+.3e}".format(self.funcs[hname].GetParameter(par) if par<deg else 0 ))
            print(hname,form)
                
    def Close(self):
        for f_ in self.files.values():
            f_.Close()

    def CreateCanvas(self, canvName='', zoom=True, nEntries=8, nEntries2=8):
        if 'canv' in self.__dict__: self.canv.Close()
        XMin, XMax = (15, 7500)
        YMin, YMax = (-2,3) if zoom else (-2,5)
        xName, yName = ('p_{T,jet} [GeV]', 'Variation [%]')
        self.canv = tdrCanvas('shapes'+self.year+canvName, XMin, XMax, YMin, YMax, xName, yName, square=kSquare, isExtraSpace=True)
        self.canv.SetLogx(True)
        self.leg  = tdrLeg(0.40,0.90-(nEntries+1)*0.040,0.90,0.90)
        if nEntries>6:
            nEntries, nCol = 6, nEntries//6
            self.leg  = tdrLeg(0.40,0.90-(nEntries+1)*0.030,0.90,0.90, 0.025)
            self.leg.SetNColumns(nCol)
        FixXAsisPartition(self.canv)
        self.lines= {}
        for y in [0]:
            self.lines[y] = rt.TLine(XMin, y, XMax, y)
            tdrDrawLine(self.lines[y], lcolor=rt.kBlack, lstyle=rt.kDashed, lwidth=1)
    
    def PlotAll(self):
        for type in self.types:
            self.Plot(pdfname=type, zoom=(type!='Rjet' and type!='nhf'), nEntries=len(self.shape_names))
            for shape in ['Hcal','ecal','trkEff']:
                self.Plot(pdfname=type, zoom=(shape!='Hcal' or type!='Rjet'), nEntries=4, extraName=shape)
        for shape in self.shape_names:
            self.Plot(pdfname=shape, zoom=(shape!='hadHcalZB106'), nEntries=len(self.types))

    def Plot(self, pdfname, nEntries, zoom=True, extraName=None):
        self.CreateCanvas(canvName=pdfname, zoom=zoom, nEntries=nEntries)
        for (type, shape), hist in self.hists.items():
            legname, color = self.info[type]
            if pdfname!=shape: legname, color = self.info[shape]
            if extraName!=None and not (extraName in type or extraName in shape): continue
            if pdfname in type or pdfname in shape:
                tdrDraw(hist, 'P', mcolor=color)
                list(hist.GetListOfFunctions())[0].SetLineColor(color)
                self.leg.AddEntry(hist, legname, 'pl')
        if extraName!=None: pdfname += extraName
        self.canv.SaveAs(os.path.join(self.outputPath, 'PF_shapes_variations_' + pdfname + self.pdfextraname + '.pdf'))



def main():
    FS = FitPFShapes()
    FS.LoadInputs()
    FS.Fit()
    FS.PlotAll()
    FS.Close()

if __name__ == '__main__':
    main()