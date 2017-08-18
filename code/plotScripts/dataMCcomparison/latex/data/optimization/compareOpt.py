#!/usr/bin/env python
'''
For comparison of two optimization results.
dongliang.zhang@cern.ch
'''
import os, sys, re
import glob
from ROOT import *
from math import sqrt
from rootUtil import useAtlasStyle, waitRootCmd, savehistory, mkupHistSimple, get_default_fig_dir
funlist=[]

sDir = get_default_fig_dir()
sTag = 'test_'
sDirectly = False
if gROOT.IsBatch(): sDirectly = True

def updateDict(d, x):
    for j in d: d[j] = d[j]*d[j]
    lines = None
    with open(x,'read') as f:
        lines = f.readlines()
    for l in lines[1:]:
        if l[0] == '#': continue
        fs = l.rstrip().split()
        if len(fs)<5:
            print fs
            continue
        Zn = float(fs[4])
        if Zn<0: Zn = 0.
        k=(fs[0],fs[1])

        if fs[0]=='400.0' and fs[1]=='25.0':
            print Zn, d.get(k,-1)
        d[k] = Zn*Zn + d.get(k,0)
        if fs[0]=='400.0' and fs[1]=='25.0':
            print fs[0], fs[1], Zn, d.get(k,0)
    for j in d: d[j] = sqrt(d[j])


class OptResults:
    ### Info pattern match
    infoP = re.compile('#(.*): (.*) \+/- (.*) \(.*')
    nPoint = 10.

    def __init__(self, txtfile=None, tag='test', ref=None):
        self.tag = tag
        self.hNBkg = None 
        self.hBkgRef = None
        self.gNsig = TGraphErrors()
        self.gZnSig = TGraph()
        self.gTBkg = None
        self.min1 = 0.
        self.ref = ref
        if txtfile is not None:
            self.load(txtfile)
    def getZnDiff(self, ref=None):
        if ref is None: ref = self.ref
        if self.gZnSig.GetN() != ref.gZnSig.GetN(): return None
        gU = TGraph()
        gD = TGraph()
        for i in range(self.gZnSig.GetN()):
            xi,yi,zi = Double(0),Double(0),Double(0)
            ref.gZnSig.GetPoint(i, xi, yi)
            self.gZnSig.GetPoint(i, xi, zi)
            if zi>yi:
                gU.SetPoint(gU.GetN(), xi, zi)
            else:
                gD.SetPoint(gD.GetN(), xi, zi)
        gU.SetMarkerStyle(26)
        gD.SetMarkerStyle(32)
        return gU, gD

    def load(self, txtfile):
        nP = self.nPoint*0.5
        fSig = self.nPoint/100.
        with open(txtfile,'r') as f1:
            lines = f1.readlines()

        nTotalBkg = None
        bkgList = []
        for line in lines:
            line = line.rstrip()
            m = self.infoP.match(line)
            if m is not None:
                print m.group(1), m.group(2), m.group(3)
                head = m.group(1)
                if head == 'Total BG':
                    nTotalBkg = float(m.group(2))
                    self.gTBkg = TGraphErrors()
                    self.gTBkg.SetPoint(0, 0, float(m.group(2)))
                    self.gTBkg.SetPointError(0, nP, float(m.group(3)))
                    self.gTBkg.SetPoint(1, self.nPoint*0.9999, float(m.group(2)))
                    self.gTBkg.SetPointError(1, nP, float(m.group(3)))
                else:
                    bkgList.append((head, float(m.group(2)), float(m.group(3))))
#                     print head
                continue
            # Now phase signal
            # 550.0 0.0 0.166 0.050 -0.211
            values = line.split() 
            if len(values) < 5: continue

            n = self.gNsig.GetN()
            self.gNsig.SetPoint(n, (n+1)*fSig, float(values[2])+nTotalBkg)
            self.gNsig.SetPointError(n, 0, float(values[3]))
            self.gZnSig.SetPoint(n, (n+1)*fSig, float(values[4]))

        nBkgC = len(bkgList)
        ### to get the same order as the reference
        keys = range(nBkgC)
        if self.ref:
            x1 = self.ref.hNBkg.GetXaxis()
            ## Get the tag list in current list
            ltemp = [x[0] for x in bkgList]
            rlist = [x1.GetBinLabel(i+1) for i in range(x1.GetNbins())]

            st = set(ltemp+rlist)
#             print len(st), len(ltemp), len(rlist)
            nKey = len(st)
            nRef = len(rlist)

            if nKey != nBkgC:
                self.hBkgRef = TH1F(self.tag+"hBkgRef","hBkg",nKey, 0., self.nPoint)

                # copy the existing ones 
                for i in range(nRef):
                    lbx = x1.GetBinLabel(i+1)
#                     print "normal:", i+1, lbx
                    self.hBkgRef.SetBinContent(i+1,self.ref.hNBkg.GetBinContent(i+1))
                    self.hBkgRef.SetBinError(i+1,self.ref.hNBkg.GetBinError(i+1))
                    self.hBkgRef.GetXaxis().SetBinLabel(i+1,x1.GetBinLabel(i+1))
                    st.remove(lbx)
                ## add the new ones
                print st
                for t in st:
#                     print 'extra:', nRef+1, t
                    self.hBkgRef.GetXaxis().SetBinLabel(nRef+1,t)
                    nRef += 1

                x1 = self.hBkgRef.GetXaxis()
                self.hNBkg = self.hBkgRef.Clone(self.tag+"hNBkg")
                self.hNBkg.Reset()
#                 self.hNBkg.Clear()

            ## find their indices in the refernece
            keys = []
            for i in range(nKey):
                lbx = x1.GetBinLabel(i+1)
#                 print lbx
                keys.append(ltemp.index(lbx) if lbx in ltemp else -1)

#             print ltemp
#             print keys
#             print [x1.GetBinLabel(i+1) for i in range(len(keys))]

                   
#             keys = [ltemp.index(x1.GetBinLabel(i+1)) for i in range(nKey)]

        if self.hNBkg == None:
            self.hNBkg = TH1F(self.tag+"hNBkg","hBkg",nBkgC, 0., self.nPoint)

        ibin = 0
        for ix in keys:
            ibin += 1
            ### some keys in reference might not be in this plot
            if ix<0: continue

            x = bkgList[ix]
            self.hNBkg.SetBinContent(ibin, x[1])
            self.hNBkg.SetBinError(ibin, x[2])
            self.hNBkg.GetXaxis().SetBinLabel(ibin, x[0])

            if self.min1 > x[1]-x[2]: self.min1 = x[1]-x[2]
#         self.hNBkg.Draw()
#         waitRootCmd()

class ComparerX:
    def __init__(self, tag='test'):
        self.tag = tag
        self.rx1 = None
        self.rx2 = None
        self.channels = ['ee_1','ee_2','emu_1','emu_2','mumu_1','mumu_2']
        self.dir0 = '/home/dzhang/links/eosOther/cloShared/save/'
        self.sDir = sDir
        self.sTag = sTag
        self.opt = 'opt'
        self.autoSave = sDirectly

    def compareChans(self, chans = ['ee_1','ee_2','emu_1','emu_2','mumu_1','mumu_2']):
        for chan in chans:
            self.showCompare(chan)
    def showCompare(self,chan):
        r1 = OptResults(self.dir0+self.rx1[0]+'/SR_SS_'+chan+'_'+self.opt+'_0.txt','r1_')
        r2 = OptResults(self.dir0+self.rx2[0]+'/SR_SS_'+chan+'_'+self.opt+'_0.txt','r2_', r1)
        lg1 = self.rx1[1]
        lg2 = self.rx2[1]

        hfirst = r2.hNBkg
        hfirst.Draw("axis")

#         r1.gTBkg.SetFillStyle(3004)
        r1.gTBkg.SetFillColor(kGray)
        r1.gTBkg.Draw("F2")
        r1.gTBkg.Draw("L")
        r2.gTBkg.Draw("F2")
        r2.gTBkg.Draw("L")
        r2.gTBkg.SetFillStyle(3005)
        r2.gTBkg.SetFillColor(46)
        r2.gTBkg.SetLineColor(2)
    #     r2.gTBkg.SetMarkerColor(2)
    #     r2.gTBkg.SetMarkerStyle(26)

        r1.gNsig.Draw("Psame")
        r2.gNsig.Draw("Psame")
        r2.gNsig.SetLineColor(2)
        r2.gNsig.SetMarkerColor(2)
        r2.gNsig.SetMarkerStyle(26)

        r1.gZnSig.Draw("Psame")
#         r2.gZnSig.Draw("Psame")
        r1.gZnSig.SetMarkerColor(kGray)
#         r2.gZnSig.SetMarkerColor(kGray+2)
#         r2.gZnSig.SetMarkerStyle(26)
        gU,gD = r2.getZnDiff()
        gU.SetMarkerColor(kGray+2)
        gD.SetMarkerColor(kGray+2)
        gU.Draw("Psame")
        gD.Draw("Psame")

        fun1 = TF1("fun1", "0", 0, 100)
        fun1.SetLineStyle(2)
    #     fun1.SetLineColor(kGray)
        fun1.Draw("same")

        hBkgRef = r2.hBkgRef if r2.hBkgRef else r1.hNBkg
        hBkgRef.Draw("histsame")
        hBkgRef.Draw("Esame")
        r2.hNBkg.Draw("Esame")
        r2.hNBkg.SetLineColor(2)
        r2.hNBkg.SetMarkerColor(2)
        r2.hNBkg.SetMarkerStyle(26)

        max1 = max(r1.gNsig.GetHistogram().GetMaximum(), r2.gNsig.GetHistogram().GetMaximum())
        min1 = min(r1.min1, r2.min1)
        dm = 0.02*(max1-min1)
        hfirst.GetYaxis().SetRangeUser(min1-dm,max1+dm)
        hfirst.GetYaxis().SetTitle("Events")
        hfirst.Draw('axis same')
        

        lg = TLegend(0.7,0.8,0.9,0.93)
        lg.SetHeader(chan)
        lg.SetFillStyle(0)
        lg.AddEntry(hBkgRef,lg1,'lp')
        lg.AddEntry(r2.hNBkg,lg2,'lp')
        lg.Draw()

        gPad.Update()
        waitRootCmd(self.sDir+self.sTag+chan,self.autoSave)

    def getZnDiff2D(self, gr1,gr2):
        grDiff2D = gr1.Clone("grDiff2D")
        xiS,yiS,ziS = grDiff2D.GetX(),grDiff2D.GetY(),grDiff2D.GetZ()
        xjS,yjS,zjS = gr2.GetX(),gr2.GetY(),gr2.GetZ()
        for i in range(len(xiS)):
            if abs(xiS[i]-xjS[i])>0.001 or abs(xjS[i]-xjS[i])>0.001:
                print "the two graph does not match"
                return
            grDiff2D.SetPoint(i, xiS[i], yiS[i], ziS[i]-zjS[i])
#         grDiff2D.Draw('colz')
#         waitRootCmd()

        return grDiff2D;

    def getSig(self, args, show, info=None,savename="test"):
        dictA = {}
        for x in args:
            updateDict(dictA, x)
        print dictA

        g2d = TGraph2D()
        i = 0
        for x in dictA:
            g2d.SetPoint(i, float(x[0]),float(x[1]),dictA[x])
            i += 1

        if info: g2d.SetTitle(info)
        if not show:
            return g2d

    #     g2d.SetPoint(g2d.GetN(), 125, 0,0)

        gStyle.SetPadRightMargin(0.16)

        g2d.Draw('colz')
        h1 = g2d.GetHistogram();
    #     h1 = g2d.GetHistogram().Clone("h1");
        h1.GetXaxis().SetTitle("m_{#tilde{#chi}^{#pm}_{1}/#tilde{#chi}^{0}_{2}} [GeV]");
        h1.GetYaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} [GeV]");
        h1.GetZaxis().SetTitle("Z_{n}");
        h1.GetXaxis().SetRangeUser(0,350)
        h1.GetZaxis().SetRangeUser(0,4)

    #     g2d.RemovePoint(g2d.GetN()-1)
    #     h1.Draw('axis')
    #     g2d.Draw('colzsame')

        fun1 = TF1("linex","x-125.5", 0, 350);
        fun1.SetLineColor(kGray)
        fun1.SetLineWidth(2)
        fun1.SetLineStyle(8)
        fun1.Draw("same");

        lg = TLegend(0.2, 0.65,0.4,0.85)
        lg.SetFillStyle(0)

        ## run 1 limit: http://www.hepdata.net/record/ins1341609?version=1&table=Table31
        fRun1 = TFile('HEPData-ins1341609-v1.root','read')
        grRun1 = fRun1.Get('Table 31/Graph1D_y1')
        grRun1.SetLineColor(2)
        grRun1.Draw("Lsame")
        grRun1E = fRun1.Get('Table 30/Graph1D_y1')
        grRun1E.SetLineColor(2)
        grRun1E.SetLineStyle(2)
        grRun1E.Draw("Lsame")
        lg.AddEntry(grRun1E, 'Run1 expected', 'l')
        lg.AddEntry(grRun1, 'Run1 limit', 'l')

        ## contours
        lines = {1:2, 2:5, 3:9} # Zn, line Style
        for lx in lines:
            ltag = str(lx)
            grL = g2d.GetContourList(lx);
            next1 = TIter(grL);
            while True:
                obj = next1()
                if obj == None: break
                obj.SetLineWidth(2);
                obj.SetLineStyle(lines[lx]);
                obj.Draw("same");
                obj.SetName(ltag);
            xt0 = gPad.GetPrimitive(ltag);
    #         if xt0==None: continue
            lg.AddEntry(xt0, 'Z_{n}='+str(lx), 'l')

        lg.Draw()
        if info:
            lt = TLatex()
            lt.DrawLatexNDC(0.2,0.9,info)

        gPad.Update()
        waitRootCmd(self.sDir+self.sTag+savename, self.autoSave)

        return g2d


    def showCompareZn(self, gr1, gr2):
        h1 = TH2F('h1',"h1;m_{#tilde{#chi}^{#pm}_{1}/#tilde{#chi}^{0}_{2}} [GeV];m_{#tilde{#chi}^{0}_{1}} [GeV]",100,125,350,100,0,120)
        h1.Draw("axis")

        fun1 = TF1("linex","x-125.5", 0, 350);
        fun1.SetLineColor(kGray)
        fun1.SetLineWidth(2)
        fun1.SetLineStyle(8)
        fun1.Draw("same");

        lg = TLegend(0.2, 0.5,0.4,0.93)
        lg.SetFillStyle(0)
        lines = {1:2, 2:5, 3:9} # Zn, line Style

        lg.AddEntry(None, gr1.GetTitle(),"")
        for lx in lines:
            ltag = "s1_"+str(lx)
            grL = gr1.GetContourList(lx);
            next1 = TIter(grL);
            while True:
                obj = next1()
                if obj == None: break
                obj.SetLineWidth(2);
                obj.SetLineStyle(lines[lx]);
                obj.Draw("same");
                obj.SetName(ltag);
            xt0 = gPad.GetPrimitive(ltag);
            lg.AddEntry(xt0, 'Z_{n}='+str(lx), 'l')

        lg.AddEntry(None, "~~~~~~~","")
        lg.AddEntry(None, gr2.GetTitle(),"")
        for lx in lines:
            ltag = "s2_"+str(lx)
            grL = gr2.GetContourList(lx);
            next1 = TIter(grL);
            while True:
                obj = next1()
                if obj == None: break
                obj.SetLineColor(2);
                obj.SetLineWidth(2);
                obj.SetLineStyle(lines[lx]);
                obj.Draw("same");
                obj.SetName(ltag);
            xt0 = gPad.GetPrimitive(ltag);
            if xt0==None: continue
            lg.AddEntry(xt0, 'Z_{n}='+str(lx), 'l')

        ## run 1 limit: http://www.hepdata.net/record/ins1341609?version=1&table=Table31
        fRun1 = TFile('HEPData-ins1341609-v1.root','read')
        grRun1 = fRun1.Get('Table 31/Graph1D_y1')
        grRun1.SetFillColor(3)
        grRun1.Draw("Fsame")
        grRun1E = fRun1.Get('Table 30/Graph1D_y1')
        grRun1E.SetFillColor(4)
        grRun1E.SetFillStyle(3004)
        grRun1E.Draw("Fsame")
    #     lg.AddEntry(grRun1E, 'Run1 expected', 'l')
    #     lg.AddEntry(grRun1, 'Run1 limit', 'l')

        lg.Draw()
        gPad.Update()
        waitRootCmd(self.sDir+self.sTag+'ZnDiff', self.autoSave)

        h1.Draw("axis")
        grDiff = self.getZnDiff2D(gr2, gr1)
        grDiff.Draw('colz')
        fun1.Draw("same");
        gPad.Update()
        waitRootCmd(self.sDir+self.sTag+'ZnDiff2D', self.autoSave)

    def compare2X(self):
        show = True
        x1 = self.rx1[0].lstrip('significance_')
        x2 = self.rx2[0].lstrip('significance_')
        print "significance_"+x1+"/SR*.txt"
        gr1 = self.getSig(glob.glob("significance_"+x1+"/SR*.txt"), show, x1,x1.replace("2.8", "2p8"))
        gr2 = self.getSig(glob.glob("significance_"+x2+"/SR*.txt"), show, x2,x2.replace("2.8", "2p8"))
        
#         self.getZnDiff2D(gr2, gr1)
        self.showCompareZn(gr1,gr2)


def test1():
    cx1 = ComparerX('tes1')
#     cx1.rx1 = ('2.8_1D/significance','2.8-1D')
#     cx1.rx2 = ('2.8_1D_loose/significance','Loose')
#     cx1.sTag = "muEta2p7_"

    ### git version
    cx1.dir0 = './'
#     cx1.rx1 = ('significance_2.8_1D_avg','pt25')
#     cx1.rx2 = ('significance_2.8_1D_pt25_avg','qF')
    cx1.rx1 = ('significance_2p8_1D_pt25_avg','pt25')
    cx1.rx2 = ('significance_2p8_1D_pt25_qF_avg','qF')
    cx1.opt = 'opt'
    cx1.sTag = "muEta2p7Opt_"

#     cx1.compareChans()
    cx1.compare2X()

def test():
    cx1 = ComparerX('tes1')
    cx1.sDir = './dCompare/'
    cx1.autoSave = True

    if len(sys.argv)<3:
        print "Need at least two arguments"
        return

    A1 = sys.argv[1]
    A2 = sys.argv[2]
    stag = ('plots_'+A1+'_VS_plots_'+A2+'_').replace('.','p')

    ### git version
    cx1.dir0 = './'
    cx1.rx1 = ('significance_'+A1,A1)
    cx1.rx2 = ('significance_'+A2,A2)
    cx1.opt = 'opt'
    cx1.sTag = stag

    cx1.compareChans()


funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    test1()
#     for fun in funlist: print fun()
