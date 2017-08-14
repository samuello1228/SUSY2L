#!/usr/bin/env python
'''
For comparison of two optimization results.
dongliang.zhang@cern.ch
'''
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory, mkupHistSimple, get_default_fig_dir
funlist=[]

sDir = get_default_fig_dir()
sTag = 'test_'
sDirectly = False
if gROOT.IsBatch(): sDirectly = True

class OptResults:
    ### Info pattern match
    infoP = re.compile('#(.*): (.*) \+/- (.*) \(.*')
    nPoint = 10.

    def __init__(self, txtfile=None, tag='test', ref=None):
        self.tag = tag
        self.hNBkg = None 
        self.gNsig = TGraphErrors()
        self.gZnSig = TGraph()
        self.gTBkg = None
        self.min1 = 0.
        self.ref = ref
        if txtfile is not None:
            self.load(txtfile)

    def load(self, txtfile):
        nP = self.nPoint*0.5
        fSig = self.nPoint/100.
        print "testing", self.tag
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
                    self.gTBkg.SetPoint(1, self.nPoint, float(m.group(2)))
                    self.gTBkg.SetPointError(1, nP, float(m.group(3)))
                else:
                    bkgList.append((head, float(m.group(2)), float(m.group(3))))
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
        self.hNBkg = TH1F(self.tag+"hBkg","hBkg",nBkgC, 0., self.nPoint)
        ibin = 1

        keys = range(len(bkgList))
        if self.ref:
            x1 = self.ref.hNBkg.GetXaxis()
            ## Get the tag list in current list
            ltemp = [x[0] for x in bkgList]
            print ltemp
            print [x1.GetBinLabel(i+1) for i in keys]

            ## find their indices in the refernece
            keys = [ltemp.index(x1.GetBinLabel(i+1)) for i in keys]
        for ix in keys:
            x = bkgList[ix]
            self.hNBkg.SetBinContent(ibin, x[1])
            self.hNBkg.SetBinError(ibin, x[2])
            self.hNBkg.GetXaxis().SetBinLabel(ibin, x[0])
            ibin += 1

            if self.min1 > x[1]-x[2]: self.min1 = x[1]-x[2]


class ComparerX:
    def __init__(self, tag='test'):
        self.tag = tag
        self.rx1 = None
        self.rx2 = None
        self.channels = ['ee_1','ee_2','emu_1','emu_2','mumu_1','mumu_2']
        self.dir0 = '/home/dzhang/links/eosOther/cloShared/save/'
        self.sDir = sDir
        self.sTag = sTag
        self.autoSave = sDirectly

    def compareChans(self, chans = ['ee_1','ee_2','emu_1','emu_2','mumu_1','mumu_2']):
        for chan in chans:
            self.showCompare(chan)
    def showCompare(self,chan):
        r1 = OptResults(self.dir0+self.rx1[0]+'/significance/SR_SS_'+chan+'_opt_0.txt','r1_')
        r2 = OptResults(self.dir0+self.rx2[0]+'/significance/SR_SS_'+chan+'_opt_0.txt', 'r2_', r1)
        lg1 = self.rx1[1]
        lg2 = self.rx2[1]

        hfirst = r1.hNBkg
        hfirst.Draw("hist")

        r1.gTBkg.SetFillStyle(3004)
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
        r2.gZnSig.Draw("Psame")
        r1.gZnSig.SetMarkerColor(kGray)
        r2.gZnSig.SetMarkerColor(kGray)
        r2.gZnSig.SetMarkerStyle(26)

        fun1 = TF1("fun1", "0", 0, 100)
        fun1.SetLineStyle(2)
    #     fun1.SetLineColor(kGray)
        fun1.Draw("same")

        r1.hNBkg.Draw("histsame")
        r1.hNBkg.Draw("Esame")
        r2.hNBkg.Draw("Esame")
        r2.hNBkg.SetLineColor(2)
        r2.hNBkg.SetMarkerColor(2)
        r2.hNBkg.SetMarkerStyle(26)


        max1 = max(r1.gNsig.GetHistogram().GetMaximum(), r2.gNsig.GetHistogram().GetMaximum())
        min1 = min(r1.min1, r2.min1)
        dm = 0.02*(max1-min1)
        hfirst.GetYaxis().SetRangeUser(min1-dm,max1+dm)
        hfirst.GetYaxis().SetTitle("Events")
        

        lg = TLegend(0.7,0.8,0.9,0.93)
        lg.SetHeader(chan)
        lg.SetFillStyle(0)
        lg.AddEntry(r1.hNBkg,lg1,'lp')
        lg.AddEntry(r2.hNBkg,lg2,'lp')
        lg.Draw()

        gPad.Update()
        waitRootCmd(self.sDir+self.sTag+chan,self.autoSave)

def test():
    cx1 = ComparerX('tes1')
    cx1.rx1 = ('2.8_1D','2.8-1D')
    cx1.rx2 = ('2.8_1D_loose','Loose')
    cx1.compareChans()

funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    for fun in funlist: print fun()
