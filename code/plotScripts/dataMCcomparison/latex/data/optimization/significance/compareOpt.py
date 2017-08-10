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

    def __init__(self, txtfile=None, tag='test'):
        self.tag = tag
        self.hNBkg = None 
        self.gNsig = TGraphErrors()
        self.gZnSig = TGraph()
        self.gTBkg = None
        self.min1 = 0.
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
                    self.gTBkg.SetPoint(0, nP, float(m.group(2)))
                    self.gTBkg.SetPointError(0, nP, float(m.group(3)))
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
        for x in bkgList:
            self.hNBkg.SetBinContent(ibin, x[1])
            self.hNBkg.SetBinError(ibin, x[2])
            self.hNBkg.GetXaxis().SetBinLabel(ibin, x[0])
            ibin += 1

            if self.min1 > x[1]-x[2]: self.min1 = x[1]-x[2]

def test():
    dir0 = '/home/dzhang/links/eosOther/cloShared/save/'
    print dir0

    r1 = OptResults(dir0+'2.8_2D/significance/SR_SS_ee_1_opt_0.txt','r1_')
    r2 = OptResults(dir0+'2.4_2D/significance/SR_SS_ee_1_opt_0.txt', 'r2_')
    lg1 = "2.8_2D"
    lg2 = "2.4_2D"

    hfirst = r1.hNBkg
    r1.hNBkg.Draw("hist")
    r2.hNBkg.Draw("Esame")
    r2.hNBkg.SetLineColor(2)
    r2.hNBkg.SetMarkerColor(2)
    r2.hNBkg.SetMarkerStyle(26)


    r1.gNsig.Draw("Psame")
    r2.gNsig.Draw("Psame")
    r2.gNsig.SetLineColor(2)
    r2.gNsig.SetMarkerColor(2)
    r2.gNsig.SetMarkerStyle(26)

    r1.gTBkg.Draw("L")
    r2.gTBkg.Draw("Lsame")
    r2.gTBkg.SetLineColor(2)
    r2.gTBkg.SetMarkerColor(2)
    r2.gTBkg.SetMarkerStyle(26)


    max1 = max(r1.gNsig.GetHistogram().GetMaximum(), r2.gNsig.GetHistogram().GetMaximum())
    min1 = min(r1.min1, r2.min1)
    dm = 0.02*(max1-min1)
    hfirst.GetYaxis().SetRangeUser(min1-dm,max1+dm)
    

    lg = TLegend(0.7,0.8,0.9,0.93)
    lg.SetFillStyle(0)
    lg.AddEntry(r1.hNBkg,lg1,'lp')
    lg.AddEntry(r2.hNBkg,lg2,'lp')
    lg.Draw()

    gPad.Update()
    waitRootCmd()
funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    for fun in funlist: print fun()
