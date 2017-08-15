#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory, mkupHistSimple, get_default_fig_dir
from math import sqrt
import glob
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

def getSig(args, show=True, info=None,savename="test"):
    dictA = {}
    for x in args:
        updateDict(dictA, x)
    print dictA

    g2d = TGraph2D()
    i = 0
    for x in dictA:
        g2d.SetPoint(i, float(x[0]),float(x[1]),dictA[x])
        i += 1

    if not show:
        return g2d

    gStyle.SetPadRightMargin(0.16)
    g2d.Draw('colz')
    h1 = g2d.GetHistogram();
    h1.GetXaxis().SetTitle("m_{#tilde{#chi}^{#pm}_{1}/#tilde{#chi}^{0}_{2}} [GeV]");
    h1.GetYaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} [GeV]");
    h1.GetZaxis().SetTitle("Z_{n}");
    h1.GetXaxis().SetRangeUser(0,350)
    h1.GetZaxis().SetRangeUser(0,4)

    lg = TLegend(0.2, 0.65,0.4,0.85)
    lg.SetFillStyle(0)
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
        lg.AddEntry(xt0, 'Z_{n}='+str(lx), 'l')
    lg.Draw()

    if info:
        lt = TLatex()
        lt.DrawLatexNDC(0.2,0.9,info)

    gPad.Update()
    waitRootCmd(sDir+savename, sDirectly)

    return g2d

def test1():
    if len(sys.argv)<2:
        print 'Usage:'
        print sys.argv[0], 'FILENAMES'
        print 'ERROR: At least 1 arguments is needed.'

    getSig(sys.argv[1:])

def test():
#     gr1 = getSig(glob.glob("SR*.txt"),True,"Jet |#eta|<2.4","JetEta_2p4")
#     gr2 = getSig(glob.glob("bugFix_Aug04/SR*.txt"),True,"Jet |#eta|<2.8","JetEta_2p8")

    show = False
    gr1 = getSig(glob.glob("SR*.txt"),show,"Jet |#eta|<2.4","JetEta_2p4")
    gr2 = getSig(glob.glob("bugFix_Aug04/SR*.txt"),show,"Jet |#eta|<2.8","JetEta_2p8")

#     c1 = TCanvas()
#     c1.Divide(2)
#     c1.cd(1)
#     gr1.Draw('colz')
#     c1.cd(2)
#     gr2.Draw('colz')
#     waitRootCmd()



#     gStyle.SetPadRightMargin(0.16)
    gr1.Draw()
    h1 = gr1.GetHistogram();
    gr2.Draw()
    h1 = gr2.GetHistogram();
    h1.Draw("axis")
    h1.GetXaxis().SetTitle("m_{#tilde{#chi}^{#pm}_{1}/#tilde{#chi}^{0}_{2}} [GeV]");
    h1.GetYaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} [GeV]");
    h1.GetZaxis().SetTitle("Z_{n}");
    h1.GetXaxis().SetRangeUser(0,350)
    h1.GetZaxis().SetRangeUser(0,4)

    fun1 = TF1("linex","x-125.5", 0, 350);
    fun1.SetLineColor(kGray)
    fun1.SetLineWidth(2)
    fun1.SetLineStyle(8)
    fun1.Draw("same");
#     waitRootCmd()

    lg = TLegend(0.2, 0.65,0.4,0.85)
    lg.SetFillStyle(0)
    lines = {1:2, 2:5, 3:9} # Zn, line Style

    lg.AddEntry(None, "|#eta|<2.4","")
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
    lg.AddEntry(None, "|#eta|<2.8","")
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

    lg.Draw()
    gPad.Update()
    waitRootCmd()

funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    for fun in funlist: print fun()
