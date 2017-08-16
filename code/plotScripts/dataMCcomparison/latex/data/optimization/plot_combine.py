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
    waitRootCmd(sDir+savename, sDirectly)

    return g2d

def compare2():
    show = False
#     gr1 = getSig(glob.glob("significance_2.8_1D/SR*.txt"),show,"Baseline","Baseline")
    gr1 = getSig(glob.glob("significance_2.8_1D_avg/SR*.txt"),show,"BaselineAvg","BaselineAvg")
    gr2 = getSig(glob.glob("significance_2.8_1D_loose_avg/SR*.txt"),show,"LooseAvg","LooseAvg")
    showCompare(gr1,gr2)

def showCompare(gr1, gr2):
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
    waitRootCmd()


def test1():
    if len(sys.argv)<2:
        print 'Usage:'
        print sys.argv[0], 'FILENAMES'
        print 'ERROR: At least 1 arguments is needed.'

    getSig(sys.argv[1:])

def test():
    gr1 = getSig(glob.glob("significance_2.8_1D_loose/SR*opt_0.txt"),True,"Muon |#eta|<2.7","MuEta_2p7")
#     gr2 = getSig(glob.glob("bugFix_Aug04/SR*.txt"),True,"Jet |#eta|<2.8","JetEta_2p8")


funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
#     for fun in funlist: print fun()
    compare2()
