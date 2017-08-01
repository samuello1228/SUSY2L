#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory, mkupHistSimple, get_default_fig_dir
from math import sqrt
funlist=[]

sDir = get_default_fig_dir()
sTag = 'test_'
sDirectly = False
if gROOT.IsBatch(): sDirectly = True


def updateDict(d, x):
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
        if Zn>0: d[(fs[0],fs[1])] = Zn*Zn + d.get((fs[0],fs[1]), 0.)
    for x in d: d[x] = sqrt(d[x])


def test():
    if len(sys.argv)<2:
        print 'Usage:'
        print sys.argv[0], 'FILENAMES'
        print 'ERROR: At least 1 arguments is needed.'

    dir0 = os.getenv('SAMPLEDIR_LAMB')
    print dir0

    gStyle.SetPadRightMargin(0.16)

    dictA = {}
    for x in sys.argv[1:]:
        updateDict(dictA, x)
    print dictA

    g2d = TGraph2D()
    i = 0
    for x in dictA:
        g2d.SetPoint(i, float(x[0]),float(x[1]),dictA[x])
        i += 1

    g2d.Draw('colz')
    h1 = g2d.GetHistogram();
    h1.GetXaxis().SetTitle("m_{#tilde{#chi}^{#pm}_{1}/#tilde{#chi}^{0}_{2}} [GeV]");
    h1.GetYaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} [GeV]");
    h1.GetZaxis().SetTitle("Z_{n}");

    grL = g2d.GetContourList(1);

    next1 = TIter(grL);
    while True:
        obj = next1()
        if obj == None: break
        obj.SetLineWidth(2);
        obj.SetLineStyle(2);
        obj.Draw("same");
        obj.SetName("s0");

    xt0 = gPad.GetPrimitive("s0");


    grL = g2d.GetContourList(2);
    next2 = TIter(grL);
    while True:
        obj = next2()
        if obj == None: break
        obj.SetLineWidth(2);
        obj.SetLineStyle(9);
        obj.Draw("same");
        obj.SetName("s1");

    xt1 = gPad.GetPrimitive("s1");


    lt = TLatex()
#     lt.DrawLatexNDC(0.2,0.9,"New")
    lt.DrawLatexNDC(0.2,0.9,"Old")

    waitRootCmd()
funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    for fun in funlist: print fun()
