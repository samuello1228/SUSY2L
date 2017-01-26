#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory

# gROOT.LoadMacro('Sample.C+')

gROOT.LoadMacro('Plot.C+')
from ROOT import Plot, Sample


funlist=[]

def useStyle(sample, style, size=None):
    for i in range(size if size else len(style)): sample.style[i]=style[i]

def test():
    dir0 = '/home/dzhang/links/outSpace/temp_test/'
    print dir0

    s1 = Sample('data','data_','Pseudo-data','Data')
    s1.tree1 = TChain('evt2l')
    s1.tree1.Add(dir0+'361659.Sherpa_CT10_lllvSFMinusjj_e4584_s2726_r7772_r7676_p2879-*.root')

    s2a = Sample('bkg_VV','VV_','VV background','VV')
    useStyle(s2a, (2,1,1,2,20,1,2))
    s2a.tree1 = TChain('evt2l')
    s2a.tree1.Add(dir0+'361659.Sherpa_CT10_lllvSFMinusjj_e4584_s2726_r7772_r7676_p2879-*.root')
    s2a.weight = 0.2

    s2b = Sample('bkg_Vg','Vg_','Vgamma background','V#gamma')
    useStyle(s2b, (4,1,1,4,20,1,4))
    s2b.tree1 = TChain('evt2l')
    s2b.tree1.Add(dir0+'361659.Sherpa_CT10_lllvSFMinusjj_e4584_s2726_r7772_r7676_p2879-*.root')
    s2b.weight = 0.78

    s3a = Sample('sig_1','sig1_','signal sample 1','Sig1')
    useStyle(s3a, (3,2,2))
    s3a.tree1 = TChain('evt2l')
    s3a.tree1.Add(dir0+'361659.Sherpa_CT10_lllvSFMinusjj_e4584_s2726_r7772_r7676_p2879-*.root')
    s3a.weight = 0.4

    s3b = Sample('sig_2','sig2_','signal sample 2','Sig2')
    useStyle(s3b, (5,9,2))
    s3b.tree1 = TChain('evt2l')
    s3b.tree1.Add(dir0+'361659.Sherpa_CT10_lllvSFMinusjj_e4584_s2726_r7772_r7676_p2879-*.root')
    s3b.weight = 0.5

    p1 = Plot()
    p1.mode = 1
    p1.sData = s1
    p1.sSM.push_back(s2a)
    p1.sSM.push_back(s2b)
    p1.sSig.push_back(s3a)
    p1.sSig.push_back(s3b)

#     p1.test()
    p1.showPlot('leps.pt', TH1F('h1','lep_pt;Leading lepton p_{T} [GeV];Events / 3 GeV',100, 0, 300))
    gPad.Update()
    waitRootCmd()

    p1.showPlot('l12.m', TH1F('h1','l12_pt;m_{ll} [GeV];Events / 2 GeV',100, 0, 200))

    gPad.Update()
    waitRootCmd()
funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    for fun in funlist: print fun()
