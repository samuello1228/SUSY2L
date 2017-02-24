#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory, get_default_fig_dir
gROOT.LoadMacro('Plot.C+')
from ROOT import Plot, Sample, SampleGroup
from sample_info import sampleInfo
from anaUtil import loadSamples 

sDir = get_default_fig_dir()
sTag = 'test_'
sDirectly = False
if gROOT.IsBatch(): sDirectly = True

def useStyle(sample, style, size=None):
    for i in range(size if size else len(style)): sample.style[i]=style[i]

def addSample(dsname,treename='evt2l'):
    tag = dsname.split()[3]
    sx = Sample(tag, tag, tag, tag)
    sx.tree1 = TChain(treename)
    sx.tree1.Add(dir0+dsname+'/*.root*')
#     sx.weight = GetWeight()

def checkCompare():
    samples,sampleG = loadSamples('dataset_v101.list', None, '/net/s3_data_home/dzhang/links/SAMPLES/R20/susyNtuple/AnalysisBase-02-04-25x/')

    sI1 = sampleInfo()
    sI1.Unit_fb = True
    sI1.loadAll()

    Lumi = 3.2;
    ### get weights
    for s in samples:
        sx = samples[s]
        if sx.name.find('ds00') == -1:
            sx.weight = sI1.getXSec('id=='+sx.info)*Lumi/sx.getStatWeight()
        print sx.name, sx.getStatWeight(), sx.tree1.GetEntries(), sx.weight

    for g in sampleG:
        print g
        if g.find('dM') == -1: continue
        gx = sampleG[g]
        gx.setUpOwnChain(TChain('evt2l'))
        sx = gx.sampleList[0]
        gx.weight = sI1.getXSec('id=='+sx.info)*Lumi/gx.getStatWeight()
        print g, sI1.getXSec('id=='+sx.info), Lumi, gx.getStatWeight()
#     return

    ### data
    dir1 = '/home/dzhang/work/bsmSearch/ewSUSY/analysis/v20_7b/SUSY2L/code/lRun/output/'
    s1 = Sample('data','data_','Pseudo-data','Data')
    s1.tree1 = TChain('evt2l')
    s1.tree1.Add(dir1+'v10.0.data/fetch/data-myOutput/*.root')

    ### set styles
    useStyle(sampleG['VV'],(3,1,1,3,20,1,3))
    useStyle(sampleG['Vgamma'],(4,1,1,4,20,1,4))
#     useStyle(sampleG['Top'],(5,1,1,5,20,1,5))
    useStyle(sampleG['ZJet'],(2,1,1,2,20,1,2))
    useStyle(sampleG['dM20'],(3,2,2))
    useStyle(sampleG['dM100'],(5,9,2))


    p1 = Plot()
    p1.mode = 1
    p1.sData = s1
#     p1.sSM.push_back(sampleG['Top'])
#     p1.sSM.push_back(sampleG['VV'])
#     p1.sSM.push_back(sampleG['Vgamma'])
    p1.sSM.push_back(sampleG['ZJet'])
    p1.sSig.push_back(sampleG['dM20'])
    p1.sSig.push_back(sampleG['dM100'])
    p1.mCut = '(evt.weight*evt.pwt*evt.ElSF*evt.MuSF)*(evt.flag==2&&l12.m>10&&(sig.trigCode&sig.trigMask)!=0&&evt.run<296939)'

#     p1.test()
#     p1.showPlot('leps.pt', TH1F('h1','lep_pt_logy;Leading lepton p_{T} [GeV];Events / 5 GeV',60, 0, 300))
#     waitRootCmd(sDir+sTag+"leps_pt", sDirectly)

    p1.showPlot('leps.eta', TH1F('h1','lep_eta_logy;lepton #eta;Events / 0.1',50, -2.5, 2.5))
    waitRootCmd(sDir+sTag+"leps_eta", sDirectly)


#     p1.showPlot('l12.m', TH1F('h1','l12_m;m_{ll} [GeV];Events / 2 GeV',100, 0, 200))
#     p1.showPlot('l12.m', TH1F('h1','l12_m;m_{ll} [GeV];Events / 2 GeV',100, 0, 20))
#     waitRootCmd(sDir+sTag+"l12_m", sDirectly)


if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    checkCompare()
