#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory, get_default_fig_dir
gROOT.LoadMacro('Sample.C+')
gROOT.LoadMacro('Plot.C+')
from ROOT import Plot, Sample, SampleGroup
from sample_info import sampleInfo
from anaUtil import loadSamples, anaSpace 

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

def readSampleFromFile():
    x = Sample(filename)
    print x.tree1.GetEntries(), x.weight, x.crossSection, x.Weight

    

def checkCompare():
    samples,sampleG = loadSamples('dataset_v101.list', None, '/net/s3_data_home/dzhang/links/SAMPLES/R20/susyNtuple/AnalysisBase-02-04-25x/')

    sI1 = sampleInfo()
    sI1.Unit_fb = True
    sI1.loadAll()

    Lumi = 36.07;
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
    s1.tree1.Add(dir1+'v11.1.data/fetch/data-myOutput/*.root')

    ### set styles
    useStyle(sampleG['VV'],(3,1,1,3,20,1,3))
    useStyle(sampleG['Wgamma'],(4,1,1,4,20,1,4))
    useStyle(sampleG['Top'],(5,1,1,5,20,1,5))
    useStyle(sampleG['ZJet'],(2,1,1,2,20,1,2))
    useStyle(sampleG['dM20'],(1,2,2))
    useStyle(sampleG['dM100'],(6,9,2))
    useStyle(sampleG['Zgamma'],(8,9,2))

    p1 = Plot()
    p1.mode = 1
    p1.sData = s1
    p1.sSM.push_back(sampleG['Wgamma'])
    p1.sSM.push_back(sampleG['VV'])
    p1.sSM.push_back(sampleG['Top'])
    p1.sSM.push_back(sampleG['ZJet'])
    p1.sSig.push_back(sampleG['dM20'])
    p1.sSig.push_back(sampleG['dM100'])
    p1.sSig.push_back(sampleG['Zgamma'])
#     p1.mCut = '(evt.weight*evt.pwt*evt.ElSF*evt.MuSF)*(evt.flag==2&&l12.m>10&&(sig.trigCode&sig.trigMask)!=0&&evt.run<296939)'

    for x in [(1, 'ee noISR', 'ee_noISR'), (2, '#mu#mu noISR','mumu_noISR'), (3, 'e#mu noISR', 'emu_noISR'), (4, 'ee ISR', 'ee_ISR'), (6, 'e#mu ISR', 'emu_ISR'), (5, '#mu#mu ISR','mumu_ISR')]:
#     for x in [(2, '#mu#mu noISR','mumu_noISR'), (1, 'ee noISR', 'ee_noISR'), (3, 'e#mu noISR', 'emu_noISR'), (4, 'ee ISR', 'ee_ISR'), (6, 'e#mu ISR', 'emu_ISR'), (5, '#mu#mu ISR','mumu_ISR')]:
#     for x in [(3, 'e#mu noISR', 'emu_noISR'),(4, 'ee ISR', 'ee_ISR'), (6, 'e#mu ISR', 'emu_ISR'), (5, '#mu#mu ISR','mumu_ISR')]:
        p1.mCut = '(evt.weight*evt.pwt*evt.ElSF*evt.MuSF)*(evt.flag=='+str(x[0])+'&&(sig.trigCode&sig.trigMask)!=0)'
        p1.showInfo = x[1]
        sTag = 'preSel_'+x[2]

        p1.showPlot('l12.m', TH1F('h1','l12_m_logy;m_{ll} [GeV];Events / 2 GeV',100, 0, 200))
        waitRootCmd(sDir+sTag+"l12_m", sDirectly)


#     p1.mCut = '(evt.weight*evt.pwt*evt.ElSF*evt.MuSF)*(evt.flag==2&&l12.m>81&&l12.m<101)'

#     p1.test()
#     p1.showPlot('leps.pt', TH1F('h1','lep_pt_logy;Leading lepton p_{T} [GeV];Events / 5 GeV',60, 0, 300))
#     waitRootCmd(sDir+sTag+"leps_pt", sDirectly)
#     p1.showPlot('sig.trigMask', TH1F('h1','sig_trigMask;Trigger Mask;Events',60, 0, 70000000))
#     waitRootCmd(sDir+sTag+"sig_trigMask", sDirectly)


#     p1.showPlot('leps.eta', TH1F('h1','lep_eta_logy;lepton #eta;Events / 0.1',50, -2.5, 2.5))
#     waitRootCmd(sDir+sTag+"leps_eta", sDirectly)


    p1.showPlot('l12.m', TH1F('h1','l12_m_logy;m_{ll} [GeV];Events / 2 GeV',100, 0, 200))
#     p1.showPlot('l12.m', TH1F('h1','l12_m;m_{ll} [GeV];Events / 2 GeV',100, 0, 20))
    waitRootCmd(sDir+sTag+"l12_m", sDirectly)


def run_test():
    ### get analysis space
    write = False
    if write:
        samples,sampleG = loadSamples('dataset_v101.list', None, '/net/s3_data_home/dzhang/links/SAMPLES/R20/susyNtuple/AnalysisBase-02-04-25x/')

        sI1 = sampleInfo()
        sI1.Unit_fb = True
        sI1.loadAll()

        Lumi = 36.07;
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

        ### data
        dir1 = '/home/dzhang/work/bsmSearch/ewSUSY/analysis/v20_7b/SUSY2L/code/lRun/output/'
        s1 = Sample('data','data_','Pseudo-data','Data')
        s1.tree1 = TChain('evt2l')
        s1.tree1.Add(dir1+'v11.1.data/fetch/data-myOutput/*.root')

        ### set styles
        useStyle(sampleG['VV'],(3,1,1,3,20,1,3))
        useStyle(sampleG['Wgamma'],(4,1,1,4,20,1,4))
        useStyle(sampleG['Top'],(5,1,1,5,20,1,5))
        useStyle(sampleG['ZJet'],(2,1,1,2,20,1,2))
        useStyle(sampleG['dM20'],(1,2,2))
        useStyle(sampleG['dM100'],(6,9,2))
        useStyle(sampleG['Zgamma'],(8,9,2))

    p1 = Plot("ana_test1.root")
    p1.mode = 1

    if write:
        p1.sData = s1
        p1.sSM.push_back(sampleG['Wgamma'])
        p1.sSM.push_back(sampleG['VV'])
        p1.sSM.push_back(sampleG['Top'])
        p1.sSM.push_back(sampleG['ZJet'])
        p1.sSig.push_back(sampleG['dM20'])
        p1.sSig.push_back(sampleG['dM100'])
        p1.sSig.push_back(sampleG['Zgamma'])
        p1.saveSamples()
        return


    p1.sData = p1.getSample("data")
    p1.sSM.push_back(p1.getSampleGroup('Wgamma'))
    p1.sSM.push_back(p1.getSampleGroup('VV'))
    p1.sSM.push_back(p1.getSampleGroup('Top'))
    p1.sSM.push_back(p1.getSampleGroup('ZJet'))
    p1.sSig.push_back(p1.getSampleGroup('dM20'))
    p1.sSig.push_back(p1.getSampleGroup('dM100'))
    p1.sSig.push_back(p1.getSampleGroup('Zgamma'))

#     for x in [(1, 'ee noISR', 'ee_noISR'), (2, '#mu#mu noISR','mumu_noISR'), (3, 'e#mu noISR', 'emu_noISR'), (4, 'ee ISR', 'ee_ISR'), (6, 'e#mu ISR', 'emu_ISR'), (5, '#mu#mu ISR','mumu_ISR')]:
    for x in [("flag2", '#mu#mu noISR','mumu_noISR')]:
        p1.useEntryList(x[0],False,"evt.flag==2")
        p1.mCut = '(evt.weight*evt.pwt*evt.ElSF*evt.MuSF)*((sig.trigCode&sig.trigMask)!=0)'
        p1.showInfo = x[1]
        sTag = 'preSel_'+x[2]
        p1.showPlot('l12.m', TH1F('h1','l12_m_logy;m_{ll} [GeV];Events / 2 GeV',100, 0, 200))
        waitRootCmd(sDir+sTag+"l12_m", sDirectly)


if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    run_test()
#     checkCompare()
