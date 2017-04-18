#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory, get_default_fig_dir
gROOT.LoadMacro('Sample.C+')
gROOT.LoadMacro('Plot.C+')
from ROOT import Plot, Sample, SampleGroup
from sample_info import sampleInfo
from anaUtil import loadSamples, anaSpace, loadSamplesL 
from array import array

gROOT.ProcessLine("Sample* noneSample(nullptr);")
from ROOT import noneSample;

sDir = get_default_fig_dir()
sTag = 'test_'
sDirectly = False
if gROOT.IsBatch(): sDirectly = True

def useStyle(sample, style, size=None):
    for i in range(size if size else len(style)): sample.style[i]=style[i]
    print sample.name,
    for i in range(len(style)): print sample.style[i]
    print

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


def buildAna(anaFileName):
    p1 = Plot(anaFileName)
    p1.mode = 1

    ### get analysis space
    write = True
    if write:
#         samples,sampleG = loadSamplesL('dataset_v130.list', None, '/home/dzhang/work/bsmSearch/ewSUSY/analysis/v20_7b/SUSY2L/code/lRun/output/v15.0.MC/data-myOutput/')
        samples,sampleG = loadSamplesL('dataset_v130.list', None, '/home/dzhang/work/bsmSearch/ewSUSY/analysis/v20_7b/SUSY2L/code/lRun/output/v16.0.MC/data-myOutput/')
#         samples1,sampleG1 = loadSamplesL('dataset_v130.list', None, '/home/dzhang/work/bsmSearch/ewSUSY/analysis/v20_7b/SUSY2L/code/lRun/output/v13.0.MC.1/data-myOutput/')

        for s in samples:
            if samples[s].getStatWeight() < 1:
                print s
#                 samples[s] = samples1[s]
#                 print samples1[s].getStatWeight()
        for s in sampleG:
            if sampleG[s].getStatWeight() < 1:
                print s
#                 sampleG[s] = sampleG1[s]
#                 print sampleG1[s].getStatWeight()

        print '////////////////////////////'
        sI1 = sampleInfo()
        sI1.Unit_fb = True
        sI1.loadAll()

        Lumi = 36.07;
        ### get weights
        for s in samples:
            sx = samples[s]
            print '-------------------'
            print sx.name, sx.tree1.GetEntries(), sx.getStatWeight()
            print '-------------------'
            if sx.name.find('ds00') == -1:
                sx.weight = sI1.getXSec('id=='+sx.tag)*Lumi/sx.getStatWeight()
            print sx.name, sx.getStatWeight(), sx.tree1.GetEntries(), sx.weight

        for g in sampleG:
            print g
            if g.find('dM') == -1: continue
            gx = sampleG[g]
            gx.setUpOwnChain(TChain('evt2l'))
            sx = gx.sampleList[0]
            gx.weight = sI1.getXSec('id=='+sx.tag)*Lumi/gx.getStatWeight()
            print g, sI1.getXSec('id=='+sx.tag), Lumi, gx.getStatWeight()

        ### data
        dir1 = '/home/dzhang/work/bsmSearch/ewSUSY/analysis/v20_7b/SUSY2L/code/lRun/output/'
        s1 = Sample('data','data_','Pseudo-data','Data')
        s1.tree1 = TChain('evt2l')
        s1.tree1.Add(dir1+'v16.0.data/fetch/data-myOutput/*.root')

        ### set styles
        useStyle(sampleG['VV'],(3,1,1,3,20,1,3))
        useStyle(sampleG['Wgamma'],(4,1,1,4,20,1,4))
        useStyle(sampleG['Top'],(5,1,1,5,20,1,5))
        useStyle(sampleG['ZJet'],(2,1,1,2,20,1,2))
        useStyle(sampleG['dM20'],(1,2,2))
        useStyle(sampleG['dM100'],(6,9,2))
#         useStyle(sampleG['Zgamma'],(8,9,2))
        useStyle(sampleG['Zgamma'],(0,9,2))

        p1.sData = s1
        p1.sSM.push_back(sampleG['Wgamma'])
        p1.sSM.push_back(sampleG['VV'])
        p1.sSM.push_back(sampleG['Top'])
        p1.sSM.push_back(sampleG['ZJet'])
        p1.sSig.push_back(sampleG['dM20'])
        p1.sSig.push_back(sampleG['dM100'])
        p1.sSig.push_back(sampleG['Zgamma'])
        p1.saveSamples()


    sels = []
    sels += [(("flag1OR7","(evt.flag==1||evt.flag==7)"), 'ee all','ee_all')]
    sels += [(("flag2OR8","(evt.flag==2||evt.flag==8)"), '#mu#mu all','mumu_all')]
    sels += [(("flag3","evt.flag==3"), 'e#mu noISR','emu_noISR')]
    sels += [(("flag9","evt.flag==9"), 'e#mu ISR','emu_ISR')]
    sels += [(("flag1","evt.flag==1"), 'ee noISR','ee_noISR')]
    sels += [(("flag2","evt.flag==2"), '#mu#mu noISR','mumu_noISR')]
    sels += [(("flag7","evt.flag==7"), 'ee ISR','ee_ISR')]
    sels += [(("flag8","evt.flag==8"), '#mu#mu ISR','mumu_ISR')]
    sels += [(("flag4","evt.flag==4"), 'ee no ISR, SS','ee_noISR_SS')]
    sels += [(("flag5","evt.flag==5"), '#mu#mu no ISR, SS','mumu_noISR_SS')]
    sels += [(("flag6","evt.flag==6"), 'e#mu no ISR, SS','emu_noISR_SS')]
    sels += [(("flag10","evt.flag==10"), 'ee  ISR, SS','ee_ISR_SS')]
    sels += [(("flag11","evt.flag==11"), '#mu#mu  ISR, SS','mumu_ISR_SS')]
    sels += [(("flag12","evt.flag==12"), 'e#mu  ISR, SS','emu_ISR_SS')]

    eTag = "_B"
    eSel = "&&(sig.trigCode&sig.trigMask)!=0&&leps[0].ElChargeID==1&&leps[1].ElChargeID==1&&l12.m>15"
    elListA = [("",""),(eTag,eSel)]
    for el in elListA:
        for x in sels:
            p1.useEntryList(x[0][0]+el[0],True,x[0][1]+el[1], True)

def test_dev(anaName="ana1_Apr04b.root"):
    p1 = Plot(anaName)
    p1.mode = 1

    h1 =  TH1F('h1','l0_pt;Leading lepton p_{T} [GeV]; Events / 5 GeV',100,0,500)

    ptBins = [0, 5, 10, 15, 20, 25, 30, 40, 50, 100, 150, 200, 500]
    h2 = TH1F('h2','l0_pt;Leading lepton p_{T} [GeV]; Events / 5 GeV',len(ptBins)-1,array('d',ptBins))
    s1 = p1.getSampleGroup('dM20')
    hx = s1.getHistFromTree("leps[0].pt", h1, "")

    useStyle(s1,(2,1,1,2,20,1,2))
    hx2 = s1.getHistFromTree("leps[0].pt", h2, "")


    hx.Draw()
    hx2.Draw("same")
    waitRootCmd()


def run_test1(anaName="ana1_Apr04b.root"):
    p1 = Plot(anaName, 6)
    p1.mode = 1

    ### get analysis space
    write = False
    if write:
        samples,sampleG = loadSamplesL('dataset_v130.list', None, '/home/dzhang/work/bsmSearch/ewSUSY/analysis/v20_7b/SUSY2L/code/lRun/output/v14.0.MC/data-myOutput/')
#         samples1,sampleG1 = loadSamplesL('dataset_v130.list', None, '/home/dzhang/work/bsmSearch/ewSUSY/analysis/v20_7b/SUSY2L/code/lRun/output/v13.0.MC.1/data-myOutput/')

        for s in samples:
            if samples[s].getStatWeight() < 1:
                print s
#                 samples[s] = samples1[s]
#                 print samples1[s].getStatWeight()
        for s in sampleG:
            if sampleG[s].getStatWeight() < 1:
                print s
#                 sampleG[s] = sampleG1[s]
#                 print sampleG1[s].getStatWeight()

        print '////////////////////////////'
        sI1 = sampleInfo()
        sI1.Unit_fb = True
        sI1.loadAll()

        Lumi = 36.07;
        ### get weights
        for s in samples:
            sx = samples[s]
            print '-------------------'
            print sx.name, sx.tree1.GetEntries(), sx.getStatWeight()
            print '-------------------'
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

        proof = ROOT.TProof.Open("lite://")
        ### data
        dir1 = '/home/dzhang/work/bsmSearch/ewSUSY/analysis/v20_7b/SUSY2L/code/lRun/output/'
        s1 = Sample('data','data_','Pseudo-data','Data')
        s1.tree1 = TChain('evt2l')
        s1.tree1.Add(dir1+'v14.0.data/fetch/data-myOutput/*.root')

        ### set styles
        useStyle(sampleG['VV'],(3,1,1,3,20,1,3))
        useStyle(sampleG['Wgamma'],(4,1,1,4,20,1,4))
        useStyle(sampleG['Top'],(5,1,1,5,20,1,5))
        useStyle(sampleG['ZJet'],(2,1,1,2,20,1,2))
        useStyle(sampleG['dM20'],(1,2,2))
        useStyle(sampleG['dM100'],(6,9,2))
#         useStyle(sampleG['Zgamma'],(8,9,2))
        useStyle(sampleG['Zgamma'],(0,9,2))

        p1.sData = s1
        p1.sSM.push_back(sampleG['Wgamma'])
        p1.sSM.push_back(sampleG['VV'])
        p1.sSM.push_back(sampleG['Top'])
        p1.sSM.push_back(sampleG['ZJet'])
        p1.sSig.push_back(sampleG['dM20'])
        p1.sSig.push_back(sampleG['dM100'])
        p1.sSig.push_back(sampleG['Zgamma'])
        p1.saveSamples()
    else:
        s_zJet = p1.getSampleGroup('ZJet')
        zmumuTags = range(364100,364113+1)
        for s in s_zJet.sampleList:
            if int(s.tag) in zmumuTags:
                s.weight *= 1.03

        p1.sData = p1.getSample("data")
        p1.sSM.push_back(p1.getSampleGroup('Wgamma'))
        p1.sSM.push_back(p1.getSampleGroup('VV'))
        p1.sSM.push_back(p1.getSampleGroup('Top'))
        p1.sSM.push_back(s_zJet)
        p1.sSig.push_back(p1.getSampleGroup('dM20'))
        p1.sSig.push_back(p1.getSampleGroup('dM100'))
    #     p1.sSig.push_back(p1.getSampleGroup('Zgamma'))

        s_zgamma = p1.getSampleGroup('Zgamma')
        p1.sSig.push_back(s_zgamma)
        useStyle(s_zgamma,(0,9,2))
#     for x in [(1, 'ee noISR', 'ee_noISR'), (2, '#mu#mu noISR','mumu_noISR'), (3, 'e#mu noISR', 'emu_noISR'), (4, 'ee ISR', 'ee_ISR'), (6, 'e#mu ISR', 'emu_ISR'), (5, '#mu#mu ISR','mumu_ISR')]:
    sels = []
#     sels += [(("flag1OR7","(evt.flag==1||evt.flag==7)"), 'ee all','ee_all')]
    sels += [(("flag2OR8","(evt.flag==2||evt.flag==8)"), '#mu#mu all','mumu_all')]
#     sels += [(("flag3","evt.flag==3"), 'e#mu noISR','emu_noISR')]
#     sels += [(("flag9","evt.flag==9"), 'e#mu ISR','emu_ISR')]
    sels += [(("flag1","evt.flag==1"), 'ee noISR','ee_noISR')]
#     sels += [(("flag2","evt.flag==2"), '#mu#mu noISR','mumu_noISR')]
    sels += [(("flag7","evt.flag==7"), 'ee ISR','ee_ISR')]
    sels += [(("flag8","evt.flag==8"), '#mu#mu ISR','mumu_ISR')]
    sels += [(("flag4","evt.flag==4"), 'ee no ISR, SS','ee_noISR_SS')]
    sels += [(("flag5","evt.flag==5"), '#mu#mu no ISR, SS','mumu_noISR_SS')]
    sels += [(("flag6","evt.flag==6"), 'e#mu no ISR, SS','emu_noISR_SS')]
    sels += [(("flag10","evt.flag==10"), 'ee  ISR, SS','ee_ISR_SS')]
    sels += [(("flag11","evt.flag==11"), '#mu#mu  ISR, SS','mumu_ISR_SS')]
    sels += [(("flag12","evt.flag==12"), 'e#mu  ISR, SS','emu_ISR_SS')]


    p1.mode = 1
#     tTag = 'gradient_m2_'
    tTag = 'cID_m2_'
    for x in sels:
#         p1.mCut = '(evt.weight*evt.pwt*evt.ElSF*evt.MuSF)*((sig.trigCode&sig.trigMask)!=0&&leps[0].ElChargeID==1&&leps[1].ElChargeID==1)'
        p1.mCut = '(evt.weight*evt.pwt*evt.ElSF*evt.MuSF)*((sig.trigCode&sig.trigMask)!=0&&leps[0].ElChargeID==1&&leps[1].ElChargeID==1&&l12.m>15)'
#         p1.mCut = '(evt.weight*evt.pwt*evt.ElSF*evt.MuSF)*((sig.trigCode&sig.trigMask)!=0&&leps[0].ElChargeID==1&&leps[1].ElChargeID==1&&(leps[0].isoPass&0x2)!=0&&(leps[1].isoPass&0x2)!=0)'
        if x[0][1].split("==")[1] in [4, 5, 6, 10, 11, 12]:
            if x[0][0] in ['flag4', 'flag10']: ### show Z control region
                p1.sData = p1.getSample("data")
                p1.sData.sCuts = "(l12.m>81&&l12.m<101)"
            else: p1.sData = noneSample
        else:
            p1.sData = p1.getSample("data")
            p1.sData.sCuts = ''

        p1.useEntryList(x[0][0],False,x[0][1], False)
        p1.showInfo = x[1]
        sTag = tTag+'preSel_'+x[2]

#         p1.showPlot('leps[0].mT', TH1F('h1','lep0_mT_logy;Leading lepton m_{T} [GeV];Events / 2 GeV',100, 0, 200))
#         waitRootCmd(sDir+sTag+"mT0", sDirectly)
#         p1.showPlot('leps[1].mT', TH1F('h1','lep1_mT_logy;Sub-leading lepton m_{T} [GeV];Events / 2 GeV',100, 0, 200))
#         waitRootCmd(sDir+sTag+"mT1", sDirectly)
#         p1.showPlot('l12.m', TH1F('h1','l12_m_logy;m_{ll} [GeV];Events / 2 GeV',100, 0, 200))
#         waitRootCmd(sDir+sTag+"mll", sDirectly)
#         p1.showPlot('leps[1].eta', TH1F('h1','l1_eta_logy;Subleading lepton #eta;Events / 0.1',60, -3, 3))
#         waitRootCmd(sDir+sTag+"_l1_eta", sDirectly)
        p1.showPlot('Length$(jets)', TH1F('h1','l1_nJets_logy;Number of jets;Events',15, 0, 15))
        waitRootCmd(sDir+sTag+"_nJet", sDirectly)
        
        jetPtBins = [i*10 for i in range(10)] + [100+i*20 for i in range(10)] + [300+i*50 for i in range(6)] + [600+i*100 for i in range(9)]
        p1.showPlot('jets[0].pt', TH1F('h1','j0_pt_logy;Leading jet p_{T} [GeV];Events / 10 GeV',len(jetPtBins)-1, array('f',jetPtBins)))
#         p1.showPlot('jets[0].pt', TH1F('h1','j0_pt_logy;Leading jet p_{T} [GeV];Events',100, 0, 200))
        waitRootCmd(sDir+sTag+"_j0_pt", sDirectly)
#         p1.showPlot('leps[1].phi', TH1F('h1','l1_eta_logy;Subleading lepton #phi [rad];Events / 0.1 rad',64, -3.2, 3.2))
#         waitRootCmd(sDir+sTag+"_l1_phi", sDirectly)
#         p1.showPlot('l12.pt', TH1F('h1','l12_pt_logy;p^{ll}_{T} [GeV];Events / 2 GeV',100, 0, 200))
#         waitRootCmd(sDir+sTag+"llpt", sDirectly)
#         p1.showPlot('sig.Met', TH1F('h1','Met_logy;E_{T}^{Miss} [GeV];Events / 2 GeV',100, 0, 200))
#         waitRootCmd(sDir+sTag+"Met", sDirectly)
    proof.Print()
if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
#     run_test()
    wa='ana_Apr18a.root'
#     run_test1("ana_Apr17a.root")
#     test_dev("ana_Apr17a.root")
#     buildAna("ana_Apr17a.root")
    buildAna(wa)
    run_test1(wa)
#     checkCompare()
