#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory, get_default_fig_dir
gROOT.LoadMacro('Plot.C+')
from ROOT import Plot, Sample, SampleGroup
from sample_info import sampleInfo
from plotMaker import plotMaker, makeChain 

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

def loadSamples(dslist, dsname=None):
    dir0 = '/net/s3_data_home/dzhang/links/SAMPLES/R20/susyNtuple/AnalysisBase-02-04-26-da7031fc/'
    ### get the samples
    samples = {}
    sampleG = {}
    with open(dslist,'read') as f1:
        grp = None
        for line in f1.readlines():
            line = line.rstrip()
            if len(line)==0 or line[0]=='#': continue
            if line.find(':')!=-1:
                tag = line[:line.find(':')]
                if dsname and tag!=dsname: continue

                grp = SampleGroup(tag,tag,tag,tag)
                print 'group',grp.name,'created!'
                continue
            if line.find('END')!=-1:
                if grp:
                    print "End of", grp.name
                    if dsname: return grp

                    sampleG[grp.name] = grp
                    grp = None
                continue
            
            ### make a sample
            fs = line.split('.')
            name = 'ds'+fs[4]
            s1 = Sample(name,fs[4],fs[4],fs[5])
            samples[s1.name] = s1
            s1.tree1 = TChain('evt2l')
            s1.tree1.Add(dir0+line+'/*.root*')
            if grp:
                grp.sampleList.push_back(s1)
    return samples, sampleG

def compare():
    oldData = loadSamples("datasets.list","Data")
    oldData.setUpOwnChain(TChain('evt2l'))

    dir1 = '/home/dzhang/work/bsmSearch/ewSUSY/analysis/v20_7b/SUSY2L/code/lRun/output/'
    s1 = Sample('newdata','newdata_','Pseudo-data2','New data')
    s1.tree1 = TChain('evt2l')
    s1.tree1.Add(dir1+'v10.0.data/fetch/data-myOutput/*.root')
    useStyle(s1, (2,1,1,2,20,1,2))

    p1 = Plot()
    p1.mode = 1
    p1.sData = oldData
    p1.sSM.push_back(s1)
    p1.mCut = ''

#     p1.test()
#     p1.showPlot('leps.pt', TH1F('h1','lep_pt_logy;Leading lepton p_{T} [GeV];Events / 3 GeV',100, 0, 300))
#     waitRootCmd(sDir+sTag+"leps_pt", sDirectly)

    p1.showPlot('l12.m', TH1F('h1','l12_pt;m_{ll} [GeV];Events / 2 GeV',100, 0, 200))
    waitRootCmd(sDir+sTag+"l12_m", sDirectly)

def checkTrigger():
    '''Check the effect of trigger requirement. Done with data only, as MC jobs are not finished yet'''
    dir1 = '/home/dzhang/work/bsmSearch/ewSUSY/analysis/v20_7b/SUSY2L/code/lRun/output/'
    ch1 = makeChain(dir1+'v10.0.data/fetch/data-myOutput/*.root','evt2l')

    pm1 = plotMaker(ch1,"")
    pm1.cmsInfo = '13 TeV, 36 fb^{-1}'

#     for channel in [(1, 'ee'), (2,'mumu'), (3, 'emu')]:
#     for channel in [(2,'mumu'), (3, 'emu')]:
    for channel in [(3, 'emu'),(4, 'ee FSR')]:
        pm1.cut0 = 'evt.flag==%d'%channel[0]
        pm1.sampleInfo = channel[1]
        pm1.sTag = "trigCheck_"+channel[1]+"_"

        cuts = [("cut0", "", "No Cut", 1)]
        cuts += [("cut1", "(sig.trigCode&sig.trigMask)!=0", "Trigger", 1)]
        pm1.compareCuts(('l12.m', TH1F('h1','l12_pt;m_{ll} [GeV];Events / 2 GeV',100, 0, 200)), cuts, '', False)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
#     compare()
    checkTrigger()
