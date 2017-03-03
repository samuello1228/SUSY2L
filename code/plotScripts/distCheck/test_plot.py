#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory, get_default_fig_dir
gROOT.LoadMacro('Sample.C+')
gROOT.LoadMacro('Plot.C+')
from ROOT import Plot, Sample, SampleGroup
from sample_info import sampleInfo
from anaUtil import loadSamples
funlist=[]

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

def test3():
    samples,sampleG = loadSamples('dataset_v101.list', None, '/net/s3_data_home/dzhang/links/SAMPLES/R20/susyNtuple/AnalysisBase-02-04-25x/')
    for x in samples:
        samples[x].WriteToFile('Analysis-Test/')
#         samples[x].m_file.Close()
# funlist.append(test3)

def getSample(fname):
    f1 = TFile(fname,'read')
    s1 = f1.Get('ds392875')
    s1.setupFromFile(fname)
    return s1


class anaSpace:
    def __init__(self,fname):
        self.m_file = TFile(fname, 'update')
        self.wSamples = []

    def getSample(self, sname):
        s1 = self.m_file.Get("Sample_"+sname)
        s1.tree1 = self.m_file.Get(sname+"_tree1")
        return s1

    def useEntryList(self, elistname):
        for s in self.wSamples:
            s.SetEntryList(elistname)

    def loadSamples(self):
        samples,sampleG = loadSamples('dataset_v101.list', None, '/net/s3_data_home/dzhang/links/SAMPLES/R20/susyNtuple/AnalysisBase-02-04-25x/')
        for x in samples:
            samples[x].writeToFile(self.m_file)

def test4():
    a1 = anaSpace("test1_ana.root")
    a1.loadSamples()

    s1 = a1.getSample('ds392875')
    print s1
    print '--------'
    print s1.tree1.GetEntries()
    print s1.GetName()
    print s1.name

    pass
funlist.append(test4)


def readTest():
    s1 = getSample("Analysis-Test/Sample_ds392875.root")
#     f1 = TFile("Analysis-Test/Sample_ds392875.root",'read')
#     s1 = f1.Get('ds392875')
#     s1.tree1 = f1.Get("tree1")
#     s1.m_file = f1
    print s1
    print '--------'
    print s1.tree1.GetEntries()
    print s1.GetName()
    print s1.name
# 
#     print f1.GetName()
# funlist.append(readTest)

def test2():
    dir0 = '/net/s3_data_home/dzhang/links/SAMPLES/R20/susyNtuple/AnalysisBase-02-04-26-da7031fc/'
    print dir0

    ### get the samples
    samples = {}
    sampleG = {}
    with open('datasets.list','read') as f1:
        grp = None
        for line in f1.readlines():
            line = line.rstrip()
            print line,'---in process'
            if len(line)==0 or line[0]=='#': continue
            if line.find(':')!=-1:
                tag = line[:line.find(':')]
                grp = SampleGroup(tag,tag,tag,tag)
                print 'group',grp.name,'created!'
                continue
            if line.find('END')!=-1:
                if grp:
                    print "End of", grp.name
                    sampleG[grp.name] = grp
                    grp = None
                continue
            
            ### make a sample
            fs = line.split('.')
            # user.clo.v8.13.364136.Sherpa_221_NNPDF30NNLO_Ztautau_MAXHTPTV140_280_BFilter_myOutput.root
            name = 'ds'+fs[4]
            s1 = Sample(name,fs[4],fs[4],fs[5])
            samples[s1.name] = s1
            s1.tree1 = TChain('evt2l')
            print dir0+line+'/*.root*'
            s1.tree1.Add(dir0+line+'/*.root*')
#             s1.tree1.Show(0)
#             return
            print s1.tree1.GetEntries()
            if grp:
                grp.sampleList.push_back(s1)


    sI1 = sampleInfo()
    sI1.Unit_fb = True
    sI1.loadAll()
#                 print sI1.getXSec('id==410000')
#                     print sI1.getXSec('id==392205')
#                         print sI1.getXSec('id==392205&&name=="127"')

    Lumi = 36.;
    ### get weights
    for s in samples:
        sx = samples[s]
        if sx.name.find('ds00') == -1:
            sx.weight = sI1.getXSec('id=='+sx.info)*Lumi/sx.getStatWeight()
        print sx.name, sx.getStatWeight(), sx.tree1.GetEntries(), sx.weight
    ### run selection

    ### set styles
    useStyle(sampleG['VV'],(3,1,1,3,20,1,3))
    useStyle(sampleG['Vgamma'],(4,1,1,4,20,1,4))
    useStyle(sampleG['Top'],(5,1,1,5,20,1,5))
    useStyle(sampleG['ZJet'],(2,1,1,2,20,1,2))
    sampleG['Data'].setUpOwnChain(TChain('evt2l'))

    p1 = Plot()
    p1.mode = 1
    p1.sData = sampleG['Data']
    p1.sSM.push_back(sampleG['Vgamma'])
    p1.sSM.push_back(sampleG['Top'])
    p1.sSM.push_back(sampleG['VV'])
    p1.sSM.push_back(sampleG['ZJet'])
    p1.mCut = 'evt.weight*(evt.flag==2&&l12.m<20)'

#     p1.test()
#     p1.showPlot('leps.pt', TH1F('h1','lep_pt_logy;Leading lepton p_{T} [GeV];Events / 3 GeV',100, 0, 300))
#     waitRootCmd(sDir+sTag+"leps_pt", sDirectly)

#     p1.showPlot('l12.m', TH1F('h1','l12_pt;m_{ll} [GeV];Events / 2 GeV',100, 0, 200))
    p1.showPlot('l12.m', TH1F('h1','l12_pt;m_{ll} [GeV];Events / 2 GeV',100, 0, 20))
    waitRootCmd(sDir+sTag+"l12_m", sDirectly)


    pass
# funlist.append(test2)

def test():
    dir0 = '/net/s3_data_home/dzhang/links/SAMPLES/R20/susyNtuple/AnalysisBase-02-04-26-da7031fc/'
    print dir0

    s1 = Sample('data','data_','Pseudo-data','Data')
    s1.tree1 = TChain('evt2l')
    s1.tree1.Add(dir0+'user.clo.v8.13.00307394.physics_Main_myOutput.root/*.root*')

    s2a = Sample('bkg_VV','VV_','VV background','VV')
    useStyle(s2a, (2,1,1,2,20,1,2))
    s2a.tree1 = TChain('evt2l')
    s2a.tree1.Add(dir0+'user.clo.v8.13.363492.Sherpa_221_NNPDF30NNLO_llvv_myOutput.root/*.root*')
    s2a.weight = 0.2
    xx = s2a.getStatWeight()
    L = 35.
    s2a.tree1.SetWeight(L/xx,'global')
    print xx
    waitRootCmd()

    s2b = Sample('bkg_Vg','Vg_','Vgamma background','V#gamma')
    useStyle(s2b, (4,1,1,4,20,1,4))
    s2b.tree1 = TChain('evt2l')
    s2b.tree1.Add(dir0+'user.clo.v8.13.301896.Sherpa_CT10_taunugammaPt35_70_myOutput.root/*.root')
    s2b.weight = 0.78

    s2c = Sample('bkg_Vg','Vg_','Vgamma background','V#gamma')
    useStyle(s2c, (4,1,1,4,20,1,4))
    s2c.tree1 = TChain('evt2l')
    s2c.tree1.Add(dir0+'user.clo.v8.13.301896.Sherpa_CT10_taunugammaPt35_70_myOutput.root/*.root')
    s2c.weight = 0.78

    s2bc = SampleGroup('bkg_Vgx','Vgx_','Vgamma background','V#gammaX')
    s2bc.sampleList.push_back(s2b)
    s2bc.sampleList.push_back(s2c)


    s3a = Sample('sig_1','sig1_','signal sample 1','Sig1')
    useStyle(s3a, (3,2,2))
    s3a.tree1 = TChain('evt2l')
    s3a.tree1.Add(dir0+'user.clo.v8.13.363491.Sherpa_221_NNPDF30NNLO_lllv_myOutput.root/*.root*')
    s3a.weight = 0.4

    s3b = Sample('sig_2','sig2_','signal sample 2','Sig2')
    useStyle(s3b, (5,9,2))
    s3b.tree1 = TChain('evt2l')
    s3b.tree1.Add(dir0+'user.clo.v8.13.361077.Sherpa_CT10_ggllvv_myOutput.root/*.root*')
    s3b.weight = 0.5

    p1 = Plot()
    p1.mode = 1
    p1.sData = s1
    p1.sSM.push_back(s2a)
#     p1.sSM.push_back(s2b)
    p1.sSM.push_back(s2bc)
    p1.sSig.push_back(s3a)
    p1.sSig.push_back(s3b)
    p1.mCut = ''

#     p1.test()
#     p1.showPlot('leps.pt', TH1F('h1','lep_pt_logy;Leading lepton p_{T} [GeV];Events / 3 GeV',100, 0, 300))
#     waitRootCmd(sDir+sTag+"leps_pt", sDirectly)

    p1.showPlot('l12.m', TH1F('h1','l12_pt;m_{ll} [GeV];Events / 2 GeV',100, 0, 200))
    waitRootCmd(sDir+sTag+"l12_m", sDirectly)

# funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    for fun in funlist: print fun()
