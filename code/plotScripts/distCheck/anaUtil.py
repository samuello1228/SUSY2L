#!/usr/bin/env python
from ROOT import *
gROOT.LoadMacro('Plot.C+')
from ROOT import Plot, Sample, SampleGroup


def useStyle(sample, style, size=None):
    for i in range(size if size else len(style)): sample.style[i]=style[i]

def loadSamplesL(dslist, dsname=None, dir0='/net/s3_data_home/dzhang/links/SAMPLES/R20/susyNtuple/AnalysisBase-02-04-26-da7031fc/'):
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
            print line
            fs = line.split('.')
            name = 'ds'+fs[0]
            s1 = Sample(name,fs[0],fs[1],fs[0])
            samples[s1.name] = s1
            s1.tree1 = TChain('evt2l')
            s1.tree1.Add(dir0+"*"+line+'*.root*')
            if grp:
                grp.sampleList.push_back(s1)
    return samples, sampleG

def loadSamples(dslist, dsname=None, dir0='/net/s3_data_home/dzhang/links/SAMPLES/R20/susyNtuple/AnalysisBase-02-04-26-da7031fc/'):
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
            print line
            fs = line.split('.')
            name = 'ds'+fs[4]
            s1 = Sample(name,fs[4],fs[4],fs[5])
            samples[s1.name] = s1
            s1.tree1 = TChain('evt2l')
            s1.tree1.Add(dir0+line+'/*.root*')
            if grp:
                grp.sampleList.push_back(s1)
    return samples, sampleG

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
