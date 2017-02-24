#!/usr/bin/env python
from ROOT import *
gROOT.LoadMacro('Plot.C+')
from ROOT import Plot, Sample, SampleGroup

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
