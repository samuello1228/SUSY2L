#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory, mkupHistSimple, get_default_fig_dir
funlist=[]

sDir = get_default_fig_dir()
sTag = 'test_'
sDirectly = False
if gROOT.IsBatch(): sDirectly = True

def getLumi(runStart=None, runEnd=None):
    ch1 = TChain('LumiMetaData')
    ch1.Add("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root")
    ch1.Add("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.lumicalc.OflLumi-13TeV-008.root")
#     ch1.Show(0)

    t1 = 0
    for i in range(ch1.GetEntries()):
        ch1.GetEntry(i)
        if runStart and ch1.RunNbr < runStart: continue
        if runEnd and ch1.RunNbr > runEnd: continue
        t1 += ch1.IntLumi
    return t1*1e-9
#     print "Lumi is", t1*1e-9, "/fb"

def test():
#     print getLumi()
#     print '2015:', getLumi(276073,284484)
    print '2016 D4-L:', getLumi(302919,311481)


#     waitRootCmd()
funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    for fun in funlist: print fun()
