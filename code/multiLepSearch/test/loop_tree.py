#!/usr/bin/env python
from ROOT import TChain, gROOT, TFile, AddressOf
gROOT.ProcessLine('.L ../macro/loader.C')
# from ROOT import SIGNATURE
# sig = SIGNATURE() 

def test():
    ch1 = TChain('evt2l')
    ch1.Add("/net/ustc_home/dzhang/work/bsmSearch/ewSUSY/analysis/evtSelect/run/v5/data-myOutput/xk.root")
#     ch1.SetBranchAddress("sig",AddressOf(sig))
    for i in range(ch1.GetEntries()):
        ch1.GetEntry(i)
        for j in ch1.jets:
            print 'jet pt:', j.pt
        for l in ch1.leps:
            print 'lepton pt:', l.pt
        print 'mT2:', ch1.mT2
#     for event in ch1:
#         print event.mT2

if __name__ == '__main__':
    test()
