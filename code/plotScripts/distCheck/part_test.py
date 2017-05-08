#!/usr/bin/env python
# import os, sys, re
# from ROOT import *
# from rootUtil import useAtlasStyle, waitRootCmd, savehistory, mkupHistSimple, get_default_fig_dir
# funlist=[]
# 
# sDir = get_default_fig_dir()
# sTag = 'test_'
# sDirectly = False
# if gROOT.IsBatch(): sDirectly = True

import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory, get_default_fig_dir
gROOT.LoadMacro('Sample.C+')
gROOT.LoadMacro('Plot.C+')
from ROOT import Plot, Sample, SampleGroup
from sample_info import sampleInfo
from anaUtil import loadSamples, anaSpace, loadSamplesL,useStyle 

gROOT.ProcessLine("Sample* noneSample(nullptr);")
from ROOT import noneSample;

sDir = get_default_fig_dir()
sTag = 'test_'
sDirectly = False
if gROOT.IsBatch(): sDirectly = True


def run_test1():
    p1 = Plot("ana1_Apr04b.root")
    p1.mode = 1
#     p1.sSM.push_back(p1.getSampleGroup('Top'))

    ZtautauL = SampleGroup("ZtautauL","ZtautauL","ZtautauL","ZtautauL")
    ZtautauH = SampleGroup("ZtautauH","ZtautauH","ZtautauH","ZtautauH")
    ZtautauM = SampleGroup("ZtautauM","ZtautauM","ZtautauM","ZtautauM")
    useStyle(ZtautauL, (4,1,1,4,20,1,4))
    useStyle(ZtautauH, (2,1,1,2,20,1,2))
    useStyle(ZtautauM, (3,1,1,3,20,1,3))

    for i in range(364210,364216): ZtautauL.m_sampleNames.push_back("ds"+str(i))
    for i in range(364128,364142): 
        if i!=364135: ZtautauH.m_sampleNames.push_back("ds"+str(i))
        else: ZtautauM.m_sampleNames.push_back("ds"+str(i))

    p1.loadSampleGroup(ZtautauL)
    p1.loadSampleGroup(ZtautauH)
    p1.loadSampleGroup(ZtautauM)

    p1.sData = noneSample
    p1.sSM.push_back(ZtautauM)
    p1.sSM.push_back(ZtautauH)
    p1.sSM.push_back(ZtautauL)

    sels = []
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


    p1.mode = 2
#     tTag = 'gradient_m2_'
    tTag = 'cID_m2_'
    for x in sels:
        p1.mCut = '(evt.weight*leps[0].ElChargeID==1&&leps[1].ElChargeID==1)'
#         p1.mCut = '(evt.weight*evt.pwt*evt.ElSF*evt.MuSF)*((sig.trigCode&sig.trigMask)!=0&&leps[0].ElChargeID==1&&leps[1].ElChargeID==1)'
#         p1.mCut = '(evt.weight*evt.pwt*evt.ElSF*evt.MuSF)*((sig.trigCode&sig.trigMask)!=0&&leps[0].ElChargeID==1&&leps[1].ElChargeID==1&&(leps[0].isoPass&0x2)!=0&&(leps[1].isoPass&0x2)!=0)'
#         if x[0][1].split("==")[1] in [4, 5, 6, 10, 11, 12]:
#             if x[0][0] in ['flag4', 'flag10']: ### show Z control region
#                 p1.sData = p1.getSample("data")
#                 p1.sData.sCuts = "(l12.m>81&&l12.m<101)"
#             else: p1.sData = noneSample
#         else:
#             p1.sData = p1.getSample("data")
#             p1.sData.sCuts = ''

        p1.useEntryList(x[0][0],False,x[0][1], False)
        p1.showInfo = x[1]
        sTag = tTag+'preSel_'+x[2]

        p1.showPlot('l12.m', TH1F('h1','l12_m_logy;m_{ll} [GeV];Events / 2 GeV',100, 0, 200))
        waitRootCmd(sDir+sTag+"mll", sDirectly)
        p1.showPlot('sig.Met', TH1F('h1','Met_logy;E_{T}^{Miss} [GeV];Events / 2 GeV',100, 0, 200))
        waitRootCmd(sDir+sTag+"Met", sDirectly)
        p1.showPlot('leps[0].mT', TH1F('h1','lep0_mT_logy;Leading lepton m_{T} [GeV];Events / 2 GeV',100, 0, 200))
        waitRootCmd(sDir+sTag+"mT0", sDirectly)
        p1.showPlot('leps[1].mT', TH1F('h1','lep1_mT_logy;Sub-leading lepton m_{T} [GeV];Events / 2 GeV',100, 0, 200))
        waitRootCmd(sDir+sTag+"mT1", sDirectly)

def test2():
    p1 = Plot("ana1_Apr04b.root")
    p1.mode = 1
#     p1.mCut = '(evt.weight*evt.pwt*evt.ElSF*evt.MuSF)*((sig.trigCode&sig.trigMask)!=0&&leps[0].ElChargeID==1&&leps[1].ElChargeID==1&&l12.m>15)'
    p1.mCut = '(evt.weight*evt.pwt*evt.ElSF*evt.MuSF)*((sig.trigCode&sig.trigMask)!=0&&leps[0].ElChargeID==1&&leps[1].ElChargeID==1&&l12.m>15&&leps[0].pt>25&&leps[1].pt>25)'
    gStyle.SetPalette(55);

    zJet = p1.getSampleGroup('ZJet')
#     zJet = p1.getSampleGroup('Top')
    p1.useEntryList(zJet, "flag2OR2OR8", "")

#     hx = zJet.getHStackFromTree("l12.pt",TH1F('h1','l12_pt_logy;p^{ll}_{T} [GeV];Events / 2 GeV',100, 0, 200), p1.mCut)
#     hx = zJet.getHStackFromTree("l12.pt",TH1F('h1','l12_pt_logy;p^{ll}_{T} [GeV];Events / 2 GeV',100, 0, 200), p1.mCut,"",False,10)
    hx = zJet.getHStackFromTree("l12.m",TH1F('h1','l12_m_logy;m_{ll} [GeV];Events / 2 GeV',100, 0, 200), p1.mCut,"",False,100)
    hx.Draw('hist')
    gPad.BuildLegend(0.75,0.75,0.95,0.95,"")
    gPad.Update()
    waitRootCmd()

def test():
    dir0 = os.getenv('SAMPLEDIR_LAMB')
    print dir0

    waitRootCmd()

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
#     run_test1()
    test2()
