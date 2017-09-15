#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory, mkupHistSimple, get_default_fig_dir
from math import pow
funlist=[]

sDir = get_default_fig_dir()
sTag = 'test_'
sDirectly = False
isDebug = False
if gROOT.IsBatch(): sDirectly = True



def dumpList():
    dir1 = '/home/dzhang/links/eosOther/cloShared/AnalysisBase-02-04-31-ebcb0e23/user.clo.v10.5.MCBG_myOutput.root/'
    dir2 = '/home/dzhang/links/eosOther/cloShared/AnalysisBase-02-04-31-ebcb0e23/user.dzhang.v10.5.MCBG_myOutput.root/'

    cut0 = "Length$(leps)==2&&(leps[0].lFlag&2)==2&&(leps[1].lFlag&2)==2&&leps[0].pt>25&&leps[1].pt>25"
    cut1Meff = 'leps[0].ID*leps[1].ID>0&&abs(leps[0].ID)/1000==11&&abs(leps[1].ID)/1000==11&&Length$(jets)==1&&Sum$((jets.jFlag&(1<<5))!=0)==0&&fabs(l12.m-91.2)>10&&leps[0].pt>65&&leps[1].pt>25&&fabs(leps[0].eta-leps[1].eta)<1.5&&(sig.HT+sig.Met)>200'
    cut1 = 'leps[0].ID*leps[1].ID>0&&abs(leps[0].ID)/1000==11&&abs(leps[1].ID)/1000==11&&Length$(jets)==1&&Sum$((jets.jFlag&(1<<5))!=0)==0&&fabs(l12.m-91.2)>10&&leps[0].pt>65&&leps[1].pt>25&&fabs(leps[0].eta-leps[1].eta)<1.5&&(sig.HT+sig.Met)>200&&(2*leps[0].pt*sig.Met*(1-cos(leps[0].MET_dPhi))>15625||2*leps[1].pt*sig.Met*(1-cos(leps[1].MET_dPhi))>15625)'

    nEvt = 0
    with open('../BGSample.txt') as f1:
        lines = f1.readlines()
    for line in lines[117:131]:

        ch1 = TChain('evt2l')
        fs = line.rstrip().split()
        ds = dir1+"*"+fs[0]+'.'+fs[1]+'*.root'
        print line.rstrip(), ds

        if fs[0]!='361070': continue

        ch1.Add(ds)

        ch1.Show(51643)
#         ch1.Scan("(sig.HT+sig.Met):sqrt(2*leps[0].pt*sig.Met*(1-cos(leps[0].MET_dPhi))):sqrt(2*leps[1].pt*sig.Met*(1-cos(leps[1].MET_dPhi))):sig.mlj","Entry$==51643")
        ch1.Scan("evt.event:sig.Met",cut0+'&&'+cut1)
        return


        ch1.GetPlayer().SetScanRedirect(True); 
#         nEvt += ch1.Scan("evt.event", cut0+'&&'+cut1)

#         ch1.GetPlayer().SetScanFileName(fs[1]+"_x.txt"); 
#         nEvt += ch1.Scan("evt.event:sqrt(2*leps[0].pt*sig.Met*(1-cos(leps[0].MET_dPhi))):sqrt(2*leps[1].pt*sig.Met*(1-cos(leps[1].MET_dPhi)))", cut0+'&&'+cut1)


        ch1.GetPlayer().SetScanFileName(fs[1]+"_meff.txt"); 
        nEvt += ch1.Scan("evt.event:sig.HT+sig.Met:sig.HT:sig.Met", cut0+'&&'+cut1Meff)

        print fs[1], ch1.GetEntries()
#         break
    print "total events:", nEvt


def getCutFlowTable(chs = None):
    chname = 'evt2l'
    if chs is None: chs = ['ee1','ee2','mumu1','mumu2','emu1','emu2']

    ### get configuration
    lines = None
    with open('../BGSample.txt') as f1:
        lines = f1.readlines()
    ### get dibosn numbers
    dir1 = '/home/dzhang/links/eosOther/cloShared/AnalysisBase-02-04-31-ebcb0e23/user.clo.v10.5.MCBG_myOutput.root/'
    dir2 = '/home/dzhang/links/eosOther/cloShared/AnalysisBase-02-04-31-ebcb0e23/user.dzhang.v10.5.MCBG_myOutput.root/'
    cut0 = "Length$(leps)==2&&(leps[0].lFlag&2)==2&&(leps[1].lFlag&2)==2&&leps[0].pt>25&&leps[1].pt>25"

    ch_cuts = {}
    ch_cuts['ee1'] = [('ee','abs(leps[0].ID)/1000==11&&abs(leps[1].ID)/1000==11') ## ee
                     ,('Njet','Length$(jets)==1')
                     ,('bVeto','Sum$((jets.jFlag&(1<<5))!=0)==0')
                     ,('mll','fabs(l12.m-91.2)>10')
                     ,('pt1','leps[0].pt>65')
#                      ,('pt2','leps[1].pt>25')
                     ,('pt2','-')
                     ,('dEta','fabs(leps[0].eta-leps[1].eta)<1.5')
                     ,('Meff','(sig.HT+sig.Met)>200')
                     ,('maxMt','(2*leps[0].pt*sig.Met*(1-cos(leps[0].MET_dPhi))>15625||2*leps[1].pt*sig.Met*(1-cos(leps[1].MET_dPhi))>15625)')
                     ,('mlj','sig.mlj<105')
                     ]
    ch_cuts['ee2'] = [('ee','abs(leps[0].ID)/1000==11&&abs(leps[1].ID)/1000==11') ## ee
                     ,('Njet','Length$(jets)>1&&Length$(jets)<4')
                     ,('bVeto','Sum$((jets.jFlag&(1<<5))!=0)==0')
                     ,('mll','fabs(l12.m-91.2)>10')
                     ,('pt1','-')
                     ,('pt2','-')
                     ,('maxMt','(2*leps[0].pt*sig.Met*(1-cos(leps[0].MET_dPhi))>{0:f}||2*leps[1].pt*sig.Met*(1-cos(leps[1].MET_dPhi))>{0:f})'.format(pow(85,2)))
                     ,('MetRel','sig.MetRel>70')
                     ,('dEta','fabs(leps[0].eta-leps[1].eta)<1.5')
                     ,('Meff','(sig.HT+sig.Met)>200')
                     ,('maxMt','(2*leps[0].pt*sig.Met*(1-cos(leps[0].MET_dPhi))>15625||2*leps[1].pt*sig.Met*(1-cos(leps[1].MET_dPhi))>15625)')
                     ,('mlj','sig.mlj<205')
                     ,('mT2','sig.mT2>15')
                     ]


    ch_cuts['mumu1'] = [('mumu','(abs(leps[0].ID)/1000+abs(leps[1].ID)/1000)==26') ## emu / mue
                     ,('Njet','Length$(jets)==1')
                     ,('bVeto','Sum$((jets.jFlag&(1<<5))!=0)==0')
                     ,('mll','-')
                     ,('pt1','-')
                     ,('pt2','-')
                     ,('maxMt','(2*leps[0].pt*sig.Met*(1-cos(leps[0].MET_dPhi))>{0:f}||2*leps[1].pt*sig.Met*(1-cos(leps[1].MET_dPhi))>{0:f})'.format(pow(100,2)))
                     ,('mlj','sig.mlj<100')
                     ]


    ch_cuts['mumu2'] = [('mumu','(abs(leps[0].ID)/1000+abs(leps[1].ID)/1000)==26') ## emu / mue
                     ,('Njet','Length$(jets)>1&&Length$(jets)<4')
                     ,('bVeto','Sum$((jets.jFlag&(1<<5))!=0)==0')
                     ,('mll','-')
                     ,('pt1','-')
                     ,('pt2','-')
                     ,('MetRel','sig.MetRel>40')
                     ,('maxMt','(2*leps[0].pt*sig.Met*(1-cos(leps[0].MET_dPhi))>{0:f}||2*leps[1].pt*sig.Met*(1-cos(leps[1].MET_dPhi))>{0:f})'.format(pow(75,2)))
                     ,('dEta','fabs(leps[0].eta-leps[1].eta)<1.5')
                     ,('Meff','(sig.HT+sig.Met)>200')
                     ,('mlj','sig.mlj<205')
                     ,('ptll','l12.pt>40')
                     ]


    ch_cuts['emu1'] = [('emu','(abs(leps[0].ID)/1000+abs(leps[1].ID)/1000)==24') ## emu / mue
                     ,('Njet','Length$(jets)==1')
                     ,('bVeto','Sum$((jets.jFlag&(1<<5))!=0)==0')
                     ,('mll','-')
                     ,('pt1','-')
                     ,('pt2','-')
                     ,('MetRel','sig.MetRel>80')
                     ,('maxMt','-')
                     ,('mlj','sig.mlj<95')
                     ]


    ch_cuts['emu2'] = [('emu','(abs(leps[0].ID)/1000+abs(leps[1].ID)/1000)==24') ## emu / mue
                     ,('Njet','Length$(jets)>1&&Length$(jets)<4')
                     ,('bVeto','Sum$((jets.jFlag&(1<<5))!=0)==0')
                     ,('pt1','-')
                     ,('pt2','-')
                     ,('mll','-')
                     ,('dEta','fabs(leps[0].eta-leps[1].eta)<1.5')
                     ,('Meff','(sig.HT+sig.Met)>200')
                     ,('MetRel','sig.MetRel>25')
                     ,('maxMt','(2*leps[0].pt*sig.Met*(1-cos(leps[0].MET_dPhi))>{0:f}||2*leps[1].pt*sig.Met*(1-cos(leps[1].MET_dPhi))>{0:f})'.format(pow(105,2)))
                     ,('mlj','sig.mlj<265')
                     ,('mT2','sig.mT2>75')
                     ]

    print ch_cuts['ee1']
#     return


    gROOT.LoadMacro('elements.C')
    from ROOT import GroupData, getBKG
    for bkg in getBKG():
        if bkg.GroupName not in ['VVV','Higgs']: continue
#         if bkg.GroupName == "VV": continue
        print bkg.lower,bkg.upper+1

        ch1 = TChain(chname)
        ### get the chain
        for line in lines[bkg.lower:bkg.upper+1]:
            fs = line.rstrip().split()
            ds = dir1+"*"+fs[0]+'.'+fs[1]+'*.root'
            ds2 = dir2+"*"+fs[0]+'.'+fs[1]+'*.root'
            print line.rstrip(), ds
            ch1.Add(ds)
            ch1.Add(ds2)
#             break

        print ch1.GetEntries(cut0)
        ### get entrylist
        ch1.Draw('>>elist',cut0,'entrylist')
        elist = gDirectory.Get('elist')
        ch1.SetEntryList(elist)
        if isDebug:
            elist.Print()
            print elist.GetN()

        chs = ['ee1']
        for chnl in chs:
            cutlist = ch_cuts[chnl]
            print chnl
            cutV = 'leps[0].ID*leps[1].ID>0'
            print 'SS', ch1.GetEntries(cutV)
            for cut in cutlist:
                if cut[1] == '-':
                    print cut[0], '-'
                    continue
                cutV += '&&'+cut[1]
                print cut[0], ch1.GetEntries(cutV)

def test():
    dir0 = os.getenv('SAMPLEDIR_LAMB')
    print dir0

    chname = 'evt2l'
    ch1 = TChain(chname)
#     ch1.Add('/home/dzhang/links/eosOther/cloShared/AnalysisBase-02-04-31-ebcb0e23/user.clo.v10.5.MCBG_myOutput.root/user.clo.mc15_13TeV.410000.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad.11771572._000094.myOutput.root')

    ### get dibosn numbers
    dir1 = '/home/dzhang/links/eosOther/cloShared/AnalysisBase-02-04-31-ebcb0e23/user.clo.v10.5.MCBG_myOutput.root/'
    lines = None
    with open('../BGSample.txt') as f1:
        lines = f1.readlines()
    for line in lines[117:131]:
        fs = line.rstrip().split()
        ds = dir1+"*"+fs[0]+'.'+fs[1]+'*.root'
        print line.rstrip(), ds
        ch1.Add(ds)
        break

#     return

    ### Scan for event check
#     ch1.Show(3999)
#     return

#     ch1.Show(0)
    cut0 = "Length$(leps)==2&&(leps[0].lFlag&2)==2&&(leps[1].lFlag&2)==2&&leps[0].pt>25&&leps[1].pt>25"
    print ch1.GetEntries(cut0)

#     ch1.Draw('>>elist',cut0,'entrylist')
#     elist = gDirectory.Get('elist')
#     ch1.SetEntryList(elist)
#     elist.Print()
#     print elist.GetN()

#     print ch1.GetEntries("")
    cuts = []

    cutV = 'leps[0].ID*leps[1].ID>0'
    cuts.append(('SS',cutV))

    cutV += '&&abs(leps[0].ID)/1000==11&&abs(leps[1].ID)/1000==11'
    cuts.append(('ee',cutV))

    cutV += '&&Length$(jets)==1'
    cuts.append(('1jet',cutV))

    cutV += '&&Sum$((jets.jFlag&(1<<5))!=0)==0'
    cuts.append(('bVeto',cutV))

#     print cut0+'&&'+cutV
#     ch1.Draw('l12.m>>h1(100,0,200)', cut0+'&&'+cutV)
#     h1 = gPad.GetPrimitive('h1')
#     h1.SetLineColor(2)
# 
    cutV += '&&fabs(l12.m-91.2)>10'
    cuts.append(('mll',cutV))

#     ch1.Draw('l12.m', cut0+'&&'+cutV,'same')
#     waitRootCmd()

    cutV += '&&leps[0].pt>65'
    cuts.append(('pt1',cutV))

    cutV += '&&leps[1].pt>25'
    cuts.append(('pt2',cutV))

    cutV += '&&fabs(leps[0].eta-leps[1].eta)<1.5'
    cuts.append(('dEta',cutV))

    cutV += '&&(sig.HT+sig.Met)>200'
    cuts.append(('Meff',cutV))


#     ch1.Scan('evt.run:evt.event:leps[0].mT:leps[1].mT:sqrt(2*leps[0].pt*sig.Met*(1-cos(leps[0].MET_dPhi))):sqrt(2*leps[1].pt*sig.Met*(1-cos(leps[1].MET_dPhi)))', cut0+'&&'+cutV)
#     return

#     cutV += '&&(leps[0].mT>125||leps[1].mT>125)'
    cutV += '&&(2*leps[0].pt*sig.Met*(1-cos(leps[0].MET_dPhi))>15625||2*leps[1].pt*sig.Met*(1-cos(leps[1].MET_dPhi))>15625)'
    cuts.append(('MaxmT',cutV))

    cutV += '&&sig.mlj<105'
    cuts.append(('mlj',cutV))

    for c in cuts:
#         print c[0], c[1]
        print c[0], ch1.GetEntries(c[1])



#     print ch1.GetEntries("Length$(leps)==2&&(leps[0].lFlag&2)==2&&(leps[1].lFlag&2)==2")
#     print ch1.GetEntries("Length$(leps)==2&&(leps[0].lFlag&2)==2&&(leps[1].lFlag&2)==2&&abs(leps[0].ID)/1000==11&&abs(leps[1].ID)/1000==11")
#     print ch1.Draw("leps.pt>>(100,0,200)","(leps.lFlag&2)==2")
#     print ch1.Draw("leps.pt>>(100,0,200)","Length$(leps)==2&&(leps[0].lFlag&2)==2&&(leps[1].lFlag&2)==2&&abs(leps[0].ID)/1000==11&&abs(leps[1].ID)/1000==11")


    waitRootCmd()
funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    getCutFlowTable()
#     for fun in funlist: print fun()
#     dumpList()
