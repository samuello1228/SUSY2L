# skimNTUP.py
# script to skim NTUP from multiLepSearch
# output a skimmed NTUP that can be fed to TMVA for training and cut optimization

import sys
import os

if len(sys.argv)!=3:
  print "Usage: python skimNTUP.py inFile outFile"
  sys.exit()

#if you have set up thr shell env "FiltreeUtilPath" i.e.
#    export FiltreeUtilPath=/path/to/FiltreeUtil
#then the prog will use that, else default to ./
path = os.environ.get( "FiltreeUtilPath", "./") 

import ROOT
ROOT.gSystem.AddDynamicPath( path )
ROOT.gROOT.Macro( path+"/rootlogon.C" )
from ROOT import FiltreeStream, TTree, TChain, TFile

#---------------------------Main-------------------------------
# The input data
t = TChain("evt2l")
t.Add( sys.argv[1])

#the output files
outFile = TFile(sys.argv[2],"RECREATE")

# The nodes in the filter chain
fsSrc = FiltreeStream("fsSrc")


#set output verbose level
fsSrc.printLv = 2

#connect the src node to the data
fsSrc.ConnectTo(t)
fsSrc.SetOutput("evt2l_skim", outFile)

#All TTree Draw cmd can be used in principle
fsSrc.AddCut("trigger" , "sig.trigCode!=0")
fsSrc.AddCut("GRL"     , "evt.cuts==1")
#------------------------------------------------------------------------------------------
#this won't crash but corrupt subsequent cutflow result when length change from 2 to 3 !!!
#fsSrc.AddCut("2Lepton1", "Length$(leps)") 
#fsSrc.AddCut("2Lepton" , "Length$(leps)==2")
#------------------------------------------------------------------------------------------
fsSrc.AddCut("2Lepton" , "@leps.size()==2")
fsSrc.AddCut("pt1"     , "leps.pt[0]>30") # unit: GeV
fsSrc.AddCut("pt2"     , "leps.pt[1]>30")
fsSrc.AddCut("mll_60"  , "l12.m>60")

#fsSrc.AddSpectator("ID0", "int(leps.ID[0]/1000)")
#fsSrc.AddSpectator("ID1", "int(leps.ID[1]/1000)")

fsSrc.SetWriteBranch("*", False)
fsSrc.SetWriteBranch("evt", True) #donno why have to keep this else crash
fsSrc.SetWriteBranch("l12", True)
fsSrc.SetWriteBranch("sig", True)
fsSrc.SetWriteBranch("leps.pt", True)
fsSrc.SetWriteBranch("leps.ID", True)

# start the filter process chain
#fsSrc.StartDataFlow(1000)
fsSrc.StartDataFlow()

#save the cut flow
tmpFile = ROOT.TFile(sys.argv[1], "READ")
hMyStreamCutFlow = fsSrc.GetCutFlowHist()
hUpStreamCutFlow = tmpFile.FindObjectAny("hCutFlow")

#print hUpStreamCutFlow


if (not hUpStreamCutFlow):
  hMyStreamCutFlow.SetNameTitle("hCombinedCutFlow", "hCombinedCutFlow")
  outFile.WriteTObject(hMyStreamCutFlow)
else:
  #combine upStream cutflow with ours
  upNCuts = hUpStreamCutFlow.GetNbinsX() 
  myNCuts = hMyStreamCutFlow.GetNbinsX()
  totNCuts = upNCuts + myNCuts -1
  hCombinedCutFlow = ROOT.TH1D("hCombinedCutFlow", "hCombinedCutFlow", totNCuts, 0, totNCuts)
  for i in range(upNCuts):
    nPass = hUpStreamCutFlow.GetBinContent(i+1)
    label = hUpStreamCutFlow.GetXaxis().GetBinLabel(i+1)
    hCombinedCutFlow.SetBinContent(i+1, nPass)
    hCombinedCutFlow.GetXaxis().SetBinLabel(i+1, label)
  for i in range(1,myNCuts):
    nPass = hMyStreamCutFlow.GetBinContent(i+1)
    label = hMyStreamCutFlow.GetXaxis().GetBinLabel(i+1)
    hCombinedCutFlow.SetBinContent(upNCuts+i, nPass)
    hCombinedCutFlow.GetXaxis().SetBinLabel(upNCuts+i, label)
  outFile.WriteTObject(hCombinedCutFlow)
  
outFile.Write()
outFile.Close()

tmpFile.Close()
