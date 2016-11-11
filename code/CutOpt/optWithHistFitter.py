import ROOT
import os
import subprocess
import xml.dom.minidom
import sys

#usage:
#python optWithHistFitter2.py inFile.root outFile.root sigCount bkgCount


inF  = ROOT.TFile(sys.argv[1], "READ")
outF = ROOT.TFile(sys.argv[2], "RECREATE")
#expected count at 3248pb-1 before applying Optimized Cuts
#startSig = 27
#startBkg = 807547
startSig = float(sys.argv[3])
startBkg = float(sys.argv[4])


heffS = inF.FindObjectAny("MVA_BDTD_effS")
heffB = inF.FindObjectAny("MVA_BDTD_effB")

hSignificance = ROOT.TH1D(heffS)
hSignificance.Reset()
hSignificance.SetNameTitle("significance", "significance")

hPval = ROOT.TH1D(heffS)
hPval.Reset()
hPval.SetNameTitle("Pvalue", "Pvalue")

outList = []

nBins = heffS.GetNbinsX()
for iBin in range(1,nBins+1, 100):
  effS = heffS.GetBinContent(iBin)
  effB = heffB.GetBinContent(iBin)
  print iBin, effS, effB

  inputNbkg  = effB*startBkg
  inputNdata = effS*startSig + inputNbkg
  
  if (inputNbkg<=0): continue  #this happens when TMVA could not find a valid cut
  if (inputNdata<=0): continue

  #You need to setup HitFitter in shell first for this to work
  histFitterDir = "/home/ytchan/HistFitterTutorial/"
  histFitterCmd = 'HistFitter.py -w -f -z %s/MyRepeatSUSY_UpperLimitAnalysis_SS.py -c "inputNdata=%f;inputNbkg=%f"' \
                  % (os.getcwd(),inputNdata,inputNbkg)
  subprocess.Popen( histFitterCmd, shell=True, cwd=histFitterDir).wait()

  hfOutFile = ROOT.TFile("%s/results/MyCutOptAnalysis_SS_Output_hypotest.root"%histFitterDir, "READ")
  fitResult = hfOutFile.Get("discovery_htr_Sig")
  sig = fitResult.Significance()
  p0  = fitResult.NullPValue()
  hfOutFile.Close()

  hSignificance.SetBinContent(iBin, sig)
  hPval.SetBinContent(iBin, p0)
  aRow = (iBin, effS, effB, inputNdata, inputNbkg, p0, sig)
  print aRow
  outList.append(aRow)

  #break
outF.WriteTObject(hSignificance)
outF.WriteTObject(hPval)
outF.Close()

for iBin, effS, effB, inputNdata, inputNbkg, p0, sig in outList:
  print iBin, effS, effB, inputNdata, inputNbkg, p0, sig
