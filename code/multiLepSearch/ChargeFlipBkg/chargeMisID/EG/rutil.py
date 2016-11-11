'''
 TUAN M. NGUYEN (PhD Student; Supervisor: Jean-Francois Arguin)
 Universite de Montreal
 tuan.nguyen.manh@cern.ch
 
 Electron charge flip
'''


import csv
import re
import math
import numpy as np
from array import array

from ROOT import TH2F, TH1F, TFile


def extractInfo(fileName):
  pat = re.compile(r'(\d.*?).csv')
  info = pat.search(fileName)
  return info.groups()[0] if info else fileName[:-4]

#######################################################
def toan(csvFile):
  with open(csvFile, 'rb') as f:
    q = list(csv.reader(f, delimiter=','))

  for row in q:
    row[:] = [float(n) for n in row]

  return np.array(q)

#######################################################
def makeRange(values):
  ticks = []
  for i in xrange(len(values)-1):
    tick = '%s-%s' % (values[i], values[i+1])
    ticks.append(tick)
  return ticks
 
#ptRates = [rates[0]]


 
#######################################################
def ratesErrorsBins(rates, errors, etas, pts):
  ratesErrors = zip(rates, errors)
  etaTicks, ptTicks = makeRange(etas), makeRange(pts)
  ptRates = [ratesErrors[i:i+len(pts)-1] for i in range(0,len(rates), len(pts)-1)]
  ptRatesLabels = [zip(ptTicks, e) for e in ptRates]
  etaPtRatesLabels = zip(makeRange(etas), ptRatesLabels)
  return etaPtRatesLabels
  
def printRatesErrorsBins(tag, etaPtRatesLabels):
  print "\n\n\n=====================================\n\n RESULT\n\n"
  print tag, "\n\n\n"
  for e in etaPtRatesLabels:
    print e[0]
    for se in e[1]:
      print se
      
##### This part was added for plotting

#######################################################
# name has the form '70.000000_110.000000_0.000000_0.000000_mc'
def beautifyName(name):
  d = '_'
  parts = name.split(d)
  nums = parts[:-1]  
  nums = map(str,map(float, nums))
  return d.join(nums+parts[-1:])

#######################################################

def hist2D(rates, errors, etas, pts, histName):
  etaArr, ptArr = array('d',etas), array('d',pts)
  nEta,nPt = len(etas) - 1, len(pts) - 1
  hist = TH2F(histName, ";#eta;P_{T}",nEta, etaArr, nPt, ptArr)
  for i, (rate,error) in enumerate(zip(rates,errors)):
    iBiny = i % nPt
    iBinx = (i - iBiny) / nPt
    hist.SetBinContent(iBinx + 1, iBiny +1, rate)
    hist.SetBinError(iBinx + 1, iBiny +1, error)
  return hist
    
def writeToFile(name, hist):
  rootFile = TFile(name,"recreate")
  hist.Write()
  rootFile.Close()

#####

if __name__ == '__main__':
  pass
  













  
  
