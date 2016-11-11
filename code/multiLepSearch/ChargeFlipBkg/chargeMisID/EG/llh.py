'''
 TUAN M. NGUYEN (PhD Student; Supervisor: Jean-Francois Arguin)
 Universite de Montreal
 tuan.nguyen.manh@cern.ch
 
 Electron charge flip
'''



#######################################################
# Minuit
#######################################################
def logLlhPair(i, j, parameters):
  nIJ = n[i,j]
  nssIJ = nss[i,j]
  epsI, epsJ = parameters[i], parameters[j]
  expectedNss = nIJ * ((1-epsI)*epsJ + (1-epsJ)*epsI)
  try:
    logLlh = nssIJ * log(expectedNss) - expectedNss
  except ValueError:
    logLlh = 0.0
  return logLlh

def negativeLog(numParameters, parameters):
  negLog = 0.0
  i = 0
  while i < numParameters[0]:
    j = 0
    while j < numParameters[0]:
      negLog -= logLlhPair(i, j, parameters)
      j += 1
    i += 1
  return negLog

# Every parameter of fcn is an array
def fcn(numParameters, derivatives, function, parameters, internalFlag):
  function[0] = negativeLog(numParameters, parameters)

def ratesAndErrors(numParameters, minuit):
  rates = np.zeros(numParameters)
  rateErrors = np.zeros(numParameters)
  rate, rateError = Double(0), Double(0)
  for i in xrange(numParameters):
    minuit.GetParameter(i, rate, rateError)
    rates[i] = float(rate)
    rateErrors[i] = float(rateError)
  return rates, rateErrors

def runMinuit(numParameters):
  minuit = TMinuit(numParameters)
  minuit.SetPrintLevel(0) 
  
  minuit.SetFCN(fcn)
  arglist = np.zeros(numParameters) + 0.01
  internalFlag, arglist[0] = Long(0), 0.5
  minuit.mnexcm("SET ERR", arglist, 1, internalFlag)
  
  initialValues = np.zeros(numParameters) + 0.01
  steps = np.zeros(numParameters) + 0.001
  
  for i in xrange(numParameters):
    name = "epsilon%s" % i
    minuit.mnparm(i, name, initialValues[i], steps[i], 0, 1, internalFlag)
  
  # arglist[0] = 2
  # minuit.mnexcm("SET STR", arglist, 1, internalFlag)
  
  arglist[0], arglist[1] = 10000, 0.1
  minuit.mnexcm("SIMPLEX", arglist, 1, internalFlag)
  minuit.mnexcm("MIGRAD", arglist, 1, internalFlag)
  
  print "FIT STATUS is " +str(minuit.GetStatus())
  return ratesAndErrors(numParameters, minuit)


#######################################################
# Run
#######################################################
from ROOT import TMinuit, Long, Double

from math import log as log
import numpy as np

import rutil
import time

start_time = time.time()
##############################

# etas = [0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.47]
etas = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.37, 1.52, 1.8, 2, 2.2, 2.47]
pts  = [20, 30, 40, 50, 60, 80, 120, 1000.0]


##############################
nEta,nPt = len(etas),len(pts)
numParameters = (nEta-1)*(nPt-1)
##############################

def doIt(n, nss, name, outFile):
  rates,errors = runMinuit(numParameters)
  ratesN, errorsN = '%s_rates.txt' % name, '%s_errors.txt' % name ##### Can also be without ".txt" at the not to save as txt
  with open(ratesN, 'wb') as ratesF, open(errorsN, 'wb') as errorsF:
    np.savetxt(ratesF, rates)   ##### if no ".txt" is given above, here use just save()
    np.savetxt(errorsF, errors) ##### if no ".txt" is given above, here use just save()
    with open(outFile, 'a') as txtOut:
      txtOut.write('%s\n%s\n' % (ratesN, errorsN))
  reLabels = rutil.ratesErrorsBins(rates, errors, etas, pts)
  rutil.printRatesErrorsBins(name, reLabels)
  
##### This part was added for TH2 plotting
  bName = rutil.beautifyName(name)
  rutil.writeToFile('%s.root' % name, rutil.hist2D(rates, errors, etas, pts, '%s_misid' % bName))
#####

  
##############################
# This stuff below cannot be put into a function, 
# because Minuit needs to see global variables n and nss!!!

if __name__ == '__main__':
  import sys
  txtFile, outFile = sys.argv[1], sys.argv[2]
  with open(txtFile, 'rb') as txt:
    lines = txt.read().splitlines()
  itLines = iter(lines)
  for ncsv in itLines:
    nsscsv = next(itLines)
    print 'Estimating rates for ... %s and %s' % (ncsv, nsscsv)
    n, nss = rutil.toan(ncsv), rutil.toan(nsscsv)
    doIt(n, nss, rutil.extractInfo(ncsv), outFile)




print("--- %s TIME ---" % (time.time() - start_time))















































