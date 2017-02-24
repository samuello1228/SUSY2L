import re

def getSumWfromNTUPList( ntupListFile ):
  """
  given a list of MC NTUP files
  get the SumW for each type of sample

  type is identified from datasetID in the fileName
  """
  outputDict = {}

  f = open(ntupListFile)
  for aLine in f:
    #remove and ignore comments in line
    aLine = aLine.split("#")[0] # content after # are considered commnets
    if len(aLine)==0: continue
  
    #read the line
    ntupName = aLine.split()[0]

    match = re.search(".[0-9]{6}.", ntupName)
    if not match:
      print "ssUtil.getSumWfromNTUPList : Cannot infer MC datasetID from filename %s , skipped" %ntupName
      continue

    sampleID = int(match.group()[1:-1])

    import ROOT
    ntupFile = ROOT.TFile(ntupName)
    sumW = ntupFile.Get("hCutFlow").GetBinContent(2) #FIXME: hardcoded bin2 for SumW

    if sampleID in outputDict:
      outputDict[sampleID] += sumW
    else:
      outputDict[sampleID] = sumW

  return outputDict

  
def loadSampleDict( sampleListFile ):
  """
  parse sample list file.
  line format in file:
  phyType phyProc generator runIDs

  eg.
  #Type   Process     Generator(Nominal)  id
  #diboson  WZ->qqll    Sherpa      361084      
  diboson WW->llvv(gg)    Sherpa      361077
  diboson Z(ee)+gamma   Sherpa      301535,301899-301901
  
  """
  outputDict = {}
  sampleListFile = open(sampleListFile)
  for aLine in sampleListFile:
    #remove and ignore comments in line
    aLine = aLine.split("#")[0] # content after # are considered commnets
    if len(aLine)==0: continue
  
    #read the line
    elements = aLine.split()
    if len(elements)!=4: 
      print "ssUtil.loadSampleDict : cannont understand line:"
      print "#" + aLine
      continue
  
    phyType   = elements[0]
    phyProc   = elements[1]
    generator = elements[2]
    runIDs    = elements[3].split(",")
  
    #interpret and expand runIDs range
    tmp = []
    for aStr in runIDs:
      if '-' in aStr:
        (startID, endID) = map(int, aStr.split('-') )
        tmp += [i for i in range(startID, endID+1)]
      else:
        tmp.append( int(aStr) )
    runIDs = tmp

    for runID in runIDs:
      outputDict[runID] = { "phyType"  : elements[0],
                "phyProc"  : elements[1],
          "generator": elements[2],
        }
  return outputDict

def guessSampleType( filename ):
   if (  "enugamma"   in filename): return "wgamma"
   if ( "munugamma"   in filename): return "wgamma"
   if ("taunugamma"   in filename): return "wgamma"
   if ("llll"         in filename): return "diboson"
   if ("lllv"         in filename): return "diboson"
   if ("llvv"         in filename): return "diboson"
   if ("ttW"          in filename): return "topX"
   if ("ttee"         in filename): return "topX"
   if ("ttmumu"       in filename): return "topX"
   if ("tttautau"     in filename): return "topX"
   if ("physics_Main" in filename): return "data"

   #match WZ and Slep samples eg MGPy8EG_A14N23LO_C1N2_WZ_300p0_250p0_3L_2L7_myOutput
   matches = re.search("[a-zA-Z0-9_]*C1N2[a-zA-Z0-9_]*", filename)
   if (matches): return matches.group(0)[0:-9] #-9 to remove _myOutput

   return None

def guessSampleDm(fileName):
  match = re.search("Slep_([0-9]+)_([0-9]+)_", fileName)
  if not match:
    print "Cannot guess dM from sample name"
    return None,None,None
 
  m1 = int(match.group(1))
  m2 = int(match.group(2))
  dm = abs(m1-m2)

  return dm,m1,m2

#--------------------------------------------------------------------------------
#Common setting that have to be sync-ed between 
#TMVA training, BDT application and HistFitter fit
#--------------------------------------------------------------------------------

#-----------------
#BDT variables
#-----------------
#format: (name, formula) 
#follow 2LSS note Ch5.2: Disciminating variables
#the formula reference variable names in obj_def.h
basicBDTVars = [
                 ( "mT2"    , "sig.mT2"    ),
                 ( "pt"     , "l12.pt"     ),
                 ( "MET"    , "sig.MetRel" ),
                 ( "Ht"     , "sig.HT"     ),
                 ( "mTl1"   , "leps.mT[0]" ),
                 ( "mTl2"   , "leps.mT[1]" ),
                 ( "ll_dPhi", "l12.dPhi"   ),
                 ( "l12m"   , "l12.m"),
               ]

#ISR region
isrBDTVars   = [
                 ( "JetMET_dPhi"  , "jets.MET_dPhi[0]"     ),
                 ( "MET_JetPt_R"  , "sig.MetRel/jets.pt[0]"),
                 ( "l1Pt_JetPt_R" , "leps.pt[0]/jets.pt[0]"),
               ]

#-----------------
#Cuts
#-----------------
trigCut   = "sig.trigCode!=0"
#tauCut    = "nTau==0"  #must hold, not filled in ssEvtSelection (as of 27Jul2016) 
#bjetCut   = "nBJet==0" #BJet info not filled in ssEvtSelection (as of 27Jul2016)
#cosmicCut = "nCosmic==0" #done in ssEvtSelection

isrCut       = "Sum$(jets.pt>20 && abs(jets.eta)<2.4) > 0" #nCentralJets>0 or ==0
nonisrCut    = "Sum$(jets.pt>20 && abs(jets.eta)<2.4) ==0" #nCentralJets>0 or ==0
zMassCut     = "!((int(abs(leps.ID[0])/1000)==11 || int(abs(leps.ID[0])/1000)==13) && int(abs(leps.ID[0])/1000) == int(abs(leps.ID[1])/1000) && fabs(l12.m - 91.1876)<=10)"
# zMassCut     = "!(int(abs(leps.ID[0])/1000) == int(abs(leps.ID[1])/1000) )"
eeCut        = "int(abs(leps.ID[0])/1000) == 11 && int(abs(leps.ID[1])/1000) == 11"
emuCut       = "(int(abs(leps.ID[0])/1000) == 11 && int(abs(leps.ID[1])/1000) == 13) || (int(abs(leps.ID[0])/1000) == 13 && int(abs(leps.ID[1])/1000) == 11)"
mumuCut      = "int(abs(leps.ID[0])/1000) == 13 && int(abs(leps.ID[1])/1000) == 13"

# For diagnostics
# isrCut = "isr"; nonisrCut = "noisr"; eeCut = "ee" ; emuCut = "emu" ; mumuCut = "mumu"

sigLepCut    = "((leps.lFlag[0] & 2) + (leps.lFlag[1] & 2) )==4" # see obj_def.h in multilepSearch, IS_SIGNAL = 2th bit in lFlag
ssCut = "((leps[0].ID>0) == (leps[1].ID>0))"

exact2LepCut = "Length$(leps.lFlag)==2"
lepptCut = "(leps.pt[0]>25.) && (leps.pt[1]>20.)"
l12mCut = "(l12.m>60.)"

sigLepSSWithDataBkgCut = "((%s)&&(%s)) || (!isMC && (qFwt+fLwt)!=0)" % (sigLepCut, ssCut)
ptMllCut = "%s && %s" % (lepptCut, l12mCut)

def getCut(ch):
  lepFlav = "1"
  whichISR = "1"

  if type(ch) is int:
    # Channels:
    #    noISR:  0=ee,  1=emu,  2=mumu,  3=combFlav,  4=SF,  
    #      ISR: 10=ee, 11=emu, 12=mumu, 13=combFlav, 14=SF, 
    #  combISR: 20=ee, 21=emu, 22=mumu, 23=combFlav, 24=SF, 
    # SFOSveto: +100
    # useISR = True if int(ch/10)==1 else False

    if int(int(ch%100)/10)==0:
      whichISR=nonisrCut
    elif int(int(ch%100)/10)==1:
      whichISR=isrCut

    if ch%10==0:
      lepFlav = eeCut
    elif ch%10==1:
      lepFlav = emuCut
    elif ch%10==2:
      lepFlav = mumuCut
    elif ch%10==4:
      lepFlav = "%s || %s" % (eeCut, mumuCut)

    if int(ch/100)==1:
      lepFlav = "(%s) && (%s || %s)" % (lepFlav, ssCut, emuCut)

    # lepFlav = "1"

  elif type(ch) is bool:
    whichISR = isrCut if ch else nonisrCut

  ptMllCut = "1"

  myCut = "&&".join(["(%s)"%cut for cut in [trigCut, whichISR, zMassCut, sigLepSSWithDataBkgCut, exact2LepCut, lepFlav, ptMllCut]])

  # return zMassCut
  return myCut
