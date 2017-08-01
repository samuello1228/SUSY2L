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
   if ("ttbarW"       in filename): return "topX"
   if ("ttee"         in filename): return "topX"
   if ("ttmumu"       in filename): return "topX"
   if ("tttautau"     in filename): return "topX"
   if ("ttH"          in filename): return "topX"
   if ("ttZ"          in filename): return "topX"
   if ("4top"         in filename): return "4top"
   if ("physics_Main" in filename): return "data"
   if ("Zee_Mll10_40" in filename): return "ZeeLowM"
   if ("Zmm_Mll10_40" in filename): return "ZmmLowM"
   if ("eegamma"      in filename): return "eegamma"
   if ("mumugamma"    in filename): return "mumugamma"
   if ("tautaugamma"  in filename): return "tautaugamma"
   if ("Ztt_Mll10_40" in filename): return "ZtautauLowM"
   if ("WqqZll"       in filename): return "diboson"
   if ("ZqqZll"       in filename): return "diboson"
   if ("Zee_Pt"       in filename): return "Zee"

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
                 ( "mT2"    , "sig.mT2 * (sig.mT2>0)"              ),
                 ( "pt"     , "l12.pt * (l12.pt>0)  "              ),
                 ( "MET"    , "sig.MetRel * (sig.MetRel>0)"        ),
                 ( "Ht"     , "sig.HT * (sig.HT>0)"                ),
                 ( "mTl1"   , "leps.mT[0] * (leps.mT[0]>0)"        ),
                 ( "mTl2"   , "leps.mT[1] * (leps.mT[1]>0)"        ),
                 ( "ll_dPhi", "l12.dPhi * (fabs(l12.dPhi) < 3.15)" ),
                 ( "l12m"   , "l12.m * (l12.m>0)"                  ),
                 # ( "l12m"   , "(int(abs(leps.ID[0]))!=int(abs(leps.ID[1])))*100+l12.m"),
               ]

#ISR region
isrBDTVars   = [
                 ( "JetMET_dPhi"  , "jets.MET_dPhi[0] * (fabs(l12.dPhi) < 3.15)"              ),
                 ( "MET_JetPt_R"  , "sig.MetRel/jets.pt[0] * (sig.MetRel > 0 && jets.pt[0]>0)"),
                 ( "l1Pt_JetPt_R" , "leps.pt[0]/jets.pt[0] * (leps.pt[0] > 0 && jets.pt[0]>0)"),
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
# zMassCut     = "!((int(abs(leps.ID[0])/1000)==11 || int(abs(leps.ID[0])/1000)==13) && int(abs(leps.ID[0])/1000) == int(abs(leps.ID[1])/1000) && fabs(l12.m - 91.1876)<=10)"
zMassCut     = "!(fabs(l12.m - 91.1876)<=10)"
# eeCut        = "int(abs(leps.ID[0])/1000) == 11 && int(abs(leps.ID[1])/1000) == 11"
eeCut        = "Sum$(abs(leps.ID))==22000"
# emuCut       = "(int(abs(leps.ID[0])/1000) == 11 && int(abs(leps.ID[1])/1000) == 13) || (int(abs(leps.ID[0])/1000) == 13 && int(abs(leps.ID[1])/1000) == 11)"
emuCut       = "Sum$(abs(leps.ID))==24000"
# mumuCut      = "int(abs(leps.ID[0])/1000) == 13 && int(abs(leps.ID[1])/1000) == 13"
mumuCut      = "Sum$(abs(leps.ID))==26000"
# cftCut       = "leps.ElChargeID[0] && leps.ElChargeID[1]" #Deprecated 

# For diagnostics
# isrCut = "isr"; nonisrCut = "noisr"; eeCut = "ee" ; emuCut = "emu" ; mumuCut = "mumu"

sigLepCut    = "((leps.lFlag[0] & 2) + (leps.lFlag[1] & 2) )==4" # see obj_def.h in multilepSearch, IS_SIGNAL = 2th bit in lFlag
ssCut = "((leps[0].ID>0) == (leps[1].ID>0))"

exact2LepCut = "Length$(leps.lFlag)==2"
lepptCut = "(leps.pt[0]>25.) && (leps.pt[1]>20.)"
l12mCut = "(l12.m>60.)"

ptMllCut = "%s && %s" % (lepptCut, l12mCut)
ptMllCut = "%s" % (lepptCut)

def getCut(ch):
  lepFlav = "1"
  whichISR = "1"

  if type(ch) is int:
    '''
    Channels:
       noISR:  0=ee,  1=emu,  2=mumu,  3=combFlav,  4=SF,  
         ISR: 10=ee, 11=emu, 12=mumu, 13=combFlav, 14=SF, 
     combISR: 20=ee, 21=emu, 22=mumu, 23=combFlav, 24=SF, 

    Allow DFOS: +100 // Only use for signal samples!! Use for dm>=100 GeV
    Allow OS:   +200 // Only use for signal samples!! Use for dm<= 50 GeV
    SS only or fake events from data otherwise

    useISR = True if int(ch/10)==1 else False
    '''

    if int(int(ch%100)/10)==0:
      whichISR=nonisrCut
    elif int(int(ch%100)/10)==1:
      whichISR=isrCut

    if ch%10==0:
      lepFlav = "%s && %s" % (eeCut, zMassCut)
    elif ch%10==1:
      lepFlav = emuCut
    elif ch%10==2:
      lepFlav = "%s && %s" % (mumuCut, zMassCut)
    elif ch%10==3:
      lepFlav = "(%s || %s)" % (emuCut, zMassCut)
    elif ch%10==4:
      lepFlav = "(%s || %s) && %s" % (eeCut, mumuCut, zMassCut)


    if int(ch/100)==1: # Allow opposite-sign different flavor
      lepFlav = "(%s) && (%s || %s)" % (lepFlav, ssCut, emuCut)
    if int(ch/100)==2: # Allow both SS and OS for all flavors
      lepFlav = lepFlav
    else: # Same-sign only or OS fake events from data
      lepFlav = "(%s) && (%s || (!isMC && (qFwt+fLwt)!=0))" % (lepFlav, ssCut)

    # lepFlav = "1"

  elif type(ch) is bool:
    whichISR = isrCut if ch else nonisrCut

  myCut = "&&".join(["(%s)"%cut for cut in [trigCut, whichISR, sigLepCut, exact2LepCut, lepFlav]])
 
  # return zMassCut
  return myCut
