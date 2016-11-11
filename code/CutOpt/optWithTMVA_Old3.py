# !!! TMVA at ROOT5.34 is used, seems not work at ROOT6

import sys 
import ROOT
import re

def setupTMVA(pathToTMVA):
  ROOT.gROOT.SetMacroPath( pathToTMVA )
  ROOT.gROOT.Macro       ( "TMVAlogon.C" )    
  ROOT.gROOT.LoadMacro   ( "TMVAGui.C" )

def setupXsecDB(pathToSUSYTools):
  ROOT.gInterpreter.AddIncludePath(pathToSUSYTools)
  ROOT.gROOT.ProcessLine( ".x $ROOTCOREDIR/scripts/load_packages.C")
  ROOT.gROOT.ProcessLine( ".L " + pathToSUSYTools + "Root/SUSYCrossSection.cxx")

def getXSECxEff( db, target, chID = 0):
  filename = target.split("/")[-1]
  runID = int(filename.split(".")[4])
  if chID==0 :
    return db.xsectTimesEff(runID)
  else:
    return db.xsectTimesEff(runID, chID)

#pathToTMVA = "/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.02.12-x86_64-slc6-gcc48-opt/tmva/test"
pathToTMVA = "/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.12-x86_64-slc6-gcc49-opt/tmva/test"
#pathToTMVA = "/home/dayabay/ROOT/root_5.34.21/tmva/test/"
pathToSUSYTools = "/home/ytchan/SUSY_RCArea_Old3/SUSYTools/"

tarLumi = 3248.28; #unit = pb-1
##----------------------------------------------------------------
##only work at ROOT6
##use this to obtain the XSECxEff for bkg & sigList
##----------------------------------------------------------------
#setupXsecDB(pathToSUSYTools)
#from ROOT.SUSY import CrossSectionDB
#xsecDB = CrossSectionDB(pathToSUSYTools + "data/mc15_13TeV")
#for (fname, entry, x) in bkgList:
#  print getXSECxEff(xsecDB, fname)
#for (fname, entry, x) in sigList:
#  print getXSECxEff(xsecDB, fname, 125)
#sys.exit()

#----------------------- Start Copy from TmMVAClassification.py ---------------------------------


# --------------------------------------------
# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser

# --------------------------------------------

# Default settings for command line arguments
DEFAULT_OUTFNAME = "cutOpt_TMVA.root"
DEFAULT_INFNAME  = "dummy.root"
DEFAULT_BKGLIST  = "bkgList.txt"
DEFAULT_SIGLIST  = "sigList.txt"
DEFAULT_ANATYPE  = "plain"
DEFAULT_TREESIG  = "NOTUSED"
DEFAULT_TREEBKG  = "NOTUSED"
DEFAULT_METHODS  = "BDTD"
#DEFAULT_METHODS  = "Cuts,CutsD,CutsPCA,CutsGA,CutsSA,Likelihood,LikelihoodD,LikelihoodPCA,LikelihoodKDE,LikelihoodMIX,PDERS,PDERSD,PDERSPCA,PDEFoam,PDEFoamBoost,KNN,LD,Fisher,FisherG,BoostedFisher,HMatrix,FDA_GA,FDA_SA,FDA_MC,FDA_MT,FDA_GAMT,FDA_MCMT,MLP,MLPBFGS,MLPBNN,CFMlpANN,TMlpANN,SVM,BDT,BDTD,BDTG,BDTB,RuleFit"

# Print usage help
def usage():
    print " "
    print "Usage: python %s [options]" % sys.argv[0]
    print "  -m | --methods    : gives methods to be run (default: all methods)"
    print "  -i | --inputfile  : name of input ROOT file (default: '%s')" % DEFAULT_INFNAME
    print "  -t | --inputtrees : input ROOT Trees for signal and background (default: '%s %s')" \
          % (DEFAULT_TREESIG, DEFAULT_TREEBKG)
    print "  -b | --bkgList    : name of bkg ROOT file list (default: '%s')" % DEFAULT_BKGLIST
    print "  -s | --sigList    : name of sig ROOT file list (default: '%s')" % DEFAULT_SIGLIST
    print "  -a | --anaType    : type of analysis: plain/doISR (default: '%s')" % DEFAULT_ANATYPE
    print "  -o | --outputfile : name of output ROOT file containing results (default: '%s')" % DEFAULT_OUTFNAME
    print "  -v | --verbose"
    print "  -? | --usage      : print this help message"
    print "  -h | --help       : print this help message"
    print " "

# Main routine
def main():

    try:
        # retrive command line options
        shortopts  = "m:i:t:b:s:a:o:vh?"
        longopts   = ["methods=", "inputfile=", "inputtrees=", "bkgList=", "sigList=", "anaType=", "outputfile=", "verbose", "help", "usage"]
        opts, args = getopt.getopt( sys.argv[1:], shortopts, longopts )

    except getopt.GetoptError:
        # print help information and exit:
        print "ERROR: unknown options in argument %s" % sys.argv[1:]
        usage()
        sys.exit(1)

    infname     = DEFAULT_INFNAME
    bkgList     = DEFAULT_BKGLIST
    sigList     = DEFAULT_SIGLIST
    anaType     = DEFAULT_ANATYPE
    treeNameSig = DEFAULT_TREESIG
    treeNameBkg = DEFAULT_TREEBKG
    outfname    = DEFAULT_OUTFNAME
    methods     = DEFAULT_METHODS
    verbose     = False
    for o, a in opts:
        if o in ("-?", "-h", "--help", "--usage"):
            usage()
            sys.exit(0)
        elif o in ("-m", "--methods"):
            methods = a
        elif o in ("-i", "--inputfile"):
            infname = a
        elif o in ("-b", "--bkgList"):
            bkgList = a
        elif o in ("-s", "--sigList"):
            sigList = a
        elif o in ("-a", "--anaType"):
            anaType = a
        elif o in ("-o", "--outputfile"):
            outfname = a
        elif o in ("-t", "--inputtrees"):
            a.strip()
            trees = a.rsplit( ' ' )
            trees.sort()
            trees.reverse()
            if len(trees)-trees.count('') != 2:
                print "ERROR: need to give two trees (each one for signal and background)"
                print trees
                sys.exit(1)
            treeNameSig = trees[0]
            treeNameBkg = trees[1]
        elif o in ("-v", "--verbose"):
            verbose = True

    # Print methods
    mlist = methods.replace(' ',',').split(',')
    print "=== TMVAClassification: use method(s)..."
    for m in mlist:
        if m.strip() != '':
            print "=== - <%s>" % m.strip()

    # Import ROOT classes
    from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut
    
    # check ROOT version, give alarm if 5.18 
    if gROOT.GetVersionCode() >= 332288 and gROOT.GetVersionCode() < 332544:
        print "*** You are running ROOT version 5.18, which has problems in PyROOT such that TMVA"
        print "*** does not run properly (function calls with enums in the argument are ignored)."
        print "*** Solution: either use CINT or a C++ compiled version (see TMVA/macros or TMVA/examples),"
        print "*** or use another ROOT version (e.g., ROOT 5.19)."
        sys.exit(1)
    
    #Edited
    setupTMVA(pathToTMVA)
    
    # Import TMVA classes from ROOT
    from ROOT import TMVA

    # Output file
    outputFile = TFile( outfname, 'RECREATE' )
    
    # Create instance of TMVA factory (see TMVA/macros/TMVAClassification.C for more factory options)
    # All TMVA output can be suppressed by removing the "!" (not) in 
    # front of the "Silent" argument in the option string
    suffix = outfname.split("/")[-1].split(".")[0]
    factory = TMVA.Factory( "TMVAClassification_%s"%suffix, outputFile, 
                            "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P:AnalysisType=Classification" )

    # Set verbosity
    factory.SetVerbose( verbose )
    
    # If you wish to modify default settings 
    # (please check "src/Config.h" to see all available global options)
    #    gConfig().GetVariablePlotting()).fTimesRMS = 8.0
    #    gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory"

    # Define the input variables that shall be used for the classifier training
    # note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    # [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

    #Edited
    #follow 2LSS note Ch5.2: Disciminating variables
    factory.AddVariable( "mT2    := sig.mT2"                      ,  'F' )
    factory.AddVariable( "pt     := l12.pt"                       ,  'F' )
    factory.AddVariable( "MET    := sig.MetRel"                   ,  'F' )
    factory.AddVariable( "Ht     := Sum$(jets.pt) + Sum$(leps.pt)",  'F' )
    factory.AddVariable( "mTl1   := leps.mT[0]"                   ,  'F' )
    factory.AddVariable( "mTl2   := leps.mT[1]"                   ,  'F' )
    factory.AddVariable( "ll_dPhi:= l12.dPhi"                     ,  'F' )
    factory.AddVariable( "l12m   := (int(abs(leps.ID[0]))!=int(abs(leps.ID[1])))*100 + l12.m", 'F')
    
    #ISR region
    if (anaType=="doISR"):
      factory.AddVariable( "JetMET_dPhi  := jets.MET_dPhi[0]"           ,  'F' )
      factory.AddVariable( "MET_JetPt_R  := sig.MetRel/jets.pt[0]"      ,  'F' )
      factory.AddVariable( "l1Pt_JetPt_R := leps.pt[0]/jets.pt[0]"      ,  'F' )

    #factory.AddSpectator( "pt1 := leps.pt[0]" , 'F' )
    #factory.AddSpectator( "pt2 := leps.pt[1]" , 'F' )
    #factory.AddSpectator( "ID1 := int(leps.ID[0])" , 'I' )
    #factory.AddSpectator( "ID2 := int(leps.ID[1])" , 'I' )
    #factory.AddSpectator( "nCentralJets := Sum$(jets.pt>20 && abs(jets.eta)<2.4)" , 'I' )

    #FIXME
    setupXsecDB(pathToSUSYTools)
    from ROOT.SUSY import CrossSectionDB
    xsecDB = CrossSectionDB(pathToSUSYTools + "data/mc15_13TeV/")

    #read in training data
    openedInFileList = []

    # Read input sig
    sigList = open(sigList,"r")
    for infname in sigList:
      inFile = TFile.Open( infname[:-1] )
      openedInFileList.append(inFile)

      hCutFlow = inFile.FindObjectAny("hCutFlow")
      mcEntry = hCutFlow.GetBinContent(1)

      #FIXME: hard coded extract runNum from filePath
      m = re.match(".*\.([0-9]{6})\..*", infname)
      runNum = int(m.groups()[0])
      xSECxEff = xsecDB.xsectTimesEff(runNum, 125) + xsecDB.xsectTimesEff(runNum, 127) #125,127 is channel no.

      # Get the trees for training
      signal  = inFile.Get( "Data_" )
    
      # Global event weights (see below for setting event-wise weights)
      #signalWeight = getXSECxEff(xsecDB, infname) * tarLumi / mcEntry
      #signalWeight = xSECxEff * tarLumi / mcEntry
      signalWeight = 1.0 * tarLumi / mcEntry #treat diff SUSY scenario with equal weight
      if signalWeight<=0 : 
        print "Encounter <=0 weight sample %s , skipped" % infname
        continue

      print "mc sig ", runNum, mcEntry, xSECxEff
      factory.AddSignalTree( signal, signalWeight )
    sigList.close()

    # Read input bkg
    bkgList = open(bkgList,"r")
    for infname in bkgList:
      inFile = TFile.Open( infname[:-1] )
      openedInFileList.append(inFile)

      if "physics" in infname:
        #its real data
	print "data bkg", infname[:-1]
        background  = inFile.Get( "CFlip_" )
	if background: factory.AddBackgroundTree( background, 1.0 )
        background  = inFile.Get( "FakeLep_" )
        if background: factory.AddBackgroundTree( background, 1.0 )
      else:
        #its MC data
        hCutFlow = inFile.FindObjectAny("hCutFlow")
        mcEntry = hCutFlow.GetBinContent(1)

        #FIXME: hard coded extract runNum from filePath
        m = re.match(".*\.([0-9]{6})\..*", infname)
        runNum = int(m.groups()[0])
        xSECxEff = xsecDB.xsectTimesEff(runNum)

        # Get  trees for training
        background  = inFile.Get( "Data_" )
      
        # Global event weights (see below for setting event-wise weights)
        backgroundWeight = xSECxEff * tarLumi / mcEntry
        if backgroundWeight<=0 : 
          print "Encounter <=0 weight sample %s , skipped" % infname

        print "mc bkg ", runNum, mcEntry, xSECxEff
        factory.AddBackgroundTree( background, backgroundWeight )
    bkgList.close()


    # event-wise weights
    #factory.SetSignalWeightExpression( "weight" )
    #factory.SetBackgroundWeightExpression( "weight" )
    factory.SetSignalWeightExpression( "ElSF*MuSF" )
    factory.SetBackgroundWeightExpression( "(CFlipWeight0*FakeLepWeight0)!=1.0 ? CFlipWeight0*FakeLepWeight0  : !TMath::IsNaN(weight)? ElSF*MuSF*weight: 0.0" )

    # Apply additional cuts on the signal and background sample. 
    # example for cut: mycut = TCut( "abs(var1)<0.5 && abs(var2-0.5)<1" )
    # trigCut   = "sig.trigCode!=0"

    #"HLT_mu24_iloose_L1MU15" for mumu emu, "HLT_e24_lhmedium_iloose_L1EM20VH" for ee
    trigCut   = "((nMu>0) && (sig.trigCode & (1<<2))) || ((nMu==0) && (sig.trigCode & (1<<26)))" 

    grlCut    = "evtInfo.passGRL==1"
    wCut      = "weight>0 && weight<1e9"
    tauCut    = "1" # "nTau==0"  FIXME nTau not properly filled in NTUP yet..
    bjetCut   = "Sum$(jets.isBJet)==0"
    cosmicCut = "Sum$(leps.isCosmic)==0"

    htCut     = "(Sum$(jets.pt) + Sum$(leps.pt))>40"
    posWCut   = "FakeLepWeight0>0"

    isrCut = "Sum$(jets.pt>20 && abs(jets.eta)<2.4) %s" % (">0" if anaType=="doISR" else "==0") #nCentralJets>0 or ==0
    zMassCut  = "!(int(abs(leps.ID[0])) == int(abs(leps.ID[1])) && fabs(l12.m - 91.1876)<=5)"

    #commonCut = "&&".join(["(%s)"%cut for cut in [trigCut , grlCut , bjetCut, cosmicCut]])
    commonCut = "&&".join(["(%s)"%cut for cut in [trigCut , grlCut , wCut , zMassCut , isrCut, tauCut, bjetCut, cosmicCut]])
    commonCut = TCut(commonCut)

    sigCut = "&&".join(["(%s)"%cut for cut in [trigCut , grlCut ,        zMassCut , isrCut, tauCut, bjetCut, cosmicCut]])
    sigCut = TCut(sigCut)

    bkgCut = "&&".join(["(%s)"%cut for cut in [trigCut , grlCut , wCut , zMassCut , isrCut, tauCut, bjetCut, cosmicCut, posWCut]])
    bkgCut = TCut(bkgCut)
    # Here, the relevant variables are copied over in new, slim trees that are
    # used for TMVA training and testing
    # "SplitMode=Random" means that the input events are randomly shuffled before
    # splitting them into training and test samples

    factory.PrepareTrainingAndTestTree( sigCut, bkgCut,
    "nTrain_Signal=0:nTrain_Background=0:nTest_Background=0:SplitMode=Random:NormMode=EqualNumEvents:!V" )
                                        #"nTrain_Signal=0:nTrain_Background=2000:SplitMode=Random:NormMode=EqualNumEvents:!V" )

    # --------------------------------------------------------------------------------------------------

    # ---- Book MVA methods
    #
    # please lookup the various method configuration options in the corresponding cxx files, eg:
    # src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
    # it is possible to preset ranges in the option string in which the cut optimisation should be done:
    # "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

    # Cut optimisation
    if "Cuts" in mlist:
        factory.BookMethod( TMVA.Types.kCuts, "Cuts",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" )

    if "CutsD" in mlist:
        factory.BookMethod( TMVA.Types.kCuts, "CutsD",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" )

    if "CutsPCA" in mlist:
        factory.BookMethod( TMVA.Types.kCuts, "CutsPCA",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" )

    if "CutsGA" in mlist:
        factory.BookMethod( TMVA.Types.kCuts, "CutsGA",
                            "H:!V:FitMethod=GA:VarProp=FSmart:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" )

    if "CutsSA" in mlist:
        factory.BookMethod( TMVA.Types.kCuts, "CutsSA",
                            "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" )

    # Likelihood ("naive Bayes estimator")
    if "Likelihood" in mlist:
        factory.BookMethod( TMVA.Types.kLikelihood, "Likelihood",
                            "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" )

    # Decorrelated likelihood
    if "LikelihoodD" in mlist:
        factory.BookMethod( TMVA.Types.kLikelihood, "LikelihoodD",
                            "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" )

    # PCA-transformed likelihood
    if "LikelihoodPCA" in mlist:
        factory.BookMethod( TMVA.Types.kLikelihood, "LikelihoodPCA",
                            "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ) 

    # Use a kernel density estimator to approximate the PDFs
    if "LikelihoodKDE" in mlist:
        factory.BookMethod( TMVA.Types.kLikelihood, "LikelihoodKDE",
                            "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ) 

    # Use a variable-dependent mix of splines and kernel density estimator
    if "LikelihoodMIX" in mlist:
        factory.BookMethod( TMVA.Types.kLikelihood, "LikelihoodMIX",
                            "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ) 

    # Test the multi-dimensional probability density estimator
    # here are the options strings for the MinMax and RMS methods, respectively:
    #      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
    #      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
    if "PDERS" in mlist:
        factory.BookMethod( TMVA.Types.kPDERS, "PDERS",
                            "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" )

    if "PDERSD" in mlist:
        factory.BookMethod( TMVA.Types.kPDERS, "PDERSD",
                            "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" )

    if "PDERSPCA" in mlist:
        factory.BookMethod( TMVA.Types.kPDERS, "PDERSPCA",
                             "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" )

   # Multi-dimensional likelihood estimator using self-adapting phase-space binning
    if "PDEFoam" in mlist:
        factory.BookMethod( TMVA.Types.kPDEFoam, "PDEFoam",
                            "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" )

    if "PDEFoamBoost" in mlist:
        factory.BookMethod( TMVA.Types.kPDEFoam, "PDEFoamBoost",
                            "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" )

    # K-Nearest Neighbour classifier (KNN)
    if "KNN" in mlist:
        factory.BookMethod( TMVA.Types.kKNN, "KNN",
                            "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" )

    # H-Matrix (chi2-squared) method
    if "HMatrix" in mlist:
        factory.BookMethod( TMVA.Types.kHMatrix, "HMatrix", "!H:!V" )

    # Linear discriminant (same as Fisher discriminant)
    if "LD" in mlist:
        factory.BookMethod( TMVA.Types.kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" )

    # Fisher discriminant (same as LD)
    if "Fisher" in mlist:
        factory.BookMethod( TMVA.Types.kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" )

    # Fisher with Gauss-transformed input variables
    if "FisherG" in mlist:
        factory.BookMethod( TMVA.Types.kFisher, "FisherG", "H:!V:VarTransform=Gauss" )

    # Composite classifier: ensemble (tree) of boosted Fisher classifiers
    if "BoostedFisher" in mlist:
        factory.BookMethod( TMVA.Types.kFisher, "BoostedFisher", 
                            "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2" )

    # Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
    if "FDA_MC" in mlist:
        factory.BookMethod( TMVA.Types.kFDA, "FDA_MC",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

    if "FDA_GA" in mlist:
        factory.BookMethod( TMVA.Types.kFDA, "FDA_GA",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

    if "FDA_SA" in mlist:
        factory.BookMethod( TMVA.Types.kFDA, "FDA_SA",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

    if "FDA_MT" in mlist:
        factory.BookMethod( TMVA.Types.kFDA, "FDA_MT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

    if "FDA_GAMT" in mlist:
        factory.BookMethod( TMVA.Types.kFDA, "FDA_GAMT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

    if "FDA_MCMT" in mlist:
        factory.BookMethod( TMVA.Types.kFDA, "FDA_MCMT",
                            "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1)(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

    # TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
    if "MLP" in mlist:
        factory.BookMethod( TMVA.Types.kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" )

    if "MLPBFGS" in mlist:
        factory.BookMethod( TMVA.Types.kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" )

    if "MLPBNN" in mlist:
        factory.BookMethod( TMVA.Types.kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ) # BFGS training with bayesian regulators

    # CF(Clermont-Ferrand)ANN
    if "CFMlpANN" in mlist:
        factory.BookMethod( TMVA.Types.kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ) # n_cycles:#nodes:#nodes:...  

    # Tmlp(Root)ANN
    if "TMlpANN" in mlist:
        factory.BookMethod( TMVA.Types.kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ) # n_cycles:#nodes:#nodes:...

    # Support Vector Machine
    if "SVM" in mlist:
        factory.BookMethod( TMVA.Types.kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" )

    # Boosted Decision Trees
    if "BDTG" in mlist:
        factory.BookMethod( TMVA.Types.kBDT, "BDTG",
                            "!H:!V:NTrees=1000:MinNodeSize=1.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" )                        

    if "BDT" in mlist:
        factory.BookMethod( TMVA.Types.kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" )

    if "BDTB" in mlist:
        factory.BookMethod( TMVA.Types.kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" )

    if "BDTD" in mlist:
        factory.BookMethod( TMVA.Types.kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=2:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" )

    # RuleFit -- TMVA implementation of Friedman's method
    if "RuleFit" in mlist:
        factory.BookMethod( TMVA.Types.kRuleFit, "RuleFit",
                            "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" )

    # --------------------------------------------------------------------------------------------------
            
    # ---- Now you can tell the factory to train, test, and evaluate the MVAs. 

    # Train MVAs
    factory.TrainAllMethods()
    
    # Test MVAs
    factory.TestAllMethods()
    
    # Evaluate MVAs
    factory.EvaluateAllMethods()    
    
    # Save the output.
    outputFile.Close()
    
    print "=== wrote root file %s\n" % outfname
    print "=== TMVAClassification is done!\n"
    
    # open the GUI for the result macros    
    #gROOT.ProcessLine( "TMVAGui(\"%s\")" % outfname )
    
    # keep the ROOT thread running
    #gApplication.Run() 

# ----------------------------------------------------------

if __name__ == "__main__":
    main()
