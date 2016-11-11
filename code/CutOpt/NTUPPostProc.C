/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TMath.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "SUSYTools/SUSYCrossSection.h"

using namespace TMVA;

void NTUPPostProc( TString inFileName = "dummy.root" ) 
{   

   //---------------------------------------------------------------

   // This loads the library
   TMVA::Tools::Instance();

   //xsec DB
   SUSY::CrossSectionDB *xsecDB(0);
   xsecDB = new SUSY::CrossSectionDB(gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc15_13TeV/"));

   std::cout << std::endl;
   std::cout << "==> Start NTUPPostProc" << std::endl;

   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

   TMVA::Reader *reader    = new TMVA::Reader( "!Color:!Silent" );    
   TMVA::Reader *readerISR = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t mT2, pt, MET, Ht, mTl1, mTl2, ll_dPhi, l12m;
   Float_t JetMET_dPhi, MET_JetPt_R, l1Pt_JetPt_R;

   reader->AddVariable( "mT2    := sig.mT2"                      , &mT2     );
   reader->AddVariable( "pt     := l12.pt"                       , &pt      );
   reader->AddVariable( "MET    := sig.MetRel"                   , &MET     );
   reader->AddVariable( "Ht     := Sum$(jets.pt) + Sum$(leps.pt)", &Ht      );
   reader->AddVariable( "mTl1   := leps.mT[0]"                   , &mTl1    );
   reader->AddVariable( "mTl2   := leps.mT[1]"                   , &mTl2    );
   reader->AddVariable( "ll_dPhi:= l12.dPhi"                     , &ll_dPhi );
   reader->AddVariable( "l12m   := (int(abs(leps.ID[0]))!=int(abs(leps.ID[1])))*100 + l12.m", &l12m);
    
   //ISR region
   readerISR->AddVariable( "mT2          := sig.mT2"                      , &mT2     );
   readerISR->AddVariable( "pt           := l12.pt"                       , &pt      );
   readerISR->AddVariable( "MET          := sig.MetRel"                   , &MET     );
   readerISR->AddVariable( "Ht           := Sum$(jets.pt) + Sum$(leps.pt)", &Ht      );
   readerISR->AddVariable( "mTl1         := leps.mT[0]"                   , &mTl1    );
   readerISR->AddVariable( "mTl2         := leps.mT[1]"                   , &mTl2    );
   readerISR->AddVariable( "ll_dPhi      := l12.dPhi"                     , &ll_dPhi );
   readerISR->AddVariable( "l12m   := (int(abs(leps.ID[0]))!=int(abs(leps.ID[1])))*100 + l12.m", &l12m);
   readerISR->AddVariable( "JetMET_dPhi  := jets.MET_dPhi[0]"             , &JetMET_dPhi  );
   readerISR->AddVariable( "MET_JetPt_R  := sig.MetRel/jets.pt[0]"        , &MET_JetPt_R  );
   readerISR->AddVariable( "l1Pt_JetPt_R := leps.pt[0]/jets.pt[0]"        , &l1Pt_JetPt_R );

   // Spectator variables declared in the training have to be added to the reader, too
   //Float_t spec1,spec2;
   //reader->AddSpectator( "spec1 := var1*2",   &spec1 );
   //reader->AddSpectator( "spec2 := var1*3",   &spec2 );


   // --- Book the MVA methods

   TString dir    = "weights/";
   TString prefix = "TMVAClassification";
   TString methodName = "";
   TString weightfile = "";

   // Book method(s)
   methodName = "BDT_dM50";
   weightfile = dir + prefix + TString("_Slep_dM50_trainResult_BDTD") + TString(".weights.xml");
   reader->BookMVA( methodName, weightfile ); 

   methodName = "BDT_dM50_ISR";
   weightfile = dir + prefix + TString("_Slep_dM50_ISR_trainResult_BDTD") + TString(".weights.xml");
   readerISR->BookMVA( methodName, weightfile ); 

   methodName = "BDT_dM100";
   weightfile = dir + prefix + TString("_Slep_newSample_trigCut_za_trainResult_BDTD") + TString(".weights.xml");
   reader->BookMVA( methodName, weightfile ); 

   methodName = "BDT_dM100_ISR";
   weightfile = dir + prefix + TString("_Slep_newSample_trigCut_za_ISR_trainResult_BDTD") + TString(".weights.xml");
   readerISR->BookMVA( methodName, weightfile ); 
   
   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   

   bool isSig  = false;
   bool isData = false; 
   
   TString sName = "unknown"; //sample name
   if (inFileName.Index(  "enugamma"  )!=-1) sName="wgamma_";
   if (inFileName.Index( "munugamma"  )!=-1) sName="wgamma_";
   if (inFileName.Index("taunugamma"  )!=-1) sName="wgamma_";
   if (inFileName.Index("llll"        )!=-1) sName="diboson_";
   if (inFileName.Index("lllv"        )!=-1) sName="diboson_";
   if (inFileName.Index("llvv"        )!=-1) sName="diboson_";
   if (inFileName.Index("ttW"         )!=-1) sName="ttbar_";
   if (inFileName.Index("ttee"        )!=-1) sName="ttbar_";
   if (inFileName.Index("ttmumu"      )!=-1) sName="ttbar_";
   if (inFileName.Index("tttautau"    )!=-1) sName="ttbar_";
   if (inFileName.Index("physics_Main")!=-1) {sName=""; isData=true;}

   int matchIdx; Ssiz_t matchLen;
   //match WZ and Slep samples eg MGPy8EG_A14N23LO_C1N2_WZ_300p0_250p0_3L_2L7_myOutput
   matchIdx = inFileName.Index( TRegexp("[a-zA-Z0-9_]*C1N2[a-zA-Z0-9_]*") , &matchLen);
   if (matchIdx!=-1 && matchLen>13){
     sName = inFileName(matchIdx, matchLen-9)+"_"; //-9 to remove _myOutput
     isSig = true;
   }

   //mc sample ID
   int mcSampleID = 0;
   matchIdx = inFileName.Index( TRegexp(".[0-9][0-9][0-9][0-9][0-9][0-9].") , &matchLen); //root's re seems broken..?
   if (matchIdx!=-1 && matchLen==8){
     mcSampleID = TString(inFileName(matchIdx+1, 6)).Atoi();
   }

   if (sName=="unknown"){
      std::cout << "ERROR: could not infer sample type from data file path name" << std::endl;
      exit(1);
   }

   if (!isData && mcSampleID==0){
      std::cout << "ERROR: could not infer mc sample ID from data file path name" << std::endl;
      exit(1);
   }

   TFile *input = new TFile( inFileName, "UPDATE");
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

   int mcEntry = 1;
   TH1F* hCutFlow = (TH1F*)input->FindObjectAny("hCutFlow");
   if (hCutFlow) mcEntry = hCutFlow->GetBinContent(1);


   TIter next(input->GetListOfKeys());
   TKey* key;
   int nKeys = input->GetNkeys();
   for (int iKey=0; iKey<nKeys; iKey++){
     key = (TKey *) next();
     if (!(TString(key->GetClassName())=="TTree")) continue;

     TString treeName = key->GetName();
     bool accept = false;

     if (sName==""){
       //real data
       if (treeName=="Data_"   ) accept = true;
       if (treeName=="CFlip_"  ) accept = true;
       if (treeName=="FakeLep_") accept = true;
     }else if (isSig){
       //mc signal
       if (treeName=="Data_"   ) accept = true;
     }else{
       //mc bkg
       if (treeName.Index("Data_")==0) accept = true;
     }

     if (!accept) continue; //only cal BDT on accepted trees

     //-----------------------------------------------------------
     // Prepare the event tree
     //-----------------------------------------------------------
     // - here the variable names have to corresponds to your tree
     // - you can use the same variables as above which is slightly faster,
     //   but of course you can use different ones and copy the values inside the event loop
     //
     std::cout << "--- Select signal sample " << treeName << std::endl;
     TTree* theTree = (TTree*)key->ReadObj();

     // FIXME: this again have to be manually sync-ed with reader var setup
     TTreeFormula fmT2    ( "mT2"     , "sig.mT2"                      , theTree );
     TTreeFormula fpt     ( "pt"      , "l12.pt"                       , theTree );
     TTreeFormula fMET    ( "MET"     , "sig.MetRel"                   , theTree );
     TTreeFormula fHt     ( "Ht"      , "Sum$(jets.pt) + Sum$(leps.pt)", theTree );
     TTreeFormula fmTl1   ( "mTl1"    , "leps.mT[0]"                   , theTree );
     TTreeFormula fmTl2   ( "mTl2"    , "leps.mT[1]"                   , theTree );
     TTreeFormula fll_dPhi( "ll_dPhi" , "l12.dPhi"                     , theTree );
     TTreeFormula fl12m   ( "l12m"    , "(int(abs(leps.ID[0]))!=int(abs(leps.ID[1])))*100 + l12.m", theTree);


     //ISR region
     TTreeFormula fJetMET_dPhi ( "JetMET_dPhi"  , "jets.MET_dPhi[0]"                     , theTree );
     TTreeFormula fMET_JetPt_R ( "MET_JetPt_R"  , "sig.MetRel/jets.pt[0]"                , theTree );
     TTreeFormula fl1Pt_JetPt_R( "l1Pt_JetPt_R" , "leps.pt[0]/jets.pt[0]"                , theTree );

     TTreeFormula fisrCut      ( "isrCut"       , "Sum$(jets.pt>20 && abs(jets.eta)<2.4)", theTree );

     TTreeFormula fweight      ( "weight"       , "weight"                               , theTree );
     TTreeFormula fweightU     ( "weightU"      , "PRW_DATASF__1up"                      , theTree );
     TTreeFormula fweightD     ( "weightD"      , "PRW_DATASF__1down"                    , theTree );

     //-----------------------------------------------------------
     // Prepare the output tree
     //-----------------------------------------------------------
     TTree* outTree = new TTree( "BDT_"+sName+theTree->GetName(), "Added extra var and AddFriend-ed original tree");
     bool isISR;
     float BDT_dM50_output(-100), BDT_dM100_output(-100);
     float genWeight(1.0), genWeightU(1.0), genWeightD(1.0);
     float weight(1.0), weightU(1.0), weightD(1.0);
     outTree->AddFriend(theTree);
     outTree->Branch("isISR", &isISR);
     outTree->Branch("BDT_dM50_output" , &BDT_dM50_output );
     outTree->Branch("BDT_dM100_output", &BDT_dM100_output);
     outTree->Branch("genWeight"       , &genWeight );
     outTree->Branch("genWeight__1up"  , &genWeightU);
     outTree->Branch("genWeight__1down", &genWeightD);

     //tmp measure to mask -NAN weight...
     outTree->Branch("weight"           , &weight);
     outTree->Branch("PRW_DATASF__1up"  , &weightU);
     outTree->Branch("PRW_DATASF__1down", &weightD);

     if (!isData){
       if (isSig){
         std::cout << "genW " 
		   << inFileName << " "
		   << xsecDB->rawxsect     (mcSampleID, 125) << " "
		   << xsecDB->kfactor      (mcSampleID, 125) << " "
		   << xsecDB->efficiency   (mcSampleID, 125) << " "
		   << xsecDB->xsectTimesEff(mcSampleID, 125) << " "
		   << xsecDB->rawxsect     (mcSampleID, 127) << " "
		   << xsecDB->kfactor      (mcSampleID, 127) << " "
		   << xsecDB->efficiency   (mcSampleID, 127) << " "
		   << xsecDB->xsectTimesEff(mcSampleID, 127) << " "
		   << mcEntry                                << " "
		   << theTree->GetEntries()                  << " "
		   << std::endl;
         genWeight  = (xsecDB->xsectTimesEff(mcSampleID, 125) + 
                       xsecDB->xsectTimesEff(mcSampleID, 127) ) * 1.0 / mcEntry; //125,127 is channel number, 1.0 is tarLumi pb-1
	 genWeightU = genWeight;
	 genWeightD = genWeight;
       }else{
         genWeight  = xsecDB->xsectTimesEff(mcSampleID)* 1.0 / mcEntry;
	 genWeightU = genWeight;
	 genWeightD = genWeight;
       }
     }

     //if (!isData){
     //  if (!isSig){
     //    std::cout << "genW " 
     //   	   << inFileName << " "
     //   	   << xsecDB->rawxsect     (mcSampleID) << " "
     //   	   << xsecDB->kfactor      (mcSampleID) << " "
     //   	   << xsecDB->efficiency   (mcSampleID) << " "
     //   	   << xsecDB->xsectTimesEff(mcSampleID) << " "
     //   	   << mcEntry                           << " "
     //   	   << theTree->GetEntries()             << " "
     //   	   << std::endl;
     //  }
     //}
     //break;

     //-----------------------------------------------------------
     // --- Event loop
     //-----------------------------------------------------------
     std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
     TStopwatch sw;
     sw.Start();
     for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
     //for (Long64_t ievt=0; ievt<-1;ievt++) {

	//if (ievt==10) break;
        if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

        theTree->GetEntry(ievt);

        // This .GetNdata() is needed else fmTl2 give 0
        // strange ROOT problem during access of array like data eg. leps.mT[1]
        // https://root.cern.ch/phpBB3/viewtopic.php?t=16775
        fmT2    .GetNdata();
        fpt     .GetNdata();
        fMET    .GetNdata();
        fHt     .GetNdata();
        fmTl1   .GetNdata();
        fmTl2   .GetNdata();
        fll_dPhi.GetNdata();
        fl12m   .GetNdata();

        fJetMET_dPhi .GetNdata();
        fMET_JetPt_R .GetNdata();
        fl1Pt_JetPt_R.GetNdata();
                                ;
        fisrCut      .GetNdata();
                                ;
        fweight      .GetNdata();
        fweightU     .GetNdata();
        fweightD     .GetNdata();


	mT2     = fmT2    .EvalInstance();
        pt      = fpt     .EvalInstance();
        MET     = fMET    .EvalInstance();
        Ht      = fHt     .EvalInstance();
        mTl1    = fmTl1   .EvalInstance();
        mTl2    = fmTl2   .EvalInstance();
        ll_dPhi = fll_dPhi.EvalInstance();
        l12m    = fl12m   .EvalInstance();

	isISR = fisrCut .EvalInstance()>0.;
	if (isISR>0.01){
          JetMET_dPhi  = fJetMET_dPhi .EvalInstance();
          MET_JetPt_R  = fMET_JetPt_R .EvalInstance();
          l1Pt_JetPt_R = fl1Pt_JetPt_R.EvalInstance();
	}

	weight  = fweight .EvalInstance();
	weightU = fweightU.EvalInstance();
	weightD = fweightD.EvalInstance();
	//catch nan and change to to large -ve value
	if (TMath::IsNaN(weight)){
	  weight  = -1000;
	  weightU = -1000;
	  weightD = -1000;
	}


        // --- Return the MVA outputs and fill into histograms
        if (isISR){
          BDT_dM50_output  = readerISR->EvaluateMVA( "BDT_dM50_ISR"  );
          BDT_dM100_output = readerISR->EvaluateMVA( "BDT_dM100_ISR" );
	}else{
          BDT_dM50_output  = reader->EvaluateMVA( "BDT_dM50"  );
          BDT_dM100_output = reader->EvaluateMVA( "BDT_dM100" );
	}

	//if (ievt==3715){
	//  theTree->Show(ievt);
	//  std::cout << "vars    " << mT2 << " " << pt << " " << MET << " " << Ht << " " << mTl1 << " " << mTl2 << " " << ll_dPhi << std::endl;
	//  std::cout << "ISRvars " << isISR << " " << JetMET_dPhi << " " << MET_JetPt_R << " " << l1Pt_JetPt_R << std::endl;
	//  std::cout << "BDTvars " << BDT_dM50_output << " " << BDT_dM100_output << std::endl;
	//}
	
	outTree->Fill();
     }

     outTree->Write();
     delete outTree;

     // Get elapsed time
     sw.Stop();
     std::cout << "--- End of event loop: "; sw.Print();

   }

   delete reader;
   delete readerISR;
    
   std::cout << "==> NTUPPostProc is done!" << std::endl << std::endl;
} 

int main( int argc, char** argv )
{
   if (argc<2){
     std::cout << "Usage:  root -l NTUPPostProc.C+('inFileName.root')" << std::endl;
   }
   NTUPPostProc(argv[1]); 
   return 0; 
}
