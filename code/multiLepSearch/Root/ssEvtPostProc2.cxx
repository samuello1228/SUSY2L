#include <TError.h>
#include <algorithm>
#include <vector>
#include <EventLoop/StatusCode.h>

#include <multiLepSearch/mCHECK.h>
#include <multiLepSearch/ssEvtPostProc2.h>
#include <multiLepSearch/obj_def.h>

#include <TH1D.h>
#include <TFile.h>
#include <TRegexp.h>
#include <TKey.h>

#include "TTreeFormula.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

using namespace std;
using namespace TMVA;

// this is needed to distribute the algorithm to the workers
ClassImp(ssEvtPostProc2)

ssEvtPostProc2 :: ssEvtPostProc2(string name):m_name(name){
  inputFile="dummy.root";
  norminalTreeName="PP1_evt2l";
  isData=true;

  //m_susyEvt = new susyEvts();
}

ssEvtPostProc2 :: ~ssEvtPostProc2(){
}


int ssEvtPostProc2 :: initialize ()
{

  // --- Create the Reader object
  reader    = new TMVA::Reader( "!Color:!Silent" );    
  readerISR = new TMVA::Reader( "!Color:!Silent" );    

  basicBDTVarValues.resize( basicBDTVarNameList.size() );
    isrBDTVarValues.resize(   isrBDTVarNameList.size() );

  // Create a set of variables and declare them to the reader
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
  char tmpStr[200];
  for (unsigned int i=0;i<basicBDTVarNameList.size();i++){
    sprintf(tmpStr, "%s := %s", basicBDTVarNameList[i].c_str(), basicBDTVarFormulaList[i].c_str() );  
    reader   ->AddVariable( tmpStr, &(basicBDTVarValues[i]) );
    readerISR->AddVariable( tmpStr, &(basicBDTVarValues[i]) );
  }
   
  //ISR variables
  for (unsigned int i=0;i<isrBDTVarNameList.size();i++){
    sprintf(tmpStr, "%s := %s", isrBDTVarNameList[i].c_str(), isrBDTVarFormulaList[i].c_str() );  
    readerISR->AddVariable( tmpStr, &(isrBDTVarValues[i]) );
  }

  // Spectator variables declared in the training have to be added to the reader, too
  //Float_t spec1,spec2;
  //reader->AddSpectator( "spec1 := var1*2",   &spec1 );
  //reader->AddSpectator( "spec2 := var1*3",   &spec2 );
  
  // --- Book the MVA methods
  // Book method(s)
  for (unsigned int i=0; i< BDTNameList.size();i++){
    if (varSetList[i] == "ISR"){
      readerISR->BookMVA( BDTNameList[i], weightFileList[i] ); 
    }else{
      reader   ->BookMVA( BDTNameList[i], weightFileList[i] ); 
    }
  }


  return 0;
}


int ssEvtPostProc2 :: runLoop ()
{
  TFile* inF = new TFile(inputFile.c_str(), "UPDATE");

  //-----------------------------------------------------------
  //Loop over objects (Trees) in the file
  //-----------------------------------------------------------
  TIter next(inF->GetListOfKeys());
  TKey* key;
  TString lastTreeName = "";
  int nKeys = inF->GetNkeys();
  for (int iKey=0; iKey<nKeys; iKey++){
    key = (TKey *) next();
    if (!(TString(key->GetClassName())=="TTree")) continue;

    TTree* rawTree = (TTree*)key->ReadObj();
    TString treeName = key->GetName();
    if (treeName == lastTreeName) continue; //FIXME only take the tree with latest cycle (assume tree with highest cycle is at top of the list)
    lastTreeName = treeName;
    if (treeName.Index(norminalTreeName)!=0) continue; //only look at norminal tree and its systematics
    
    //-----------------------------------------------------------
    // TTreeFormulars to read in BDT variables from the rawTree
    //-----------------------------------------------------------
    vector<TTreeFormula*> basicBDTVarFormulas;
    for (unsigned int i=0;i<basicBDTVarNameList.size();i++){
      basicBDTVarFormulas.push_back( new TTreeFormula( basicBDTVarNameList[i].c_str(), basicBDTVarFormulaList[i].c_str(), rawTree ) );
    }
   
    vector<TTreeFormula*> isrBDTVarFormulas;
    for (unsigned int i=0;i<isrBDTVarNameList.size();i++){
      isrBDTVarFormulas.push_back( new TTreeFormula( isrBDTVarNameList[i].c_str(), isrBDTVarFormulaList[i].c_str(), rawTree ) );
    }

    TTreeFormula fisrCut      ( "isrCut"       , "Sum$(jets.pt>20 && abs(jets.eta)<2.4)", rawTree );

    //-----------------------------------------------------------
    // Prepare the output tree
    //-----------------------------------------------------------
    //string tmpStr = outTreePrefix + "_BDT_" + rawTree->GetName();
    string tmpStr("BDT_"); tmpStr += rawTree->GetName();
    TTree* outTree = new TTree( tmpStr.c_str(), "Added extra BDT var and AddFriend-ed original tree");

    bool isISR;
    vector<float> bdtOutputs; bdtOutputs.reserve( BDTNameList.size() );

    outTree->AddFriend(rawTree);
    outTree->Branch("mcEvtW", &mcEvtW);
    outTree->Branch("isISR", &isISR);
    for (unsigned int i=0; i<BDTNameList.size(); i++){
      outTree->Branch( BDTNameList[i].c_str() , &(bdtOutputs[i]) );
    }

    //-----------------------------------------------------------
    // Loop over events in the tree
    //-----------------------------------------------------------
    std::cout << "--- Processing: " << rawTree->GetName() << " " << rawTree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();
    for (Long64_t ievt=0; ievt<rawTree->GetEntries();ievt++) {
    //for (Long64_t ievt=0; ievt<-1;ievt++)

       //if (ievt==10) break;
       if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

       rawTree->GetEntry(ievt);

       // This .GetNdata() is needed else some formula give 0
       // strange ROOT problem during access of array like data eg. leps.mT[1]
       // https://root.cern.ch/phpBB3/viewtopic.php?t=16775
       for (unsigned int i=0;i<basicBDTVarFormulas.size();i++){
         basicBDTVarFormulas[i]->GetNdata();
	 basicBDTVarValues[i] = basicBDTVarFormulas[i]->EvalInstance();
       }
                               
       fisrCut      .GetNdata();
       isISR = fisrCut .EvalInstance()>0.;

       if (isISR>0.01){
         for (unsigned int i=0;i<isrBDTVarFormulas.size();i++){
           isrBDTVarFormulas[i]->GetNdata();
	   isrBDTVarValues[i] = isrBDTVarFormulas[i]->EvalInstance();
         }
       }

       // --- Return the MVA outputs and fill into histograms
       for (unsigned int i=0; i<BDTNameList.size(); i++){
         bdtOutputs[i] = 0.0;
         if      (( isISR) && varSetList[i]=="ISR"){ bdtOutputs[i] = readerISR->EvaluateMVA( BDTNameList[i] );}
	 else if ((!isISR) && varSetList[i]!="ISR"){ bdtOutputs[i] = reader   ->EvaluateMVA( BDTNameList[i] );}
       }

       
       outTree->Fill();
    }

    for (unsigned int i=0;i<basicBDTVarNameList.size();i++){basicBDTVarFormulas[i]->Delete();}
    for (unsigned int i=0;i<  isrBDTVarNameList.size();i++){  isrBDTVarFormulas[i]->Delete();}

    outTree->Write();
    delete outTree;

    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();
  }

  inF->Flush();
  inF->Close();
  delete inF;

  
  return 0;
}


int ssEvtPostProc2 :: finalize ()
{
  delete reader;
  delete readerISR;

  return 0;
}


int ssEvtPostProc2 :: execute ()
{
  initialize();
  runLoop();
  finalize();
  std::cout << "done" << endl;
  return 0;
}



