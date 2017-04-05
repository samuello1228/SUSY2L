/* SSfromMC.
 * Gabriel Gallardo 27 Mar 2017
 * Based on getSSPrediction.C
 *
 * Get SS Prediction in data from reweighting MC
 * 
 * Execution
 * - Be sure to have checked out and compiled SUSYTools via rcSetup
 * - `root -l -b- q SSfromMC.C+`
 */
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include <TString.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TLine.h>

#include "ChargeFlipBkg/common/evt2l.C"
#include "ChargeFlipBkg/ChargeFlipTool/ChargeFlipTool.cpp"

#include "SUSYTools/SUSYCrossSection.h"


// ========= CONFIGURATION =========== //
TString defaultDir = "/afs/cern.ch/user/g/ggallard/private/SUSY2L/code/multiLepSearch/ChargeFlipBkg/";
TString defaultMClist = defaultDir + "common/inFileList-ZeePowheg.txt";
// TString defaultMClist = defaultDir + "common/ZP.txt";
TString defaultDataList = defaultDir + "common/inFileList-data.txt";
// TString defaultDataList = defaultDir + "common/shortD.txt";
TString defaultOut= defaultDir + "QiD-on/FromMC";
TString dPtfile=defaultDir + "QiD-on/Powheg/ptcorr/dEhistos.root";
TString qSFfile=defaultDir + "ChargeCorrectionSF.Medium_FixedCutTightIso_CFTMedium.root";
TString qSFnameOS="SFCentral_RunNumber296939_311481_OS";
TString qSFnameSS="SFCentral_RunNumber296939_311481_SS";

bool onlySignal=true;
bool passQID=true;
bool ZWindowOnly=true;
bool applyQF=true;

// ========= INFRASTRUCTURE =========== //
class Histos;
evt2l* dataTree=0;
evt2l* mcTree=0;
TFile* fOut=0;
std::vector<Histos> histVector;
ChargeFlipTool* cfTool=0;
double mcEvtW=0;

enum inType{DATA, MC};
enum sign{OS=0, SS};

class Histos
{
  protected:
   Histos(){ SS = ExpSS = OS = RatioOS = RatioSS = 0;}
  public:
   TString name;
   TH1D* SS;
   TH1D* ExpSS;
   TH1D* OS;
   TH1D* ExpOS;
   TH1D* RatioSS;
   TH1D* RatioOS;
   // TCanvas* c;
   TLegend* l;
   double* legendXY;

   Histos(const TString name, const TString title, const TString xLabel, const TString yLabel, const int nBins, const double xMin, const double xMax){
      this->name = name;
      SS = new TH1D(name+"SS", "Observed "+title+" distribution of same sign electron pairs;"+xLabel+";"+yLabel, nBins, xMin, xMax);
      ExpSS = new TH1D(name+"ExpSS", "Estimated "+title+" distribution of same sign electron pairs;"+xLabel+";"+yLabel, nBins, xMin, xMax);
      OS = new TH1D(name+"OS", "Observed "+title+" of opposite sign electron pairs;"+xLabel+";"+yLabel, nBins, xMin, xMax);
      ExpOS = new TH1D(name+"ExpOS", "Estimated "+title+" of opposite sign electron pairs;"+xLabel+";"+yLabel, nBins, xMin, xMax);

      SS->SetDirectory(fOut);
      ExpSS->SetDirectory(fOut);
      OS->SetDirectory(fOut);
      ExpOS->SetDirectory(fOut);

      legendXY = new double[4];
      legendXY[0] = 0.5; legendXY[1] = 0.9; legendXY[2] = 0.9; legendXY[3] = 0.7;

      SS->Sumw2();
      ExpSS->Sumw2();
      OS->Sumw2();
      ExpOS->Sumw2();
      histVector.push_back(*this);
   }
   ~Histos(){}

   void Fill(const double value, const double w, const inType t, const sign s){
      if      (t==MC    && s==sign::SS) ExpSS->Fill(value,w);
      else if (t==MC    && s==sign::OS) ExpOS->Fill(value,w);
      else if (t==DATA  && s==sign::SS) SS   ->Fill(value,w);
      else if (t==DATA  && s==sign::OS) OS   ->Fill(value,w); 
      return;
   }

   void Draw(const sign s){
      // Hide info box of plot
      gStyle->SetOptStat(0);
      cout << "Drawing " << name << ((s==sign::SS)? " SS" : " OS") << endl;

      // Split pad
      TCanvas* c = new TCanvas("c"+name, "c"+name, 1); 
      c->Divide(1,2);
      c->cd(1)->SetPad(0, 0.25, 1, 1);
      c->cd(2)->SetPad(0, 0, 1, 0.25);
      c->cd(1)->SetLogy();

      TH1* obs = (s==sign::SS) ? SS : OS;
      TH1* exp = (s==sign::SS) ? ExpSS : ExpOS;

      // ====  UPPER PLOT ====== //
      // Set range of Y axis()
      exp->GetYaxis()->SetRangeUser(max(1.0,obs->GetMinimum()*0.5), max(obs->GetMaximum(),exp->GetMaximum())*5);

      // Set colors
      exp->SetLineColor(kAzure+1); exp->SetFillColor(kAzure+1); exp->SetMarkerColor(kAzure+1); 
      obs->SetLineColor(kBlack); obs->SetMarkerColor(kBlack); obs->SetMarkerStyle(kFullCircle);

      exp->Draw("bar");
      obs->Draw("same e");

      l = new TLegend(legendXY[0], legendXY[1], legendXY[2], legendXY[3]);
      l->AddEntry(obs, "Data", "lp");
      l->AddEntry(exp, "MC", "f");
      l->Draw("same");
      cout << "Upper plot done" <<endl; 
      // ====  LOWER PLOT ====== //
      c->cd(2)->SetGridy();

      TH1D* r = (TH1D*) obs->Clone((TString)(name+"Ratio")+((s==sign::SS)?"SS":"OS"));
      cout << "1" << endl;
      r->Divide(exp);
      r->SetTitle(";;Data/MC");
      TAxis* RatioX = r->GetXaxis();
      TAxis* RatioY = r->GetYaxis();

      RatioX->SetLabelOffset(999);
      RatioX->SetLabelSize(0);
      
      RatioY->SetTitleSize(0.1);
      RatioY->SetTitleOffset(0.4);
      RatioY->SetLabelSize(0.08);
      RatioY->SetRangeUser(0.5,2);

      r->Draw("p e");
      r->SetDirectory(fOut);

      TLine line;
      line.SetLineStyle(2);
      line.DrawLine(RatioX->GetXmin(), 1, RatioX->GetXmax(), 1); 
      cout << "2" << endl;

      c->cd();
      TString printName=(TString) obs->GetName()+".pdf";
      TString printTitle=(TString) "Title:"+name;
      c->Print(printName, printTitle);
      // c->Print("Hello.pdf", "Title:Hello");
      cout << "3" << endl;

      (s==sign::SS) ? RatioSS = r : RatioOS = r;
      cout << "Lower plot done" <<endl<<endl;
      delete l; delete c;
   }

   void Write(){
      if (SS) SS->Write();
      if (ExpSS) ExpSS->Write();
      if (OS) OS->Write();
      if (ExpOS) ExpOS->Write();
      if (RatioOS) RatioOS->Write();
      if (RatioSS) RatioSS->Write();
   }

   void SetLegendXY(const double x1, const double y1, const double x2, const double y2){
      legendXY[0] = x1;
      legendXY[1] = y1;
      legendXY[2] = x2;
      legendXY[3] = y2;
   }
};

// ========= FUNCTIONS =========== //
void SSfromMC(const TString outputDir, const TString mcList, const TString dataList, const TString opt);
bool initialize(const TString outputDir, const TString mcList, const TString dataList, const TString opt);
TChain* loadData(const TString);
void processMC(const TString mcList);
bool passCut(const evt2l* tree, const inType t);
int getOrigElecI(const evt2l* tree, int i);

inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
   if ( (x != n) && (x % (n/100+1) != 0) ) return;

   float ratio    =  x/(float)n;
   unsigned int c =  ratio * w;

   std::cout << setw(3) << (int)(ratio*100) << "% [";
   for (unsigned int x=0; x<c; x++) std::cout << "=";
   for (unsigned int x=c; x<w; x++) std::cout << " ";
   std::cout << "]\r" << flush;
}

void SSfromMC(const TString outputDir=defaultOut, 
                     const TString mcList=defaultMClist, 
                     const TString dataList=defaultDataList,
                     const TString opt="")
{
   
   if(!initialize(outputDir, mcList, dataList, opt)){
      std::cout << "Exited from initialize()" << std::endl;
      return;
   }

   ///////////////////////////////////
   // Setup output file
   ///////////////////////////////////
   if(!gSystem->cd(outputDir)){
      gSystem->mkdir(outputDir);
      if(!gSystem->cd(outputDir)){
         std::cout << "Could not create output directory " << outputDir << std::endl;
         return;
      }
   }

   fOut = TFile::Open(outputDir + "/SSpredictions.root", "recreate");
   if(!fOut){
      std::cout << "Could not create output file " << outputDir << "/SSpredictions.root" << std::endl;
      return;
   }

   // Initialize histograms here ----------------------------------------
   // gDirectory->cd(fOut);
   Histos hMass("hMass", "invariant mass", "m_{ll}", "Events/2 GeV", 100, 60, 260);
   Histos hLeadingPt("hLeadingPt", "leading p_{T}", "p_{T}", "Events/4 GeV", 50, 20, 200);
   Histos hSubleadingPt("hSubleadingPt", "subleading p_{T}", "p_{T}", "Events/4 GeV", 50, 20, 200);
   Histos hLeadingEta("hLeadingEta", "leading |#eta|", "Leading |#eta|", "", 40, -2.47, 2.47);
   Histos hSubleadingEta("hSubleadingEta", "subleading |#eta|", "Subleading |#eta|", "", 40, -2.47, 2.47);
   Histos hLeadingPhi("hLeadingPhi", "leading #phi", "Leading |#phi|", "", 40, -3.15, 3.15);
   Histos hSubleadingPhi("hSubleadingPhi", "subleading #phi", "Subleading |#phi|", "", 40, -3.15, 3.15);

   hMass.SetLegendXY(0.7, 0.2, 0.9, 0.35);

   // Load xsecDB
   gROOT->Macro("$ROOTCOREDIR/scripts/load_packages.C");
   SUSY::CrossSectionDB *xsecDB = new SUSY::CrossSectionDB(gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc15_13TeV/"));
   if (!xsecDB){
      cout << "CrossSectionDB could not be loaded" << endl;
      return;
   }

   // Load charge scale factors
   TFile *fQSF = TFile::Open(qSFfile);
   TH2* hQSF_OS = (TH2*) fQSF->Get(qSFnameOS);
   TH2* hQSF_SS = (TH2*) fQSF->Get(qSFnameSS);

   //////////////////////////
   // Initialize input trees
   //////////////////////////
   // MC
   std::ifstream mcIn(mcList);
   if(!mcIn.is_open()) { cout << "Could not open " << mcList << endl; return; }
   else cout << "Reading MC from " << mcList << endl;
   // Data
   TChain *dataCh = loadData(dataList);
   if(!dataCh){ std::cout << "Could not load data from " << dataList << std::endl; return; }
   dataTree = new evt2l(dataCh);
   if(!dataTree) { std::cout << "Could not initialize evt2l object from " << dataList << std::endl; return; }

   //////////////////////////
   // MC Loop
   //////////////////////////
   // Define mcEvtW ---------------
   vector<string> mcFiles; double sumW = 0;
   for (std::string line; getline( mcIn, line ); )
   {
      TFile f(line.c_str());
      if(!f.IsOpen()){ cout << "Could not open " << line << endl; continue; }

      TH1* h = (TH1*) f.Get("hCutFlow");
      if (!h){ cout << "Could not get hCutFlow histogram for " << line << endl; return; }

      sumW += h->GetBinContent(2);
      f.Close(0);
      mcFiles.push_back(line);
   }
   double xSecxEff = xsecDB->xsectTimesEff(361106);
   double mcEvtW = xSecxEff * 33257.2 / sumW;

   // Run main MC loop
   for(auto line : mcFiles){

      // Open file ----------------
      TFile f(line.c_str());
      if(!f.IsOpen()){ cout << "Could not open " << line << endl; continue; }

      TTree* t = (TTree*) f.Get("evt2l");
      if(!t){ cout << "No tree found in "<< line << endl; continue; }

      evt2l* tree = new evt2l(t);
      if(!tree){ cout << "Could not load evt2l of " << line << endl; continue; }

      // Loop in one file --------------
      int nEntries = tree->fChain->GetEntries();
      cout << line << ":\n" ;
      for(long long i=0; i<nEntries; i++){
         loadbar(i+1, nEntries);
         tree->GetEntry(i);

         if(!passCut(tree, MC)) continue;

         sign s = ((tree->leps_ID[0]>0) == (tree->leps_ID[1]>0))? SS : OS;

         // Reject if no original electron from Z
         int tp1 = getOrigElecI(tree, 0); int tp2 = getOrigElecI(tree, 1);
         if(tp1<0 || tp2<0) continue;

         // Get electron charge correction SF // positron has +ve ID, -ve pdgID
         double w1 = (tree->leps_ID[0]>0 == (tree->truths_pdgId[tp1] < 0))?  
         1 :hQSF_SS->GetBinContent(hQSF_SS->FindBin(tree->leps_pt[0]>=150? 149: tree->leps_pt[0], fabs(tree->leps_eta[0])))
            / hQSF_OS->GetBinContent(hQSF_OS->FindBin(tree->leps_pt[0]>=150? 149: tree->leps_pt[0], fabs(tree->leps_eta[0])));
            // hQSF_OS->GetBinContent(hQSF_OS->FindBin(tree->leps_pt[0]>=150? 149: tree->leps_pt[0], fabs(tree->leps_eta[0])))
            // : hQSF_SS->GetBinContent(hQSF_SS->FindBin(tree->leps_pt[0]>=150? 149: tree->leps_pt[0], fabs(tree->leps_eta[0]))) ;

         double w2 = (tree->leps_ID[1]>0 == (tree->truths_pdgId[tp2] < 0))?  
         1 :hQSF_SS->GetBinContent(hQSF_SS->FindBin(tree->leps_pt[1]>=150? 149: tree->leps_pt[1], fabs(tree->leps_eta[1])))
            / hQSF_OS->GetBinContent(hQSF_OS->FindBin(tree->leps_pt[1]>=150? 149: tree->leps_pt[1], fabs(tree->leps_eta[1])));
            // hQSF_OS->GetBinContent(hQSF_OS->FindBin(tree->leps_pt[1]>=150? 149: tree->leps_pt[1], fabs(tree->leps_eta[1])))
            // : hQSF_SS->GetBinContent(hQSF_SS->FindBin(tree->leps_pt[1]>=150? 149: tree->leps_pt[1], fabs(tree->leps_eta[1]))) ;


         // double w1, w2; w1 = w2 = 1;
         double w = tree->evt_weight * tree->evt_pwt *tree->evt_ElSF * tree->evt_MuSF * mcEvtW *w1 *w2;

         hMass.Fill(tree->l12_m, w*w1*w2, MC, s);
         hLeadingPt.Fill(tree->leps_pt[0], w*w1, MC, s);
         hSubleadingPt.Fill(tree->leps_pt[1], w*w2, MC, s);
         hLeadingEta.Fill(tree->leps_eta[0], w*w1, MC, s);
         hSubleadingEta.Fill(tree->leps_eta[1], w*w2, MC, s);
         hLeadingPhi.Fill(tree->leps_phi[0], w*w1, MC, s);
         hSubleadingPhi.Fill(tree->leps_phi[1], w*w2, MC, s);
      }

      // // Cleanup
      // delete h; delete tree; delete t;
      // f.Close();
      cout << endl;
   }
   
   //////////////////////////
   // Data loop
   //////////////////////////
   cout << endl << "Data tree:" << endl;
   long long nEntries = dataTree->fChain->GetEntries();
   for(long long i=0; i<nEntries; i++){
      loadbar(i+1, nEntries);
      dataTree->GetEntry(i);

      if(!passCut(dataTree, DATA)) continue;

      sign s = ((dataTree->leps_ID[0]>0) == (dataTree->leps_ID[1]>0))? SS : OS;
      
      double w = 1; // dataTree->evt_weight * dataTree->evt_pwt *dataTree->evt_ElSF * dataTree->evt_MuSF;

      hMass.Fill(dataTree->l12_m, w, DATA, s);
      hLeadingPt.Fill(dataTree->leps_pt[0], w, DATA, s);
      hSubleadingPt.Fill(dataTree->leps_pt[1], w, DATA, s);
      hLeadingEta.Fill(dataTree->leps_eta[0], w, DATA, s);
      hSubleadingEta.Fill(dataTree->leps_eta[1], w, DATA, s);
      hLeadingPhi.Fill(dataTree->leps_phi[0], w, DATA, s);
      hSubleadingPhi.Fill(dataTree->leps_phi[1], w, DATA, s);
   }
   cout << endl << endl;

   // Cleanup
   delete dataTree; delete dataCh;
   delete xsecDB;
   //////////////////////////
   // Draw and write histograms
   //////////////////////////
   for(auto h : histVector){
      h.Draw(sign::OS);
      // h.Draw(sign::SS);
      // h.Write();
   }
   fOut->Write();
   return;
}

bool passCut(const evt2l* tree, const inType t)
{
   bool baseCut = (tree->sig_trigCode>0 && tree->leps_ == 2
                  && int(abs(tree->leps_ID[0]/1000))==11 && int(abs(tree->leps_ID[1]/1000))==11); 
   if(!baseCut) return false;

   // Select signal events
   if(onlySignal && !(((tree->leps_lFlag[0] & 2)/2) && ((tree->leps_lFlag[1] & 2)/2))) return false;

   // For SS Data events, select ZWindow only
   if(t==DATA && ZWindowOnly 
      && ((tree->leps_ID[0]>0) == (tree->leps_ID[1]>0)) && fabs(tree->l12_m - 91)>10) 
      return false;

   // Select events which pass electron ChargeFlipTagger
   if(passQID && !(tree->leps_ElChargeID[0] && tree->leps_ElChargeID[1])) return false;

   // In MC, plot only events with truth-matched electrons 
   // if (t==MC &&
   //    (tree->leps_truthI[0] < 0 || tree->leps_truthI[0] < 0 
   //      || abs(tree->truths_pdgId[tree->leps_truthI[0]])!=11 
   //      || abs(tree->truths_pdgId[tree->leps_truthI[0]])!=11) )
   //    return false;

   return true;
}

bool initialize(const TString outputDir, const TString mcList, const TString dataList, const TString opt)
{
   ///////////////////////////////////
   // Read options
   ///////////////////////////////////
   if(opt!="") // Finished if no options to read
   {
      TString key, value;
      int nChar = opt.Length();
      for(int i=0; i<nChar; i++){
         key = value = "";
         while(opt[i]!='=' && opt[i]!='\0'){
            key += opt[i];
            i++;
         } i++;
         while(opt[i]!=',' && opt[i]!='\0'){
            value += opt[i];
            i++;
         }
   
         if(key=="el"){
            if(value=='s') onlySignal = true;
            else if (value=='l') onlySignal = false;
            else {
               std::cout << "Option unrecognized. " << std::endl;
               return false;
            }
         } 
         
         else {
            std::cout << "Option unrecognized. " << std::endl;
            return false;
         }
      }
   }

   if(onlySignal){
      std::cout << "Selecting Signal events" << std::endl;
   } else {
      std::cout << "Selecting LooseBaseline events" << std::endl;
   }

   // Everything okay!
   std::cout << "initialize() complete" << std::endl;
   return true;
}

TChain* loadData(const TString fileList){
   //return a TChain linked to the data files
   TChain* tc = new TChain("evt2l");

   if (fileList){
      std::ifstream inFiles(fileList);
      if (inFiles.is_open()){
         cout << "Reading from " << fileList << endl;
         for( std::string line; getline( inFiles, line ); ){tc->Add(line.c_str());}
         inFiles.close();
      }
      else {
         cout << "Could not open " << fileList << endl;
         return 0;
      }
   }
   return tc;
}

int getOrigElecI(const evt2l* tree, int i){
   stringstream tempOut;
   // Find matching truth particle
   int tp = tree->leps_truthI[i];
   if (tp < 0) return -1; // tp = -1 if no truth particle;
   int matchedPDGID = tree->truths_pdgId[tp];
   tempOut << "PDGIDs: " << matchedPDGID << '\t';
   if (abs(matchedPDGID) != 11) {
      tempOut << " Not an electron" << endl;
      return -1; // Reject if not matched to an electron
   }

  // Find original electron from Z
   int mother = tree->truths_motherI[tp];

   // while mother exists and is not photon or electron
   while(mother>=0 && tp >= 0 && (abs(tree->truths_pdgId[mother]) == 11 || tree->truths_pdgId[mother] == 22)  && abs(tree->truths_pdgId[mother]) < 100){ 
      tp = mother;
      mother = tree->truths_motherI[tp];
      tempOut << tree->truths_pdgId[tp] << '\t';
   }
   tempOut << tree->truths_pdgId[mother] << '\t';  

   // Reject if the particle found is invalid
   if (tp <0 || mother <0){ 
      tempOut << "Particle invalid" << endl;
      return -1; 
   }

   // Reject if particle found is not electron
   if (abs(tree->truths_pdgId[tp])==11){
      return tp;
   } else if (tree->truths_pdgId[tp]==22){
      return -2;
   }
   else{ 
      return -1;
   }
}
