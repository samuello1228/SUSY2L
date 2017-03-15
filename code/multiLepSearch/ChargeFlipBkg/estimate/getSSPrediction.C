/* getSSPrediction.
 * Gabriel Gallardo 8 Aug 2016
 * 
 * Get observed SS and predicted SS and save to histograms
 * Complete rewrite. Based on extractSSpredction and ChargeFlipTool
 * 
 */
#include <vector>
#include <iostream>
#include <fstream>

#include <TString.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TSystem.h>

#include "../common/evt2l.C"
#include "../ChargeFlipTool/ChargeFlipTool.cpp"

// ========= CONFIGURATION =========== //
TString defaultOut="outputDirYL-PtCorr";
TString defaultNTupleList="../common/inFileList-DataYL.txt";
// TString defaultMisIdfile="../QID-on/rates_wSys.root";
TString defaultMisIdfile="../QID-on/80.000000_100.000000_0.000000_0.000000_DATA.root";
// TString defaultMisIdfile="../common/chargeMisID_Zee_data_signal_wSys.root";
TString misIDhistname="80.0_100.0_0.0_0.0_DATA_misid";
//TString misIDhistname="80.0_100.0_20.0_20.0_DATA_misid";
// TString misIDhistname="hFlipProb_MCtruth";
// TString misIDhistname="hFlipProb_MCLH";
// TString misIDhistname="hFlipProb_data";
TString defaultdPtfile="../QID-on/ptcorr/dEhistos.root";
bool onlySignal=true;
bool applyPtCorrection=true;
bool passQID=true;
bool ZWindowOnly = true;
// ========= INFRASTRUCTURE =========== //
class Histos;
evt2l* evt2lTree=0;
TFile* fOut=0;
std::vector<Histos> histVector;
ChargeFlipTool* cfTool=0;

class Histos
{
  protected:
   Histos(){ SS = ExpSS = OS = Ratio = 0;}
  public:
   TString name;
   TH1D* SS;
   TH1D* ExpSS;
   TH1D* OS;
   TH1D* Ratio;
   TCanvas* c;
   TLegend* l;

   Histos(const TString name, const TString title, const TString xLabel, const TString yLabel, const int nBins, const double xMin, const double xMax){
      this->name = name;
      SS = new TH1D(name+"SS", "Observed "+title+" distribution of same sign electron pairs;"+xLabel+";"+yLabel, nBins, xMin, xMax);
      ExpSS = new TH1D(name+"ExpSS", "Estimated "+title+" distribution of same sign electron pairs;"+xLabel+";"+yLabel, nBins, xMin, xMax);
      OS = new TH1D(name+"OS", "Observed "+title+" of opposite sign electron pairs;"+xLabel+";"+yLabel, nBins, xMin, xMax);
      c = new TCanvas("c"+name, "c"+name, 1); 
      l = new TLegend(0.5, 0.9, 0.9, 0.7);

      SS->SetDirectory(fOut);
      ExpSS->SetDirectory(fOut);
      OS->SetDirectory(fOut);

      SS->Sumw2();
      ExpSS->Sumw2();
      OS->Sumw2();
      histVector.push_back(*this);
   }
   ~Histos(){}

   void Fill(const double value, const double w, const double chargeFlipWeight){
      if(chargeFlipWeight<0){
         SS->Fill(value, w);
      } else{
         ExpSS->Fill(value, w*chargeFlipWeight);
         OS->Fill(value, w);
      }      
   }

   void Draw(){
      // Hide info box of plot
      gStyle->SetOptStat(0);

      // Split pad
      c->Divide(1,2);
      c->cd(1)->SetPad(0, 0.25, 1, 1);
      c->cd(2)->SetPad(0, 0, 1, 0.25);
      c->cd(1)->SetLogy();

      // ====  UPPER PLOT ====== //
      // Set range of Y axis()
      // if(SS->GetMinimum()==0) SS->GetYaxis()->SetRangeUser(1, max(SS->GetMaximum(),ExpSS->GetMaximum())*5);
      SS->GetYaxis()->SetRangeUser(max(1.0,SS->GetMinimum()*0.5), max(SS->GetMaximum(),ExpSS->GetMaximum())*5);

      // Set colors
      SS->SetLineColor(kAzure+1); SS->SetFillColor(kAzure+1); SS->SetMarkerColor(kAzure+1); 
      ExpSS->SetLineColor(kRed); ExpSS->SetMarkerColor(kRed); 

      SS->Draw("bar");
      ExpSS->Draw("same e");

      l->AddEntry(SS, "Observed SS events", "f");
      l->AddEntry(ExpSS, "Predicted SS events from OS", "lp");
      l->Draw("same");

      // ====  LOWER PLOT ====== //
      c->cd(2)->SetGridy();

      Ratio = (TH1D*) ExpSS->Clone(name+"Ratio");
      Ratio->Divide(SS);
      Ratio->SetTitle("");
      TAxis* RatioX = Ratio->GetXaxis();
      TAxis* RatioY = Ratio->GetYaxis();

      RatioX->SetLabelOffset(999);
      RatioX->SetLabelSize(0);
      RatioX->SetTitle("");
      
      RatioY->SetTitle("Exp/Obs");
      RatioY->SetTitleSize(0.1);
      RatioY->SetTitleOffset(0.4);
      RatioY->SetLabelSize(0.08);
      RatioY->SetRangeUser(0.5,2);

      Ratio->Draw("p e");
      Ratio->SetDirectory(fOut);
      
      TLine line;
      line.SetLineStyle(2);
      line.DrawLine(RatioX->GetXmin(), 1, RatioX->GetXmax(), 1); 

      c->Print(name+".pdf", "Title:"+name);
   }

   void DrawOS(){
      TCanvas* can = new TCanvas("c"+name+"OS", "c"+name+"OS", 600, 600);
      OS->Draw();
      TString OSname = OS->GetName();
      can->Print(OSname+".pdf", "Title:"+OSname);
   }

   void Write(){
      SS->Write();
      ExpSS->Write();
      OS->Write();
      Ratio->Write();
   }

   void SetLegendXY(const double x1, const double y1, const double x2, const double y2){
      if(l) delete l;
      l = new TLegend(x1, y1, x2, y2);
   }
};

// ========= FUNCTIONS =========== //
void getSSPrediction(const TString outputDir, const TString ntupleList, const TString misIDfile, const TString dPtfile, const TString opt);
bool initialize(const TString outputDir, const TString ntupleList, const TString misIDfile, const TString dPtfile, const TString opt);
TChain* loadData(const TString);

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

void getSSPrediction(const TString outputDir=defaultOut, 
                     const TString ntupleList=defaultNTupleList, 
                     const TString misIDfile=defaultMisIdfile,
                     const TString dPtfile=defaultdPtfile, 
                     const TString opt="")
{
   
   if(!initialize(outputDir, ntupleList, misIDfile, dPtfile, opt)){
      std::cout << "Exited from initialize()" << std::endl;
      return;
   }

   // Initialize histograms here ----------------------------------------
   Histos hMass("hMass", "invariant mass", "m_{ll}", "Events/2 GeV", ZWindowOnly?20:100, ZWindowOnly?70:60, ZWindowOnly?110:260);
   Histos hLeadingPt("hLeadingPt", "leading p_{T}", "p_{T}", "Events/4 GeV", 50, 20, 200);
   Histos hSubleadingPt("hSubleadingPt", "subleading p_{T}", "p_{T}", "Events/4 GeV", 50, 20, 200);
   Histos hLeadingEta("hLeadingEta", "leading |#eta|", "Leading |#eta|", "", 40, -2.47, 2.47);
   Histos hSubleadingEta("hSubleadingEta", "subleading |#eta|", "Subleading |#eta|", "", 40, -2.47, 2.47);
   Histos hLeadingPhi("hLeadingPhi", "leading #phi", "Leading |#phi|", "", 40, -3.15, 3.15);
   Histos hSubleadingPhi("hSubleadingPhi", "subleading #phi", "Subleading |#phi|", "", 40, -3.15, 3.15);

  hMass.SetLegendXY(0.7, 0.2, 0.9, 0.35);

   long long nEntries = evt2lTree->fChain->GetEntries();
   for(long long i=0; i<nEntries; i++){
      loadbar(i+1, nEntries);
      evt2lTree->GetEntry(i);

      bool baseCut = (evt2lTree->sig_trigCode>0 && evt2lTree->leps_ == 2
                     && int(abs(evt2lTree->leps_ID[0]/1000))==11 && int(abs(evt2lTree->leps_ID[1]/1000))==11); 
      if(!baseCut) continue;

      // Select signal events
      if(onlySignal && !(((evt2lTree->leps_lFlag[0] & 2)/2) && ((evt2lTree->leps_lFlag[1] & 2)/2))) continue;

      // Select Zee events
      if(ZWindowOnly && fabs(evt2lTree->l12_m - 91)>10) continue;

      // Select events which pass electron ChargeFlipTagger
      if(passQID && !(evt2lTree->leps_ElChargeID[0] && evt2lTree->leps_ElChargeID[1])) continue;

      bool SSevent = ((evt2lTree->leps_ID[0]>0) == (evt2lTree->leps_ID[1]>0));

      double pt1 = evt2lTree->leps_pt[0];
      double pt2 = evt2lTree->leps_pt[1];

      // ptCorrection for OS events
      // Correct pt of electron with highest probability to flip.
      auto flip1 = cfTool->getChargeFlipRate(evt2lTree->leps_eta[0], evt2lTree->leps_pt[0]);
      auto flip2 = cfTool->getChargeFlipRate(evt2lTree->leps_eta[1], evt2lTree->leps_pt[1]);
      if(applyPtCorrection && !SSevent){
         if(flip1.first > flip2.first){
            pt1 = (cfTool->getCorrectedPt(evt2lTree->leps_eta[0], evt2lTree->leps_pt[0])).first;
         } else {
            pt2 = (cfTool->getCorrectedPt(evt2lTree->leps_eta[1], evt2lTree->leps_pt[1])).first;
         }
      }
      TLorentzVector p1, p2;
      p1.SetPtEtaPhiM(pt1, evt2lTree->leps_eta[0], evt2lTree->leps_phi[0], 0.000511);
      p2.SetPtEtaPhiM(pt2, evt2lTree->leps_eta[1], evt2lTree->leps_phi[1], 0.000511);
      double mll = (p1+p2).M();

      // Fill histograms here ----------------------------------------
      double chargeFlipWeight=0;
      if(!SSevent){
         double pSS = flip1.first + flip2.first - 2*flip1.first*flip2.first;
         chargeFlipWeight = pSS/(1-pSS);
      } else chargeFlipWeight = -1;

      double w = 1;//evt2lTree->evt_weight * evt2lTree->evt_pwt;

      hMass.Fill(mll, w, chargeFlipWeight);
      hLeadingPt.Fill(pt1, w, chargeFlipWeight);
      hSubleadingPt.Fill(pt2, w, chargeFlipWeight);
      hLeadingEta.Fill(evt2lTree->leps_eta[0], w, chargeFlipWeight);
      hSubleadingEta.Fill(evt2lTree->leps_eta[1], w, chargeFlipWeight);
      hLeadingPhi.Fill(evt2lTree->leps_phi[0], w, chargeFlipWeight);
      hSubleadingPhi.Fill(evt2lTree->leps_phi[1], w, chargeFlipWeight);
   }

   for(auto h : histVector){
      h.Draw();
      h.DrawOS();
      h.Write();
   }
}

bool initialize(const TString outputDir, const TString ntupleList, const TString misIDfile, const TString dPtfile, const TString opt)
{
   ///////////////////////////////////
   // Load data
   ///////////////////////////////////
   TChain *ch = loadData(ntupleList);
   if(!ch){
      std::cout << "Could not load data from " << ntupleList << std::endl;
      return false;
   }
   evt2lTree = new evt2l(ch);
   std::cout << evt2lTree->fChain->GetEntries() << " events in input\n";

   ///////////////////////////////////
   // Load ChargeFlipTool 
   ///////////////////////////////////
   // cfTool = new ChargeFlipTool(misIDfile, dPtfile);
   cfTool = new ChargeFlipTool();
   if(!cfTool->loadMisIDfile(misIDfile, misIDhistname))
   {
      cout << "Failed to load misID file\n";
      return false;
   }
   if(applyPtCorrection && !cfTool->loadDPTfile(dPtfile))
   {
      cout << "Failed to load dPt file\n";
      return false;
   }

   ///////////////////////////////////
   // Setup output file
   ///////////////////////////////////
   if(!gSystem->cd(outputDir)){
      gSystem->mkdir(outputDir);
      if(!gSystem->cd(outputDir)){
         std::cout << "Could not create output directory " << outputDir << std::endl;
         return false;
      }
   }

   fOut = TFile::Open("SSpredictions.root", "recreate");
   if(!fOut){
      std::cout << "Could not create output file " << outputDir << "/SSpredictions.root" << std::endl;
      return false;
   }

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
         } else {
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
         for( std::string line; getline( inFiles, line ); ){tc->Add(line.c_str());}
         inFiles.close();
      }
   }
   return tc;
}
