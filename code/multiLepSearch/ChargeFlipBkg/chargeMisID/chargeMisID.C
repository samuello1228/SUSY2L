/* chargeMisID.C
 * calculate electron charge misID from UM/HKU SUSY ntuples for Run 2
 * 
 * Written by Gabriel Gallardo February 2016, with contribution by YatLong Chan Nov 2015
 *
 * Execution instructions: 
 * - Specify isMC before running! Running on isMC=true with data will cause segmentation fault.
 * - Set config variables appropriately
 * $ root -l -b
 * root[0] . x chargeMisID.C("outputDir/", "inputFileList.txt", "opt") ## if outputDir does not exist, it will be created
 *
 * Execution recommendations:
 * - Use the "opt" argument when scripting. See documentation of options in README.md
 * 
 */

// #include "evt2l.C"
// #include "obj_def.h"
#include "../common/common.C"
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TText.h>
#include <TF1.h>
#include <TPie.h>
#include <TMath.h>
#include <math.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <TColor.h>
#include <TSystem.h>
#include <TString.h>

// ============= CONFIG VARIABLES ================== //

Bool_t isMC = true;          // true for MC, false for data
Bool_t onlySignal = true;       // true for signal, false for looseBaseline
const Bool_t debugMC = true;
const Bool_t MCtruthPt = isMC && false; // true for finding MC LH based on pt of matched truth electron.
Bool_t newDerivation = true;
Bool_t applyPRW = true;
Bool_t findPtCorr = isMC && true;

// ---- Configure output options --- ///
const Bool_t  print = true;         // Set to true to output to pdf
const Bool_t  toFile = true;        // Set to true to output to output.root
const Bool_t  drawEta = false;      // Draw integrated rates over eta
const Bool_t  drawPt  = false;      // Draw integrated rates over pt
const Bool_t  drawEtaPt = true;     // Draw rates over eta in diff pt bins
const Bool_t  drawPtEta = true;     // Draw rates over pt in diff eta bins
const Bool_t  logPt = true;

// const Bool_t onlyZMCevents = false;

// ---- Configure Pt Bins --- //
// const Double_t tPtEdges[]  = {20, 30, 40, 50, 60, 80, 120}; // 6 bins
// const Double_t tPtEdges1[] = {20, 30, 35, 40, 45, 50, 55, 60, 80, 120}; // 9 bins, for finer binning in 30 -60
// const Double_t tPtEdges2[]  = {20, 30, 40, 50, 60, 80, 120, 300}; // 7 bins, for extended high pt range
// const Double_t tPtEdges3[] = {20.0, 60.0, 90.0, 130.0, 1000.0}; // 4 bins, EGamma chargeMisID tool
const Double_t tPtEdges4[]  = {20, 30, 40, 50, 60, 80, 120, 1000.0}; // 7 bins, for extended high pt range
const Double_t   *PtEdges  = tPtEdges4;
const Int_t      nPtBins   = 7;

// ---- Configure Eta Bins --- //
const bool ABS_ETA = true;
// const Double_t tEtaEdges1[]  = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.37}; // 6 bins, for testing
// const Double_t tEtaEdges2[]  = {1.5, 1.6, 1.7, 1.8}; //3 bins, for testing
// const Double_t tEtaEdges3[]  = {-2.5, -2, -1.8, -1.52, -1.37, -1, -0.5, 0, 0.5, 1, 1.37, 1.52, 1.8, 2, 2.5}; // 14 bins, for testing eta (vs. |eta|)
const Double_t tEtaEdges4[] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.37, 1.52, 1.8, 2, 2.2, 2.47}; // 12 bins, for fine binning
const Double_t tEtaEdges[]   = {0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.47}; // 7 bins
const Double_t            *EtaEdges = tEtaEdges4; 
const Int_t                nEtaBins = 12;

// ------ Configure event cut ----- //
double Z_MASS_UPPER_LIMIT = 100;
double Z_MASS_LOWER_LIMIT = 75;
// const double dPhiCut = 2;

// ------ Process only some entries ---- //
long long maxEntries = 0;

// ----- Configure sideband subtraction for likelihood method --- //
// Width of sidebands in GeV
// Set to zero for now
double SIDEBAND_L = 0; 
double SIDEBAND_R = 0;

// ----- Configure tag conditions for TP --- //
bool doTP = false;

// ================== INFRASTRUCTURE =======================// 

// "enum" for accessing elements in histogram arrays
const int ERROR = 0;
const int N_PROBES = 1;
const int N_SSPROBES = 2;
const int N_OSPROBES = 3;

// const int DEFAULT = 0;
// const int LOOSE = 1;
// const int MEDIUM = 2;
// const int TIGHT = 3;
// const int MC = 10;

// const int AUTO = 0;
// const int MAN = 1;

// --- For error(eta, pt) (TP) --- 
TH2D           *hEP_TP[4];

// --- For error(eta, pt) (LH) --- 
Double_t *aLHnOS;
Double_t *aLHnSS;
TH2D* hFlipProb = NULL;
ROOT::Math::Minimizer* mini;

// --- For sideband subtraction of LH method --- 
Double_t *aSideLnOS;
Double_t *aSideLnSS;
Double_t *aSideRnOS;
Double_t *aSideRnSS;

// For error(eta, pt) (MC truth)
TH2      *hEP_MC[4];

// For dPt
TFile*   ptOut       = NULL;
TH3*     hDPT        = NULL;
TH3*     hDPTflipped = NULL;
TH3*     hDPTok      = NULL;

// output 
TFile             *outputFile = 0;
std::ofstream     decayTree;

// ==================== FUNCTION PROTOTYPES ==================== //
void chargeMisID(const char *out = "outputDir", TString fileList = "inFileList.txt", const char *opt = "");
inline void loadbar(unsigned int x, unsigned int n, unsigned int w);

TChain* loadData(TString fileList);
bool initialize(const string opt);

void processEvents(evt2l& evt2lTree);
void processMC(evt2l& evt2lTree);
void fillCountArr(double* arr, double pt0, double eta0, double pt1, double eta1, double weight);
int getTruthElecI(evt2l& evt2lTree, int i);

void EMuTP(evt2l& evt2lTree);
void FillHist(double eta, double pt, TH2* h[4], bool iseeSS, double w=1);
void DivideHist(TH1* h[]);

void subtractBkg(double* central, const double* left, const double* right);
double likelihood(const double *x);
double getCountArr(const double* arr, int pt0Bin, int eta0Bin, int pt1Bin, int eta1Bin);

void draw();
void write_matrix();

// ==================== FUNCTION DEFINITIONS ==================== //
// Progress bar 
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

// ------- MAIN --------- //
void chargeMisID(const char * out, TString fileList, const char* opt)
{

   TChain* tc = loadData(fileList);
   if(!tc){
      std::cout << "NTuples could not be loaded" << endl;
      return;
   }
   evt2l evt2lTree(tc);
   std::cout << evt2lTree.fChain->GetEntries() << " events in input \n";

   if(!gSystem->cd(out)){
      gSystem->mkdir(out);
      gSystem->cd(out);
   }
   if(!initialize(opt)){
      std::cout << "Error in initialize()." << std::endl;
      return;
   }

   // ------------- LIKELIHOOD METHOD ON ZEE-------------------- //
   if(!doTP){
      processEvents(evt2lTree); // Likelihood loop over all entries in TChain/TTree
      std::cout << std::endl << "processEvents() done. \n";
      std::cout << std::endl << "Calculating misID... \n";

      if(SIDEBAND_L!=0 || SIDEBAND_R !=0){
         subtractBkg(aLHnSS, aSideLnSS, aSideRnSS);
         subtractBkg(aLHnOS, aSideLnOS, aSideRnOS);
      }

      //ROOT::Math::Minimizer* myMinimizer =
      mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

      mini->SetMaxFunctionCalls(1000000); 

      int nOptVar = nPtBins*nEtaBins;
      ROOT::Math::Functor f(&likelihood, nOptVar);
      mini->SetFunction(f);
      mini->SetErrorDef(0.5); // Define scale for errors. 0.5 for likelihood minimization

      //set intital value and step size 
      double* step     = new double[nOptVar];
      double* flipProb = new double[nOptVar]; //the var we want to find out via minuit

      std::fill_n(     step, nOptVar, 1e-5); // increment
      std::fill_n( flipProb, nOptVar, 1e-4); // initial flip prob

      // Set the free variables to be minimized!   
      for (int  ptBin = 0;  ptBin<  nPtBins;  ptBin++){
      for (int etaBin = 0; etaBin< nEtaBins; etaBin++){
         int gBin = ptBin*nEtaBins + etaBin;
         char varName[50];
         sprintf(varName, "pt%02d_eta%02d", ptBin, etaBin);
         mini->SetLowerLimitedVariable(gBin,varName, flipProb[gBin], step[gBin], 0.0);
      }}

      // std::cout << "pt0Bin\teta0Bin\tpt1Bin\teta1Bin\tnOS\tnSS\n";
      // for (int  pt0Bin = 0;  pt0Bin<  nPtBins;  pt0Bin++){
      // for (int  pt1Bin = 0;  pt1Bin<  nPtBins;  pt1Bin++){
      // for (int eta0Bin = 0; eta0Bin< nEtaBins; eta0Bin++){
      // for (int eta1Bin = 0; eta1Bin< nEtaBins; eta1Bin++){
      //    double nOS  = getCountArr(aLHnOS, pt0Bin, eta0Bin, pt1Bin, eta1Bin);
      //    double nSS  = getCountArr(aLHnSS, pt0Bin, eta0Bin, pt1Bin, eta1Bin);
      //    if(nOS != 0 || nSS != 0)
      //    std::cout << pt0Bin << "\t" << eta0Bin <<"\t" << pt1Bin << "\t" <<eta1Bin << "\t" <<nOS << "\t" <<nSS << std::endl;
      // }}}}

      // -ve matrix elements to zero
      int nTotBins = nOptVar*nOptVar;
      for(int i=0; i<nTotBins; i++){
         if(aLHnSS[i]<0) aLHnSS[i]=0;
         if(aLHnOS[i]<0) aLHnOS[i]=0;
      }

      mini->Minimize();

      //save output
      const double *xs = mini->X();
      const double *err = mini->Errors();
      for (int etaBin = 0; etaBin< nEtaBins; etaBin++){
      for (int  ptBin = 0;  ptBin<  nPtBins;  ptBin++){
         int gBin = ptBin*nEtaBins + etaBin;
         // if (abs(xs[gBin]-1e-4)<1e-10) continue; //no change from init val
         hFlipProb->SetBinContent(etaBin+1, ptBin+1, xs[gBin]);
         hFlipProb->SetBinError(etaBin+1, ptBin+1, err[gBin]);
      }}
      std::cout << "Min func val : " << mini->MinValue()  << std::endl;

      write_matrix();

   }

   // -------------- TAG-AND-PROBE METHOD ----------------- //
   
   if(doTP){ 
      EMuTP(evt2lTree);
   }

   // missing systematic errors!

   // ---- MC TRUTH --- //
   if(isMC) DivideHist((TH1**) hEP_MC);

   // ---- DRAW ----- //
   draw();

   // CLEANUP
   if(toFile) {
      outputFile->Write();
      std::cout << "Files saved to " << out << std::endl;
   }
   if(findPtCorr){
      ptOut->Write();
      std::cout << "Pt correction histograms saved" << std::endl;
   }
   return;
}

TChain* loadData(TString fileList){
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

bool initialize(const string opt){

   // ------ Configure options ----- //
   if(opt != ""){
      int pos=-1;
      string var;
      do{
         pos++;
         var = opt.substr(pos, 3);
         pos += 3;
         if(var=="mr=") {
            Z_MASS_UPPER_LIMIT = std::stoi(opt.substr(pos, 3));
            pos++;
         }
         else if (var=="ml=") Z_MASS_LOWER_LIMIT = std::stoi(opt.substr(pos, 2));
         else if (var=="sb=") SIDEBAND_R = SIDEBAND_L = std::stoi(opt.substr(pos, 2));
         else if (var=="sl=") SIDEBAND_L = std::stoi(opt.substr(pos, 2));
         else if (var=="sr=") SIDEBAND_R = std::stoi(opt.substr(pos, 2));
         else if (var=="el=") {
            if(opt[pos]=='l') onlySignal = false;
            else if (opt[pos]=='s') onlySignal = true;
            else return false;
            pos--;
         }
         else if (var=="mc="){
            if(opt[pos]=='y') isMC = true;
            else if (opt[pos]=='n') isMC = false;
            else return false;
            pos--;

         }
         else if (var=="mx="){
            maxEntries = std::stoi(opt.substr(pos,3));
            pos++; 
         }
         else if (var=="nd=") {
            if(opt[pos]=='y') newDerivation = true;
            else if (opt[pos]=='n') newDerivation = false;
            else return false;
            pos--;
         }
         else if (var=="pt=") {
            if(opt[pos]=='y') findPtCorr = isMC && true;
            else if (opt[pos]=='n') findPtCorr = false;
            else return false;
            pos--;
         }
         else if (var=="tt=") {
            if(opt[pos]=='y') doTP = true;
            else if (opt[pos]=='n') doTP = false;
            else return false;
            pos--;
         }
         else {
            std::cout << "options in wrong format" << std::endl;
            return false;
         }
         pos += 2;
      } while (opt[pos]==',');
   }

   if(isMC) SIDEBAND_R = SIDEBAND_L = 0;

   std::cout << std::endl;
   if(isMC) std::cout << "MC "; else std::cout << "Data ";
   if(onlySignal) std::cout << "Signal ";
   else std::cout << "LooseBaseline ";
   std::cout << "electrons" << std::endl;
   if(doTP){
   	std::cout << "Electron muon tag-and-probe for ttbar" << std::endl;
   } else {
	   std::cout << Z_MASS_LOWER_LIMIT << " < Zmass < " << Z_MASS_UPPER_LIMIT << std::endl;
	   std::cout << "Sidebands: L=" << SIDEBAND_L << ", R=" << SIDEBAND_R << std::endl;
   }
   std::cout << std::endl;

   // ----------- MISID(PT, ETA) HISTOGRAMS ---------------- //
   if(doTP){
      hEP_TP[ERROR] = new TH2D("hTPFlip", "Flip prob (TP);#eta;p_{T}", nEtaBins, EtaEdges, nPtBins, PtEdges);
      hEP_TP[N_PROBES] = new TH2D("hEPTPNProbes", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
      hEP_TP[N_SSPROBES] = new TH2D("hEPTPNSSProbes", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
      hEP_TP[N_OSPROBES] = new TH2D("hEPTPNOSProbes", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   }

   // ------------- LIKELIHOOD METHOD ---------------- //
   int nTotBins = nEtaBins * nEtaBins * nPtBins * nPtBins;
   aLHnSS = new Double_t[nTotBins];
   aLHnOS = new Double_t[nTotBins];
   std::fill_n(aLHnOS, nTotBins, 0);
   std::fill_n(aLHnSS, nTotBins, 0);

   aSideLnOS = new Double_t[nTotBins]; std::fill_n(aSideLnOS, nTotBins, 0);
   aSideLnSS = new Double_t[nTotBins]; std::fill_n(aSideLnSS, nTotBins, 0);
   aSideRnOS = new Double_t[nTotBins]; std::fill_n(aSideRnOS, nTotBins, 0);
   aSideRnSS = new Double_t[nTotBins]; std::fill_n(aSideRnSS, nTotBins, 0);

   hFlipProb = new TH2D("hFlipProb", "Flip prob (likelihood method);#eta;p_{T}", nEtaBins, EtaEdges, nPtBins, PtEdges);

   // -------------- MC TRUTH -------------- //
   if (isMC){
      hEP_MC[ERROR] = new TH2D("hMCFlip", "Flip prob (MC);#eta;p_{T}", nEtaBins, EtaEdges, nPtBins, PtEdges);
      hEP_MC[N_PROBES] = new TH2D("hEPMCNel", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
      hEP_MC[N_SSPROBES] = new TH2D("hEPMCNelNoFlip", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
      hEP_MC[N_OSPROBES] = new TH2D("hEPMCNelFlipped", "", nEtaBins, EtaEdges, nPtBins, PtEdges);

      // srand(time(NULL)); // Random number generator for FSR photon
   }
   
   // ------------- OUTPUT ---------------------- //
   if(toFile) {
      outputFile = new TFile("output.root", "RECREATE");
      hFlipProb->SetDirectory(outputFile);
      if(doTP) hEP_TP[ERROR]->SetDirectory(outputFile);
      if(isMC) {
         hEP_MC[ERROR]->SetDirectory(outputFile); 
      }
   }
   if(debugMC) decayTree.open("decayTree.txt");

   // -------------- dPT ------------------------ //
   if(findPtCorr){
      ptOut = new TFile("dPThistos.root", "recreate");

      const int ndPTedges = 101;
      double dPTbins[ndPTedges];
      for(int i=0; i<ndPTedges; i++) dPTbins[i]=i-50;

      hDPT = new TH3D("hDPT", "#Delta p_{T} for all electrons;Reco |#eta|;Reco p_{T};#Delta p_{T} (GeV)", nEtaBins, EtaEdges, nPtBins, PtEdges, ndPTedges-1, dPTbins);
      hDPT->SetDirectory(ptOut);

      hDPTok = (TH3D*) hDPT->Clone("hDPTok");
      hDPTok->SetTitle("#Delta p_{T} for charge ok electrons");
      hDPTok->SetDirectory(ptOut);

      hDPTflipped = (TH3D*) hDPT->Clone("hDPTflipped");
      hDPTflipped->SetTitle("#Delta p_{T} for charge flipped electrons");
      hDPTok->SetDirectory(ptOut);
   }

   return true;
}

void printMC(const evt2l& evt2lTree){
   std::cout << "Printing decay info…" << std::endl
            << "#\tCharge\tEta\tPt\tDecay" << std::endl;
   for (int i=0; i<2; i++){
      char charge = (evt2lTree.leps_ID[i] > 0) ? '+' : '-';
      printf("%i\t%c\t%1.2f\t%3.2f\t", i, charge, evt2lTree.leps_eta[i], evt2lTree.leps_pt[i]);

      int part = evt2lTree.leps_truthI[i];
      int part1 = 0;
      int part2 = 0;
      while (part >= 0){
         std::cout << evt2lTree.truths_pdgId[part] <<' ';
         if (part1==part) break;
         if (part2==part) break;
         part2 = part1;
         part1 = part;
         part = evt2lTree.truths_motherI[part];
      }
      std::cout << std::endl;
   }
}

void processMCalt(evt2l& evt2lTree){
   for (int i=0; i<2; i++){
      if (evt2lTree.leps_truthI[i] < 0) continue;
      int charge = 0;
      if (evt2lTree.leps_truthType[i] == 4){
         int p = rand() % 100;
         if(p<50) charge = 1;
         else charge = -1;
      }
      else if (evt2lTree.leps_truthType[i] == 2){
         charge = (evt2lTree.leps_ID[0] > 0) - (evt2lTree.leps_ID[0] < 0);
      }
      if (charge == 0) continue;

      bool hasFlipped = ((charge*evt2lTree.leps_ID[i]) < 0);
      FillHist(evt2lTree.leps_eta[i], evt2lTree.leps_pt[i], hEP_MC, hasFlipped);

   }
}

int getTruthElecI(evt2l& evt2lTree, int i){
   // Find matching truth particle
   int tp = evt2lTree.leps_truthI[i];
   if (tp < 0) return -1; // tp = -1 if no truth particle;
   int matchedPDGID = evt2lTree.truths_pdgId[tp];
   if(debugMC) decayTree << "Particle " << i+1 << ": " << matchedPDGID << " ";
   // if(debugMC) printf("Particle %i: %i ", i, evt2lTree.truths_pdgId[tp]);
   if (abs(matchedPDGID) != 11) return -1;

   int mother = evt2lTree.truths_motherI[tp];

   if(doTP){
		while(mother>=0 && tp>=0 && abs(evt2lTree.truths_pdgId[mother])!=24 && abs(evt2lTree.truths_pdgId[mother]) < 100)
		{	// While mother and tp are valid, and mother is not W
         if(debugMC) decayTree << evt2lTree.truths_pdgId[mother] << " ";
			tp = mother;
			mother = evt2lTree.truths_motherI[tp];
		}
      if(debugMC) decayTree << evt2lTree.truths_pdgId[mother] << " ";

      if(tp<0)
      {  // Truth particle invalid
         if(debugMC) decayTree << "tp <0" << std::endl;
         return -1;
      }

      if(mother<0)
      {  // Mother invalid
         if(debugMC) decayTree << "mother <0" << std::endl;
         return -1;
      }

      if(abs(evt2lTree.truths_pdgId[mother])!=24)
      {  // Mother of particle not W
         if(debugMC) decayTree << "Mother pdgID = " << evt2lTree.truths_pdgId[mother] << std::endl;
         return -1;
      }

		if(abs(evt2lTree.truths_motherI[mother])<0)
		{	// If mother of W is invalid, return
         if(debugMC) decayTree << "Mother of W invalid" << std::endl;
			return -1;
		}

      if(abs(evt2lTree.truths_pdgId[evt2lTree.truths_motherI[mother]])!=6)
      {  // Mother of W not top
         if(debugMC) decayTree << "Mother of W = " << evt2lTree.truths_pdgId[evt2lTree.truths_motherI[mother]] << " != top" << std::endl;
         if(abs(evt2lTree.truths_pdgId[evt2lTree.truths_motherI[mother]])<=6 
            || abs(evt2lTree.truths_pdgId[evt2lTree.truths_motherI[mother]])==21
            || abs(evt2lTree.truths_pdgId[evt2lTree.truths_motherI[mother]])==9) return tp;
         return -1;
      }

		if (abs(evt2lTree.truths_pdgId[tp])!=11)
      {  // Particle found not electron
         if(debugMC) decayTree << "Particle from t -> W decay is " << evt2lTree.truths_pdgId[tp] << std::endl;
         return -1;
      }

      return tp;
   } else {
	   // Find original electron from Z
	   while(mother>=0 && tp >= 0 && evt2lTree.truths_pdgId[mother] != 23 && abs(evt2lTree.truths_pdgId[mother]) < 100){ // while mother exists and is not Z
	      if(debugMC) decayTree << evt2lTree.truths_pdgId[mother] << " ";
	      // if(debugMC) printf("%i ", evt2lTree.truths_pdgId[mother]);
	      tp = mother;
	      mother = evt2lTree.truths_motherI[tp];
	      // while (evt2lTree.truths_pdgId[tp] == evt2lTree.truths_pdgId[mother]){
	      //    tp = mother;
	      //    mother = evt2lTree.truths_motherI[tp];
	      // }
	   }
	   if(debugMC) decayTree << evt2lTree.truths_pdgId[mother] << " ";

	   // Reject if the particle found is not from Z
	   if (tp <0 || mother <0 || evt2lTree.truths_pdgId[mother]!=23){ 
	      if(debugMC){
	         if (tp < 0) decayTree << "tp <0" << std::endl;
	         if (mother < 0) decayTree << "mother <0" << std::endl;
	         else if (evt2lTree.truths_pdgId[mother]!=23) decayTree << "Mother pdgID = " << evt2lTree.truths_pdgId[mother] << std::endl;
	      }
	      return -1;
	   }

	   // Reject if particle found is not electron
	   if (abs(evt2lTree.truths_pdgId[tp])==11) return tp;
	   else {
	      if(debugMC) decayTree << "Orig from Z is " << evt2lTree.truths_pdgId[tp] << std::endl;
	      return -1;
	   }
	}
}

void processMC(evt2l& evt2lTree){
   // printMC(evt2lTree);
   for(int i=0; i<2; i++){
      // Skip if not electron, for emuTP
   	if(int(abs(evt2lTree.leps_ID[i]/1000))!=11) continue; 
      int tp = getTruthElecI(evt2lTree, i);
      if (tp<0) continue;

      int charge = (evt2lTree.truths_pdgId[tp] < 0) ? 1 : -1;
      bool hasFlipped = ((charge*evt2lTree.leps_ID[i]) < 0);

      double pt;
      if(MCtruthPt) pt = evt2lTree.truths_pt[tp];
      else pt = evt2lTree.leps_pt[i];

      double w = 1.0 * evt2lTree.evt_weight*evt2lTree.evt_ElSF*evt2lTree.evt_MuSF;
      if (applyPRW) w *= evt2lTree.evt_pwt;

      FillHist(evt2lTree.leps_eta[i], pt, hEP_MC, hasFlipped, w);

      if(findPtCorr){
         double deltaPt = evt2lTree.leps_pt[i] - evt2lTree.truths_pt[tp];
         hDPT->Fill(evt2lTree.leps_eta[i], evt2lTree.leps_pt[i], deltaPt, w);
         if(hasFlipped) hDPTflipped->Fill(evt2lTree.leps_eta[i], evt2lTree.leps_pt[i], deltaPt, w);
         else hDPTok->Fill(evt2lTree.leps_eta[i], evt2lTree.leps_pt[i], deltaPt, w);
      }
   }
   if(debugMC) decayTree << "\n";
   return; 
}

void fillCountArr(double* arr, double pt0, double eta0, double pt1, double eta1, double weight=1){
  //util to find appropiate bin and fill count to our countOS/SS array
  //std::cout << pt0 << " " << eta0 << " " << pt1 << " " << eta1 << std::endl;
   if(ABS_ETA){
      eta0 = abs(eta0);
      eta1 = abs(eta1);
   }

   if (pt0<PtEdges[0] || pt0>PtEdges[nPtBins]) return;
   if (pt1<PtEdges[0] || pt1>PtEdges[nPtBins]) return;
   if (eta0<EtaEdges[0] || eta0>EtaEdges[nEtaBins]) return;
   if (eta1<EtaEdges[0] || eta1>EtaEdges[nEtaBins]) return;

   int  pt0Bin = 0;
   int  pt1Bin = 0;
   int eta0Bin = 0;
   int eta1Bin = 0;

   for ( pt0Bin=0; pt0Bin< nPtBins; pt0Bin++){ if ( PtEdges[ pt0Bin+1]> pt0) break;}
   for ( pt1Bin=0; pt1Bin< nPtBins; pt1Bin++){ if ( PtEdges[ pt1Bin+1]> pt1) break;}
   for (eta0Bin=0;eta0Bin<nEtaBins;eta0Bin++){ if (EtaEdges[eta0Bin+1]>eta0) break;}
   for (eta1Bin=0;eta1Bin<nEtaBins;eta1Bin++){ if (EtaEdges[eta1Bin+1]>eta1) break;}


   int idx=0;
                     idx+=  pt0Bin; //i.e. 1234 = ((1*10+2)*10+3)*10+4
   idx*= nEtaBins;   idx+= eta0Bin;
   idx*=  nPtBins;   idx+=  pt1Bin;
   idx*= nEtaBins;   idx+= eta1Bin;

   arr[idx]+= weight;

   // std::cout << idx << std::endl;
}

void processEvents(evt2l& evt2lTree){
   // fill count of SS, OS event to pt,eta binned array  

   // Checks
   TFile* fChecks = TFile::Open("checks.root", "recreate");
   std::vector<TH1D*> hEta;
   std::vector<TH1D*> hPt;
   std::vector<TH1D*> hMass;
   TH1I* hCutflow = new TH1I("hCutflow", "Number of events passed", 8, 0, 8);
   hCutflow->GetXaxis()->SetBinLabel(1, "Total");
   hCutflow->GetXaxis()->SetBinLabel(2, "Trig+GRL+2e");
   hCutflow->GetXaxis()->SetBinLabel(3, "Pass LooseBaseline");
   hCutflow->GetXaxis()->SetBinLabel(4, "Pass Signal");
   hCutflow->GetXaxis()->SetBinLabel(5, "Zmass+SB");
   hCutflow->GetXaxis()->SetBinLabel(6, "Zmass");
   hCutflow->GetXaxis()->SetBinLabel(7, "Left SB");
   hCutflow->GetXaxis()->SetBinLabel(8, "Right SB");

   hEta.push_back(new TH1D("hEtaAll", "Eta distribution of all electrons", 200, -2.47, 2.47));
   hEta.push_back(new TH1D("hEtaPreselected", "Eta distribution of preselected electrons", 200, -2.47, 2.47));
   hEta.push_back(new TH1D("hEtaLoose", "Eta distribution of electrons in LooseBaseline pairs", 200, -2.47, 2.47));
   hEta.push_back(new TH1D("hEtaSignal", "Eta distribution of electrons in Signal pairs", 200, -2.47, 2.47));
   hEta.push_back(new TH1D("hEtaSignalZSB", "Eta distribution of electrons in Signal pairs within Zmass+SB window", 200, -2.47, 2.47));
   hEta.push_back(new TH1D("hEtaSignalZ", "Eta distribution of electrons in Signal pairs within Zmass window", 200, -2.47, 2.47));
   hEta.push_back(new TH1D("hEtaSignalLSB", "Eta distribution of electrons in Signal pairs within left SB", 200, -2.47, 2.47));
   hEta.push_back(new TH1D("hEtaSignalRSB", "Eta distribution of electrons in Signal pairs within right SB", 200, -2.47, 2.47));

   hPt.push_back(new TH1D("hPtAll", "Pt distribution of all electrons", 200, 20, 200));
   hPt.push_back(new TH1D("hPtPreselected", "Pt distribution of preselected electrons", 200, 20, 200));
   hPt.push_back(new TH1D("hPtLoose", "Pt distribution of electrons in LooseBaseline pairs", 200, 20, 200));
   hPt.push_back(new TH1D("hPtSignal", "Pt distribution of electrons in Signal pairs", 200, 20, 200));
   hPt.push_back(new TH1D("hPtSignalZSB", "Pt distribution of electrons in Signal pairs within Zmass+SB window", 200, 20, 200));
   hPt.push_back(new TH1D("hPtSignalZ", "Pt distribution of electrons in Signal pairs within Zmass window", 200, 20, 200));
   hPt.push_back(new TH1D("hPtSignalLSB", "Pt distribution of electrons in Signal pairs within left SB", 200, 20, 200));
   hPt.push_back(new TH1D("hPtSignalRSB", "Pt distribution of electrons in Signal pairs within right SB", 200, 20, 200));

   hMass.push_back(new TH1D("hMassAll", "Mass distribution of all electrons", 200, 20, 200));
   hMass.push_back(new TH1D("hMassPreselected", "Mass distribution of preselected electrons", 200, 20, 200));
   hMass.push_back(new TH1D("hMassLoose", "Mass distribution of electrons in LooseBaseline pairs", 200, 20, 200));
   hMass.push_back(new TH1D("hMassSignal", "Mass distribution of electrons in Signal pairs", 200, 20, 200));
   hMass.push_back(new TH1D("hMassSignalZSB", "Mass distribution of electrons in Signal pairs within Zmass+SB window", 200, 20, 200));
   hMass.push_back(new TH1D("hMassSignalZ", "Mass distribution of electrons in Signal pairs within Zmass window", 200, 20, 200));
   hMass.push_back(new TH1D("hMassSignalLSB", "Mass distribution of electrons in Signal pairs within left SB", 200, 20, 200));
   hMass.push_back(new TH1D("hMassSignalRSB", "Mass distribution of electrons in Signal pairs within right SB", 200, 20, 200));


   #define CUTFLOW(i)                                                                            \
   hCutflow->Fill(i);                                                                            \
   hEta[i]->Fill(evt2lTree.leps_eta[0], weight); hEta[i]->Fill(evt2lTree.leps_eta[1], weight);   \
   hPt[i]->Fill(evt2lTree.leps_pt[0], weight); hPt[i]->Fill(evt2lTree.leps_pt[1], weight);       \
   hMass[i]->Fill(mll, weight);                                                                  

   long long nEntries = evt2lTree.fChain->GetEntries();
   for (long long i =0; i<nEntries; i++){
      if(maxEntries != 0) if (i>maxEntries) return;
      loadbar(i+1, nEntries); // loadbar
      evt2lTree.GetEntry(i);

      double weight = evt2lTree.evt_weight*evt2lTree.evt_ElSF*evt2lTree.evt_MuSF;
      if(applyPRW) weight *= evt2lTree.evt_pwt;
      double mll = COMMON::getMll(evt2lTree);

      CUTFLOW(0);
      // ----- EVENT SELECTION -------//
      // At least one trigger, exactly two leptons, and leptons 0 and 1 are both electrons
      bool baseCut = (evt2lTree.sig_trigCode>0 && evt2lTree.leps_ == 2
                     && int(abs(evt2lTree.leps_ID[0]/1000))==11 && int(abs(evt2lTree.leps_ID[1]/1000))==11); 
      if(!newDerivation) baseCut = baseCut && evt2lTree.evt_cuts==1;
      if(!baseCut) continue;

      CUTFLOW(1);

      if(!newDerivation){
         bool isLooseBaseLine = ((evt2lTree.leps_lFlag[0] & 1<<0) && (evt2lTree.leps_lFlag[1] & 1<<0));
         if (!isLooseBaseLine) continue;
      }

      CUTFLOW(2);

      bool isSignal[2]; 
      if(newDerivation){
         isSignal[0] = (evt2lTree.leps_lFlag[0] & 2)/2;
         isSignal[1] = (evt2lTree.leps_lFlag[1] & 2)/2;
      }
      else {
         isSignal[0] = (evt2lTree.leps_lFlag[0] & 1<<1);
         isSignal[1] = (evt2lTree.leps_lFlag[1] & 1<<1);
      }
      
      bool iseeOS = ((evt2lTree.leps_ID[0]>0) != (evt2lTree.leps_ID[1]>0));
      bool iseeSS = !iseeOS;

      if(onlySignal && !(isSignal[0] && isSignal[1])) continue;

      CUTFLOW(3);

      // Z->ee events
      bool zeeCut = false;
      bool zeeLeft = false;
      bool zeeRight = false;

      if (mll >= Z_MASS_LOWER_LIMIT && mll <= Z_MASS_UPPER_LIMIT) zeeCut = true;
      else if (SIDEBAND_L!=0 && mll > (Z_MASS_LOWER_LIMIT - SIDEBAND_L) && mll < Z_MASS_LOWER_LIMIT) zeeLeft = true;
      else if (SIDEBAND_R!=0 && mll > Z_MASS_UPPER_LIMIT && mll < (Z_MASS_UPPER_LIMIT+SIDEBAND_R)) zeeRight = true;
      
      if (!(zeeCut || zeeLeft || zeeRight)) continue;

      CUTFLOW(4);

      double pt1, pt2;
      if(MCtruthPt){
         int truth1 = getTruthElecI(evt2lTree, 0);
         int truth2 = getTruthElecI(evt2lTree, 1);
         if(truth1 < 0 || truth2 < 0) continue;
         pt1 = evt2lTree.truths_pt[truth1];
         pt2 = evt2lTree.truths_pt[truth2];
      }
      else{
         pt1 = evt2lTree.leps_pt[0];
         pt2 = evt2lTree.leps_pt[1];
      }

      // ================= LIKELIHOOD =================== //
      if (iseeOS){ 
         if(zeeLeft) {
            fillCountArr(aSideLnOS, pt1, evt2lTree.leps_eta[0], pt2, evt2lTree.leps_eta[1], weight);
            CUTFLOW(6);
         }
         else if(zeeRight){
            fillCountArr(aSideRnSS, pt1, evt2lTree.leps_eta[0], pt2, evt2lTree.leps_eta[1], weight);
            CUTFLOW(7);
         }
         else{
            fillCountArr(aLHnOS, pt1, evt2lTree.leps_eta[0], pt2, evt2lTree.leps_eta[1], weight);
            CUTFLOW(5);
         }
      }

      if (iseeSS){
         if(zeeLeft){
            fillCountArr(aSideLnSS, pt1, evt2lTree.leps_eta[0], pt2, evt2lTree.leps_eta[1], weight);
            CUTFLOW(6);
            
         }
         else if(zeeRight){
            fillCountArr(aSideRnSS, pt1, evt2lTree.leps_eta[0], pt2, evt2lTree.leps_eta[1], weight);
            CUTFLOW(7);
         }
         else{
            fillCountArr(aLHnSS, pt1, evt2lTree.leps_eta[0], pt2, evt2lTree.leps_eta[1], weight);
            CUTFLOW(5);
         }
      }

      // ============ MC TRUTH ============//
      if(isMC && zeeCut) processMC(evt2lTree);
      // std::cout << "\n";
      // if (i>2500) break;
   }

   #undef CUTFLOW

   for(auto h : hPt) {h->SetDirectory(fChecks); h->Write();}
   for(auto h : hEta) {h->SetDirectory(fChecks); h->Write();}
   for(auto h : hMass) {h->SetDirectory(fChecks); h->Write();}
   hCutflow->SetDirectory(fChecks); hCutflow->Write();
   fChecks->Close();
}

void EMuTP(evt2l& evt2lTree){
   // Checks
   TFile* fChecks = TFile::Open("checks.root", "recreate");
   std::vector<TH1D*> hEta;
   std::vector<TH1D*> hPt;
   std::vector<TH1D*> hMass;
   TH1I* hCutflow = new TH1I("hCutflow", "Number of events passed", 6, 0, 6);
   hCutflow->GetXaxis()->SetBinLabel(1, "Total");
   hCutflow->GetXaxis()->SetBinLabel(2, "Trig+GRL+emu");
   hCutflow->GetXaxis()->SetBinLabel(3, "nBJets==1OR2");
   hCutflow->GetXaxis()->SetBinLabel(4, "Signal Mu");
   hCutflow->GetXaxis()->SetBinLabel(5, "Signal E");

   hEta.push_back(new TH1D("hEtaAll", "Eta distribution of all leptons", 200, -2.47, 2.47));
   hEta.push_back(new TH1D("hEtaPreselected", "Eta distribution of preselected leptons", 200, -2.47, 2.47));
   hEta.push_back(new TH1D("hEta1or2b", "Eta distribution of leptons in pairs with nBJets==1OR2", 200, -2.47, 2.47));
   hEta.push_back(new TH1D("hEtaESigMu", "Eta distribution of leptons in emu pairs with signal mu", 200, -2.47, 2.47));
   hEta.push_back(new TH1D("hEtaSigESigMu", "Eta distribution of leptons in emu pairs with signal emu", 200, -2.47, 2.47));

   hPt.push_back(new TH1D("hPtAll", "Pt distribution of all leptons", 200, 20, 200));
   hPt.push_back(new TH1D("hPtPreselected", "Pt distribution of preselected leptons", 200, 20, 200));
   hPt.push_back(new TH1D("hPt1or2b", "Pt distribution of leptons in pairs with nBJets==1OR2", 200, 20, 200));
   hPt.push_back(new TH1D("hPtESigMu", "Pt distribution of leptons in emu pairs with signal mu", 200, 20, 200));
   hPt.push_back(new TH1D("hPtSigESigMu", "Pt distribution of leptons in emu pairs with signal emu", 200, 20, 200));

   hMass.push_back(new TH1D("hMassAll", "Mass distribution of all leptons", 200, 20, 200));
   hMass.push_back(new TH1D("hMassPreselected", "Mass distribution of preselected leptons", 200, 20, 200));
   hMass.push_back(new TH1D("hMass1or2b", "Mass distribution of leptons in pairs with nBJets==1OR2", 200, 20, 200));
   hMass.push_back(new TH1D("hMassESigMu", "Mass distribution of leptons in emu pairs with signal mu", 200, 20, 200));
   hMass.push_back(new TH1D("hMassSigESigMu", "Mass distribution of leptons in emu pairs with signal emu", 200, 20, 200));

   #define CUTFLOW(i)                                                                            \
   hCutflow->Fill(i);                                                                            \
   hEta[i]->Fill(evt2lTree.leps_eta[0], weight); hEta[i]->Fill(evt2lTree.leps_eta[1], weight);   \
   hPt[i]->Fill(evt2lTree.leps_pt[0], weight); hPt[i]->Fill(evt2lTree.leps_pt[1], weight);       \
   hMass[i]->Fill(mll, weight);                                                                  

   long long nEntries = evt2lTree.fChain->GetEntries();
   for (long long i =0; i<nEntries; i++){
      if(maxEntries != 0) if (i>maxEntries) return;
      loadbar(i+1, nEntries); // loadbar
      evt2lTree.GetEntry(i);

      double weight = evt2lTree.evt_weight*evt2lTree.evt_ElSF*evt2lTree.evt_MuSF;
      if(applyPRW) weight *= evt2lTree.evt_pwt;
      double mll = COMMON::getMll(evt2lTree);

      CUTFLOW(0);
      // ----- EVENT SELECTION -------//
      // At least one trigger, exactly two leptons, and leptons 0 and 1 form emu pair
      bool baseCut = (evt2lTree.sig_trigCode>0 && evt2lTree.leps_ == 2);
      bool isEMu = (int(abs(evt2lTree.leps_ID[0]/1000))==13 && int(abs(evt2lTree.leps_ID[1]/1000))==11)
                  || (int(abs(evt2lTree.leps_ID[0]/1000))==11 && int(abs(evt2lTree.leps_ID[1]/1000))==13);
      if(!(baseCut && isEMu)) continue;

      CUTFLOW(1);

      int nBJet = 0;
      for(int k=0;k<evt2lTree.jets_;k++){
         if((evt2lTree.jets_jFlag[k] & 1<<5)) nBJet++;
      }

      if(nBJet!=1 && nBJet!=2) continue;
      CUTFLOW(2);

      for(int i=0; i<2; i++){
         if(!(int(abs(evt2lTree.leps_ID[i]/1000))==13 && int(abs(evt2lTree.leps_ID[1-i]/1000))==11)) continue;

         if(!((evt2lTree.leps_lFlag[i] & 2)/2)) continue;
         CUTFLOW(3);

         if(onlySignal && !((evt2lTree.leps_lFlag[1-i] & 2)/2)) continue;
         CUTFLOW(4);

         bool isEMuSS = ((evt2lTree.leps_ID[i])*(evt2lTree.leps_ID[1-i])>0);

         FillHist(evt2lTree.leps_eta[1-i], evt2lTree.leps_pt[1-i], (TH2**) hEP_TP, isEMuSS, weight);
      }

      if(isMC) processMC(evt2lTree);
   }

   DivideHist((TH1**) hEP_TP);
}

void FillHist(double eta, double pt, TH2* h[4], bool hasFlipped, double w){
   if(ABS_ETA) eta = fabs(eta);
   h[N_PROBES]->Fill(eta, pt, w);
   if(hasFlipped) h[N_SSPROBES]->Fill(eta, pt, w);
   else h[N_OSPROBES]->Fill(eta, pt, w);

   return;
}

void DivideHist(TH1* h[]){
   h[N_PROBES]->Sumw2();
   h[N_SSPROBES]->Sumw2();
   if(h[N_OSPROBES] != 0) h[N_OSPROBES]->Sumw2();

   h[ERROR]->Add(h[N_SSPROBES]);
   h[ERROR]->Divide(h[N_PROBES]);
}

double getCountArr(const double* arr, int pt0Bin, int eta0Bin, int pt1Bin, int eta1Bin){
  //util to get count from our aLHnOS/SS array given intuitive index

   if (pt0Bin<0 || pt0Bin>=nPtBins) return -1;
   if (pt1Bin<0 || pt1Bin>=nPtBins) return -1;
   if (eta0Bin<0 || eta0Bin>=nEtaBins) return -1;
   if (eta1Bin<0 || eta1Bin>=nEtaBins) return -1;

   int idx=0;
                     idx+=  pt0Bin; //i.e. 1234 = ((1*10+2)*10+3)*10+4
   idx*= nEtaBins;   idx+= eta0Bin;
   idx*=  nPtBins;   idx+=  pt1Bin;
   idx*= nEtaBins;   idx+= eta1Bin;

   return arr[idx];
}

void subtractBkg(double* central, const double* left, const double* right){
   int nTotBins = nEtaBins * nEtaBins * nPtBins * nPtBins;
   double weight = 1.0/(SIDEBAND_L+SIDEBAND_R);
   for(int i = 0; i<nTotBins; i++){
      central[i] = central[i] - weight*(SIDEBAND_L*left[i] + SIDEBAND_R*right[i]);
      central[i] = (central[i]<0)? 0 : central[i];
   }
}

double likelihood(const double *x){
  //cost function to minimize by optimizing charge flip prob at different pt, eta
   double cost = 0.;

   for (int  pt0Bin = 0;  pt0Bin<  nPtBins;  pt0Bin++){
   for (int  pt1Bin = 0;  pt1Bin<  nPtBins;  pt1Bin++){
   for (int eta0Bin = 0; eta0Bin< nEtaBins; eta0Bin++){
   for (int eta1Bin = 0; eta1Bin< nEtaBins; eta1Bin++){
      double nOS  = getCountArr(aLHnOS, pt0Bin, eta0Bin, pt1Bin, eta1Bin);
      double nSS  = getCountArr(aLHnSS, pt0Bin, eta0Bin, pt1Bin, eta1Bin);
      double nAll = nOS + nSS;
      double eff  = x[pt0Bin*nEtaBins + eta0Bin] + x[pt1Bin*nEtaBins + eta1Bin];
      if (nAll==0) continue;
      cost -= log(nAll*eff)*nSS - nAll*eff;
   }}}}
   return cost;
}

void draw(){

   // ===== Set styles ======
   gStyle->SetOptStat(0);

   TString measurementLabel = "From Z#rightarrow ee LH";
   if(doTP){
      hFlipProb = hEP_TP[ERROR];
      measurementLabel = "From e#mu TP";
   }

   hFlipProb->SetMarkerStyle(20);
   hFlipProb->SetMarkerColor(kGreen+4);
   hFlipProb->SetLineColor(kGreen+4);

   if(isMC){
      hEP_MC[ERROR]->SetMarkerStyle(20);
      hEP_MC[ERROR]->SetMarkerColor(kRed); 
      hEP_MC[ERROR]->SetLineColor(kRed);
   }

   // ======= Draw histograms ======= 

   if(drawEtaPt){
      TCanvas *cEtaPt[nPtBins];
      TLegend *lEtaPt[nPtBins];
   
      string sPdfFilename = "EtaPt.pdf";
      TCanvas c("c", "c", 600, 600);
   
      if (print) c.Print((sPdfFilename+"[").c_str());

      for(int ptBin=0; ptBin<nPtBins; ptBin++){
         // Set names of Canvases
         std::ostringstream lowPt, highPt;
         lowPt << std::fixed << std::setprecision(0) << PtEdges[ptBin];
         highPt << std::fixed << std::setprecision(0) << PtEdges[ptBin+1];

         
         string sCanTitle = "Error rate as a function of eta in " + lowPt.str() + "<= Pt < " + highPt.str();
         string sCanPad = "cEtaPt" + lowPt.str();
         string sHistName = "hEP" + lowPt.str();
         string sHistName_MC = sHistName + "_MC";
      
         cEtaPt[ptBin] = new TCanvas(sCanPad.c_str(), sCanTitle.c_str(), 600, 600);
         cEtaPt[ptBin]->SetLogy();

         lEtaPt[ptBin] = new TLegend(0.2, 0.8, 0.65, 0.88);
         string sLegendHeader = lowPt.str() + " #leq p_{T} < " + highPt.str();
         lEtaPt[ptBin]->SetHeader(sLegendHeader.c_str());

         // ----------------- Measurement plot --------------- //
         TH1* hEtaPt_LH = 0;
         hEtaPt_LH = hFlipProb->ProjectionX((sHistName+"LH").c_str(), ptBin+1, ptBin+1, "e");

         hEtaPt_LH->GetYaxis()->SetRangeUser(0.0001, 0.1);
         hEtaPt_LH->Draw("e1");
         if(ABS_ETA) hEtaPt_LH->GetXaxis()->SetTitle("|#eta|");

         lEtaPt[ptBin]->AddEntry(hEtaPt_LH, measurementLabel, "lpf");

         // --------------- MC plots --------------- //

         if(isMC){
            TH1* hEtaPt_MC = hEP_MC[ERROR]->ProjectionX((sHistName+"MC").c_str(), ptBin+1, ptBin+1, "e");
            hEtaPt_MC->Draw("e1 same");
            lEtaPt[ptBin]->AddEntry(hEtaPt_MC, "From MC truth", "lpf");
         }

         lEtaPt[ptBin]->Draw();

         if(print){
            string sPdfTitle = "Title:" + sLegendHeader;
            cEtaPt[ptBin]->Print(sPdfFilename.c_str(), sPdfTitle.c_str());
         }
      }

      if(print) c.Print((sPdfFilename+"]").c_str());
   }

   if(drawPtEta){
      TCanvas *cPtEta[nEtaBins];
      TLegend *lPtEta[nEtaBins];
   
      string sPdfFilename = "PtEta.pdf";
      TCanvas c("c", "c", 600, 600);
   
      if (print) c.Print((sPdfFilename+"[").c_str());

      for(int etaBin=0; etaBin<nEtaBins; etaBin++){
         // Set names of Canvases
         std::ostringstream lowEta, highEta, lowEta100;
         lowEta << std::fixed << std::setprecision(2) << EtaEdges[etaBin];
         lowEta100 << std::fixed << std::setprecision(0) << EtaEdges[etaBin]*100;
         highEta << std::fixed << std::setprecision(2) << EtaEdges[etaBin+1];

         
         string sCanTitle = "Error rate as a function of pt in " + lowEta.str() + "<=";
         if(ABS_ETA) sCanTitle += "|Eta| < " + highEta.str();
         else sCanTitle += "Eta  < " + highEta.str();
         string sCanPad = "cPtEta" + lowEta100.str();
         string sHistName = "hPE" + lowEta100.str();
         string sHistName_MC = sHistName + "_MC";
      
         cPtEta[etaBin] = new TCanvas(sCanPad.c_str(), sCanTitle.c_str(), 600, 600);
         cPtEta[etaBin]->SetLogy();
         if(logPt)cPtEta[etaBin]->SetLogx();

         lPtEta[etaBin] = new TLegend(0.2, 0.8, 0.65, 0.88);
         string sLegendHeader = lowEta.str() + " #leq #eta < " + highEta.str();
         lPtEta[etaBin]->SetHeader(sLegendHeader.c_str());

         // ----------------- Measurement plot --------------- //
         TH1* hPtEta_LH = 0;
         hPtEta_LH = hFlipProb->ProjectionY((sHistName+"LH").c_str(), etaBin+1, etaBin+1, "e");

         hPtEta_LH->GetYaxis()->SetRangeUser(0.0001, 0.1);
         hPtEta_LH->Draw("e1");

         lPtEta[etaBin]->AddEntry(hPtEta_LH, measurementLabel, "lpf");

         // -----------------  MC Truth ---------------- //
         if(isMC){
            TH1* hPtEta_MC = hEP_MC[ERROR]->ProjectionY((sHistName+"MC").c_str(), etaBin+1, etaBin+1, "e");
            hPtEta_MC->Draw("e1 same");
            lPtEta[etaBin]->AddEntry(hPtEta_MC, "From MC truth", "lpf");
         }

         lPtEta[etaBin]->Draw();

         if(print){
            string sPdfTitle = "Title:" + sLegendHeader;
            cPtEta[etaBin]->Print(sPdfFilename.c_str(), sPdfTitle.c_str());
         }
      }

      if(print) c.Print((sPdfFilename+"]").c_str());
   }

   if(drawPt){
      std::cout << "Now this part" << std::endl;
      TCanvas* cPt = new TCanvas("cPt", "cPt", 600, 600); cPt->SetLogy();
      if(logPt) cPt->SetLogx();
      TLegend* lPt = new TLegend(0.2, 0.8, 0.65, 0.88);

      TH1* hPt_LH = hFlipProb->ProjectionY("hPt_LH");
      hPt_LH->GetYaxis()->SetRangeUser(0.0005, 0.5);
      hPt_LH->Draw("e1");
      lPt->AddEntry(hPt_LH, measurementLabel, "lpf");


      if(isMC){
         TH1* hPt_MC = hEP_MC[ERROR]->ProjectionY("hPt_MC");
         hPt_MC->Draw("e1 same");
         lPt->AddEntry(hPt_MC, "From MC", "lpf");
      }

      lPt->Draw();
      std::cout <<"Elephants ";
      if(print){
         cPt->Print("Pt.pdf", "Title:Pt");
      }
      std::cout <<"in love";
   }

   if(drawEta){
      TCanvas* cEta = new TCanvas("cEta", "cEta", 600, 600); cEta->SetLogy();
      TLegend* lEta = new TLegend(0.2, 0.8, 0.65, 0.88);

      TH1* hEta_LH = hFlipProb->ProjectionX("hEta_LH");
      hEta_LH->GetYaxis()->SetRangeUser(0.0005, 0.5);
      hEta_LH->Draw("e1");
      lEta->AddEntry(hEta_LH, measurementLabel, "lpf");

      if(isMC){
         TH1* hEta_MC = hEP_MC[ERROR]->ProjectionX("hEta_MC");
         hEta_MC->Draw("e1 same");
         lEta->AddEntry(hEta_MC, "From MC", "lpf");
      }

      lEta->Draw();
      if(print){
         cEta->Print("Eta.pdf", "Title:Eta");
      }
   }

   gStyle->SetPalette(kLightTerrain);

   if(doTP){
      TCanvas cTP("cTP", "cTP", 1280, 720);
      TH2D* hTP = (TH2D*) hEP_TP[ERROR]->Clone("hTP"); 
      hTP->SetMarkerColor(kBlack);
      if(logPt) cTP.SetLogy();
      hTP->Draw("colz text e");
      cTP.Print("hTP.pdf", "Title:hTP");
      // if(toFile) {hTP->SetDirectory(outputFile); hTP->Write();}
   } else {
      TCanvas cLH("cLH", "cLH", 1280, 720); 
      TH2D* hLH = (TH2D*) hFlipProb->Clone("hLH"); 
      hLH->SetMarkerColor(kBlack);
      if(logPt) cLH.SetLogy();
      hLH->Draw("colz text e"); 
      cLH.Print("hLH.pdf", "Title:hLH");
   }
   if(isMC) {
      TCanvas cMC("cMC", "cMC", 1280, 720);
      TH2D* hMC = (TH2D*) hEP_MC[ERROR]->Clone("hMC"); 
      hMC->SetMarkerColor(kBlack);
      if(logPt) cMC.SetLogy();
      hMC->Draw("colz text e");
      cMC.Print("hMC.pdf", "Title:hMC");
      // if(toFile) {hMC->SetDirectory(outputFile); hMC->Write();}
   }
}

void write_matrix()
{
   std::ofstream fileN("n.csv");
   std::ofstream fileNSS("nss.csv");
   fileN   << std::fixed << std::setprecision(1); 
   fileNSS << std::fixed << std::setprecision(1);
   for(int eta0Bin=0; eta0Bin<nEtaBins; eta0Bin++){
   for(int pt0Bin=0 ; pt0Bin <nPtBins ; pt0Bin++ ){
   for(int eta1Bin=0; eta1Bin<nEtaBins; eta1Bin++){
   for(int pt1Bin=0 ; pt1Bin <nPtBins ; pt1Bin++ ){
      int idx=0;
                        idx+=  pt0Bin; //i.e. 1234 = ((1*10+2)*10+3)*10+4
      idx*= nEtaBins;   idx+= eta0Bin;
      idx*=  nPtBins;   idx+=  pt1Bin;
      idx*= nEtaBins;   idx+= eta1Bin;

      if(eta1Bin==0 && pt1Bin==0){ fileN << aLHnOS[idx] + aLHnSS[idx]; fileNSS << aLHnSS[idx]; }
      else {fileN << "," << aLHnSS[idx] + aLHnOS[idx]; fileNSS << "," << aLHnSS[idx];}

      if(eta1Bin==nEtaBins-1 && pt1Bin==nPtBins-1) {fileN << std::endl; fileNSS << std::endl;}

   }}}}
   fileN.close(); fileNSS.close();

}