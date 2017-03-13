/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/*
 TUAN M. NGUYEN (PhD Student; Supervisor: Jean-Francois Arguin)
 Universite de Montreal
 tuan.nguyen.manh@cern.ch
 
 Electron charge flip
 
 */
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/



#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>

#include "TChain.h"

// for testPlots()
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"

using namespace std;

typedef vector<vector<float> > matrix;
typedef vector<float> row;

enum inputType {EGAMMA, SUSY};
inputType in=SUSY;

enum sel{LOOSEBASELINE, SIGNAL};
sel selection=SIGNAL;

bool passQID=true;

//##################################
//# Choose the bool for MC or DATA! #

// bool m_DT = false; //for MC
bool m_DT = true; //for DATA

bool applyPRW = !m_DT && true;

/***************************************************************************/
/***************************************************************************/
/* Binning information */

// row VETAS = {0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.47}; 
row VETAS = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.37, 1.52, 1.8, 2, 2.2, 2.47};
row VPTS  = {20, 30, 40, 50, 60, 80, 120, 1000.0}; 

// PTSCALE makes sure GeV is used consistently in the code. If the input data
// is in unit of MeV then PTSCALE = 1000.0. If the input data is in unit of
// GeV then PTSCALE should be set to 1.0, so that no further conversion is done.
// The quantities affected are pt and mll.
float PTSCALE  = 1.0;

//##################################

string TAG = "";
//float m_weight = 1.;

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

/***************************************************************************/
/***************************************************************************/
/*
 The input ntuple is required to have the following variables (the meaning is
 obvious from the name)
 
 float weight
 float mll
 int ch1, float pt1, float eta1
 int ch2, float pt2, float eta2
 
 *** IMPORTANT ***
 -----------------
 THE DATA TYPES OF YOUR ROOT FILES MUST BE OF TYPES AS BELOW. OTHERWISE
 STRANGE THINGS MAY HAPPEN
 -----------------
 *****************
 
*/

class Input {
public:
  
  float MCEvtWeight;
  float MCPileupWeight;
  float Zcand_M;

  unsigned long int trigCode;
  
  // Electron 1
  int elCand1_charge;
  float elCand1_pt;
  float elCand1_cl_eta;
  float elCand1_phi;
  int elCand1_ID;
  int elCand1_origCharge;
  bool elCand1_qID;
 
  // Electron 2
  int elCand2_charge;
  float elCand2_pt;
  float elCand2_cl_eta;
  float elCand2_phi;
  int elCand2_ID;
  int elCand2_origCharge;
  bool elCand2_qID; 

  // == Variables for cuts (EGamma)
  int isTagTag;

  int elCand1_isLooseAndBLayerLL2015_v8;
  int elCand2_isLooseAndBLayerLL2015_v8;
  int elCand1_isMediumLL2015_v8;
  int elCand2_isMediumLL2015_v8;
  int elCand1_isTightLL2015_v8;
  int elCand2_isTightLL2015_v8;

  int elCand1_passd0z0; // |d0sig| < 5
  int elCand2_passd0z0; // |z0 sin(theta)| < 0.5mm

  int elCand1_isolGradientLoose;
  int elCand2_isolGradientLoose;

  // == Variables for cuts (SUSY)
  int elCand1_flag;
  int elCand2_flag;

  // == Variables for cuts (both)

  // == Truth pt
  // float elCand1_truthPt;
  // float elCand2_truthPt;


  // int elCand1_expectHitInBLayer;
  // int elCand2_expectHitInBLayer;


};

void connect_input(Input &ip, TChain &events) {
  events.SetBranchStatus("*",0);
  #define CONNECT(b) ip.b = 0; events.SetBranchStatus(#b,1); events.SetBranchAddress(#b,&(ip.b)); 

  CONNECT(MCEvtWeight)
  CONNECT(MCPileupWeight); 
  CONNECT(Zcand_M);
  
  // Electron 1
  CONNECT(elCand1_pt);
  CONNECT(elCand1_cl_eta);
  CONNECT(elCand1_charge);
  CONNECT(elCand1_phi);
  
  // Electron 2
  CONNECT(elCand2_pt);
  CONNECT(elCand2_cl_eta);
  CONNECT(elCand2_charge);
  CONNECT(elCand2_phi);


  // Variables for cuts
  if(in==EGAMMA){  
    CONNECT(isTagTag);
    CONNECT(elCand1_isLooseAndBLayerLL2015_v8);
    CONNECT(elCand2_isLooseAndBLayerLL2015_v8);
    CONNECT(elCand1_isMediumLL2015_v8);
    CONNECT(elCand2_isMediumLL2015_v8);
    CONNECT(elCand1_isTightLL2015_v8);
    CONNECT(elCand2_isTightLL2015_v8);

    CONNECT(elCand1_passd0z0);
    CONNECT(elCand2_passd0z0);

    CONNECT(elCand1_isolGradientLoose);
    CONNECT(elCand2_isolGradientLoose);

    // CONNECT(elCand1_expectHitInBLayer);
    // CONNECT(elCand2_expectHitInBLayer);
  }

  if(in==SUSY){
    CONNECT(elCand1_flag);
    CONNECT(elCand2_flag);
    CONNECT(trigCode);
    CONNECT(elCand1_ID);
    CONNECT(elCand2_ID);
    CONNECT(elCand1_origCharge);
    CONNECT(elCand2_origCharge);
    CONNECT(elCand1_qID);
    CONNECT(elCand2_qID);
  }

  #undef CONNECT
}


/***************************************************************************/
/* Matrix operations */
/***************************************************************************/
float row_min(row i){
  return *min_element(i.begin(), i.end());
}

float row_max(row i) {
  return *max_element(i.begin(), i.end());
}


matrix s_mul(float s, matrix A) {
  matrix M(A.size(), row(A[0].size(), 0.0));
  
  for (size_t i(0); i < M.size(); i++){
    for (size_t j(0); j < M[i].size(); j++){
      M[i][j] += (A[i][j] * s);
    }
  }
  return M;

}
matrix m_add(matrix A, matrix B) {
  matrix M(A.size(), row(A[0].size(), 0.0));
  
  for (size_t i(0); i < M.size(); i++){
    for (size_t j(0); j < M[i].size(); j++){
      M[i][j] = (A[i][j] + B[i][j]);
    }
  }
  return M;
}

matrix m_sub(matrix A, matrix B) {
  return m_add(A, s_mul(-1.0, B));
}

/* Return a matrix whose negative values become zero */
matrix negative_to_zero(matrix M){
  matrix N(M.size(), row(M[0].size(), 0.0));
  
  for (size_t i(0); i < M.size(); i++){
    for (size_t j(0); j < M[i].size(); j++){
      if (M[i][j] >= 0)
        N[i][j] = M[i][j];
    }
  }
  return N;
}

/* ofn is outfile name, "x.txt" */
void write_matrix_to_file(string ofn, matrix M){
  ofstream text;
  text.open(ofn);
  for (size_t i(0); i < M.size(); i++){
    text << M[i][0];
    for (size_t j(1); j < M[i].size(); j++){
      text << "," << M[i][j];
    }
    text << endl;
  }
}

/***************************************************************************/
/***************************************************************************/
/* Check if an event has pt or eta of either electron in the pair out of range */
bool is_out_eta_pt(float elCand1_eta, float elCand1_pt,
                   float elCand2_eta, float elCand2_pt,
                   float min_eta, float max_eta,
                   float min_pt, float max_pt)
{
  return
  (elCand1_eta < min_eta) || (elCand1_eta > max_eta) ||
  (elCand2_eta < min_eta) || (elCand2_eta > max_eta) ||
  (elCand1_pt < min_pt) || (elCand1_pt > max_pt) ||
  (elCand2_pt < min_pt) || (elCand2_pt > max_pt);
}

bool is_out_Zcand_M(float Zcand_M, float lZcand_M, float rZcand_M, float bl, float br){
  return
  (Zcand_M < (lZcand_M - bl)) ||
  (Zcand_M > (rZcand_M + br));
}


/* 
 Assume eta, pt are already within range
 The return -1 is thus only there to remove C++ warning
 */

int eta_bin(float eta){
  for (size_t i(1); i < VETAS.size(); i++){
    if (eta <= VETAS[i])
      return i-1;
  }
  return -1;
}

int pt_bin(float pt){
  for (size_t i(1); i < VPTS.size(); i++){
    if (pt <= VPTS[i])
      return (i-1);
  }
  return -1;
}

int bin_id(int etaB, int ptB, int npt){
  return (etaB * (npt-1)) + ptB;
}

/************************************************/
/************************************************/
/* Assumed that there is no out-of-range mll  */

void bin_Zcand_M(Input &ip,
            float Zcand_M,
            int bid1, int bid2,
            float lZcand_M, float rZcand_M,
            matrix &central,
            matrix &left,
            matrix &right)
{

  float weight = 1;
  if(!m_DT){ 
     weight *= ip.MCEvtWeight;
    if (applyPRW) weight *= ip.MCPileupWeight;
  }
  
  if ((Zcand_M > lZcand_M) && (Zcand_M < rZcand_M))
    central[bid1][bid2] += weight;
  else if (Zcand_M < lZcand_M)
    left[bid1][bid2] += weight;
  else if (Zcand_M > rZcand_M)
    right[bid1][bid2] += weight;
  else return;
}

/************************************************/
/************************************************/

matrix subtract_background(float bl,
                           float br,
                           matrix central,
                           matrix left,
                           matrix right)
{
  
  matrix M = m_sub(central, s_mul(1.0 / (bl+br),
                                  m_add(s_mul(bl, left),
                                        s_mul(br, right))));
  return negative_to_zero(M);
}

/***************************************************************************/
/***************************************************************************/


void print_setting(float lZcand_M,
                   float rZcand_M,
                   string tag,
                   bool bkg=false,
                   float bl=0.0,
                   float br=0.0)
{
  cout << "========== " << tag << " ==========" << endl;
  cout << "Setting -- LMLL: " << lZcand_M << " -- RMLL: " << rZcand_M;
  if (bkg)
    cout << " -- BL: " << bl << " -- BR: " << br << endl;
  else
    cout << endl;
}

string file_format(string name,
                   float lZcand_M, float rZcand_M,
                   float bl, float br,
                   string tag,
                   string extension)
{
  
  string
  slZcand_M = to_string(lZcand_M),
  srZcand_M = to_string(rZcand_M),
  sbl = to_string(bl),
  sbr = to_string(br);
  
  string s = "_";
  return
  name + s + slZcand_M + s + srZcand_M + s + sbl + s + sbr +
  "_" + tag + "." + extension;
}

// ofn is outfile name, "x.txt"
void write_matrix(string ofn, matrix M){
  ofstream text;
  text.open(ofn);
  text << std::fixed << std::setprecision(1);
  for (size_t i(0); i < M.size(); i++){
    text << M[i][0];
    for (size_t j(1); j < M[i].size(); j++){
      text << "," << M[i][j];
    }
    text << endl;
  }
}

// Print the n and nss matrices and write them into csv files
void print_write_info(string header, string info, matrix central, string name){
  cout << endl << endl;
  cout << "------------------------------------" << endl;
  cout << "------------------------------------" << endl;
  cout << header << endl;
  cout << "----------" << endl;
  write_matrix(name, central);
  cout << info << " has been written to file " << name << endl;
  cout << "----------" << endl;
}

void write_to_file(string outFile,
                   string n_name,
                   string nss_name,
                   string des,
                   string tag,
                   matrix n,
                   matrix nss,
                   float lZcand_M,
                   float rZcand_M,
                   float bl=0.0,
                   float br=0.0) 
{
  ofstream binFile;
  binFile.open(outFile, ios_base::app);
  
  string n_n = file_format(n_name, lZcand_M, rZcand_M, bl, br, tag, "csv");
  binFile << n_n << endl;
  print_write_info(des,
                   "N",
                   n,
                   n_n);
  
  string nss_n = file_format(nss_name, lZcand_M, rZcand_M, bl, br, tag, "csv");
  binFile << nss_n << endl;
  print_write_info(des,
                   "NSS",
                   nss,
                   nss_n);
}


/***************************************************************************/
/***************************************************************************/
/*
 The output are csv files showing how N and NSS are distributed in the eta-pt
 bins.
 
 lmll: The left mll value
 rmll: The right mll value 
 bkg: if true, background substraction is done; if false, no background subtraction
      is done
 bl: the width of the left mll band when doing background subtraction
 br: the width of the right mlll band when doing background subtraction
  */

void binData(TChain &events,
             Input &ip,
             row &etas,
             row &pts,
             float lZcand_M,
             float rZcand_M,
             string tag,
             string outFile,
             bool bkg=false,
             float bl=0.0,
             float br=0.0)
{
  // Checks
  char checkName[200];
  sprintf(checkName, "checks_%d_%d_%d_%d.root", (int) lZcand_M, (int) rZcand_M, (int) bl, (int) br);
  TFile f(checkName, "recreate");
  std::vector<TH1D*> hEta;
  std::vector<TH1D*> hPt;
  std::vector<TH1D*> hMass;
  TH1I* hCutflow = new TH1I("hCutflow", "Number of events passed", 9, 0, 9);
  hCutflow->GetXaxis()->SetBinLabel(1, "Total");
  hCutflow->GetXaxis()->SetBinLabel(2, "Trig+GRL+2e");
  hCutflow->GetXaxis()->SetBinLabel(3, "Pass LooseBaseline");
  hCutflow->GetXaxis()->SetBinLabel(4, "Pass Signal");
  hCutflow->GetXaxis()->SetBinLabel(5, "Pass qID");
  hCutflow->GetXaxis()->SetBinLabel(6, "Zmass+SB");
  hCutflow->GetXaxis()->SetBinLabel(7, "Zmass");
  hCutflow->GetXaxis()->SetBinLabel(8, "Left SB");
  hCutflow->GetXaxis()->SetBinLabel(9, "Right SB");

  hEta.push_back(new TH1D("hEtaAll", "Eta distribution of all electrons", 200, -2.47, 2.47));
  hEta.push_back(new TH1D("hEtaPreselected", "Eta distribution of preselected electrons", 200, -2.47, 2.47));
  hEta.push_back(new TH1D("hEtaLoose", "Eta distribution of electrons in LooseBaseline pairs", 200, -2.47, 2.47));
  hEta.push_back(new TH1D("hEtaSignal", "Eta distribution of electrons in Signal pairs", 200, -2.47, 2.47));
  hEta.push_back(new TH1D("hEtaSignalQID", "Eta distribution of electrons in Signal pairs passing QID", 200, -2.47, 2.47));
  hEta.push_back(new TH1D("hEtaSignalZSB", "Eta distribution of electrons in Signal pairs within Zmass+SB window", 200, -2.47, 2.47));
  hEta.push_back(new TH1D("hEtaSignalZ", "Eta distribution of electrons in Signal pairs within Zmass window", 200, -2.47, 2.47));
  hEta.push_back(new TH1D("hEtaSignalLSB", "Eta distribution of electrons in Signal pairs within left SB", 200, -2.47, 2.47));
  hEta.push_back(new TH1D("hEtaSignalRSB", "Eta distribution of electrons in Signal pairs within right SB", 200, -2.47, 2.47));

  hPt.push_back(new TH1D("hPtAll", "Pt distribution of all electrons", 200, 20, 200));
  hPt.push_back(new TH1D("hPtPreselected", "Pt distribution of preselected electrons", 200, 20, 200));
  hPt.push_back(new TH1D("hPtLoose", "Pt distribution of electrons in LooseBaseline pairs", 200, 20, 200));
  hPt.push_back(new TH1D("hPtSignal", "Pt distribution of electrons in Signal pairs", 200, 20, 200));
  hPt.push_back(new TH1D("hPtSignalQID", "Pt distribution of electrons in Signal pairs passing QID", 200, 20, 200));
  hPt.push_back(new TH1D("hPtSignalZSB", "Pt distribution of electrons in Signal pairs within Zmass+SB window", 200, 20, 200));
  hPt.push_back(new TH1D("hPtSignalZ", "Pt distribution of electrons in Signal pairs within Zmass window", 200, 20, 200));
  hPt.push_back(new TH1D("hPtSignalLSB", "Pt distribution of electrons in Signal pairs within left SB", 200, 20, 200));
  hPt.push_back(new TH1D("hPtSignalRSB", "Pt distribution of electrons in Signal pairs within right SB", 200, 20, 200));

  hMass.push_back(new TH1D("hMassAll", "Mass distribution of all electrons", 200, 20, 200));
  hMass.push_back(new TH1D("hMassPreselected", "Mass distribution of preselected electrons", 200, 20, 200));
  hMass.push_back(new TH1D("hMassLoose", "Mass distribution of electrons in LooseBaseline pairs", 200, 20, 200));
  hMass.push_back(new TH1D("hMassSignal", "Mass distribution of electrons in Signal pairs", 200, 20, 200));
  hMass.push_back(new TH1D("hMassSignalQID", "Mass distribution of electrons in Signal pairs passing QID", 200, 20, 200));
  hMass.push_back(new TH1D("hMassSignalZSB", "Mass distribution of electrons in Signal pairs within Zmass+SB window", 200, 20, 200));
  hMass.push_back(new TH1D("hMassSignalZ", "Mass distribution of electrons in Signal pairs within Zmass window", 200, 20, 200));
  hMass.push_back(new TH1D("hMassSignalLSB", "Mass distribution of electrons in Signal pairs within left SB", 200, 20, 200));
  hMass.push_back(new TH1D("hMassSignalRSB", "Mass distribution of electrons in Signal pairs within right SB", 200, 20, 200));
  int NETA = etas.size(), NPT = pts.size();
  int SIZE = (NETA - 1) * (NPT - 1);
  
  float MINETA = row_min(etas), MAXETA = row_max(etas);
  float MINPT = row_min(pts), MAXPT = row_max(pts);
  
  matrix nC(SIZE, row(SIZE, 0.0));
  matrix nL(SIZE, row(SIZE, 0.0));
  matrix nR(SIZE, row(SIZE, 0.0));
  
  matrix nssC(SIZE, row(SIZE, 0.0));
  matrix nssL(SIZE, row(SIZE, 0.0));
  matrix nssR(SIZE, row(SIZE, 0.0));
  
  float elCand1_eta, elCand1_pt, elCand2_eta, elCand2_pt, Zcand_M;
  int bid1, bid2;
  
  cout << endl << endl;
  cout << "*****************************************************************************" << endl;
  cout << "*****************************************************************************" << endl;
  cout << "*****************************************************************************" << endl;
  cout << "RUNNING... " << endl << endl << endl;
  
  #define CUTFLOW(i)                                                    \
  hCutflow->Fill(i);                                                    \
  hEta[i]->Fill(ip.elCand1_cl_eta); hEta[i]->Fill(ip.elCand2_cl_eta);   \
  hPt[i]->Fill(elCand1_pt); hPt[i]->Fill(elCand2_pt);                   \
  hMass[i]->Fill(Zcand_M);   

  int nEntries(events.GetEntries());
  for (int i(0); i < nEntries; i++) {
    loadbar(i+1, nEntries);
    events.GetEntry(i);
    elCand1_eta = fabs(ip.elCand1_cl_eta);
    elCand2_eta = fabs(ip.elCand2_cl_eta);
    elCand1_pt = ip.elCand1_pt / PTSCALE;
    elCand2_pt = ip.elCand2_pt / PTSCALE;
    Zcand_M = ip.Zcand_M / PTSCALE;

    // My cuts =========================================================================================
    // Signal electrons.

    if (in==EGAMMA){
      // LooseBaseline
      if (elCand1_pt < 20 || elCand2_pt < 20) continue;
      if (elCand1_eta > 2.47 || elCand2_eta > 2.47) continue;

      if(!(ip.elCand1_isLooseAndBLayerLL2015_v8 && ip.elCand2_isLooseAndBLayerLL2015_v8)) continue;

      if (selection == SIGNAL){
        if (elCand1_pt > 300 && !ip.elCand1_isMediumLL2015_v8) continue;
        if (elCand1_pt <= 300 && !ip.elCand1_isTightLL2015_v8) continue;

        if (elCand2_pt > 300 && !ip.elCand2_isMediumLL2015_v8) continue;
        if (elCand2_pt <= 300 && !ip.elCand2_isTightLL2015_v8) continue;

        if(!ip.elCand1_isolGradientLoose) continue;
        if(!ip.elCand2_isolGradientLoose) continue;

        // |d0/sig_d0| < 5, |z0sin| < 0.5
        if(!ip.elCand1_passd0z0 || !ip.elCand2_passd0z0) continue;
      }
    }

    if (in==SUSY){
      CUTFLOW(0);

      if(ip.trigCode<=0) continue;

      // Select electrons only
      if (int(fabs(ip.elCand1_ID/1000))!=11 || int(fabs(ip.elCand2_ID/1000))!=11) continue;

      CUTFLOW(1);

      // Select looseBaseline electrons
      // if(!(ip.elCand1_flag & 1<<0)) continue;
      // if(!(ip.elCand2_flag & 1<<0)) continue;

      CUTFLOW(2);

      // if(selection == SIGNAL){ // Select signal electrons
      //   if(!(ip.elCand1_flag & 1 << 1)) continue;
      //   if(!(ip.elCand2_flag & 1 << 1)) continue;
      //   CUTFLOW(3);
      // }
      if(selection == SIGNAL){
        if(!((ip.elCand1_flag & 2)/2)) continue;
        if(!((ip.elCand2_flag & 2)/2)) continue;
      }
      CUTFLOW(3);

      if(passQID && !(ip.elCand1_qID && ip.elCand2_qID)) continue;
      CUTFLOW(4);
    }

    //=========================================================================================

    if (is_out_Zcand_M(Zcand_M, lZcand_M, rZcand_M, bl, br)) continue;

    CUTFLOW(5);

    if (is_out_eta_pt(elCand1_eta, elCand1_pt, elCand2_eta, elCand2_pt, MINETA, MAXETA, MINPT, MAXPT))
      continue;

    if(Zcand_M>=lZcand_M && Zcand_M<=rZcand_M) CUTFLOW(6);
    if(Zcand_M>=lZcand_M-bl && Zcand_M<lZcand_M) CUTFLOW(7);
    if(Zcand_M>rZcand_M && Zcand_M<=rZcand_M+br) CUTFLOW(8);

    bid1 = bin_id(eta_bin(elCand1_eta), pt_bin(elCand1_pt), NPT);
    bid2 = bin_id(eta_bin(elCand2_eta), pt_bin(elCand2_pt), NPT);
       
    bin_Zcand_M(ip, Zcand_M, bid1, bid2, lZcand_M, rZcand_M, nC, nL, nR);
    if (ip.elCand1_charge == ip.elCand2_charge)
      bin_Zcand_M(ip, Zcand_M, bid1, bid2, lZcand_M, rZcand_M, nssC, nssL, nssR);
  }

  #undef CUTFLOW

  cout << endl << endl << "=======================================> RESULTS: " << endl << endl;
  print_setting(lZcand_M, rZcand_M, "N and NSS", bkg, bl, br);
  
  if (!bkg){
    write_to_file(outFile, "n_nbg", "nss_nbg", "NO BKG SUBTRACTION", tag,
                  nC, nssC, lZcand_M, rZcand_M, bl, br);
  } else {
    matrix n_s = subtract_background(bl, br, nC, nL, nR);
    matrix nss_s = subtract_background(bl, br, nssC, nssL, nssR);

    write_to_file(outFile, "n_bg", "nss_bg", "BKG SUBTRACTION", tag,
                  n_s, nss_s, lZcand_M, rZcand_M, bl, br);
    
  }

  for(auto h : hPt) h->Write();
  for(auto h : hEta) h->Write();
  for(auto h : hMass) h->Write();
  hCutflow->Write();
  f.Close();
}

/***************************************************************************/
/***************************************************************************/
/*
 The following functions allow convenient binning 
*/


/* Binning with background subtraction on and off */
void binDataBkgOnOff(TChain &events,
                     Input &ip,
                     row &etas,
                     row &pts,
                     float lZcand_M,
                     float rZcand_M,
                     float bl,
                     float br,
                     string tag,
                     string outFile) 
{
  binData(events, ip, etas, pts, lZcand_M, rZcand_M, tag, outFile, false);
  binData(events, ip, etas, pts, lZcand_M, rZcand_M, tag, outFile, true, bl, br);
}

/* Binning where the width of the sideband may be varied 

  Typically, when doing a sideband variation, the width should be 
  chosen at least two times, one smaller than the original width,
  one greater than the orignial width. In the following, the prefix
  i (as in ibl) means increased, and d (as in dbl) means decreased. 
  
*/


/*
  ibl: width of the left band 
  ibr: width of the right band
  
  If another variation of the width of the sideband is needed,
  call the function with more set to true, and supply:
  
  dbl: width of the left band
  dbr: width of the right band
  
   
*/
void binDataSideBandVar(TChain &events,
                        Input &ip,
                        row &etas,
                        row &pts,
                        float lZcand_M,
                        float rZcand_M,
                        float bl,
                        float br,
                        float ibl,
                        float ibr,
                        string tag,
                        string outFile,
                        bool more=true,
                        float dbl=0.0,
                        float dbr=0.0)
{
  binData(events, ip, etas, pts, lZcand_M, rZcand_M, tag, outFile, true, bl, br);
  binData(events, ip, etas, pts, lZcand_M, rZcand_M, tag, outFile, true, ibl, ibr);
  if (more)
    binData(events, ip, etas, pts, lZcand_M, rZcand_M, tag, outFile, true, dbl, dbr);
  
}



/*
  Binning where the central mll band is varied. 
  
  nlmll: the new value of the left side mll
  nrmll: the new value of the right side mll

*/

void binDataCentralBandVar(TChain &events,
                           Input &ip,
                           row &etas,
                           row &pts,
                           float lZcand_M,
                           float rZcand_M,
                           float nlZcand_M,
                           float nrZcand_M,
                           string tag,
                           string outFile,
                           bool bkg=false,
                           float bl=0.0,
                           float br=0.0)
{
  binData(events, ip, etas, pts, lZcand_M, rZcand_M, tag, outFile, bkg, bl, br);
  binData(events, ip, etas, pts, nlZcand_M, nrZcand_M, tag, outFile, bkg, bl, br);
  
}

/*
  Binning where all of the following can be done at once (refer to the functions
  above)
  
    background on/off
    central band variation
    sideband variation   

*/

void binDataBkgOOSBVar(TChain &events,
                       Input &ip,
                       row &etas,
                       row &pts,
                       float lZcand_M,
                       float rZcand_M,
                       float bl,
                       float br,
                       float ibl,
                       float ibr,
                       string tag,
                       string outFile,
                       bool more=false,
                       float dbl=0.0,
                       float dbr=0.0)
{
  // Background on/off
  binDataBkgOnOff(events, ip, etas, pts, lZcand_M, rZcand_M, bl, br, tag, outFile);
  // Sideband variation
  binData(events, ip, etas, pts, lZcand_M, rZcand_M, tag, outFile, true, ibl, ibr);
  if (more)
    binData(events, ip, etas, pts, lZcand_M, rZcand_M, tag, outFile, true, dbl, dbr);
}



void binDataBkgOOSBVarCBVar(TChain &events,
                            Input &ip,
                            row &etas,
                            row &pts,
                            float lZcand_M,
                            float rZcand_M,
                            float bl,
                            float br,
                            float nlZcand_M,
                            float nrZcand_M,
                            float ibl,
                            float ibr,
                            string tag,
                            string outFile,
                            bool more=false,
                            float dbl=0.0,
                            float dbr=0.0)
{
  // Background on/off
  binDataBkgOnOff(events, ip, etas, pts, lZcand_M, rZcand_M, bl, br, tag, outFile);
  // Sideband variation
  binData(events, ip, etas, pts, lZcand_M, rZcand_M, tag, outFile, true, ibl, ibr);
  if (more)
    binData(events, ip, etas, pts, lZcand_M, rZcand_M, tag, outFile, true, dbl, dbr);
  // Central band variation
  binData(events, ip, etas, pts, nlZcand_M, nrZcand_M, tag, outFile, true, bl, br);
  
}

void testPlots(TChain &events, Input &ip){
  TH1D* ptOne = new TH1D("pt", "pt", 200, 0, 200);
  TH1D* ptTwo = new TH1D("ptTwo", "pt", 200, 0, 200);
  for (int i(0); i < events.GetEntries(); i++){
    ptOne->Fill(ip.elCand1_pt);
    ptTwo->Fill(ip.elCand2_pt);
  }

  TCanvas *c = new TCanvas("can", "canvas");
  ptOne->Draw();
  ptTwo->Draw("same");
  c->Print("test.pdf");
  return;

}

void binMC(TChain &events, Input &ip, row &etas, row &pts, float lZcand_M, float rZcand_M)
{
  TH2* hTruth = new TH2D("hMCTruthRate", "MC Truth misID rate", etas.size()-1, &etas[0], pts.size()-1, &pts[0]);
  TH2* hN = (TH2*) hTruth->Clone("hMC_N");
  TH2* hNflipped = (TH2*) hTruth->Clone("hMC_Nflipped");


  char fName[100];
  sprintf(fName, "MCTruth_%d_%d.root", (int) lZcand_M, (int) rZcand_M);
  TFile* fOut = new TFile(fName, "recreate");
  hTruth->SetDirectory(fOut);
  hN->SetDirectory(fOut);
  hNflipped->SetDirectory(fOut);


  float elCand1_eta, elCand1_pt, elCand2_eta, elCand2_pt, Zcand_M;


  int nEntries(events.GetEntries());
  for(int i(0); i< nEntries; i++)
  {
    loadbar(i+1, nEntries);
    events.GetEntry(i);

    elCand1_eta = fabs(ip.elCand1_cl_eta);
    elCand2_eta = fabs(ip.elCand2_cl_eta);
    elCand1_pt = ip.elCand1_pt / PTSCALE;
    elCand2_pt = ip.elCand2_pt / PTSCALE;
    Zcand_M = ip.Zcand_M / PTSCALE;

    if(ip.trigCode<=0) continue;
    if (int(fabs(ip.elCand1_ID/1000))!=11 || int(fabs(ip.elCand2_ID/1000))!=11) continue;
    if(!((ip.elCand1_flag & 2)/2)) continue;
    if(!((ip.elCand2_flag & 2)/2)) continue;

    if (is_out_Zcand_M(Zcand_M, lZcand_M, rZcand_M, 0, 0)) continue;

    if (is_out_eta_pt(elCand1_eta, elCand1_pt, elCand2_eta, elCand2_pt, -2.47, 2.47, 20, 1000))
      continue;
      
    if(passQID && !(ip.elCand1_qID && ip.elCand2_qID)) continue;
    double weight = 1;
    weight *= ip.MCEvtWeight;
    if(applyPRW) weight *= ip.MCPileupWeight;

    if(ip.elCand1_origCharge!=0)
    {
      hN->Fill(elCand1_eta, elCand1_pt, weight);
      if(ip.elCand1_origCharge != ip.elCand1_charge) hNflipped->Fill(elCand1_eta, elCand1_pt, weight);
    }
    
    if(ip.elCand2_origCharge!=0)
    {
      hN->Fill(elCand2_eta, elCand2_pt, weight);
      if(ip.elCand2_origCharge != ip.elCand2_charge) hNflipped->Fill(elCand2_eta, elCand2_pt, weight);
    }

  }

  hTruth->Sumw2();
  hTruth->Add(hNflipped);
  hTruth->Divide(hN);

  fOut->Write();
  fOut->Close();
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/* g++ -O3 -Wall -Wextra -std=c++11 -o flip_rates flip_rates.cxx `root-config --cflags --glibs` */
/* 
   The main function requires three arguments:
   1. The name of the tree (example: eventTreee)
   2. Name of the output file (example: bins.txt)
      Note that this file is written in append mode 
   3. Name of the data sets, separated by commas (example: 1.root,2.root,3.root)
*/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/



int main(int argc, char *argv[]){
  if (argc!=3 && argc!=4) {
    cout << "Wrong number of arguments. " << argc-1 << "provided, 2 or 3 required." << endl;
    return 1;
  }

  if (m_DT) TAG = "DATA";
  else TAG = "MC";


  TChain events(argv[1]);
  string outFile = string(argv[2]);
  string dataName = string(argv[3]);
  
  istringstream ss(dataName);
  string token;
  while(getline(ss, token, ',')) {
    cout << endl << endl;
    cout << "*******************************************" << endl;
    cout << "Reading in data ..." << endl;
    cout << token << endl;
    events.Add(token.c_str());
    cout << "*******************************************" << endl;
  }

  Input ip;
  connect_input(ip, events);

  // testPlots(events, ip);
  // return 0;
  
  // binDataBkgOOSBVarCBVar(events, ip, VETAS, VPTS, 75, 100, 25, 25, 80, 110, 20, 20, "dt", outFile, true, 30, 30);
  
  // binDataCentralBandVar(events, ip, VETAS, VPTS, 75, 100, 80, 110, "mc_truth", outFile);
  // binDataCentralBandVar(events, ip, VETAS, VPTS, 75, 100, 80, 110, "mc", outFile);
  
  
  /////////////////////////////
  //binDataBkgOOSBVar(events, ip, VETAS, VPTS, 75, 100, 25, 25, 20, 20, TAG, outFile, true, 30, 30);
  // binDataBkgOOSBVar(events, ip, VETAS, VPTS, 70, 110, 25, 25, 20, 20, "dt", outFile, true, 30, 30);
  // binDataBkgOOSBVar(events, ip, VETAS, VPTS, 75, 105, 25, 25, 20, 20, "dt", outFile, true, 30, 30);
  // binDataBkgOOSBVar(events, ip, VETAS, VPTS, 80, 100, 25, 25, 20, 20, "dt", outFile, true, 30, 30);
  
  //binDataCentralBandVar(events, ip, VETAS, VPTS, 75, 100, 70, 110, "mc", outFile);
  //binDataCentralBandVar(events, ip, VETAS, VPTS, 75, 105, 80, 100, "mc", outFile);
  /////////////////////////////
  
  
  // if (m_DT) binData(events, ip, VETAS, VPTS, 80, 100, TAG, outFile);//, true, 20, 20);
  // binData(events, ip, VETAS, VPTS, 75, 100, TAG, outFile);

  //if (m_DT) binDataCentralBandVar(events, ip, VETAS, VPTS, 80, 110, 85, 105, TAG, outFile, true, 20, 20);
  //binDataCentralBandVar(events, ip, VETAS, VPTS, 80, 110, 85, 105, TAG, outFile);

  // binData(events, ip, VETAS, VPTS, 80, 100, TAG, outFile);

  // if (m_DT) binDataBkgOnOff(events, ip, VETAS, VPTS, 80, 100, 20, 20, TAG, outFile);
  // else binData(events, ip, VETAS, VPTS, 80, 100, TAG, outFile);

  if (m_DT)
  {
    // Nominal 80-100, 20
    binData(events, ip, VETAS, VPTS, 80, 100, TAG, outFile, true, 20, 20);

    // No SB
    binData(events, ip, VETAS, VPTS, 80, 100, TAG, outFile);

    // Wide SB 25 GeV
    binData(events, ip, VETAS, VPTS, 80, 100, TAG, outFile, true, 25, 25);

    // Narrow SB
    binData(events, ip, VETAS, VPTS, 80, 100, TAG, outFile, true, 15, 15);

    // Down-shifted mass window
    binData(events, ip, VETAS, VPTS, 75, 100, TAG, outFile, true, 20, 20);

    // Wide mass window
    binData(events, ip, VETAS, VPTS, 75, 105, TAG, outFile, true, 20, 20);

  }
  else
  { 
    binData(events, ip, VETAS, VPTS, 80, 100, TAG, outFile, false, 0, 0);
    binMC(events, ip, VETAS, VPTS, 80, 100);
  }

  return 0;
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

