#define TagProbePlot_cxx

#include "TagProbePlot.h"
//#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TText.h>
#include <TF1.h>
#include <TStyle.h>
#include <TPie.h>

// Define Pt bins. Remember to update nPtBins in .h accordingly
static const Double_t tPtEdges[] = {20, 30, 40, 50, 60, 80, 120}; // 6 bins
static const Double_t tPtEdges1[] = {20, 30, 35, 40, 45, 50, 55, 60, 80, 120}; // 9 bins, for finer binning in 30 -60
const Double_t* TagProbePlot::PtEdges = tPtEdges1;

// Define Eta bins. Remember to update nEtaBins in .h accordingly
static const Double_t tEtaEdges1[] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.37}; // 6 bins, for testing
static const Double_t tEtaEdges2[] = {1.5, 1.6, 1.7, 1.8}; //3 bins, for testing
static const Double_t tEtaEdges3[] = {-2.5, -2, -1.8, -1.52, -1.37, -1, -0.5, 0, 0.5, 1, 1.37, 1.52, 1.8, 2, 2.5}; // 14 bins, for testing eta (vs. |eta|)
static const Double_t tEtaEdges[] = {0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.5}; // 7 bins
const Double_t* TagProbePlot::EtaEdges = tEtaEdges;

static inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
    if ( (x != n) && (x % (n/100+1) != 0) ) return;
 
    float ratio    =  x/(float)n;
    unsigned int c =  ratio * w;
 
    cout << setw(3) << (int)(ratio*100) << "% [";
    for (unsigned int x=0; x<c; x++) cout << "=";
    for (unsigned int x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
}


//Opens
void TagProbePlot::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   //TString option = GetOption();
   if (toFile) outputFile = new TFile("output.root", "RECREATE");
}


void TagProbePlot::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   //TString option = GetOption();

}


Bool_t TagProbePlot::Process(Long64_t entry)
{

   fChain->GetEntry(entry);
   nEvents++;
   loadbar(nEvents, nTotEvents);
      
   //=========================================================//
   // Run on all events (MC and data) that pass the event cut //
   //=========================================================//
   //*
   if(PassedEventCut()){ 
      ProcessData();
      
      //===============================================//
      // Comparisons run on saved MC truth information //
      //===============================================//
      //*
      if(fMCevent){
         ProcessMC();
      }
   }
   
  
   //*/
   return kTRUE;
}

// Does nothing
void TagProbePlot::SlaveTerminate()
{

}

// Draws histograms
void TagProbePlot::Terminate()
{
   std::cout << "Terminate" << std::endl;
   hPt[DEFAULT][ERROR]->SetMinimum(0.001);
   hPt[DEFAULT][ERROR]->SetMaximum(0.03);
   
   hEta[DEFAULT][ERROR]->SetMinimum(0.0001);
   hEta[DEFAULT][ERROR]->SetMaximum(0.03);
   // Draw eta(Pt) histograms 
   //DivideHist(hPtError, hPt, hPtnSSProbes);
   //hPtError->SetLineColor(kBlue);
   
   //TF1 funcTwo("two", "2"); // Divide by two?

   // Print only name, no stats on histogram
   gStyle->SetOptStat(0);
   
   // Prep eta, pt, invMass, dPhi, d0Sig distributions
   DivideHist((TH1**)hPtMass[TIGHT]);
   DivideHist((TH1**)hEtaMass[TIGHT]);
   DivideHist((TH1**)hPtMassNoZ[TIGHT]);
   DivideHist((TH1**)hEtaMassNoZ[TIGHT]);
   
   // Prep tag-and-probe histograms
   DivideHist(hPt[LOOSE]);
   //hPtErrorLoose->Divide(&funcTwo, 1);
   hPt[LOOSE][ERROR]->SetMarkerStyle(20);
   hPt[LOOSE][ERROR]->SetMarkerColor(kBlue);
   hPt[LOOSE][ERROR]->SetLineColor(kBlue);
   
   DivideHist(hPt[MEDIUM]);
   //hPtErrorMedium->Divide(&funcTwo, 1);
   hPt[MEDIUM][ERROR]->SetMarkerStyle(21);
   hPt[MEDIUM][ERROR]->SetMarkerColor(kBlue);
   hPt[MEDIUM][ERROR]->SetLineColor(kBlue);
   
   DivideHist(hPt[TIGHT]);
   //hPtErrorTight->Divide(&funcTwo, 1);
   hPt[TIGHT][ERROR]->SetMarkerStyle(22);
   hPt[TIGHT][ERROR]->SetMarkerColor(kBlue);
   hPt[TIGHT][ERROR]->SetLineColor(kBlue);
   
   DivideHist(hMass[MEDIUM]);
   
   DivideHist(hEta[LOOSE]);
   //hEta[LOOSE][ERROR]->Divide(&funcTwo, 1);
   hEta[LOOSE][ERROR]->SetMarkerStyle(20);
   hEta[LOOSE][ERROR]->SetMarkerColor(kBlue);
   hEta[LOOSE][ERROR]->SetLineColor(kBlue);
   
   DivideHist(hEta[MEDIUM]);
   //hEta[MEDIUM][ERROR]->Divide(&funcTwo, 1);
   hEta[MEDIUM][ERROR]->SetMarkerStyle(21);
   hEta[MEDIUM][ERROR]->SetMarkerColor(kBlue);
   hEta[MEDIUM][ERROR]->SetLineColor(kBlue);
   
   DivideHist(hEta[TIGHT]);
   //hEta[TIGHT][ERROR]->Divide(&funcTwo, 1);
   hEta[TIGHT][ERROR]->SetMarkerStyle(22);
   hEta[TIGHT][ERROR]->SetMarkerColor(kBlue);
   hEta[TIGHT][ERROR]->SetLineColor(kBlue);

   /*
   for(int i=0; i<=nPtBins; i++){
      DivideHist(hEtaPtLoose[i]);
      DivideHist(hEtaPtMedium[i]);
      DivideHist(hEtaPtTight[i]);
      
      hEtaPtLoose[i][ERROR]->SetMarkerStyle(22);
      hEtaPtLoose[i][ERROR]->SetMarkerColor(kBlue);
      hEtaPtLoose[i][ERROR]->SetLineColor(kBlue);
      hEtaPtMedium[i][ERROR]->SetMarkerStyle(22);
      hEtaPtMedium[i][ERROR]->SetMarkerColor(kBlue);
      hEtaPtMedium[i][ERROR]->SetLineColor(kBlue);
      hEtaPtTight[i][ERROR]->SetMarkerStyle(22);
      hEtaPtTight[i][ERROR]->SetMarkerColor(kBlue);
      hEtaPtTight[i][ERROR]->SetLineColor(kBlue);
   }
   //*/
   
   PRINT_STAT("pt,eta,mass", "done");
   
   DivideHist((TH1**) hEP[LOOSE]);
   DivideHist((TH1**) hEP[MEDIUM]);
   DivideHist((TH1**) hEP[TIGHT]); 
   
   PRINT_STAT("hep", "done");  
   
   // Prep MC histograms
   DivideHist(hPt_MC[DEFAULT]);
   hPt_MC[DEFAULT][ERROR]->SetMarkerColor(kRed);
   hPt_MC[DEFAULT][ERROR]->SetLineColor(kRed);

   DivideHist(hPt_MC[LOOSE]);
   hPt_MC[LOOSE][ERROR]->SetMarkerStyle(20);
   hPt_MC[LOOSE][ERROR]->SetMarkerColor(kRed);
   hPt_MC[LOOSE][ERROR]->SetLineColor(kRed);

   DivideHist(hPt_MC[MEDIUM]);
   hPt_MC[MEDIUM][ERROR]->SetMarkerStyle(21);
   hPt_MC[MEDIUM][ERROR]->SetMarkerColor(kRed);
   hPt_MC[MEDIUM][ERROR]->SetLineColor(kRed);

   DivideHist(hPt_MC[TIGHT]);
   hPt_MC[TIGHT][ERROR]->SetMarkerStyle(22);
   hPt_MC[TIGHT][ERROR]->SetMarkerColor(kRed);
   hPt_MC[TIGHT][ERROR]->SetLineColor(kRed);

   DivideHist(hMass_MC[LOOSE]);
   hMass_MC[LOOSE][ERROR]->SetMarkerColor(kRed);
   hMass_MC[LOOSE][ERROR]->SetLineColor(kRed);
   
   DivideHist(hMass_MC[MEDIUM]);
   hMass_MC[MEDIUM][ERROR]->SetMarkerColor(kRed);
   hMass_MC[MEDIUM][ERROR]->SetLineColor(kRed);

   DivideHist(hMass_MC[TIGHT]);
   hMass_MC[TIGHT][ERROR]->SetMarkerColor(kRed);
   hMass_MC[TIGHT][ERROR]->SetLineColor(kRed);
   
   DivideHist(hEta_MC[DEFAULT]);
   //hEta_MC[DEFAULT][ERROR]->SetMarkerStyle(20);
   hEta_MC[DEFAULT][ERROR]->SetMarkerColor(kRed);
   hEta_MC[DEFAULT][ERROR]->SetLineColor(kRed);
   
   DivideHist(hEta_MC[LOOSE]);
   hEta_MC[LOOSE][ERROR]->SetMarkerStyle(20);
   hEta_MC[LOOSE][ERROR]->SetMarkerColor(kRed);
   hEta_MC[LOOSE][ERROR]->SetLineColor(kRed);
   
   DivideHist(hEta_MC[MEDIUM]);
   hEta_MC[MEDIUM][ERROR]->SetMarkerStyle(21);
   hEta_MC[MEDIUM][ERROR]->SetMarkerColor(kRed);
   hEta_MC[MEDIUM][ERROR]->SetLineColor(kRed);
   
   DivideHist(hEta_MC[TIGHT]);
   hEta_MC[TIGHT][ERROR]->SetMarkerStyle(22);
   hEta_MC[TIGHT][ERROR]->SetMarkerColor(kRed);
   hEta_MC[TIGHT][ERROR]->SetLineColor(kRed);
   
   DivideHist((TH1**) hEP_MC[LOOSE]);
   hEP_MC[LOOSE][ERROR]->SetMarkerStyle(20);
   hEP_MC[LOOSE][ERROR]->SetMarkerColor(kRed);
   hEP_MC[LOOSE][ERROR]->SetLineColor(kRed);
   
   DivideHist((TH1**) hEP_MC[MEDIUM]);
   hEP_MC[MEDIUM][ERROR]->SetMarkerStyle(21);
   hEP_MC[MEDIUM][ERROR]->SetMarkerColor(kRed);
   hEP_MC[MEDIUM][ERROR]->SetLineColor(kRed);
   
   DivideHist((TH1**) hEP_MC[TIGHT]);
   hEP_MC[TIGHT][ERROR]->SetMarkerStyle(22);
   hEP_MC[TIGHT][ERROR]->SetMarkerColor(kRed);
   hEP_MC[TIGHT][ERROR]->SetLineColor(kRed);
   
   /*
   for(int i=0; i<nPtBins; i++){
      DivideHist(hEtaPt_MC[i]);
      hEtaPt_MC[i][ERROR]->SetMarkerStyle(20);
      hEtaPt_MC[i][ERROR]->SetMarkerColor(kRed);
      hEtaPt_MC[i][ERROR]->SetLineColor(kRed);
   }
   //*/

   std::cout << "whoops" << std::endl;
   // == FOR MISID(PT) == //
   if(drawPt){
      TCanvas *cPT = new TCanvas("cPT", "Comparison of misID rate (Pt) from TP and MC truth");
      cPT->Divide(1,2);
      TPad* cPTBottom= (TPad*) cPT->cd(2);
      cPTBottom->Divide(2, 1);
   
      //TCanvas *cPrint = new TCanvas();
   
      /// For loose cuts
      //TCanvas *cLoose = new TCanvas("errorPTLoose", "Comparison of error rate (Pt) from tag-and-probe and MC truth (loose cuts)");
      cPTBottom->cd(2)->SetLogy();
      // hPt[DEFAULT][ERROR]->Draw();
      hPt[LOOSE][ERROR]->SetMaximum(hPt_MC[LOOSE][ERROR]->GetMaximum()*2);
      hPt[LOOSE][ERROR]->SetMinimum(hPt_MC[LOOSE][ERROR]->GetMinimum()/2);
      hPt[LOOSE][ERROR]->Draw("E");
      hPt_MC[LOOSE][ERROR]->Draw("E SAME");

      // Legend for error(pt)
      TLegend *lLoose = new TLegend(0.25,0.75,0.45,0.88);
      lLoose->SetHeader("For loose cuts:");
      lLoose->AddEntry(hPt[LOOSE][ERROR], "From tag-and-probe", "lpf");
      lLoose->AddEntry(hPt_MC[LOOSE][ERROR], "From MC truth", "lpf");
      lLoose->Draw();
   
      /// For medium cuts
      //TCanvas *cMedium = new TCanvas("errorPTMedium", "Comparison of error rate (Pt) from tag-and-probe and MC truth (medium cuts)");
      cPTBottom->cd(1)->SetLogy();
      // hPt[DEFAULT][ERROR]->Draw();
      hPt[MEDIUM][ERROR]->SetMaximum(hPt_MC[MEDIUM][ERROR]->GetMaximum()*2);
      hPt[MEDIUM][ERROR]->SetMinimum(hPt_MC[MEDIUM][ERROR]->GetMinimum()/2);
      hPt[MEDIUM][ERROR]->Draw("E");
      hPt_MC[MEDIUM][ERROR]->Draw("e SAME");

      // Legend for error(pt)
      TLegend *lMedium = new TLegend(0.25,0.75,0.45,0.88);
      lMedium->SetHeader("For medium cuts:");
      //lPT->AddEntry(hPtError, "From tag-and-probe", "lpf");
      lMedium->AddEntry(hPt[MEDIUM][ERROR], "From tag-and-probe", "lpf");
      lMedium->AddEntry(hPt_MC[MEDIUM][ERROR], "From MC truth", "lpf");
      lMedium->Draw();
 
      /// For tight cuts
      //TCanvas *cTight = new TCanvas("errorPTTight", "Comparison of error rate (Pt) from tag-and-probe and MC truth (tight cuts)");
      cPT->cd(1)->SetLogy();
      // hPt[DEFAULT][ERROR]->Draw();
      hPt[TIGHT][ERROR]->SetMaximum(hPt_MC[TIGHT][ERROR]->GetMaximum()*2);
      hPt[TIGHT][ERROR]->SetMinimum(hPt_MC[TIGHT][ERROR]->GetMinimum()/2);
      hPt[TIGHT][ERROR]->Draw("e SAME");
      hPt_MC[TIGHT][ERROR]->Draw("e SAME");

      // Legend for error(pt)
      TLegend *lTight = new TLegend(0.25,0.75,0.45,0.88);
      lTight->SetHeader("For tight cuts");
      //lPT->AddEntry(hPtError, "From tag-and-probe", "lpf");
      lTight->AddEntry(hPt[TIGHT][ERROR], "From tag-and-probe", "lpf");
      lTight->AddEntry(hPt_MC[TIGHT][ERROR], "From MC truth", "lpf");
      lTight->Draw();
      if(print) cPT->Print("Pt.pdf");
   }
   
   // === FOR MISID(MASS) === //
   /*
   TCanvas *cMassE = new TCanvas("cMass", "Comparison of error rates across the invariant mass");
   hMass[MEDIUM][ERROR]->Draw();
   hMass_MC[MEDIUM][ERROR]->Draw("e same");
   
   TLegend *lMass = new TLegend(0.25, 0.75, 0.45, 0.88);
   lMass->SetHeader("Mass");
   lMass->AddEntry(hMass[MEDIUM][ERROR], "From tag-and-probe (medium cuts)", "lpf");
   lMass->AddEntry(hMass_MC[MEDIUM][ERROR], "From MC truth", "lpf");
   lMass->Draw();
   if(print) cMassE->Print("Mass.pdf");
   //*/
   
   // ===== FOR MISID(ETA) ===== // 
   if(drawEta){
      TCanvas *cEta = new TCanvas("cEta", "Error rates across eta");
      cEta->Divide(1,2);
      TPad* cEtaBottom = (TPad*) cEta->cd(2);
      cEtaBottom->Divide(2,1);
    
       /// For loose cuts
      //TCanvas *cLoose = new TCanvas("errorPTLoose", "Comparison of error rate (Pt) from tag-and-probe and MC truth (loose cuts)");
      cEtaBottom->cd(2)->SetLogy();
      // hEta[DEFAULT][ERROR]->Draw();
      hEta[LOOSE][ERROR]->SetMaximum(hEta_MC[LOOSE][ERROR]->GetMaximum()*2);
      hEta[LOOSE][ERROR]->SetMinimum(hEta_MC[LOOSE][ERROR]->GetMinimum()/2);
      hEta[LOOSE][ERROR]->Draw("e");
      hEta_MC[LOOSE][ERROR]->Draw("e SAME");
      if(ABS_ETA) hEta[LOOSE][ERROR]->GetXaxis()->SetTitle("|#eta|");
   
      // Legend for error(Eta)
      TLegend *lEtaLoose = new TLegend(0.25,0.75,0.45,0.88);
      lEtaLoose->SetHeader("For loose cuts:");
      lEtaLoose->AddEntry(hEta[LOOSE][ERROR], "From tag-and-probe", "lpf");
      lEtaLoose->AddEntry(hEta_MC[LOOSE][ERROR], "From MC truth", "lpf");
      lEtaLoose->Draw();

      /// For medium cuts
      //TCanvas *cMedium = new TCanvas("errorEtaMedium", "Comparison of error rate (Eta) from tag-and-probe and MC truth (medium cuts)");
      cEtaBottom->cd(1)->SetLogy();
      // hEta[DEFAULT][ERROR]->Draw();
      hEta[MEDIUM][ERROR]->SetMaximum(hEta_MC[MEDIUM][ERROR]->GetMaximum()*2);
      hEta[MEDIUM][ERROR]->SetMinimum(hEta_MC[MEDIUM][ERROR]->GetMinimum()/2);
      hEta[MEDIUM][ERROR]->Draw("e");
      hEta_MC[MEDIUM][ERROR]->Draw("e SAME");
      if(ABS_ETA) hEta[MEDIUM][ERROR]->GetXaxis()->SetTitle("|#eta|");

      // Legend for error(Eta)
      TLegend *lEtaMedium = new TLegend(0.25,0.75,0.45,0.88);
      lEtaMedium->SetHeader("For medium cuts:");
      //lEta->AddEntry(hEta[DEFAULT][ERROR], "From tag-and-probe", "lpf");
      lEtaMedium->AddEntry(hEta[MEDIUM][ERROR], "From tag-and-probe", "lpf");
      lEtaMedium->AddEntry(hEta_MC[MEDIUM][ERROR], "From MC truth", "lpf");
      lEtaMedium->Draw();
 
      /// For tight cuts
      //TCanvas *cTight = new TCanvas("errorEtaTight", "Comparison of error rate (Eta) from tag-and-probe and MC truth (tight cuts)");
      cEta->cd(1)->SetLogy();
      // hEta[DEFAULT][ERROR]->Draw();
      hEta[TIGHT][ERROR]->SetMaximum(hEta_MC[TIGHT][ERROR]->GetMaximum()*2);
      hEta[TIGHT][ERROR]->SetMinimum(hEta_MC[TIGHT][ERROR]->GetMinimum()/2);
      hEta[TIGHT][ERROR]->Draw("e");
      hEta_MC[TIGHT][ERROR]->Draw("e SAME");
      if(ABS_ETA) hEta[MEDIUM][ERROR]->GetXaxis()->SetTitle("|#eta|");
  

      // Legend for error(Eta)
      TLegend *lEtaTight = new TLegend(0.25,0.75,0.45,0.88);
      lEtaTight->SetHeader("For tight cuts");
      //lEta->AddEntry(hEta[DEFAULT][ERROR], "From tag-and-probe", "lpf");
      lEtaTight->AddEntry(hEta[TIGHT][ERROR], "From tag-and-probe", "lpf");
      lEtaTight->AddEntry(hEta_MC[TIGHT][ERROR], "From MC truth", "lpf");
      lEtaTight->Draw();
      if(print) cEta->Print("eta.pdf");
   }
   
   // === FOR MISID(ETA, PT) ===
   //TCanvas *cEtaPt = new TCanvas("cEtaPt", "Error rates as a function of eta across different Pt ranges", 600, 600);
   /*
   for(int i = 0; i<nPtBins; i++){
      hEtaPt[i][ERROR]->SetMaximum(1);
      hEtaPt[i][ERROR]->SetMinimum(0.00001);
   }
   //*/
   //int toDraw[] = {0, 1, 2, 3, 4, 5, 6, 7};
   //cEtaPt->Divide(1, 8);
   
   if(drawEtaPt){
      TCanvas *cEtaPt[4][nPtBins];
      TLegend *lEtaPt[4][nPtBins];
   
      for(int j=1; j<=3; j++){
         string sCut;
         switch(j){ // comment out continue to include that particular cut
         case LOOSE:    sCut = "Loose";
                        continue;
                        break;
         case MEDIUM:   sCut = "Medium";
                        // continue;
                        break;
         case TIGHT:    sCut = "Tight";
                        continue;
                        break;
         }
         string sPdfFilename = "EtaPt" + sCut + ".pdf";
         TCanvas c("c", "c", 600, 600);
      
         if (print) c.Print((sPdfFilename+"[").c_str());
         for(int i=0; i<nPtBins; i++){
            // Set names of Canvases
            std::ostringstream lowPt, highPt;
            lowPt << std::fixed << std::setprecision(0) << PtEdges[i];
            highPt << std::fixed << std::setprecision(0) << PtEdges[i+1];

            
            string sCanTitle = "Error rate as a function of eta in " + lowPt.str() + "<= Pt < " + highPt.str() + " (" + sCut + " cut)";
            string sCanPad = "cEtaPt" + lowPt.str() + sCut;
            string sHistName = "hEP" + lowPt.str() + sCut;
            string sHistName_MC = sHistName + "_MC";
         
            cEtaPt[j][i] = new TCanvas(sCanPad.c_str(), sCanTitle.c_str(), 600, 600);
            cEtaPt[j][i]->SetLogy();
            //hEP[j][ERROR]->Draw();
            TH1* hEtaPt = hEP[j][ERROR]->ProjectionX(sHistName.c_str(), i+1, i+1, "e");
            hEtaPt->GetYaxis()->SetRangeUser(0.0001, 0.1);
            hEtaPt->SetMarkerStyle(20);
            hEtaPt->SetMarkerColor(kBlue);
            hEtaPt->SetLineColor(kBlue);
            hEtaPt->Draw("e");
            TH1* hEtaPt_MC = hEP_MC[j][ERROR]->ProjectionX(sHistName_MC.c_str(), i+1, i+1, "e");
            hEtaPt_MC->SetMarkerStyle(20);
            hEtaPt_MC->SetMarkerColor(kRed);
            hEtaPt_MC->SetLineColor(kRed);
            hEtaPt_MC->Draw("e same");
            if(ABS_ETA) hEtaPt->GetXaxis()->SetTitle("|#eta|");
         
            lEtaPt[j][i] = new TLegend(0.2, 0.8, 0.65, 0.88);
            string sLegendHeader = lowPt.str() + " #leq p_{T} < " + highPt.str() + " (" + sCut + " cut)";
            lEtaPt[j][i]->SetHeader(sLegendHeader.c_str());
            lEtaPt[j][i]->AddEntry(hEtaPt, "From tag-and-probe", "lpf");
            lEtaPt[j][i]->AddEntry(hEtaPt_MC, "From MC truth", "lpf");
            lEtaPt[j][i]->Draw();
         
            if(print){
               string sPdfTitle = "Title:" + sLegendHeader;
               cEtaPt[j][i]->Print(sPdfFilename.c_str(), sPdfTitle.c_str());
            }
         
         }
         if(print) c.Print((sPdfFilename+"]").c_str());
      }
   }
   
   if(drawPtEta){
      TCanvas *cPtEta[4][nEtaBins];
      TLegend *lPtEta[4][nEtaBins];
   
      for(int j=1; j<=3; j++){
         string sCut;
         switch(j){ // comment out continue to include that particular cut
         case LOOSE:    sCut = "Loose";
                        continue;
                        break;
         case MEDIUM:   sCut = "Medium";
                        // continue;
                        break;
         case TIGHT:    sCut = "Tight";
                        continue;
                        break;
         }
         string sPdfFilename = "PtEta" + sCut + ".pdf";
         TCanvas c("c", "c", 600,600);
      
         if (print) c.Print((sPdfFilename+"[").c_str());
         for(int i=0; i<nEtaBins; i++){
            // Set names of Canvases
            std::ostringstream lowEta, highEta;
            lowEta << std::fixed << std::setprecision(2) << EtaEdges[i];
            highEta << std::fixed << std::setprecision(2) << EtaEdges[i+1];	

         
            string sCanTitle;
            if (ABS_ETA) sCanTitle = "Error rate as a function of Pt in " + lowEta.str() + "<= |eta| < " + highEta.str() + " (" + sCut + " cut)";
            else sCanTitle = "Error rate as a function of Pt in " + lowEta.str() + "<= eta < " + highEta.str() + " (" + sCut + " cut)";
            string sCanPad = "cPtEta" + lowEta.str() + sCut;
            string sHistName = "hPE" + lowEta.str() + sCut;
            string sHistName_MC = sHistName + "_MC";
         
            cPtEta[j][i] = new TCanvas(sCanPad.c_str(), sCanTitle.c_str(), 600, 600);
            cPtEta[j][i]->SetLogy();
            //hEP[j][ERROR]->Draw();
            TH1* hPtEta = hEP[j][ERROR]->ProjectionY(sHistName.c_str(), i+1, i+1, "e");
            hPtEta->GetYaxis()->SetRangeUser(0.0001, 0.1);
            hPtEta->SetMarkerStyle(20);
            hPtEta->SetMarkerColor(kBlue);
            hPtEta->SetLineColor(kBlue);
            hPtEta->Draw("e");
            TH1* hPtEta_MC = hEP_MC[j][ERROR]->ProjectionY(sHistName_MC.c_str(), i+1, i+1, "e");
            hPtEta_MC->SetMarkerStyle(20);
            hPtEta_MC->SetMarkerColor(kRed);
            hPtEta_MC->SetLineColor(kRed);
            hPtEta_MC->Draw("e same");
         
            lPtEta[j][i] = new TLegend(0.2, 0.8, 0.65, 0.88);
            string sLegendHeader;
            if (ABS_ETA) sLegendHeader = lowEta.str() + " #leq |#eta| < " + highEta.str() + " (" + sCut + " cut)";
            else sLegendHeader = lowEta.str() + " #leq #eta < " + highEta.str() + " (" + sCut + " cut)";
            lPtEta[j][i]->SetHeader(sLegendHeader.c_str());
            lPtEta[j][i]->AddEntry(hPtEta, "From tag-and-probe", "lpf");
            lPtEta[j][i]->AddEntry(hPtEta_MC, "From MC truth", "lpf");
            lPtEta[j][i]->Draw();
         
            if(print){
               string sPdfTitle = "Title:" + sLegendHeader;
               cPtEta[j][i]->Print(sPdfFilename.c_str(), sPdfTitle.c_str());
            }
         
         }
         if(print) c.Print((sPdfFilename+"]").c_str());
      }
   }
   //*/
   
   // ==== FOR ORIGIN ====
   /*
   TCanvas *cPairOrig = new TCanvas("cPairOrig", "Origins of electrons SS pairs");
   hPairOrigin->Draw("lego");
   if(print) cPairOrig->Print("MC SS Pair origin.pdf");
   //*/
   
   
   // === TO FILES === 
   ToCSV();
   if(toFile) outputFile->Write();
   // ToFile();
   
   // === TERMINAL OUTPUT ===
   PRINT_STAT("Percentage of data that are real electrons", nRealE/(double) nTotalE);
   PRINT_STAT("Total number of events analyzed", nEvents);
   PRINT_STAT("No. events that passed the event cut", nPassedEvents);
   PRINT_STAT("Number of SS events", nSSEvents);
   PRINT_STAT("Number of SS events (MC)", nSSEvents_MC);
   PRINT_STAT( "","" );
   PRINT_STAT("Total number of electrons", nTotalE);
   PRINT_STAT("MisIDed from Z", hParticleOriginSS->GetBinContent(24));
   PRINT_STAT("Number from Z ", hParticleOrigin->GetBinContent(24));
   for(int i = 1; i<=200; i++){
      if (i==24) continue;
      if (hParticleOrigin->GetBinContent(i) == 0) continue;
      PRINT_STAT("MisIDed from " << i-1 , hParticleOriginSS->GetBinContent(i)); 
      PRINT_STAT("Number from " << i-1 << " ", hParticleOrigin->GetBinContent(i));  
   }
   PRINT_STAT( "" , "" );
   PRINT_STAT("Electrons from Z", nFromZ);
   PRINT_STAT("Status", "Done");
   file_eParents.close();
   file_misID.close();
   file_nProbes.close();
   // if(toFile) outputFile->Close();
}



/** ProcessData()
 * Runs on data
 * Separated from Process() for readability
 */
inline void TagProbePlot::ProcessData(){
   //============//
   // Statistics //
   //============//
   nPassedEvents++;
   if(onlyZMCevents && !EventIsZMC()) return;

   //==============================//
   // Begin tag-and-probe analysis //
   //==============================//
   if(ElectronIsTag(1)){
      nTpPairs++;
      if(PassedLHCut(2, LOOSE)) {
         FillEPHist(fE2_fEta, fE2_fPt, hEP[LOOSE], 0, ABS_ETA);
         FillHist(fE2_fPt, hPt[LOOSE], 0);
         FillHist(fE2_fEta, hEta[LOOSE], 0, ABS_ETA);
         FillHist(fInvMass, hMass[LOOSE], 0);
      }
      if(PassedLHCut(2, MEDIUM)){ 
         FillEPHist(fE2_fEta, fE2_fPt, hEP[MEDIUM], 0, ABS_ETA);
         FillHist(fE2_fPt, hPt[MEDIUM], 0);
         FillHist(fE2_fEta, hEta[MEDIUM], 0, ABS_ETA);
         FillHist(fInvMass, hMass[MEDIUM], 0);
      }
      if(PassedLHCut(2, TIGHT)) {
         FillEPHist(fE2_fEta, fE2_fPt, hEP[TIGHT], 0, ABS_ETA);
         FillHist(fE2_fPt, hPt[TIGHT], 0);
         FillHist(fE2_fEta, hEta[TIGHT], 0, ABS_ETA);
         FillHist(fInvMass, hMass[TIGHT], 0);
         
         FillEPHist(fE2_fEta, fInvMass, hEtaMass[TIGHT], 0, ABS_ETA);
         FillEPHist(fE2_fPt, fInvMass, hPtMass[TIGHT], 0);
         
         if(fMCevent && !EventIsZMC()){
            FillEPHist(fE2_fEta, fInvMass, hEtaMassNoZ[TIGHT], 0, ABS_ETA);
            FillEPHist(fE2_fPt, fInvMass, hPtMassNoZ[TIGHT], 0);
         
         }
         
      }

   }

   if(ElectronIsTag(2)) {
      nTpPairs++;
      if(PassedLHCut(1, LOOSE)) {
         FillEPHist(fE1_fEta, fE1_fPt, hEP[LOOSE], 0, ABS_ETA);
         FillHist(fE1_fPt, hPt[LOOSE], 0);
         FillHist(fE1_fEta, hEta[LOOSE], 0, ABS_ETA);
         FillHist(fInvMass, hMass[LOOSE], 0);
      }
      if(PassedLHCut(1, MEDIUM)){ 
         FillEPHist(fE1_fEta, fE1_fPt, hEP[MEDIUM], 0, ABS_ETA);
         FillHist(fE1_fPt, hPt[MEDIUM], 0);
         FillHist(fE1_fEta, hEta[MEDIUM], 0, ABS_ETA);
         FillHist(fInvMass, hMass[MEDIUM], 0);
      }
      if(PassedLHCut(1, TIGHT)) {
         FillEPHist(fE1_fEta, fE1_fPt, hEP[TIGHT], 0, ABS_ETA);
         FillHist(fE1_fPt, hPt[TIGHT], 0);
         FillHist(fE1_fEta, hEta[TIGHT], 0, ABS_ETA);
         FillHist(fInvMass, hMass[TIGHT], 0);
         
         FillEPHist(fE1_fEta, fInvMass, hEtaMass[TIGHT], 0, ABS_ETA);
         FillEPHist(fE1_fPt, fInvMass, hPtMass[TIGHT], 0);
         
         if(fMCevent && !EventIsZMC()){
            FillEPHist(fE1_fEta, fInvMass, hEtaMassNoZ[TIGHT], 0, ABS_ETA);
            FillEPHist(fE1_fPt, fInvMass, hPtMassNoZ[TIGHT], 0);
         
         }
      }
   }
}

/** ProcessMC()
 * Run on MC truth. Separated from Process() for readability
 */
inline void TagProbePlot::ProcessMC(){
   nMCevents++;
   if(EventIsZMC()){
      //===============================//
      // Measure misID on all electrons//
      //===============================//
      
      if(ElectronIsTag(2, false)){         
         if(PassedLHCut(1, LOOSE)) {
            FillEPHist(fE1_fEta, ftE1_fPt, hEP_MC[LOOSE], 1, ABS_ETA);
            FillHist(fE1_fPt, hPt_MC[LOOSE], 1);
            FillHist(fE1_fEta, hEta_MC[LOOSE], 1, ABS_ETA);
            FillHist(fInvMass, hMass_MC[LOOSE], 1);
         }
         if(PassedLHCut(1, MEDIUM)){ 
            FillEPHist(fE1_fEta, ftE1_fPt, hEP_MC[MEDIUM], 1, ABS_ETA);
            FillHist(fE1_fPt, hPt_MC[MEDIUM], 1);
            FillHist(fE1_fEta, hEta_MC[MEDIUM], 1, ABS_ETA);
            FillHist(fInvMass, hMass_MC[MEDIUM], 1);
         }
         if(PassedLHCut(1, TIGHT)) {
            FillEPHist(fE1_fEta, ftE1_fPt, hEP_MC[TIGHT], 1, ABS_ETA);
            FillHist(fE1_fPt, hPt_MC[TIGHT], 1);
            FillHist(fE1_fEta, hEta_MC[TIGHT], 1, ABS_ETA);
            FillHist(fInvMass, hMass_MC[TIGHT], 1);
         }
      }
      
      if(ElectronIsTag(1, false)){         
         if(PassedLHCut(2, LOOSE)) {
            FillEPHist(fE2_fEta, ftE2_fPt, hEP_MC[LOOSE], 2, ABS_ETA);
            FillHist(fE2_fPt, hPt_MC[LOOSE], 2);
            FillHist(fE2_fEta, hEta_MC[LOOSE], 2, ABS_ETA);
            FillHist(fInvMass, hMass_MC[LOOSE], 2);
         }
         if(PassedLHCut(2, MEDIUM)){ 
            FillEPHist(fE2_fEta, ftE2_fPt, hEP_MC[MEDIUM], 2, ABS_ETA);
            FillHist(fE2_fPt, hPt_MC[MEDIUM], 2);
            FillHist(fE2_fEta, hEta_MC[MEDIUM], 2, ABS_ETA);
            FillHist(fInvMass, hMass_MC[MEDIUM], 2);
         }
         if(PassedLHCut(2, TIGHT)) {
            FillEPHist(fE2_fEta, ftE2_fPt, hEP_MC[TIGHT], 2, ABS_ETA);
            FillHist(fE2_fPt, hPt_MC[TIGHT], 2);
            FillHist(fE2_fEta, hEta_MC[TIGHT], 2, ABS_ETA);
            FillHist(fInvMass, hMass_MC[TIGHT], 2);
         }
      }
      
      //===============//
      // Stat counters //
      //===============//
      if(ftE1_fFromZ) nFromZ++;
      if(ftE2_fFromZ) nFromZ++;

      nTotalE += 2;
      if(ftE1_fIsElectron) nRealE++;
      if(ftE2_fIsElectron) nRealE++;
      if(EventIsZMC()) nRealZ++;
      if(ftSSevent) nSSEvents_MC++;

      hParticleOrigin->Fill(ftE2_fParticleOrigin);
      if(ftE2_fCharge != fE2_fCharge) hParticleOriginSS->Fill(ftE2_fParticleOrigin);
      hParticleOrigin->Fill(ftE1_fParticleOrigin);
      if(ftE1_fCharge != fE1_fCharge) hParticleOriginSS->Fill(ftE1_fParticleOrigin);
      /*if (ftSSevent){
         
      }*/
   }
   else {
      PrintParticleOrigin();
      hPairOrigin->Fill(ftE1_fParticleOrigin, ftE2_fParticleOrigin);
   }
}

inline bool TagProbePlot::MisId(int code){
   if(code == 0) return (fE1_fCharge == fE2_fCharge);
   
   else if (code == 1){
      if( fabs(ftE1_fCharge) != 1){
         throw 1;
      }

      if(ftE1_fTruthType == 4){
         if(ftE1_fTruthBkgType == 2) return (fE1_fCharge != ftE1_fTruthBkgCharge);
         else throw 1;
      }
      else return (fE1_fCharge != ftE1_fCharge);
   }
   else if (code == 2){
      if( fabs(ftE2_fCharge) != 1){
         throw 1;
      }

      if(ftE2_fTruthType == 4) {
         if (ftE2_fTruthBkgType == 2) return (fE2_fCharge != ftE2_fTruthBkgCharge);
         else throw 1;
      }
      else return (fE2_fCharge != ftE2_fCharge);
   }
   else throw 2;
}

/** FillHist()
 * Fills the total histogram and the SS histogram with var.
 *
 * var   - the variable to be filled
 * h     - set of histograms to be filled
 * code  - the condition code for misID:
         - 0: data;     1: MC electron 1;       2: MC electron 2
 */
void TagProbePlot::FillHist(Double_t var, TH1* h[], int code, bool absolute){
   Bool_t misID;
   try{
      misID = MisId(code);
      if(absolute) var = fabs(var);
      h[N_PROBES]->Fill(var);
      if (misID) h[N_SSPROBES]->Fill(var); 
      else h[N_OSPROBES]->Fill(var);
   }
   catch(int e){
   
   }
   return;
}

void TagProbePlot::FillEPHist(Double_t eta, Double_t pt, TH2* h[4], int code, bool absEta){
   Bool_t misID;
   
   try{
      misID = MisId(code);
      if(absEta) eta = fabs(eta);
      h[N_PROBES]->Fill(eta, pt);
      if(misID) h[N_SSPROBES]->Fill(eta, pt);
      else h[N_OSPROBES]->Fill(eta, pt);
   }
   catch(int e){
   
   }
   return;
}

/** DivideHist()
 * Produces the "error" histogram by hError = hSS/hTot
 */
void TagProbePlot::DivideHist(TH1* h[]){
   h[N_PROBES]->Sumw2();
   h[N_SSPROBES]->Sumw2();
   if(h[N_OSPROBES] != 0) h[N_OSPROBES]->Sumw2();

   h[ERROR]->Add(h[N_SSPROBES]);
   h[ERROR]->Divide(h[N_PROBES]);
}


/** PassedEventCut()
 * Function stand-in for event cuts. Put all conditions here.
 *
 * Current implementation returns true iff the event passes the following cut:
 *  - invariant mass is within the Z mass window
 *  - dPhi > 2
 */
inline Bool_t TagProbePlot::PassedEventCut(){
   if (!(fInvMass < Z_MASS_UPPER_LIMIT && fInvMass > Z_MASS_LOWER_LIMIT)) return false;
   else if (fabs(fdPhi) < 2) return false;
   
   else return true;
}


/** PassedLHCut()
 * Extracts the LH cut information
 * elec - electron number
 * cut - 1 = light cut; 2 = medium cut; 3 = tight cut; (only on reconstructed e)
 */
Bool_t TagProbePlot::PassedLHCut(int elec, int cut){
   //if(ftSSevent) return false;
   int LHcut;
   if(elec == 1) LHcut = fE1_fPassCut;
   else if (elec == 2) LHcut = fE2_fPassCut;
   else return false;
   
   
   if (cut == LOOSE) return ((LHcut / 2) % 2 == 1); // passed light cut
   if (cut == MEDIUM) return ((LHcut / 4) % 2 == 1); // passed medium cut
   if (cut == TIGHT) return ((LHcut / 8) % 2 == 1); // passed tight cut
   
   return false;
   
   
}


/**EventIsZMC()
 * Returns whether the both electrons are from Z
 */
inline Bool_t TagProbePlot::EventIsZMC(){
   return ftZEEevent;
//    if(ftE1_fPt != 0 && ftE2_fPt != 0) return true;
   
   
   return false;
}


/**ElectronIsTag()
 * Applies conditions to choose a tag electron
 * 
 * elec - electron number
 * code - criteria: 0: default (check fTag flag)
 *                  1: user defined
 */
Bool_t TagProbePlot::ElectronIsTag(int elec, bool isData, int mode){
   if(mode % 10 == AUTO){
      if (elec == 1) return fE1_fTag;
      else if (elec == 2) return fE2_fTag;
      else return false;
   }

   double eta;
   double pt;
   double d0;
   double d0error;
   double BHits;
   double SiHits;

   if (elec == 1){
      eta = fE1_fEta;
      pt = fE1_fPt;
      d0 = fE1_fD0;
      d0error = fE1_fD0Error;
      BHits = fE1_fBHits;
      SiHits = fE1_fSiHits;
   }
   else if (elec == 2){
      eta = fE2_fEta;
      pt = fE2_fPt;
      d0 = fE2_fD0;
      d0error = fE2_fD0Error;
      BHits = fE2_fBHits;
      SiHits = fE2_fSiHits;
   }
   else return false;

   if(isData) hTagCuts->Fill(0);
   
   if(!PassedLHCut(elec,TIGHT)) return false;
   if(isData) hTagCuts->Fill(1);

   if (!(fabs(eta) <= 1.37)) return false;
   if(isData) hTagCuts->Fill(2);

   if ((d0error == 0) || !(fabs(d0/d0error) < 1.5)) return false;
   if (isData) hTagCuts->Fill(3);

   if (BHits < 1) return false;
   if(isData) hTagCuts->Fill(4);

   if (SiHits < 9) return false;
   if(isData) hTagCuts->Fill(5);

   return true;
}


void TagProbePlot::PrintParticleOrigin(){
   file_eParents << ftE1_fParticleOrigin << "\t" << ftE2_fParticleOrigin <<endl;
}


void TagProbePlot::ToCSV(){
   file_nProbes << "Number of SS probes / Number of probes \n,\n";
   file_misID << "Misidentification rate (%) \n,\n";

   file_nProbes << "Tight cuts" << std::endl;
   file_misID << "Tight cuts" << std::endl; 
   ToCSV(TIGHT);
   
   file_nProbes << std::endl << std::endl << "Medium cuts" << std::endl;
   file_misID << std::endl << std::endl << "Medium cuts" << std::endl; 
   ToCSV(MEDIUM);
   
   file_nProbes << std::endl << std::endl << "Loose cuts" << std::endl;
   file_misID << std::endl << std::endl << "Loose cuts" << std::endl; 
   ToCSV(LOOSE);
   
   file_nProbes << std::endl << std::endl << "MC tight" << std::endl;
   file_misID << std::endl << std::endl << "MC tight" << std::endl; 
   ToCSV(MC+TIGHT);
   
   file_nProbes << std::endl << std::endl << "MC medium" << std::endl;
   file_misID << std::endl << std::endl << "MC medium" << std::endl; 
   ToCSV(MC+MEDIUM);
   
   file_nProbes << std::endl << std::endl << "MC loose" << std::endl;
   file_misID << std::endl << std::endl << "MC loose" << std::endl; 
   ToCSV(MC+LOOSE);
}


void TagProbePlot::ToCSV(int cut){
   TH1 **hE = 0;
   TH1 **hP = 0;
   TH2 **hEtaPt = 0;
   switch(cut){
   case LOOSE:
      hE = hEta[LOOSE];
      hP = hPt[LOOSE];
      hEtaPt = hEP[LOOSE];
      break;
   case MEDIUM:
      hE = hEta[MEDIUM];
      hP = hPt[MEDIUM];
      hEtaPt = hEP[MEDIUM];
      break;
   case TIGHT:
      hE = hEta[TIGHT];
      hP = hPt[TIGHT];
      hEtaPt = hEP[TIGHT];
      break;
   case (MC+LOOSE):
      hE = hEta_MC[LOOSE];
      hP = hPt_MC[LOOSE];
      hEtaPt = hEP_MC[LOOSE];
      break;
   case (MC+MEDIUM):
      hE = hEta_MC[MEDIUM];
      hP = hPt_MC[MEDIUM];
      hEtaPt = hEP_MC[MEDIUM];
      break;
   case (MC+TIGHT):
      hE = hEta_MC[TIGHT];
      hP = hPt_MC[TIGHT];
      hEtaPt = hEP_MC[TIGHT];
      break;
   default:
      return;
   }
   
   // nProbes
   ///////////////
   file_nProbes << "Eta bins,Pt bins" << std::endl;
   file_nProbes <<",";
   
   for(int i = 0; i<= nPtBins; i++){
      file_nProbes << PtEdges[i] << ",";  
   }
   file_nProbes << "Total," << std::endl;
   //*
   for(int i = 1; i <= nEtaBins; i++){
      file_nProbes << EtaEdges[i-1] << ",";
      for(int j = 1; j<= nPtBins; j++){
         file_nProbes << hEtaPt[N_SSPROBES]->GetBinContent(i,j) 
                     << " / " 
                     << hEtaPt[N_PROBES]->GetBinContent(i, j) 
                     << ",";
      }
      file_nProbes << ",";
      file_nProbes << hE[N_SSPROBES]->GetBinContent(i) 
                  << " / " 
                  << hE[N_PROBES]->GetBinContent(i) 
                  << std::endl;
   }
   file_nProbes << EtaEdges[nEtaBins] << std::endl;
   file_nProbes << "Total,";
   for(int i=1; i<= nPtBins; i++){
      file_nProbes << hP[N_SSPROBES]->GetBinContent(i) 
                  << " / " 
                  << hP[N_PROBES]->GetBinContent(i) 
                  << "," ;
   }
   file_nProbes << std::endl << "," << std::endl << "," << std::endl;
   
   // misID
   /////////////////////
   file_misID << "Eta bins,Pt bins" << std::endl;
   file_misID << ",";
   file_misID.unsetf(std::ios_base::floatfield);
   for(int i = 0; i<= nPtBins; i++){
      file_misID << PtEdges[i] << ",";  
   }
   file_misID << "Total" << std::endl;
   
   for(int i = 1; i <= nEtaBins; i++){
      file_misID.unsetf(std::ios_base::floatfield);
      file_misID << EtaEdges[i-1] << ",";
      for(int j = 1; j<= nPtBins; j++){
         file_misID.precision(5);
         file_misID << std::fixed << hEtaPt[ERROR]->GetBinContent(i, j) * 100 << ",";
      }
      file_misID << "," << hE[ERROR]->GetBinContent(i) << std::endl;
   }
   file_misID << EtaEdges[nEtaBins] << std::endl;
   file_misID << "Total,";
   for(int i=1; i<= nPtBins; i++){
      file_misID << hP[ERROR]->GetBinContent(i) << ",";
   }
   file_misID << std::endl << "," << std::endl << "," << std::endl;
}

void TagProbePlot::ToFile(){
   //outputFile = new TFile("output.root", "RECREATE", "TpHistograms");
   // outputFile->Write();
   // hEP[MEDIUM][N_PROBES]->Write();
   // hEP[MEDIUM][N_SSPROBES]->Write();
}