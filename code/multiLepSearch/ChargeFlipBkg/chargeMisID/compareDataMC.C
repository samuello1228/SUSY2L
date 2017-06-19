/* compareDataMC.cxx
 * Written by Gabriel Gallardo Feb 2016
 *
 * A script to compare the chargeMisID rates obtained by MC truth and likelihood method on data
 * Rates are obtained by the use of chargeMisID.cxx
 *
 * Execution instructions:
 * - Choose which graphs should be drawn (drawEtaPt, drawPtEta, drawPt, drawEta)
 * $ root -l -b 
 * root[0] .x compareDataMC.C("path/dataNoSub/histograms.root", ""path/dataSubbed/histograms.root", "path/MC/histograms.root", "outputDir")
 * 
 * Excecution recommendations:
 * - Run in batch mode! Otherwise a lot of historgrams will pop out at once. 
 *    Suggest to print to pdf, then open the plots in a pdf viewer.
 * - compilation not necessary. (i.e. no need for .x compareDataMC.cxx+)
 */

 // === Config variables === //
   
const bool drawEtaPt = true;    // draw rates over eta in different pt bins
const bool drawPtEta = true;    // draw rates of pt in different eta bins
const bool drawPt = false;      // draw projection of rates over pt 
const bool drawEta = false;     // draw projection of rates over eta 
const bool logPt = true;

const bool print = true;        // output to pdf
const bool toFile = true;       // output to root file

const bool drawLegend = true;
const bool drawData = true;
const bool drawSub = false;
const bool drawMC = true;
const bool drawMCLH = false;
const bool drawCompare = true && (drawData && drawMC);

// ===== Set names in legend ==== // 
const char sData[]="MC LH";
const char sDataSubbed[]="";
const char sMC[]="MC truth";
const char sMCLH[]="";  
const char sCompare[]="LH/truth";

void splitPad(TCanvas* c){
   c->Divide(1,2);
   c->cd(1)->SetPad(0, 0.25, 1, 1);
   c->cd(2)->SetPad(0, 0, 1, 0.25);
}

void drawRatio(TCanvas* c, TH1* hDenom, TH1* hNum1, TH1* hNum2=0){
   TH1* hCompare = (TH1*) hNum1->Clone("hCompare");
   hCompare->Divide(hDenom);
   c->cd(2);
   hCompare->GetXaxis()->SetLabelOffset(999);
   hCompare->GetXaxis()->SetLabelSize(0);
   hCompare->GetXaxis()->SetTitle("");
   TAxis* hCompareY = hCompare->GetYaxis();
   hCompareY->SetTitle(sCompare);
   hCompareY->SetTitleSize(0.1);
   hCompareY->SetTitleOffset(0.4);
   hCompareY->SetLabelSize(0.08);
   hCompare->SetMarkerStyle(20);
   hCompareY->SetRangeUser(0.5,2);
   hCompare->Draw();

   if(hNum2){
      TH1* hCompare1 = (TH1*) hNum2->Clone("hCompare1");
      hCompare1->Divide(hDenom);
      hCompare1->Draw("same");
   }

   TLine l;
   TLine *l1 = l.DrawLine(hCompare->GetXaxis()->GetXmin(), 1., hCompare->GetXaxis()->GetXmax(), 1.);
   l1->SetLineStyle(2);
   l1->SetLineWidth(2);
}

bool draw(TH2* hData, TH2* hDataSubbed, TH2* hMC, TH2* hMCLH){

   // ================================================== //
   const int nPtBins = hData->GetYaxis()->GetNbins(); 
   const int nEtaBins = hData->GetXaxis()->GetNbins();

   const double* PtEdges = hData->GetYaxis()->GetXbins()->GetArray();
   const double* EtaEdges = hData->GetXaxis()->GetXbins()->GetArray();

   const bool ABS_ETA = (EtaEdges[0] == 0);

   TFile *outputFile = 0;
   if(toFile) outputFile = new TFile("comparisons.root", "RECREATE");

   // ===== Set styles ======
   gStyle->SetOptStat(0);

   hData->SetMarkerStyle(20); 
   hData->SetMarkerColor(kGreen+4);
   hData->SetLineColor(kGreen+4);

   if(drawSub){
      hDataSubbed->SetMarkerStyle(22);
      hDataSubbed->SetMarkerColor(kBlue);
      hDataSubbed->SetLineColor(kBlue);
   }

   if(drawMCLH){
      hMCLH->SetMarkerStyle(21);
      hMCLH->SetMarkerColor(kMagenta+2); 
      hMCLH->SetLineColor(kMagenta+2);
   }

   if(drawMC){
      hMC->SetMarkerStyle(23);
      hMC->SetMarkerColor(kRed); 
      hMC->SetLineColor(kRed);  
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
      
         if(drawCompare) splitPad(cEtaPt[ptBin]);
         cEtaPt[ptBin]->cd(1);
         cEtaPt[ptBin]->cd(1)->SetLogy();

         lEtaPt[ptBin] = new TLegend(0.6, 0.35, 0.9, 0.15);
         string sLegendHeader = lowPt.str() + " #leq p_{T} < " + highPt.str() + " (GeV)";
         lEtaPt[ptBin]->SetHeader(sLegendHeader.c_str());
         // lEtaPt[ptBin]->SetTextSize(0.07);
         // lEtaPt[ptBin]->SetBorderSize(0);

         // ----------------- Data plot --------------- //
         TH1* hEtaPt_LH = 0;
         hEtaPt_LH = hData->ProjectionX((sHistName+"LH").c_str(), ptBin+1, ptBin+1, "e");

         if(ptBin <= 2)hEtaPt_LH->GetYaxis()->SetRangeUser(1e-4, 0.01);
         else hEtaPt_LH->GetYaxis()->SetRangeUser(5e-4, 0.05);
         hEtaPt_LH->SetTitle("");
         hEtaPt_LH->GetYaxis()->SetTitle("misID rate");
         hEtaPt_LH->GetYaxis()->SetTitleOffset(1.4);
         if(drawData) hEtaPt_LH->Draw("e1");
         if(ABS_ETA) hEtaPt_LH->GetXaxis()->SetTitle("|#eta|");

         if(drawLegend && drawData) lEtaPt[ptBin]->AddEntry(hEtaPt_LH, sData, "lp");

         TH1* hEtaPt_LHSubbed = 0;
         if (drawSub){
            hEtaPt_LHSubbed = hDataSubbed->ProjectionX((sHistName+"LHSubbed").c_str(), ptBin+1, ptBin+1, "e");
            hEtaPt_LHSubbed->Draw("e1 same");
            if(drawLegend) lEtaPt[ptBin]->AddEntry(hEtaPt_LHSubbed, sDataSubbed, "lpf");
         }
         
         // --------------- MC plots --------------- //
         TH1* hEtaPt_MCLH = 0;
         if (drawMCLH){
            hEtaPt_MCLH = hMCLH->ProjectionX((sHistName+"MCLH").c_str(), ptBin+1, ptBin+1, "e");
            hEtaPt_MCLH->Draw("e1 same");
            if(drawLegend)  lEtaPt[ptBin]->AddEntry(hEtaPt_MCLH, sMCLH, "lpf");
         }

         TH1* hEtaPt_MC = 0;
         if(drawMC){
            hEtaPt_MC = hMC->ProjectionX((sHistName+"MC").c_str(), ptBin+1, ptBin+1, "e");
            hEtaPt_MC->Draw("e1 same");
            if(drawLegend)  lEtaPt[ptBin]->AddEntry(hEtaPt_MC, sMC, "lpf");
         }

         lEtaPt[ptBin]->Draw();

         if(drawCompare) drawRatio(cEtaPt[ptBin], hEtaPt_MC, hEtaPt_LH, hEtaPt_LHSubbed);

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

         string sCanTitle = "Charge misID rate as a function of p_{T} in " + lowEta.str() + "<=";
         if(ABS_ETA) sCanTitle += "|#eta| < " + highEta.str();
         else sCanTitle += "Eta  < " + highEta.str();
         string sCanPad = "cPtEta" + lowEta100.str();
         string sHistName = "hPE" + lowEta100.str();
         string sHistName_MC = sHistName + "_MC";
      
         cPtEta[etaBin] = new TCanvas(sCanPad.c_str(), sCanTitle.c_str(), 600, 600);

         if(drawCompare) splitPad(cPtEta[etaBin]);
         cPtEta[etaBin]->cd(1);
         cPtEta[etaBin]->cd(1)->SetLogy();
         if(logPt) cPtEta[etaBin]->cd(1)->SetLogx();

         lPtEta[etaBin] = new TLegend(0.6, 0.35, 0.9, 0.15);
         // if (etaBin <=2) lPtEta[etaBin] = new TLegend(0.15, 0.83, 0.5, 0.7);
         // else lPtEta[etaBin] = new TLegend(0.15, 0.43, 0.5, 0.3);
         string sLegendHeader = lowEta.str() + " #leq |#eta| < " + highEta.str();
         lPtEta[etaBin]->SetHeader(sLegendHeader.c_str());
         // lPtEta[etaBin]->SetTextSize(0.07);
         // lPtEta[etaBin]->SetBorderSize(0);

         // ----------------- Data plot --------------- //
         TH1* hPtEta_LH = 0;
         hPtEta_LH = hData->ProjectionY((sHistName+"LH").c_str(), etaBin+1, etaBin+1, "e");

         if(etaBin <=4) hPtEta_LH->GetYaxis()->SetRangeUser(1e-4, 1e-2);
         else hPtEta_LH->GetYaxis()->SetRangeUser(5e-4, 5e-2);
         hPtEta_LH->SetTitle("");
         hPtEta_LH->GetXaxis()->SetTitle("p_{T} (GeV)");
         hPtEta_LH->GetYaxis()->SetTitleOffset(1.4);
         hPtEta_LH->GetYaxis()->SetTitle("misID rate");
         hPtEta_LH->GetXaxis()->SetMoreLogLabels();
         if(drawData) hPtEta_LH->Draw("e1");

         if(drawLegend && drawData) lPtEta[etaBin]->AddEntry(hPtEta_LH, sData, "lp");

         TH1* hPtEta_LHSubbed = 0;
         if(drawSub){
            hPtEta_LHSubbed = hDataSubbed->ProjectionY((sHistName+"LHSubbed").c_str(), etaBin+1, etaBin+1, "e");
            hPtEta_LHSubbed->Draw("e1 same");
            if(drawLegend) lPtEta[etaBin]->AddEntry(hPtEta_LHSubbed, sDataSubbed, "lpf");
         }

         // -----------------  MC Truth ---------------- //
         TH1* hPtEta_MCLH = 0;
         if (drawMCLH){
            hPtEta_MCLH = hMCLH->ProjectionY((sHistName+"MCLH").c_str(), etaBin+1, etaBin+1, "e");
            hPtEta_MCLH->Draw("e1 same");
            if(drawLegend) lPtEta[etaBin]->AddEntry(hPtEta_MCLH, sMCLH, "lpf");
         }
         
         TH1* hPtEta_MC = 0;
         if (drawMC){
            hPtEta_MC = hMC->ProjectionY((sHistName+"MC").c_str(), etaBin+1, etaBin+1, "e");
            hPtEta_MC->Draw("e1 same");
            if(drawLegend) lPtEta[etaBin]->AddEntry(hPtEta_MC, sMC, "lpf");
         }

         lPtEta[etaBin]->Draw();

         if(drawCompare) {cPtEta[etaBin]->cd(2)->SetLogx(); drawRatio(cPtEta[etaBin], hPtEta_MC, hPtEta_LH, hPtEta_LHSubbed);}

         if(print){
            string sPdfTitle = "Title:" + sLegendHeader;
            cPtEta[etaBin]->Print(sPdfFilename.c_str(), sPdfTitle.c_str());
         }
      }

      if(print) c.Print((sPdfFilename+"]").c_str());
   }

   if(drawPt){
      TCanvas* cPt = new TCanvas("cPt", "cPt", 600, 600); cPt->SetLogy();
      TLegend* lPt = new TLegend(0.6, 0.35, 0.9, 0.15);

      if(drawCompare) splitPad(cPt);
      cPt->cd(1);
      cPt->cd(1)->SetLogy();
      if(logPt) cPt->cd(1)->SetLogx();

      TH1* hPt_LH = hData->ProjectionY("hPt_LH");
      hPt_LH->GetYaxis()->SetRangeUser(0.01, 0.5);
      hPt_LH->SetTitle("");
      hPt_LH->GetYaxis()->SetTitle("misID rate");
      hPt_LH->GetYaxis()->SetTitleOffset(1.1);
      hPt_LH->GetXaxis()->SetTitle("p_{T} (GeV)");
      hPt_LH->GetXaxis()->SetMoreLogLabels();
      if(drawData) hPt_LH->Draw("e1");
      if(drawData) lPt->AddEntry(hPt_LH, sData, "lpf");

      TH1* hPt_LHSubbed = 0;
      TH1* hPt_MCLH = 0;
      TH1* hPt_MC = 0;

      if(drawSub){
         hPt_LHSubbed = hDataSubbed->ProjectionY("hPt_LHSubbed");
         hPt_LHSubbed->Draw("e1 same");
         lPt->AddEntry(hPt_LHSubbed, sDataSubbed, "lpf");
      }

      if (drawMCLH){
         hPt_MCLH = hMCLH->ProjectionY("hPt_MCLH");
         hPt_MCLH->Draw("e1 same");
         lPt->AddEntry(hPt_MCLH, sMCLH, "lpf");
      }

      if (drawMC){
         hPt_MC = hMC->ProjectionY("hPt_MC");
         hPt_MC->Draw("e1 same");
         lPt->AddEntry(hPt_MC, sMC, "lpf");
      }

      lPt->Draw();

      if(drawCompare) {cPt->cd(2)->SetLogx(); drawRatio(cPt, hPt_MC, hPt_LH, hPt_LHSubbed);}

      if(print){
         cPt->Print("Pt.pdf", "Title:Pt");
      }
   }

   if(drawEta){
      TCanvas* cEta = new TCanvas("cEta", "cEta", 600, 600); cEta->SetLogy();
      TLegend* lEta = new TLegend(0.6, 0.35, 0.9, 0.15);

      if(drawCompare) splitPad(cEta);
      cEta->cd(1);
      cEta->cd(1)->SetLogy();

      TH1* hEta_LH = hData->ProjectionX("hEta_LH");
      hEta_LH->SetTitle("");
      hEta_LH->GetYaxis()->SetTitle("misID rate");
      hEta_LH->GetYaxis()->SetRangeUser(0.001, 0.4);
      if(ABS_ETA) hEta_LH->GetXaxis()->SetTitle("|#eta|");
      if(drawData) hEta_LH->Draw("e1");
      
      if(drawData) lEta->AddEntry(hEta_LH, sData, "lpf");

      TH1* hEta_LHSubbed = 0;
      if(drawSub){
         hEta_LHSubbed = hDataSubbed->ProjectionX("hEta_LHSubbed");
         hEta_LHSubbed->Draw("e1 same");
         lEta->AddEntry(hEta_LHSubbed, sDataSubbed, "lpf");
      }

      TH1* hEta_MCLH = 0;
      if(drawMCLH){
         hEta_MCLH = hMCLH->ProjectionX("hEta_MCLH");
         hEta_MCLH->Draw("e1 same");
         lEta->AddEntry(hEta_MCLH, sMCLH, "lpf");
      }

      TH1* hEta_MC = 0;
      if(drawMC){
         hEta_MC = hMC->ProjectionX("hEta_MC");
         hEta_MC->Draw("e1 same");
         lEta->AddEntry(hEta_MC, sMC, "lpf");
      }
      
      lEta->Draw();

      if(drawCompare) drawRatio(cEta, hEta_MC, hEta_LH, hEta_LHSubbed);

      if(print){
         cEta->Print("Eta.pdf", "Title:Eta");
      }
   }

   gStyle->SetPalette(kLightTerrain);

   string filename;

   TCanvas *cLH = new TCanvas("cLH", "cLH", 1280, 720); 
   cLH->SetLogy();
   TH2D* hLHclone = (TH2D*) hData->Clone("hLH"); 
   hLHclone->SetMarkerColor(kBlack);
   hLHclone->SetTitle("");
   hLHclone->GetYaxis()->SetTitle("p_{T} (GeV)");
   hLHclone->GetYaxis()->SetMoreLogLabels();
   hLHclone->SetMarkerSize(1.08);
   if(ABS_ETA) hLHclone->GetXaxis()->SetTitle("|#eta|");
   hLHclone->Draw("colz text 1.3e");  
   filename = sData;
   filename += ".pdf";
   cLH->Print(filename.c_str());

   if(drawSub){
      TCanvas *cSub = new TCanvas("cSub", "cSub", 1280, 720);
      cSub->SetLogy();
      hDataSubbed->SetMarkerColor(kBlack);
      hDataSubbed->SetTitle("");
      hDataSubbed->GetYaxis()->SetTitle("p_{T} (GeV)");
      hDataSubbed->GetYaxis()->SetMoreLogLabels();
      hDataSubbed->SetMarkerSize(1.08);
      if(ABS_ETA) hDataSubbed->GetXaxis()->SetTitle("|#eta|");
      hDataSubbed->Draw("colz text 1.3e");
      filename = sDataSubbed;
      filename += ".pdf";
      cSub->Print(filename.c_str());
   }

   if(drawMC){
      TCanvas *cMC = new TCanvas("cMC", "cMC", 1280, 720);
      cMC->SetLogy();
      TH2D* hMCclone = (TH2D*) hMC->Clone("hMC"); 
      hMCclone->SetMarkerColor(kBlack);
      if(ABS_ETA) hMCclone->GetXaxis()->SetTitle("|#eta|");
      hMCclone->SetTitle("");
      hMCclone->GetYaxis()->SetTitle("p_{T} (GeV)");
      hMCclone->SetMarkerSize(1.08);
      hMCclone->GetYaxis()->SetMoreLogLabels();
      hMCclone->Draw("colz text 1.3e");
      filename = sMC;
      filename += ".pdf";
      cMC->Print(filename.c_str());
   }

   if(drawMCLH){
      TCanvas *cMCLH = new TCanvas("cMCLH", "cMCLH", 1280, 720);
      cMCLH->SetLogy();
      TH2D* hMCLHclone = (TH2D*) hMCLH->Clone("hMCLH"); 
      hMCLHclone->SetMarkerColor(kBlack);
      hMCLHclone->SetTitle("");
      hMCLHclone->GetYaxis()->SetTitle("p_{T} (GeV)");
      hMCLHclone->SetMarkerSize(1.08);
      hMCLHclone->GetYaxis()->SetMoreLogLabels();
      if(ABS_ETA) hMCLHclone->GetXaxis()->SetTitle("|#eta|");
      hMCLHclone->Draw("colz text 1.3e");
      filename = sMCLH;
      filename += ".pdf";
      cMCLH->Print(filename.c_str());
   }
   if(toFile) outputFile->Write();

   return true;
}


bool compareDataMC(const char* dataFile, const char* dataSubbed, const char* mcFile, const char* out="."){

   TFile *fData = 0;
   fData = TFile::Open(dataFile);
   if(!fData){
      cout << "Failed to open " << dataFile << endl;
      return false;
   }

   // TFile *fDataSubbed = 0;
   // fDataSubbed = TFile::Open(dataSubbed);
   // if(!fDataSubbed){
   //    cout << "Failed to open " << dataSubbed << endl;
   //    return false;
   // } 

   TFile *fMC = 0;
   fMC = TFile::Open(mcFile);
   if(!fMC){
      cout << "Failed to open " << mcFile << endl;
      return false;
   }

   TH2* hData = 0; 
   TH2* hDataSubbed = 0;
   TH2* hMC = 0;
   TH2* hMCLH = 0;

   //====================
   hData = (TH2*) fData->Get("80.0_100.0_0.0_0.0_MC_misid");
   if(!hData){
      cout << "Failed to get hLH histogram from " << dataFile << endl;
      return false;
   }
   hData->SetName("hMCLH");

   //====================
   // hDataSubbed = (TH2*) fDataSubbed->Get("80.0_100.0_20.0_20.0_DATA_misid");
   // if(!hDataSubbed){
   //    cout << "Failed to get hLH histogram from " << dataSubbed << endl;
   //    return false;
   // }
   // hDataSubbed->SetName("hEGamma");

   // ====================
   hMC = (TH2*) fMC->Get("hMCTruthRate");
   if(!hMC){
      cout << "Failed to get hMC histogram from " << mcFile << endl;
      return false;
   }
   hMC->SetName("hMCTruthRate");

   // ====================
   // hMCLH = (TH2*) fMC->Get("80.0_100.0_20.0_20.0_DATA_misid");
   // if(!hMCLH){
   //    cout << "Failed to get hMCLH" << endl;
   //    return false;
   // }
   // hMCLH->SetName("hEGammaSubbed");


   if(!gSystem->cd(out)){
      gSystem->mkdir(out);
      gSystem->cd(out);
   }

   cout << "Calling draw" << endl;
   return draw(hData, hDataSubbed, hMC, hMCLH);   
}
