#include "Sample.C"
#include <iostream>
#include <vector>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TProof.h>
#include <TString.h>
#include <TLatex.h>
using std::cout;
using std::endl;

class Plot{

 public:
  Plot(TString f1="", int nProof=0){if(f1!="") m_file = new TFile(f1,"update");
    if(nProof>0){TProof *proof = TProof::Open(TString::Format("workers=%d", nProof));}
  }
  Sample* sData;
  vector<Sample*> sSM;
  vector<Sample*> sSig;
  vector<Sample*> sExtra;
  int mode;
  TCanvas* pCanvas{nullptr};
  TPad* pPad1{nullptr};
  TPad* pPad2{nullptr};
  TString mCut;
  TLatex lt;
  TFile* m_file{nullptr};
  std::string showInfo{""};

  void showPlot(TString var, TH1F* h1){

    TString canName("cav1");

    if(pCanvas) {delete pCanvas; pCanvas = nullptr;}

    if(mode>0){
      pCanvas = new TCanvas(canName,canName, 700, 600);
      float yMinP1=0.305;
      float bottomMarginP1=0.005;

      pPad1 = new TPad(Form("%s_pPad1",canName.Data()),Form("%s_pPad1",canName.Data()),0.,yMinP1,.99,1);
      pPad1->SetBottomMargin(bottomMarginP1);
      pPad1->SetFillColor(kWhite);
      pPad1->SetTickx();
      pPad1->SetTicky();
      pPad1->Draw();

      pPad2 = new TPad(Form("%s_pPad2",canName.Data()),Form("%s_pPad2",canName.Data()),0.,0.01,.99,0.295);
      pPad2->SetTopMargin(0.005);
      pPad2->SetBottomMargin(0.3);
      pPad2->SetFillColor(kWhite);
      pPad2->Draw();

    }else if(mode == 0){
//       pCanvas = new TCanvas(canName,canName, 800, 600);
      pCanvas = new TCanvas(canName);
    }else{
      cout << "Unknown mode. Do nothing" << endl;
      return;
    }

    if(pPad1) pPad1->cd();
    else pCanvas->cd();

    auto lg = new TLegend(0.6,0.7,0.9,0.88);
    lg->SetFillStyle(0);
    /// get histograms
    bool updateData(false);
    TH1F* hfirst(nullptr);
    
    TH1F* hData(nullptr);
    if(sData){
      auto sgData = dynamic_cast< SampleGroup* >(sData);
      if(sgData && sgData->status!=0) hData = sgData->getHistFromTree(var, h1, mCut);
      else hData = sData->getHistFromTree(var, h1, mCut);

      lg->AddEntry(hData, sData->leg, "p");
      hData->Draw("E");
      hfirst = hData;
     }

    /// Stack SM
    TH1F* htotal(nullptr);
    if(sSM.size()>0){

      THStack *hs = new THStack("hs","");
      for(auto& s: sSM){

        TH1F* hx(nullptr);
        SampleGroup* sg = dynamic_cast< SampleGroup* >(s);
        if(sg) hx = sg->getHistFromTree(var, h1, mCut);
        else hx = s->getHistFromTree(var, h1, mCut);

        hs->Add(hx);
        lg->AddEntry(hx, s->leg, "f");

        if(!htotal) htotal = (TH1F*)hx->Clone("htotal");
        else htotal->Add(hx);
       }
      /// plot Data and sSig
      if(hData) {
        hs->Draw("samehist");
        updateData = true;
       }else{
        hs->Draw("hist");
        hfirst = htotal;
      }

      htotal->SetFillStyle(0);
      htotal->SetLineColor(sSM[sSM.size()-1]->style[0]);
      htotal->Draw("samehist");
      lg->AddEntry(htotal, "Total SM", "l");
    }

    /// Stack SM
    if(sSig.size()>0){
      for(auto& s: sSig){
        TH1F* hx(nullptr);
        auto sg = dynamic_cast< SampleGroup* >(s);
        if(sg) hx = sg->getHistFromTree(var, h1, mCut);
        else hx = s->getHistFromTree(var, h1, mCut);

        hx->Draw("histsame");
        lg->AddEntry(hx, s->leg, "l");
       }
    }

    if(updateData){
      hData->Draw("sameE");
      hData->Draw("sameaxis");
     }

    string x(hfirst->GetTitle());
    if(x.find("_logy")!=std::string::npos) pPad1->SetLogy();
    lg->Draw();

    if(showInfo!=""){
      lt.DrawLatexNDC(0.2,0.9,showInfo.c_str());
     }

    /// show ratio if needed
    if(pPad2 && htotal){
      pPad2->cd();
      auto total_copy = (TH1F*)htotal->Clone("total_copy1");
      TString rTitle("");

      TString sameOpt="";
      TH1F* data_copy(nullptr);
      if(hData){
        data_copy = (TH1F*)hData->Clone("data_copy1");
        if(mode == 1){
          data_copy->Divide(htotal);
          rTitle = "Data / SM";

          /// MC
          for(int i=0; i<total_copy->GetNbinsX()+2; i++){
            auto x = total_copy->GetBinContent(i);
            total_copy->SetBinContent(i, 1);
            if(x==0) continue;
            total_copy->SetBinError(i, total_copy->GetBinError(i)/x);
           }
         }else if(mode==2){
          data_copy->Add(htotal, -1);
          rTitle = "Data - SM";

          /// MC
          for(int i=0; i<total_copy->GetNbinsX()+2; i++){
            total_copy->SetBinContent(i, 0);
//             total_copy->SetBinError(i, total_copy->GetBinError(i));
           }
         }

        data_copy->Draw("E");
        sameOpt = "same";
       }
      total_copy->SetLineColor(41);
      total_copy->SetMarkerSize(0);
      total_copy->SetFillStyle(1001);
      total_copy->SetFillColor(41);
      total_copy->Draw("E2"+sameOpt);

      auto xxs = total_copy->GetXaxis();
      auto yxs = total_copy->GetYaxis();
      if(data_copy){
        data_copy->Draw("Esame");
        xxs = data_copy->GetXaxis();
        yxs = data_copy->GetYaxis();
       }

      yxs->SetNoExponent();
      yxs->SetTitle(rTitle);
      yxs->SetLabelSize(0.13);
      yxs->SetNdivisions(503);         
      xxs->SetLabelSize(0.13);
      yxs->SetTitleSize(0.14);
      xxs->SetTitleSize(0.14);
//       yxs->SetTitleOffset(0.35);
      yxs->SetTitleOffset(0.45);
      xxs->SetTitleOffset(1.);
      yxs->SetLabelOffset(0.01);
      xxs->SetLabelOffset(0.03);
      xxs->SetTickLength(0.06);
      yxs->CenterTitle(); 

      pPad2->SetGridy();
      // ratio plot cosmetics
//       int firstbin = xxs->GetFirst();
//       int lastbin =  xxs->GetLast();
//       double xmax =  xxs->GetBinUpEdge(lastbin) ;
//       double xmin =  xxs->GetBinLowEdge(firstbin) ;
// 
//       TLine* l = new TLine(xmin,1.,xmax,1.);
//       TLine* l2 = new TLine(xmin,0.5,xmax,0.5);
//       TLine* l3 = new TLine(xmin,1.5,xmax,1.5);
//       TLine* l4 = new TLine(xmin,2.,xmax,2.);
//       TLine* l5 = new TLine(xmin,2.5,xmax,2.5);
//       l->SetLineWidth(1);
//       l->SetLineStyle(2);
//       l2->SetLineStyle(3);
//       l3->SetLineStyle(3);
//       l4->SetLineStyle(3);
//       l5->SetLineStyle(3);
// 
//       l->Draw("same");
//       l2->Draw("same");
//       l3->Draw("same");
//       l4->Draw("same");
//       l5->Draw("same");



      float scale = (0.295-0.001)/(1-0.305);
      auto yxs1 = hfirst->GetYaxis();
      yxs1->SetLabelSize(0.13*scale);
      yxs1->SetNdivisions(506);         
      yxs1->SetTitleSize(0.14*scale);
      yxs1->SetTitleOffset(0.45/scale);
      yxs1->SetLabelOffset(0.01*scale);

      pPad2->Update();
    }

    pCanvas->cd();
    pCanvas->Update();
  }

  Sample* getSample(TString name){
    Sample* s1 = (Sample*)m_file->Get("Sample_"+name);
    if(s1) s1->tree1 = (TChain*)m_file->Get(name+"_tree1");
    else cout << "Sample " << name << " is not found!!!" << endl;
    return s1;
   };

  SampleGroup* getSampleGroup(TString name){
    SampleGroup* s1 = (SampleGroup*)m_file->Get("Sample_"+name);
    s1->tree1 = (TChain*)m_file->Get(name+"_tree1");
//     for(size_t i=0; i<s1->m_sampleNames.size(); i++) s1->sampleList[i] = getSample(s1->m_sampleNames[i]);
    for(size_t i=0; i<s1->m_sampleNames.size(); i++){
      std::cout << s1->m_sampleNames[i] << std::endl;
      s1->sampleList[i] = getSample(s1->m_sampleNames[i]);
    }
    return s1;
   };

   void loadSampleGroup(SampleGroup* s1){
//     for(size_t i=0; i<s1->m_sampleNames.size(); i++) s1->sampleList[i] = getSample(s1->m_sampleNames[i]);
    s1->sampleList.resize(s1->m_sampleNames.size(), nullptr);
    for(size_t i=0; i<s1->m_sampleNames.size(); i++){
      std::cout << s1->m_sampleNames[i] << std::endl;
      s1->sampleList[i] = getSample(s1->m_sampleNames[i]);
     }
    return;
   };

  void saveSamples(){
    sData->writeToFile(m_file);
    for(auto x: sSM) saveSample(x);
    for(auto x: sSig) saveSample(x);
    for(auto x: sExtra) saveSample(x);
   }
  void saveSample(Sample* s){
    auto s1 = dynamic_cast< SampleGroup* >(s);
    if(s1) s1->writeToFile(m_file);
    else s->writeToFile(m_file);
   }

  void useEntryList(TString elist, bool includeExtra=false, TString cutx="", bool remake=false){
    Info("useEntryList", "Setting up entrylist: %s", elist.Data());
    if(remake && m_file){
      Info("useEntryList", "Deleting exist entrylists");
      m_file->Delete("EL_"+elist+"_*;*");
     }

    if(sData) useEntryList(sData, elist,cutx);
    for(auto x: sSM) useEntryList(x, elist, cutx);
    for(auto x: sSig) useEntryList(x, elist, cutx);
    if(includeExtra){for(auto x: sExtra) useEntryList(x, elist, cutx);}
   }

  void useEntryList(Sample* s1, TString elist, TString cutx){
    auto sx = dynamic_cast< SampleGroup* >(s1);
    if(sx){
      for(auto s2: sx->sampleList) useEntryList(s2, elist, cutx);
     }

    if(s1->tree1){
      /// try to get the entry list
      TString elistname = "EL_"+elist+"_"+s1->name;
      auto elist1 = (TEntryList*)m_file->Get(elistname);
      if(!elist1 && cutx!=""){
        /// make the entry list
        Info("useEntryList", "entrylist %s not exist for sample %s, creating it with cut: %s", elist.Data(), s1->name.c_str(), cutx.Data());
        s1->tree1->SetEntryList(0);
        s1->tree1->Draw(">>"+elistname, cutx, "entrylist");
        elist1 = (TEntryList*)gDirectory->Get(elistname);
        if(m_file){
          m_file->cd();
          elist1->Write(elistname);
         }
       }
      s1->tree1->SetEntryList((TEntryList*)m_file->Get("EL_"+elist+"_"+s1->name));
    }
   }

  void test(){
    std::cout << sData->name << std::endl;
    TH1F* h1 = new TH1F("h1","h1",100,0,1);
    h1->Fill(0.3);
    h1->Draw();
   };
};
