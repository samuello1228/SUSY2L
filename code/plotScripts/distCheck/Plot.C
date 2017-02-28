#include "Sample.C"
#include <iostream>
#include <vector>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
using std::cout;
using std::endl;

class Plot{

 public:
  Sample* sData;
  vector<Sample*> sSM;
  vector<Sample*> sSig;
  int mode;
  TCanvas* pCanvas;
  TPad* pPad1;
  TPad* pPad2;
  TString mCut;
  TLatex lt;
  std::string showInfo{""};

  void showPlot(TString var, TH1F* h1){

    TString canName("cav1");

    if(pCanvas) {delete pCanvas; pCanvas = nullptr;}

    if(mode == 1){
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
    TH1F* hData(nullptr);
    auto sgData = dynamic_cast< SampleGroup* >(sData);
    if(sgData && sgData->status!=0) hData = sgData->getHistFromTree(var, h1, mCut);
    else hData = sData->getHistFromTree(var, h1, mCut);

    lg->AddEntry(hData, sData->leg, "p");
    hData->Draw("E");

    bool updateData(false);

    /// Stack SM
    TH1F* htotal(nullptr);
    if(sSM.size()>0){
      lg->AddEntry(htotal, "Total SM", "l");

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
      hs->Draw("samehist");

      htotal->SetFillStyle(0);
      htotal->Draw("samehist");

      updateData = true;
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

    string x(hData->GetTitle());
    if(x.find("_logy")!=std::string::npos) pPad1->SetLogy();
    lg->Draw();

    if(showInfo!=""){
      lt.DrawLatexNDC(0.2,0.9,showInfo.c_str());
     }

    /// show ratio if needed
    if(pPad2 && htotal){
      pPad2->cd();
      auto data_copy = (TH1F*)hData->Clone("data_copy1");
      data_copy->Divide(htotal);
      data_copy->Draw("E");

      auto xxs = data_copy->GetXaxis();
      auto yxs = data_copy->GetYaxis();
      yxs->SetTitle("Data / SM");
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
      auto yxs1 = hData->GetYaxis();
      yxs1->SetLabelSize(0.13*scale);
      yxs1->SetNdivisions(506);         
      yxs1->SetTitleSize(0.14*scale);
      yxs1->SetTitleOffset(0.45/scale);
      yxs1->SetLabelOffset(0.01*scale);

    }

    pCanvas->cd();
    pCanvas->Update();
  }

  void test(){
    std::cout << sData->name << std::endl;
    TH1F* h1 = new TH1F("h1","h1",100,0,1);
    h1->Fill(0.3);
    h1->Draw();
   };
};
