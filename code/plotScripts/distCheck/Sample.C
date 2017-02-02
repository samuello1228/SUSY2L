#include<TChain.h>
#include<TH1F.h>
#include<string>
#include<iostream>
#include<vector>
#include<TFile.h>
using std::cout;
using std::endl;
using std::vector;
using std::string;

class Sample{
 public:
   Sample(string _name="sample1"):name(_name),tag("s1_"){};
   Sample(string _name, string _tag, string _info, string _leg):name(_name),tag(_tag),info(_info),leg(_leg),mode(0),style{-1,-1,-1,-1,-1,-1,-1,-1,-1}{
   };

   std::string name;
   std::string tag;
   std::string info;
   TString leg;
   int mode;

   int style[9]; // lc, lw, lz, mc, ms, mz, fc, fs
   vector<TH1F*> hists;
   TChain* tree1;

   float weight = -1;
   TString wtExp;

   TH1F* getHistFromTree(TString var, TH1F* h1, TString cut, TString opt="", bool dress=true){
     TString hname(tag+h1->GetName());
     if(mode==0) gDirectory->Delete(hname);

     auto hx = (TH1F*)h1->Clone(hname);
     hx->Sumw2();
     tree1->Draw(var+">>"+hname,cut,opt+"goff");
     std::cout << var+">>"+hname << " " << hx->GetEntries() << std::endl;
    
     if(dress){
       if(style[0]>=0) hx->SetLineColor(style[0]);
       if(style[1]>=0) hx->SetLineStyle(style[1]);
       if(style[2]>=0) hx->SetLineWidth(style[2]);
       if(style[3]>=0) hx->SetMarkerColor(style[3]);
       if(style[4]>=0) hx->SetMarkerStyle(style[4]);
       if(style[5]>=0) hx->SetMarkerSize(style[5]);
       if(style[6]>=0) hx->SetFillColor(style[6]);
       if(style[7]>=0) hx->SetFillStyle(style[7]);
      }

     if(weight>0) hx->Scale(weight);

     hists.push_back(hx); 
     return hx;
    }
};

int main(){
  auto s1 = new Sample("test");
  std::cout << s1->info << std::endl;

  return 0;
}
