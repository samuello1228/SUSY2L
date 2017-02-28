#include<TChain.h>
#include <TError.h>
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
   Sample(string _name, string _tag, string _info, string _leg):name(_name),tag(_tag),info(_info),leg(_leg),mode(0),style{-1,-1,-1,-1,-1,-1,-1,-1,-1},tree1(nullptr){
   };
   virtual ~Sample(){};

   std::string name;
   std::string tag;
   std::string info;
   TString sCuts{""};
   TString leg;
   int mode;

   int style[9]; // lc, lw, lz, mc, ms, mz, fc, fs
   TChain* tree1;
   vector<TH1F*> hists;

   float weight = -1;
   TString wtExp;

   float getStatWeight(TString key="aSumW"){
    float total(0);
    for (auto fs: *(tree1->GetListOfFiles())){
      TFile f(fs->GetTitle(),"read");
      TH1D* h1 = (TH1D*)f.Get("hCutFlow");
      total += h1->GetBinContent(h1->GetXaxis()->FindBin(key));
     }
    return total;
   }

   TH1F* getHistFromTree(TString var, TH1F* h1, TString cut, TString opt="", bool dress=true){
     if(sCuts != ""){
       if(cut != "") cut += "&&";
       cut += sCuts;
      }

     TString hname(tag+h1->GetName());
     if(mode==0) gDirectory->Delete(hname);

     auto hx = (TH1F*)h1->Clone(hname);
     hx->Sumw2();
     tree1->Draw(var+">>"+hname,cut,opt+"goff");
     std::cout << var+">>"+hname << " " << hx->GetEntries() << std::endl;
    
     if(dress) dressHist(hx);
     if(weight>0) hx->Scale(weight);

     hists.push_back(hx); 
     return hx;
    }

 protected:
   void dressHist(TH1* hx){
       if(style[0]>=0) hx->SetLineColor(style[0]);
       if(style[1]>=0) hx->SetLineStyle(style[1]);
       if(style[2]>=0) hx->SetLineWidth(style[2]);
       if(style[3]>=0) hx->SetMarkerColor(style[3]);
       if(style[4]>=0) hx->SetMarkerStyle(style[4]);
       if(style[5]>=0) hx->SetMarkerSize(style[5]);
       if(style[6]>=0) hx->SetFillColor(style[6]);
       if(style[7]>=0) hx->SetFillStyle(style[7]);
    }
};

class SampleGroup: public Sample{
 public:
   SampleGroup(string _name="sample1"):Sample(_name),status(-1){};
   SampleGroup(string _name, string _tag, string _info, string _leg):Sample(_name, _tag, _info, _leg),status(-1){
   };
   int status;
   std::vector< Sample* > sampleList;

   TH1F* getHistFromTree(TString var, TH1F* h1, TString cut, TString opt="", bool dress=true){
     TString hname(tag+h1->GetName());
     if(mode==0) gDirectory->Delete(hname);

     auto hx = (TH1F*)h1->Clone(hname);
     for(auto& s: sampleList){
       hx->Add(s->getHistFromTree(var, h1, cut, opt, false));
      }
    
     if(dress) dressHist(hx);
     if(weight>0) hx->Scale(weight);

     hists.push_back(hx); 
     return hx;
    }

   void setUpOwnChain(TChain* ch1=nullptr,std::string chainName="evt2l"){
     Info("setUpOwnChain", "sample group: %s", name.c_str());
     if(ch1) tree1 = ch1;
     else if(tree1){
       tree1 = new TChain(chainName.c_str());
       Info("setUpOwnChain", "creating new chain:%s", chainName.c_str());
      }

     for(auto& s: sampleList){
       tree1->Add(s->tree1);
     Info("setUpOwnChain", "adding %lld, now %lld (+%s)", s->tree1->GetEntries(), tree1->GetEntries(), s->name.c_str());
      }
     status = 0;
    }
 };

int main(){
  auto s1 = new Sample("test");
  std::cout << s1->info << std::endl;

  return 0;
}
