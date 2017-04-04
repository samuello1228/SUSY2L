#ifndef Sample_C
#define Sample_C
#include<TChain.h>
#include <TError.h>
#include<TH1F.h>
#include<TEntryList.h>
#include<string>
#include<iostream>
#include<vector>
#include<TFile.h>
using std::cout;
using std::endl;
using std::vector;
using std::string;

class Sample: public TObject{
 public:
   Sample(string _name="sample1"):name(_name),tag("s1_"){};
   Sample(string _name, string _tag, string _info, string _leg):name(_name),tag(_tag),info(_info),leg(_leg),mode(0),style{-1,-1,-1,-1,-1,-1,-1,-1,-1},tree1(nullptr){
   };
   virtual ~Sample(){
     delete m_file;
    };

   string name;
   string tag;
   string info;
   TString sCuts{""};
   TString leg;
   int mode;

   int style[9]; // lc, lw, lz, mc, ms, mz, fc, fs
   TChain* tree1; //!
   vector<TH1F*> hists;
   TFile* m_file{nullptr}; //!

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

   void useEntryList(TString elistName, TString cut=""){
     /// find the elistName in list, or root file. If not found, create it and save it.
     auto elist = (TEntryList*)m_file->Get(elistName);
     if(!elist && cut!=""){
       tree1->Draw(">>"+elistName,cut,"entrylist");
       m_file->cd();
       elist->Write(elistName);
     }
     tree1->SetEntryList(elist);
    }

   void setupFromFile(TString& filename){
     m_file = new TFile(filename, "update");
     tree1 = (TChain*) m_file->Get("tree1");
    }

   void writeToFile(TFile* f1){
     Info("Sample", "%s is wittern out", name.c_str());
     f1->cd();
     Info("Sample", "get into the files now [%s]", name.c_str());
     this->Write(("Sample_"+name).c_str());
     if(tree1) tree1->Write((name+"_tree1").c_str());
     Info("Sample", "Done [%s]", name.c_str());
    }

   void writeToFile(string dir1="", string filename=""){
     if(!m_file){
       if (filename=="") filename = "Sample_"+name+".root";
       m_file = new TFile((dir1+filename).c_str(),"recreate");
      }
     m_file->cd();
     this->Write(name.c_str());
     tree1->Write("tree1");
    }

 protected:
   void dressHist(TH1* hx){
     cout << name << " " << style[0] << endl;
       if(style[0]>=0) hx->SetLineColor(style[0]);
       if(style[1]>=0) hx->SetLineStyle(style[1]);
       if(style[2]>=0) hx->SetLineWidth(style[2]);
       if(style[3]>=0) hx->SetMarkerColor(style[3]);
       if(style[4]>=0) hx->SetMarkerStyle(style[4]);
       if(style[5]>=0) hx->SetMarkerSize(style[5]);
       if(style[6]>=0) hx->SetFillColor(style[6]);
       if(style[7]>=0) hx->SetFillStyle(style[7]);
    }

   ClassDef(Sample, 1)
};

class SampleGroup: public Sample{
 public:
   SampleGroup(string _name="sample1"):Sample(_name),status(-1){};
   SampleGroup(string _name, string _tag, string _info, string _leg):Sample(_name, _tag, _info, _leg),status(-1){
   };
   int status;
   std::vector< Sample* > sampleList;
   std::vector< string > m_sampleNames;

   float getStatWeight(TString key="aSumW"){
    float total(0);
    for (auto s: sampleList){
      total += s->getStatWeight(key);
     }
    return total;
   }

   TH1F* getHistFromTree(TString var, TH1F* h1, TString cut, TString opt="", bool dress=true){
     if(sCuts != ""){
       if(cut != "") cut += "&&";
       cut += sCuts;
      }
     cout << cut << endl;

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

   void setUpOwnChain(TChain* ch1=nullptr, string chainName="evt2l"){
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

   void writeToFile(TFile* f1){
     Info("Sample", "%s is wittern out", name.c_str());
     f1->cd();
     m_sampleNames.clear();
     m_sampleNames.reserve(sampleList.size()); 
     for(auto s: sampleList){
       s->writeToFile(f1);
       m_sampleNames.push_back(s->name);
      }
     this->Write(("Sample_"+name).c_str());
     if(tree1) tree1->Write((name+"_tree1").c_str());
    }


 private:

   ClassDef(SampleGroup, 2)
 };

int main(){
  auto s1 = new Sample("test");
  std::cout << s1->info << std::endl;

  return 0;
}
#endif
