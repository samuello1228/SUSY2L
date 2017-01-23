#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TChain.h>
#include <TH1F.h>

Int_t ID1;
Int_t ID2;
Double_t pt1;
Double_t pt2;
Double_t eta1;
Double_t eta2;
Double_t mll;
Double_t ptll;
Double_t MET;
Double_t mTtwo;
Double_t mt1;
Double_t mt2;
Double_t mtm;
Double_t HT;
Double_t R2;
Double_t l12_dPhi;
Double_t l12_MET_dPhi;
Int_t nJet;
Double_t jetpt;
Double_t jeteta;
Double_t jetphi;
Int_t nBJet;
Double_t bjetpt;
Double_t bjeteta;
Double_t bjetphi;
Int_t nCJet;
Double_t cjetpt;
Double_t cjeteta;
Double_t cjetphi;
Int_t nFJet;
Double_t fjetpt;
Double_t fjeteta;
Double_t fjetphi;
Double_t weight;
Double_t qFwt;
Double_t fLwt;
Double_t averageMu;
Int_t nVtx;

struct nEvent
{
    TString name;
    int n;
    double nw;
    double nAOD;
};
void skimming2(TString const& SamplePath,TString const& tag,TString const& SampleName,bool isPP1,std::vector<nEvent>& nSS)
{
    //get the "evt2l"
    TChain *tree1 = new TChain("evt2l");
    TChain *tree1P = nullptr;
    if(isPP1) tree1P = new TChain("PP1_evt2l");
    {
        TString fileName = SamplePath;
        fileName += "user.clo.";
        //fileName += "user.ychan.";
        fileName += tag;
        fileName += ".";
        fileName += SampleName;
        fileName += "_myOutput.root/*.root*";
        
        fileName = "/Users/samuel/Atlas/ntuple/test.root";
        
        cout<<fileName.Data()<<endl;
        tree1->Add(fileName.Data());
        if(isPP1) tree1P->Add(fileName.Data());
    }
    
    //channels
    std::vector<TString> channel;
    {
        TString ISR[2] = {"nonISR","ISR"};
        TString sign[2] = {"OS","SS"};
        TString lepton[3] = {"ee","mumu","emu"};
        for(int i=0;i<2;i++)
        {
            for(int j=0;j<2;j++)
            {
                for(int k=0;k<3;k++)
                {
                    TString element = "";
                    element += ISR[i];
                    element += "_";
                    element += sign[j];
                    element += "_";
                    element += lepton[k];
                    channel.push_back(element);
                }
            }
        }
    }
    
    //histograms
    TH1F* h2[channel.size()];
    for(unsigned int j=0;j<channel.size();j++)
    {
        TString hName2 = "hist_";
        hName2 += channel[j];
        
        h2[j] = new TH1F(hName2.Data(), "cut flow", 20, 0, 20);
        h2[j]->GetXaxis()->SetBinLabel(1,"AOD");
        h2[j]->GetXaxis()->SetBinLabel(2,"ntuple");
        h2[j]->GetXaxis()->SetBinLabel(3,"trigger");
        h2[j]->GetXaxis()->SetBinLabel(4,"=2SigLep");
        h2[j]->GetXaxis()->SetBinLabel(5,"fake");
        h2[j]->GetXaxis()->SetBinLabel(6,"pt1");
        h2[j]->GetXaxis()->SetBinLabel(7,"pt2");
        h2[j]->GetXaxis()->SetBinLabel(8,channel[j].Data());
    }
    
    //fill histograms
    {
        TObjArray* fileList = tree1->GetListOfFiles();
        for(int i=0;i<fileList->GetEntries();i++)
        {
            cout<<fileList->At(i)->GetTitle()<<endl;
            TFile *f1 = TFile::Open(fileList->At(i)->GetTitle());
            TH1F *h1 = (TH1F*) f1->Get("hCutFlow");
            for(unsigned int j=1;j<30;j++)
            {
                //cout<<j<<" "<<h1->GetBinContent(j)<<endl;;
            }
            for(unsigned int j=0;j<channel.size();j++)
            {
                h2[j]->Fill("AOD",h1->GetBinContent(2));
                h2[j]->Fill("ntuple",tree1->GetEntries());
            }
            delete f1;
        }
    }
    
    {
        int expN;
        TString Cut = "1";
        
        Cut += " && sig.trigCode!=0";
        expN = tree1->Draw("l12.m",Cut.Data());
        for(unsigned int j=0;j<channel.size();j++)
        {
            h2[j]->Fill("trigger",expN);
        }
        
        Cut += " && Length$(leps)==2 && (leps[0].lFlag & 2) && (leps[1].lFlag & 2)";
        expN = tree1->Draw("l12.m",Cut.Data());
        for(unsigned int j=0;j<channel.size();j++)
        {
            h2[j]->Fill("=2SigLep",expN);
        }
        
        Cut += " && ( ( int(abs(leps[0].ID)/1000) == 11 && leps[0].pt > 25 ) || ( int(abs(leps[0].ID)/1000) == 13 && leps[0].pt > 20 ) )";
        expN = tree1->Draw("l12.m",Cut.Data());
        for(unsigned int j=0;j<channel.size();j++)
        {
            h2[j]->Fill("pt1",expN);
        }
        
        Cut += " && ( ( int(abs(leps[1].ID)/1000) == 11 && leps[1].pt > 15 ) || ( int(abs(leps[1].ID)/1000) == 13 && leps[1].pt > 10 ) )";
        expN = tree1->Draw("l12.m",Cut.Data());
        for(unsigned int j=0;j<channel.size();j++)
        {
            h2[j]->Fill("pt2",expN);
        }
        
        for(unsigned int j=0;j<channel.size();j++)
        {
            TString ChannelCut = " && evt.flag == ";
            ChannelCut += TString::Itoa(j+1,10);
            ChannelCut = Cut + ChannelCut;
            //cout<<ChannelCut<<endl;
            expN = tree1->Draw("l12.m",ChannelCut.Data());
            h2[j]->Fill(channel[j].Data(),expN);
        }
        
        {
            TString SSCut = " && ( (evt.flag >=4 && evt.flag <=6) || (evt.flag >=10 && evt.flag <=12) )";
            SSCut = Cut + SSCut;
            nEvent element;
            element.name = SampleName;
            element.nAOD = h2[0]->GetBinContent(1);
            element.n = tree1->Draw("l12.m",SSCut.Data());;
            element.nw = 0;
            nSS.push_back(element);
        }
    }
    
    for(int k=1;k<=7;k++)
    {
        cout<<long(h2[1]->GetBinContent(k))<<endl;
    }
    cout<<endl;
    for(unsigned int j=0;j<channel.size();j++)
    {
        cout<<int(h2[j]->GetBinContent(8))<<endl;
    }
    cout<<endl;
    
    for(unsigned int j=0;j<channel.size();j++)
    {
        delete h2[j];
    }
    
    delete tree1;
    if(isPP1) delete tree1P;
}

void GetSampleName(std::vector<TString>& SampleName, TString const type, int const skip)
{
    ifstream fin;
    {
        TString fileName = type;
        fileName += "Sample.txt";
        fin.open(fileName.Data());
    }
    
    while(!fin.eof())
    {
        TString SampleNameTemp2;
        fin>>SampleNameTemp2;
        if(fin.eof()) break;
        SampleNameTemp2 += ".";
        
        TString SampleNameTemp;
        fin>>SampleNameTemp;
        SampleNameTemp2 += SampleNameTemp;
        SampleName.push_back(SampleNameTemp2);
        
        for(int i=1;i<=skip;i++)
        {
            double temp;
            fin>>temp;
        }
    }
    fin.close();
}

void analysis2()
{
    //TString SamplePath = "root://eosatlas//eos/atlas/user/c/clo/ntuple/";
    //TString SamplePath = "/srv/SUSY/ntuple/";
    //TString SamplePath = "/srv/SUSY/ychan/v8d6/";
    //TString SamplePath = "/afs/cern.ch/work/y/ychan/public/SUSY_NTUP/v7d11/";
    //TString SamplePath = "/afs/cern.ch/work/y/ychan/public/SUSY_NTUP/v8d6/";
    TString SamplePath = "/Users/samuel/Atlas/ntuple/";
    //TString SamplePath = "/Users/samuel/Atlas/ntuple/ychan/";
    
    //SamplePath += "AnalysisBase-02-04-17-414981/";
    //SamplePath += "AnalysisBase-02-04-17-419618/";
    //SamplePath += "AnalysisBase-02-04-17-419618-wt/";
    SamplePath += "AnalysisBase-02-04-18-f8c85e6b/";
    
    //TString tag = "v7.8";
    //TString tag = "v8.0";
    TString tag = "v8.4";
    //TString tag = "v8.6";
    //TString tag = "v7.11";
    
    std::vector<nEvent> nSS;
    
    //Data
    //if(true)
    if(false)
    {
        SamplePath += "data/";
        //tag += "b.Data";
        std::vector<TString> DataSampleName;
        DataSampleName.reserve(20);
        GetSampleName(DataSampleName,"Data",1);
        //for(unsigned int i=0;i<=1;i++)
        for(unsigned int i=0;i<DataSampleName.size();i++)
        {
            skimming2(SamplePath,tag,DataSampleName[i],true,nSS);
        }
    }
    
    //Background
    //if(true)
    if(false)
    {
        //SamplePath += "bkg/";
        //tag += ".MCBkg";
        std::vector<TString> BGSampleName;
        BGSampleName.reserve(20);
        GetSampleName(BGSampleName,"BG",4);
        //for(unsigned int i=76;i<=77;i++)
        for(unsigned int i=0;i<BGSampleName.size();i++)
        {
            skimming2(SamplePath,tag,BGSampleName[i],false,nSS);
        }
    }
    
    //Signal
    if(true)
    //if(false)
    {
        //SamplePath += "sig/";
        //tag += ".MCSig";
        std::vector<TString> SigSampleName;
        SigSampleName.reserve(20);
        GetSampleName(SigSampleName,"Sig",4);
        for(unsigned int i=0;i<=0;i++)
        //for(unsigned int i=0;i<SigSampleName.size();i++)
        {
            skimming2(SamplePath,tag,SigSampleName[i],false,nSS);
        }
    }
    
    for(unsigned int i=0;i<nSS.size();i++)
    {
        cout<<nSS[i].name<<" "<<nSS[i].n<<" "<<nSS[i].nw<<" "<<nSS[i].nAOD<<endl;
    }
}
