#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>

#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TString.h>
#include "AtlasLabels.C"
#include "AtlasStyle.C"

#include <TList.h>
#include <Fit/Fitter.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Math/WrappedMultiTF1.h>
#include <HFitInterface.h>

#include "RooStats/NumberCountingUtils.h"
#include "RooStats/RooStatsUtils.h"

bool dorw = 1;
bool simple = 1;
bool combined = 0;

bool doOptimize = 0;

//for Zpt reweighting
Double_t ptll;
Double_t rw;

TBranch* b_ptll;

//for expN
Double_t weight;
TBranch* b_weight;

//for charge filp BG
Double_t pt1;
Double_t pt2;
Double_t eta1;
Double_t eta2;
Double_t cfw;

TBranch* b_pt1;
TBranch* b_pt2;
TBranch* b_eta1;
TBranch* b_eta2;

void initializeTree1(std::vector< std::vector<TChain*> >& tree, std::vector<TString>& SampleID, std::vector<TString>& channel)
{
    for(unsigned int i=0;i<SampleID.size();i++)
    {
        for(unsigned int ChannelIndex=0;ChannelIndex<channel.size();ChannelIndex++)
        {
            TChain* treeTemp = new TChain("tree");
            
            TString fileName = "skimming/skimming.";
            fileName += SampleID[i];
            fileName += "_";
            fileName += channel[ChannelIndex];
            fileName += ".root";
            //cout<<fileName<<endl;
            
            treeTemp->Add(fileName.Data());
            tree[ChannelIndex].push_back(treeTemp);
            
            tree[ChannelIndex][i]->SetBranchAddress("ptll", &ptll, &b_ptll);
            tree[ChannelIndex][i]->SetBranchAddress("weight", &weight, &b_weight);
            
            tree[ChannelIndex][i]->SetBranchAddress("pt1", &pt1, &b_pt1);
            tree[ChannelIndex][i]->SetBranchAddress("pt2", &pt2, &b_pt2);
            tree[ChannelIndex][i]->SetBranchAddress("eta1", &eta1, &b_eta1);
            tree[ChannelIndex][i]->SetBranchAddress("eta2", &eta2, &b_eta2);
        }
    }
}

void initializeTree2(std::vector<TChain*>& tree2,std::vector<unsigned int>& SetOfChannel, std::vector<TString>& SampleID, std::vector<TString>& channel)
{
    for(unsigned int i=0;i<SampleID.size();i++)
    {
        TChain* treeTemp = new TChain("tree");
        for(unsigned int j=0;j<SetOfChannel.size();j++)
        {
            TString fileName = "skimming/skimming.";
            fileName += SampleID[i];
            fileName += "_";
            fileName += channel[SetOfChannel[j]];
            fileName += ".root";
            
            treeTemp->Add(fileName.Data());
        }
        
        treeTemp->SetBranchAddress("weight", &weight, &b_weight);
        tree2.push_back(treeTemp);
    }
}

struct GroupData
{
    TString GroupName;
    TString LegendName;
    TString LatexName;
    unsigned int lower;
    unsigned int upper;
};

struct Group
{
    GroupData* info;
    TH1F* h2;
};

bool compare2(Group Group1,Group Group2)
{
    return Group1.h2->GetSumOfWeights() < Group2.h2->GetSumOfWeights();
}

// definition of shared parameter
// mumu channel
int iparM[3] = { 0, // pol(2)
                 1, //
                 2  //
};

// ee channel
int iparE[4] = { 3, // pol0
                 0, // pol2
                 1, //
                 2  //
};

struct GlobalChi2
{
    const  ROOT::Math::IMultiGenFunction*  fChi2_1;
    const  ROOT::Math::IMultiGenFunction*  fChi2_2;
    
    GlobalChi2(ROOT::Math::IMultiGenFunction & f1, ROOT::Math::IMultiGenFunction & f2) :
    fChi2_1(&f1), fChi2_2(&f2) {}
    
    double operator() (const double* par) const
    {
        double p1[3];
        for(int i=0;i<3;++i)
        {
            p1[i] = par[iparM[i]];
        }
        
        double p2[4];
        for(int i=0;i<4;++i)
        {
            p2[i] = par[iparE[i]];
        }
        return (*fChi2_1)(p1) + (*fChi2_2)(p2);
    }
};

void analysis1()
{
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
    
    //For Data
    std::vector<TString> DataSampleID;
    DataSampleID.reserve(20);
    double sumDataL = 0; //cross section in pb-1
    
    {
        //read DataSample.txt
        ifstream fin;
        fin.open("DataSample.txt");
        while(!fin.eof())
        {
            TString SampleIDTemp;
            fin>>SampleIDTemp;
            if(fin.eof()) break;
            
            TString SampleNameTemp;
            fin>>SampleNameTemp;
            
            SampleIDTemp += ".";
            SampleIDTemp += SampleNameTemp;
            DataSampleID.push_back(SampleIDTemp);
            
            double DataLTemp;
            fin>>DataLTemp;
            sumDataL += DataLTemp;
        }
        fin.close();
    }
    cout<<"Total Luminosity: "<<sumDataL<<endl;
    
    /*
    //tree
    std::vector< std::vector<TChain*> > tree1Data;
    for(unsigned int ChannelIndex=0;ChannelIndex<channel.size();ChannelIndex++)
    {
        std::vector<TChain*> element;
        tree1Data.push_back(element);
    }
    initializeTree1(tree1Data,DataSampleID,channel);
    */
    
    
    //For BGMC
    std::vector<TString> BGMCSampleID;
    std::vector<double> BGMCXS; //cross section in pb
    BGMCSampleID.reserve(20);
    BGMCXS.reserve(20);
    
    {
        //read BGSample.txt
        ifstream fin;
        fin.open("BGSample.txt");
        while(!fin.eof())
        {
            TString SampleIDTemp;
            fin>>SampleIDTemp;
            if(fin.eof()) break;
            
            TString SampleNameTemp;
            fin>>SampleNameTemp;
            
            SampleIDTemp += ".";
            SampleIDTemp += SampleNameTemp;
            BGMCSampleID.push_back(SampleIDTemp);
            
            double BGMCXSTemp;
            double BGMCXSTemp2;
            fin>>BGMCXSTemp2;
            fin>>BGMCXSTemp;
            BGMCXSTemp2 *= BGMCXSTemp;
            fin>>BGMCXSTemp;
            BGMCXSTemp2 *= BGMCXSTemp;
            
            BGMCXS.push_back(BGMCXSTemp2);
            
            fin>>BGMCXSTemp;
        }
        fin.close();
    }
    cout<<"Total number of BG files: "<<BGMCSampleID.size()<<endl;
    
    //Get number of events in AOD
    std::vector<unsigned int> BGMCnAOD;
    for(unsigned int i=0;i<BGMCSampleID.size();i++)
    {
        TString NameTemp = "skimming/skimming.";
        NameTemp += BGMCSampleID[i];
        cout<<NameTemp<<": ";
        NameTemp += "_";
        NameTemp += channel[0];
        NameTemp += ".root";
        
        TFile* file = new TFile(NameTemp.Data(),"READ");
        TH1F *h1 = (TH1F*) file->Get("hist");
        unsigned int nAOD = h1->GetBinContent(1);
        cout<<nAOD<<endl;
        BGMCnAOD.push_back(nAOD);
        
        delete file;
    }
    
    //tree
    std::vector< std::vector<TChain*> > tree1BGMC;
    for(unsigned int ChannelIndex=0;ChannelIndex<channel.size();ChannelIndex++)
    {
        std::vector<TChain*> element;
        tree1BGMC.push_back(element);
    }
    initializeTree1(tree1BGMC,BGMCSampleID,channel);
    
    
    //For Signal MC
    struct SigInfo
    {
        int MassDiff;
        int ID;
        int colour;
    };
    std::vector<SigInfo> SigMassSplitting;
    {
        SigInfo element;
        
        element.MassDiff = 50;    element.ID = 18;  element.colour = 9;  SigMassSplitting.push_back(element);
        element.MassDiff = 100;   element.ID = 19;  element.colour = 11;  SigMassSplitting.push_back(element);
        element.MassDiff = 200;   element.ID = 20;  element.colour = 12;  SigMassSplitting.push_back(element);
    }
    
    std::vector<TString> SigSampleID;
    std::vector<double> SigXS; //cross section in pb
    std::vector<int> SigMass1;
    std::vector<int> SigMass2;
    SigSampleID.reserve(20);
    SigXS.reserve(20);
    SigMass1.reserve(20);
    SigMass2.reserve(20);
    
    {
        //read SigSample.txt
        ifstream fin;
        fin.open("SigSample.txt");
        while(!fin.eof())
        {
            TString SampleIDTemp;
            //for 125
            fin>>SampleIDTemp;
            if(fin.eof()) break;
            
            TString SampleNameTemp;
            fin>>SampleNameTemp;
            
            SampleIDTemp += ".";
            SampleIDTemp += SampleNameTemp;
            SigSampleID.push_back(SampleIDTemp);
            
            int SigMass;
            fin>>SigMass;
            SigMass1.push_back(SigMass);
            fin>>SigMass;
            SigMass2.push_back(SigMass);
            
            double SigXSTemp2;
            fin>>SigXSTemp2;
            
            double SigXSTemp;
            fin>>SigXSTemp;
            SigXSTemp2 *= SigXSTemp;

            /*
            fin>>SigXSTemp2;
            fin>>SigXSTemp;
            SigXSTemp2 *= SigXSTemp;
            fin>>SigXSTemp;
            SigXSTemp2 *= SigXSTemp;
            
            fin>>SigXSTemp;
            fin>>SigXSTemp;
            
            //next line for 127
            fin>>SigXSTemp;
            
            double SigXSTemp3;
            fin>>SigXSTemp3;
            fin>>SigXSTemp;
            SigXSTemp3 *= SigXSTemp;
            fin>>SigXSTemp;
            SigXSTemp3 *= SigXSTemp;
            
            fin>>SigXSTemp;
            fin>>SigXSTemp;
            
            //125 + 127
            SigXSTemp2 += SigXSTemp3;
            */
            
            //cout<<SigXSTemp2<<endl;
            SigXS.push_back(SigXSTemp2);
        }
        fin.close();
    }
    
    for(unsigned int i=0;i<SigMassSplitting.size();i++)
    {
        cout<<"Mass splitting: "<<SigMassSplitting[i].MassDiff<<endl;
        cout<<"index:"<<SigMassSplitting[i].ID<<endl;
        cout<<"Name: "<<SigSampleID[SigMassSplitting[i].ID].Data()<<endl;
        cout<<"MassDiff: "<<SigMass1[SigMassSplitting[i].ID]-SigMass2[SigMassSplitting[i].ID]<<endl;
        cout<<"XS: "<<SigXS[SigMassSplitting[i].ID]<<endl<<endl;
        
        cout<<"All samples with the same mass splitting "<<SigMassSplitting[i].MassDiff<<" GeV:"<<endl;
        for(unsigned int j=0;j<SigSampleID.size();j++)
        {
            if(SigMass1[j]-SigMass2[j] == SigMassSplitting[i].MassDiff) cout<<"index:"<<j<<" "<<SigSampleID[j]<<endl;
        }
        cout<<endl;
    }
    
    //Get number of events in AOD
    std::vector<unsigned int> SignAOD;
    for(unsigned int i=0;i<SigSampleID.size();i++)
    {
        TString NameTemp = "skimming/skimming.";
        NameTemp += SigSampleID[i];
        cout<<NameTemp<<": ";
        NameTemp += "_";
        NameTemp += channel[0];
        NameTemp += ".root";
        
        TFile* file = new TFile(NameTemp.Data(),"READ");
        TH1F *h1 = (TH1F*) file->Get("hist");
        unsigned int nAOD = h1->GetBinContent(1);
        cout<<nAOD<<endl;
        SignAOD.push_back(nAOD);
        
        delete file;
    }
    
    /*
    //tree
    std::vector< std::vector<TChain*> > tree1Sig;
    for(unsigned int ChannelIndex=0;ChannelIndex<channel.size();ChannelIndex++)
    {
        std::vector<TChain*> element;
        tree1Sig.push_back(element);
    }
    initializeTree1(tree1Sig,SigSampleID,channel);
    */
    
    //Group for MC background
    std::vector<GroupData> BGMCGroupData;
    {
        GroupData element;
        element.GroupName = "Zee"; element.LegendName = "Z#rightarrow ee"; element.LatexName = "Z$\\rightarrow ee$";
        element.lower = 48;  element.upper = 71; BGMCGroupData.push_back(element);
        
        element.GroupName = "Zmumu"; element.LegendName = "Z#rightarrow #mu#mu"; element.LatexName = "Z$\\rightarrow\\mu\\mu$";
        element.lower = 24;  element.upper = 47; BGMCGroupData.push_back(element);
        
        element.GroupName = "Ztautau"; element.LegendName = "Z#rightarrow #tau#tau"; element.LatexName = "Z$\\rightarrow\\tau\\tau$";
        element.lower = 0;   element.upper = 23; BGMCGroupData.push_back(element);
        
        element.GroupName = "ttbar"; element.LegendName = "t#bar{t}"; element.LatexName = "$t\\bar{t}$";
        element.lower = 72;  element.upper = 72; BGMCGroupData.push_back(element);
        
        element.GroupName = "Wt"; element.LegendName = "Wt"; element.LatexName = "Wt";
        element.lower = 73;  element.upper = 74; BGMCGroupData.push_back(element);
        
        element.GroupName = "VV"; element.LegendName = "VV"; element.LatexName = "VV";
        element.lower = 75;  element.upper = 88; BGMCGroupData.push_back(element);
        
        element.GroupName = "Vgamma"; element.LegendName = "V + #gamma"; element.LatexName = "V$+\\gamma$";
        element.lower = 89;  element.upper = 108;BGMCGroupData.push_back(element);
    }
    
    //Group for data-driven background
    std::vector<GroupData> BGDataGroupData;
    {
        GroupData element;
        element.GroupName = "charge flip"; element.LegendName = element.GroupName; element.LatexName = element.GroupName;
        element.lower = 0;  element.upper = 0; BGDataGroupData.push_back(element);
        
        element.GroupName = "fake lepton"; element.LegendName = element.GroupName; element.LatexName = element.GroupName;
        element.lower = 0;  element.upper = 0; BGDataGroupData.push_back(element);
    }
    
    //plotting
    //Variables for plotting
    struct VarData
    {
        TString VarName;
        TString VarTitle;
        TString unit;
        TString latexName;
        
        int bin;
        double xmin;
        double xmax;
        
        bool log;
        double ymin;
        double ymax;
    };
    
    std::vector<VarData> Var;
    {
        VarData element;
        
        element.VarName = "pt1";        element.VarTitle = "pt of the leading lepton";          element.unit = "[GeV]";
        element.bin=40;         element.xmin=25;                element.xmax=250;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\text{p}_{\\text{T}}$ of the leading lepton";
        Var.push_back(element);
        
        element.VarName = "pt2";        element.VarTitle = "pt of the subleading lepton";       element.unit = "[GeV]";
        element.bin=40;         element.xmin=20;                element.xmax=250;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\text{p}_{\\text{T}}$ of the subleading lepton";
        Var.push_back(element);
        
        element.VarName = "eta1";       element.VarTitle = "eta of the leading lepton";         element.unit = "";
        element.bin=40;         element.xmin=-3;                element.xmax=3;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\eta$ of the leading lepton";
        Var.push_back(element);
        
        element.VarName = "eta2";       element.VarTitle = "eta of the subleading lepton";      element.unit = "";
        element.bin=40;         element.xmin=-3;                element.xmax=3;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\eta$ of the subleading lepton";
        Var.push_back(element);
        
        element.VarName = "mll";        element.VarTitle = "Dilepton invariant mass";           element.unit = "[GeV]";
        element.bin=40;         element.xmin=60;                element.xmax=250;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = element.VarTitle;
        Var.push_back(element);
        
        element.VarName = "ptll";       element.VarTitle = "Dilepton pt";                       element.unit = "[GeV]";
        element.bin=40;         element.xmin=0;                 element.xmax=200;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "Dilepton $\\text{p}_{\\text{T}}$";
        Var.push_back(element);
        
        element.VarName = "MET";        element.VarTitle = "MET";                               element.unit = "[GeV]";
        element.bin=40;         element.xmin=0;                 element.xmax=200;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = element.VarTitle;
        Var.push_back(element);
        
        element.VarName = "mTtwo";      element.VarTitle = "mT2";                               element.unit = "[GeV]";
        element.bin=40;         element.xmin=0;                 element.xmax=150;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = element.VarTitle;
        Var.push_back(element);
        
        element.VarName = "averageMu";  element.VarTitle = "averageMu";                         element.unit = "";
        element.bin=35;         element.xmin=0;                 element.xmax=35;
        element.log=0;          element.ymin=0;                 element.ymax=0;
        element.latexName = "Average number of interactions per bunch crossing";
        Var.push_back(element);
        
        element.VarName = "nVtx";       element.VarTitle = "Number of vertices";                element.unit = "";
        element.bin=30;         element.xmin=0;                 element.xmax=30;
        element.log=0;          element.ymin=0;                 element.ymax=0;
        element.latexName = element.VarTitle;
        Var.push_back(element);
        
        element.VarName = "mt1";        element.VarTitle = "mT of the leading lepton";          element.unit = "[GeV]";
        element.bin=40;         element.xmin=0;                 element.xmax=300;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\text{m}_{\\text{T}}$ of the leading lepton";
        Var.push_back(element);
        
        element.VarName = "mt2";        element.VarTitle = "mT of the subleading lepton";       element.unit = "[GeV]";
        element.bin=40;         element.xmin=0;                 element.xmax=300;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\text{m}_{\\text{T}}$ of the subleading lepton";
        Var.push_back(element);
        
        element.VarName = "nJet";       element.VarTitle = "Number of jets";                    element.unit = "";
        element.bin=15;         element.xmin=0;                 element.xmax=15;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = element.VarTitle;
        Var.push_back(element);
        
        element.VarName = "jetpt";      element.VarTitle = "pT of the leading jet";             element.unit = "[GeV]";
        element.bin=40;         element.xmin=20;                element.xmax=300;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\text{p}_{\\text{T}}$ of the leading jet";
        Var.push_back(element);
        
        element.VarName = "jeteta";     element.VarTitle = "eta of the leading jet";            element.unit = "";
        element.bin=40;         element.xmin=-3;                element.xmax=3;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\eta$ of the leading jet";
        Var.push_back(element);
        
        element.VarName = "jetphi";     element.VarTitle = "phi of the leading jet";            element.unit = "";
        element.bin=40;         element.xmin=-TMath::Pi();      element.xmax=TMath::Pi();
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\phi$ of the leading jet";
        Var.push_back(element);
        
        element.VarName = "nBJet";      element.VarTitle = "Number of b-jets";                  element.unit = "";
        element.bin=8;          element.xmin=0;                 element.xmax=8;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = element.VarTitle;
        Var.push_back(element);
        
        element.VarName = "bjetpt";     element.VarTitle = "pT of the leading b-jet";           element.unit = "[GeV]";
        element.bin=40;         element.xmin=20;                element.xmax=300;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\text{p}_{\\text{T}}$ of the leading b-jet";
        Var.push_back(element);
        
        element.VarName = "bjeteta";    element.VarTitle = "eta of the leading b-jet";          element.unit = "";
        element.bin=40;         element.xmin=-3;                element.xmax=3;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\eta$ of the leading b-jet";
        Var.push_back(element);
        
        element.VarName = "bjetphi";    element.VarTitle = "phi of the leading b-jet";          element.unit = "";
        element.bin=40;         element.xmin=-TMath::Pi();      element.xmax=TMath::Pi();
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\phi$ of the leading b-jet";
        Var.push_back(element);
        
        element.VarName = "nCJet";      element.VarTitle = "Number of central light jets";      element.unit = "";
        element.bin=15;         element.xmin=0;                 element.xmax=15;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = element.VarTitle;
        Var.push_back(element);
        
        element.VarName = "cjetpt";     element.VarTitle = "pT of the leading central light jets";                              element.unit = "[GeV]";
        element.bin=40;         element.xmin=20;                element.xmax=300;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\text{p}_{\\text{T}}$ of the leading central light jet";
        Var.push_back(element);
        
        element.VarName = "cjeteta";    element.VarTitle = "eta of the leading central light jets";                             element.unit = "";
        element.bin=40;         element.xmin=-3;                element.xmax=3;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\eta$ of the leading central light jet";
        Var.push_back(element);
        
        element.VarName = "cjetphi";    element.VarTitle = "phi of the leading central light jets";                             element.unit = "";
        element.bin=40;         element.xmin=-TMath::Pi();      element.xmax=TMath::Pi();
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\phi$ of the leading central light jet";
        Var.push_back(element);
        
        element.VarName = "nFJet";      element.VarTitle = "Number of forward jets";            element.unit = "";
        element.bin=8;          element.xmin=0;                 element.xmax=8;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = element.VarTitle;
        Var.push_back(element);
        
        element.VarName = "fjetpt";     element.VarTitle = "pT of the leading forward jet";     element.unit = "[GeV]";
        element.bin=40;         element.xmin=30;                element.xmax=300;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\text{p}_{\\text{T}}$ of the leading forward jet";
        Var.push_back(element);
        
        element.VarName = "fjeteta";    element.VarTitle = "eta of the leading forward jet";    element.unit = "";
        element.bin=40;         element.xmin=-3;                element.xmax=3;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\eta$ of the leading forward jet";
        Var.push_back(element);
        
        element.VarName = "fjetphi";    element.VarTitle = "phi of the leading forward jet";    element.unit = "";
        element.bin=40;         element.xmin=-TMath::Pi();      element.xmax=TMath::Pi();
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\phi$ of the leading forward jet";
        Var.push_back(element);
        
        element.VarName = "HT";         element.VarTitle = "HT";                                element.unit = "[GeV]";
        element.bin=40;         element.xmin=65;                 element.xmax=500;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = element.VarTitle;
        Var.push_back(element);
        
        element.VarName = "R2";         element.VarTitle = "R2";                                element.unit = "";
        element.bin=40;         element.xmin=0;                 element.xmax=1;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = element.VarTitle;
        Var.push_back(element);
        
        element.VarName = "l12_dPhi";    element.VarTitle = "phi difference between the two leptons";                           element.unit = "";
        element.bin=40;         element.xmin=-TMath::Pi();      element.xmax=TMath::Pi();
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\phi$ difference between the two leptons";
        Var.push_back(element);
        
        element.VarName = "l12_MET_dPhi";element.VarTitle = "phi difference between l12 and MET";        element.unit = "";
        element.bin=40;         element.xmin=-TMath::Pi();      element.xmax=TMath::Pi();
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\phi$ difference between l12 and MET";
        Var.push_back(element);
    }
    SetAtlasStyle();
    
    //if(dorw)
    if(false)
    {
        //Z pt reweighting
        int VarIndex=5;
        std::vector<unsigned int> ChannelIndex_ll[2];
        ChannelIndex_ll[0].push_back(0);
        ChannelIndex_ll[0].push_back(6);
        ChannelIndex_ll[1].push_back(1);
        ChannelIndex_ll[1].push_back(7);

        //calculate ratio plot
        TH1F* h2Ratio_rw[2];
        for(int LeptonIndex=0;LeptonIndex<=1;LeptonIndex++)
        {
            std::vector<TChain*> tree2Data;
            initializeTree2(tree2Data,ChannelIndex_ll[LeptonIndex],DataSampleID,channel);
            std::vector<TChain*> tree2BGMC;
            initializeTree2(tree2BGMC,ChannelIndex_ll[LeptonIndex],BGMCSampleID,channel);

            //initialize histograms
            TString title;
            title += Var[VarIndex].VarTitle;
            
            TString xaxis;
            xaxis += Var[VarIndex].VarTitle;
            xaxis += " ";
            xaxis += Var[VarIndex].unit;
            
            //h2Data
            TH1F* h2Data[DataSampleID.size()];
            std::vector<TString> hName2Data;
            for(unsigned int j=0;j<DataSampleID.size();j++)
            {
                TString NameTemp = "Data_";
                NameTemp += TString::Itoa(j,10);
                h2Data[j] = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                hName2Data.push_back(NameTemp);
            }
            
            //h2DataSum
            TH1F* h2DataSum;
            {
                TString NameTemp = "DataSum";
                h2DataSum = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
            }
            
            //h2BGMC
            TH1F* h2BGMC[BGMCSampleID.size()];
            std::vector<TString> hName2BGMC;
            for(unsigned int j=0;j<BGMCSampleID.size();j++)
            {
                TString NameTemp = "BG_";
                NameTemp += TString::Itoa(j,10);
                h2BGMC[j] = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                hName2BGMC.push_back(NameTemp);
            }
            
            //h2BGGruop
            Group BGGroup[BGMCGroupData.size()];
            for(unsigned int j=0;j<BGMCGroupData.size();j++)
            {
                BGGroup[j].info = &(BGMCGroupData[j]);
                TString NameTemp = BGMCGroupData[j].GroupName;
                BGGroup[j].h2 = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
            }
            
            //h2BGSum
            TH1F* h2BGSum;
            {
                TString NameTemp = "BGSum";
                h2BGSum = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
            }
            
            //ratio plot for reweighting
            {
                TString NameTemp = "Ratio_";
                NameTemp += TString::Itoa(LeptonIndex,10);
                h2Ratio_rw[LeptonIndex] = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                h2Ratio_rw[LeptonIndex]->GetXaxis()->SetTitle(xaxis.Data());
                h2Ratio_rw[LeptonIndex]->GetYaxis()->SetTitle("(Data - NonZjets)/Zjets");
                h2Ratio_rw[LeptonIndex]->SetMarkerColor(LeptonIndex+1);
                h2Ratio_rw[LeptonIndex]->SetMarkerSize(1);
                h2Ratio_rw[LeptonIndex]->SetLineColor(LeptonIndex+1);
                h2Ratio_rw[LeptonIndex]->Sumw2();
            }
            
            //fill histograms form trees
            for(unsigned int j=0;j<DataSampleID.size();j++)
            {
                TString temp;
                temp += Var[VarIndex].VarName;
                temp += ">>";
                temp += hName2Data[j];
                tree2Data[j]->Draw(temp.Data(),"fLwt==0");
            }
            
            for(unsigned int j=0;j<BGMCSampleID.size();j++)
            {
                TString temp;
                temp += Var[VarIndex].VarName;
                temp += ">>";
                temp += hName2BGMC[j];
                tree2BGMC[j]->Draw(temp.Data(),"weight");
            }
            
            //add data
            for(unsigned int j=0;j<DataSampleID.size();j++)
            {
                h2DataSum->Add(h2Data[j]);
            }
            h2Ratio_rw[LeptonIndex]->Add(h2DataSum);
            
            //normalization
            for(unsigned int j=0;j<BGMCSampleID.size();j++)
            {
                h2BGMC[j]->Scale(BGMCXS[j]/BGMCnAOD[j] *sumDataL);
            }
            
            //Non Z+jets
            for(unsigned int j=3;j<BGMCGroupData.size();j++)
            {
                for(unsigned int k=BGMCGroupData[j].lower;k<=BGMCGroupData[j].upper;k++)
                {
                    BGGroup[j].h2->Add(h2BGMC[k]);
                }
                h2BGSum->Add(BGGroup[j].h2);
            }
            h2BGSum->Scale(-1);
            h2Ratio_rw[LeptonIndex]->Add(h2BGSum);
            
            //Z+jets
            h2BGSum->Scale(0);
            for(int j=0;j<=2;j++)
            {
                for(unsigned int k=BGMCGroupData[j].lower;k<=BGMCGroupData[j].upper;k++)
                {
                    BGGroup[j].h2->Add(h2BGMC[k]);
                }
                h2BGSum->Add(BGGroup[j].h2);
            }
            h2Ratio_rw[LeptonIndex]->Divide(h2BGSum);
            
            //delete
            //h2Data
            for(unsigned int j=0;j<DataSampleID.size();j++)
            {
                delete h2Data[j];
            }
            
            //h2DataSum
            delete h2DataSum;
            
            //h2BGMC
            for(unsigned int j=0;j<BGMCSampleID.size();j++)
            {
                delete h2BGMC[j];
            }
            
            //h2BGGruop
            for(unsigned int j=0;j<BGMCGroupData.size();j++)
            {
                delete BGGroup[j].h2;
            }
            
            //h2BGSum
            delete h2BGSum;
            
            //delete tree2
            for(unsigned int i=0;i<DataSampleID.size();i++)
            {
                delete tree2Data[i];
            }
            for(unsigned int i=0;i<BGMCSampleID.size();i++)
            {
                delete tree2BGMC[i];
            }
        }
        
        TF1* fun[2];
        //simple fit
        if(simple)
        {
            for(int LeptonIndex=0;LeptonIndex<=1;LeptonIndex++)
            {
                TString NameTemp = "fun_";
                NameTemp += TString::Itoa(LeptonIndex,10);
                fun[LeptonIndex] = new TF1(NameTemp.Data(),"pol2",Var[VarIndex].xmin,Var[VarIndex].xmax);
                fun[LeptonIndex]->SetLineColor(LeptonIndex+1);
                fun[LeptonIndex]->SetLineStyle(2);
                h2Ratio_rw[LeptonIndex]->Fit(NameTemp.Data(),"R");
            }
        }
        
        //combined fit
        if(combined)
        {
            fun[1] = new TF1("fun_mumu_OS","pol2(0)",Var[VarIndex].xmin,Var[VarIndex].xmax);
            fun[0] = new TF1("fun_ee_OS","pol0(0)*(pol2(1))",Var[VarIndex].xmin,Var[VarIndex].xmax);
            
            ROOT::Math::WrappedMultiTF1 wfM(*(fun[1]),1);
            ROOT::Math::WrappedMultiTF1 wfE(*(fun[0]),1);
            
            ROOT::Fit::DataOptions opt;
            ROOT::Fit::DataRange rangeB;
            rangeB.SetRange(Var[VarIndex].xmin,Var[VarIndex].xmax);
            
            ROOT::Fit::BinData dataM(opt,rangeB);
            ROOT::Fit::FillData(dataM,h2Ratio_rw[1]);
            ROOT::Fit::BinData dataE(opt,rangeB);
            ROOT::Fit::FillData(dataE,h2Ratio_rw[0]);
            
            ROOT::Fit::Chi2Function chi2_M(dataM,wfM);
            ROOT::Fit::Chi2Function chi2_E(dataE,wfE);
            GlobalChi2 globalChi2(chi2_M,chi2_E);
            
            const int Npar = 4;
            double par0[Npar] = {1,0,0,1};
            ROOT::Fit::Fitter fitter;
            fitter.Config().SetParamsSettings(Npar,par0);
            fitter.Config().ParSettings(0).SetLimits(0.8,1.1);
            fitter.Config().ParSettings(1).SetLimits(-0.5,0.5);
            fitter.Config().ParSettings(2).SetLimits(0,0.1);
            fitter.Config().ParSettings(3).SetLimits(0.8,1.1);
            
            fitter.Config().MinimizerOptions().SetPrintLevel(0);
            fitter.Config().SetMinimizer("Minuit2","Migrad");
            
            fitter.FitFCN(Npar,globalChi2,0,dataM.Size()+dataE.Size(),true);
            ROOT::Fit::FitResult result = fitter.Result();
            result.Print(std::cout);
            gStyle->SetOptFit(1111);
            
            fun[1]->SetFitResult(result,iparM);
            fun[1]->SetRange(rangeB().first,rangeB().second);
            fun[1]->SetLineColor(2);
            fun[1]->SetLineStyle(2);
            h2Ratio_rw[1]->GetListOfFunctions()->Add(fun[1]);
            
            fun[0]->SetFitResult(result,iparE);
            fun[0]->SetRange(rangeB().first,rangeB().second);
            fun[0]->SetLineColor(1);
            fun[0]->SetLineStyle(2);
            h2Ratio_rw[0]->GetListOfFunctions()->Add(fun[0]);
        }
        
        TLegend* leg;
        {
            Double_t xl1=0.38, yl1=0.8, xl2=xl1+0.24, yl2=yl1+0.1;
            leg = new TLegend(xl1,yl1,xl2,yl2);
            leg->SetFillStyle(0);
            leg->SetTextFont(32);
            leg->SetBorderSize(0);
            leg->AddEntry(h2Ratio_rw[0],"ee channel","p");
            leg->AddEntry(h2Ratio_rw[1],"mumu channel","p");
        }
        
        TCanvas* c2 = new TCanvas();
        c2->cd();
        gStyle->SetOptStat(0);
        h2Ratio_rw[0]->Draw();
        h2Ratio_rw[1]->Draw("same");
        leg->Draw();
        
        {
            //export histograms in eps format
            TString NameTemp = "plot/";
            NameTemp += Var[VarIndex].VarName;
            NameTemp += "_rw";
            NameTemp += ".eps";
            c2->Print(NameTemp,"eps");
        }
        
        //delete h2Ratio_rw
        delete h2Ratio_rw[0];
        delete h2Ratio_rw[1];
        //delete legend
        delete leg;
        //delete canvas
        delete c2;
        
        //add weight in the tree
        for(int LeptonIndex=0;LeptonIndex<=1;LeptonIndex++)
        {
            for(unsigned int i=0;i<ChannelIndex_ll[LeptonIndex].size();i++)
            {
                TString FileName = "skimming/rw_";
                const int ChannelIndex = ChannelIndex_ll[LeptonIndex][i];
                FileName += TString::Itoa(ChannelIndex,10);
                FileName += ".root";
                TFile* frw = new TFile(FileName.Data(),"RECREATE");
                
                for(int j=0;j<=2;j++)
                {
                    for(unsigned int k=BGMCGroupData[j].lower;k<=BGMCGroupData[j].upper;k++)
                    {
                        TString NameTemp = "tree_rw_";
                        NameTemp += TString::Itoa(k,10);
                        
                        TTree* treeTemp = new TTree(NameTemp.Data(),NameTemp.Data());
                        treeTemp->Branch("rw",&rw,"rw/D");
                        for(int m=0;m<tree1BGMC[ChannelIndex][k]->GetEntries();m++)
                        {
                            tree1BGMC[ChannelIndex][k]->GetEntry(m);
                            if(simple)
                            {
                                rw=fun[LeptonIndex]->Eval(ptll);
                            }
                            if(combined)
                            {
                                rw=fun[1]->Eval(ptll);
                            }
                            treeTemp->Fill();
                        }
                        
                        frw->cd();
                        treeTemp->Write();
                    }
                }
                delete frw;
            }
        }
        
        //delete functions
        delete fun[0];
        delete fun[1];
    }
    
    {
        //Charge Flip root files
        TFile* f_cf_data;
        TH2D* h_cf_data;
        {
            TString NameTemp = "root_files/chargeMisID_Zee_data_signal.root";
            f_cf_data = new TFile(NameTemp.Data(),"READ");
            h_cf_data = (TH2D*) f_cf_data->Get("hFlipProb");
        }
        
        TFile* f_cf_mc;
        TH2D* h_cf_mc;
        {
            TString NameTemp = "root_files/chargeMisID_Zee_MC_signal.root";
            f_cf_mc = new TFile(NameTemp.Data(),"READ");
            h_cf_mc = (TH2D*) f_cf_mc->Get("hFlipProb");
        }
        
        //add charge flip weight in the tree
        for(int ChannelIndex=3;ChannelIndex<=9;ChannelIndex+=6)
        {
            TString FileName = "skimming/cfw_";
            FileName += TString::Itoa(ChannelIndex,10);
            FileName += ".root";
            TFile* fcfw = new TFile(FileName.Data(),"RECREATE");
            for(unsigned int k=BGMCGroupData[0].lower;k<=BGMCGroupData[0].upper;k++)
            {
                TString NameTemp = "tree_cfw_";
                NameTemp += TString::Itoa(k,10);
                
                TTree* treeTemp = new TTree(NameTemp.Data(),NameTemp.Data());
                treeTemp->Branch("cfw",&cfw,"cfw/D");
                for(int m=0;m<tree1BGMC[ChannelIndex][k]->GetEntries();m++)
                {
                    tree1BGMC[ChannelIndex][k]->GetEntry(m);
                    
                    double pt1s = pt1;
                    double pt2s = pt2;
                    double eta1s = fabs(eta1);
                    double eta2s = fabs(eta2);
                    
                    if(pt1s>=1000) pt1s = 999;
                    if(pt2s>=1000) pt2s = 999;
                    if(eta1s>=2.47) eta1s = 2.46;
                    if(eta2s>=2.47) eta2s = 2.46;
                    
                    double data1 = h_cf_data->GetBinContent(h_cf_data->FindBin(eta1s,pt1s));
                    double data2 = h_cf_data->GetBinContent(h_cf_data->FindBin(eta2s,pt2s));
                    
                    double mc1 = h_cf_mc->GetBinContent(h_cf_mc->FindBin(eta1s,pt1s));
                    double mc2 = h_cf_mc->GetBinContent(h_cf_mc->FindBin(eta2s,pt2s));
                    
                    cfw = ( data1*(1-data2) + data2*(1-data1) )/( mc1*(1-mc2) + mc2*(1-mc1) );
                    treeTemp->Fill();
                }
                
                fcfw->cd();
                treeTemp->Write();
            }
            delete fcfw;
        }
        delete f_cf_data;
        delete f_cf_mc;
    }
    
    //delete tree1
    for(unsigned int ChannelIndex=0;ChannelIndex<channel.size();ChannelIndex++)
    {
        for(unsigned int i=0;i<DataSampleID.size();i++)
        {
            //delete tree1Data[ChannelIndex][i];
        }
        for(unsigned int i=0;i<BGMCSampleID.size();i++)
        {
            delete tree1BGMC[ChannelIndex][i];
        }
        for(unsigned int i=0;i<SigSampleID.size();i++)
        {
            //delete tree1Sig[ChannelIndex][i];
        }
    }
    
    struct RegionData
    {
        TString RegionName;
        std::vector<unsigned int> setOfChannel;
        bool showData;
        bool showSignificance;
        bool isSS;
        bool isSS_ee;
        std::vector<unsigned int> qFChannel;
        TString Cut;
        std::vector<TString> setOfBGMC;
        std::vector<TString> setOfBGData;
    };
    
    std::vector<RegionData> RegionInfo;
    {
        RegionData element;
        
        for(unsigned int ChannelIndex=0;ChannelIndex<channel.size();ChannelIndex++)
        {
            element.RegionName = channel[ChannelIndex];
            element.setOfChannel.clear();
            element.qFChannel.clear();
            element.setOfChannel.push_back(ChannelIndex);
            unsigned int OS = ChannelIndex;
            if(OS>=6) OS -= 6;
            element.isSS = OS>=3;
            element.isSS_ee = OS==3;
            if(element.isSS_ee) element.qFChannel.push_back(ChannelIndex-3);
            element.showData = !element.isSS;
            element.showSignificance = element.isSS;
            element.Cut = "";
            
            RegionInfo.push_back(element);
        }
        
        
        //Control Region
        element.isSS = true;
        element.showData = true;
        element.showSignificance = false;
        element.Cut = "&& mll>81.18 && mll<101.18";
        
        //SS_ee
        element.isSS_ee = true;
        
        element.RegionName = "CR_nonISR_SS_ee";
        element.setOfChannel.clear();
        element.qFChannel.clear();
        element.setOfChannel.push_back(3);
        element.qFChannel.push_back(0);
        RegionInfo.push_back(element);
        
        element.RegionName = "CR_ISR_SS_ee";
        element.setOfChannel.clear();
        element.qFChannel.clear();
        element.setOfChannel.push_back(9);
        element.qFChannel.push_back(6);
        RegionInfo.push_back(element);
        
        element.RegionName = "CR_SS_ee";
        element.setOfChannel.clear();
        element.qFChannel.clear();
        element.setOfChannel.push_back(3);
        element.setOfChannel.push_back(9);
        element.qFChannel.push_back(0);
        element.qFChannel.push_back(6);
        RegionInfo.push_back(element);
        
        //SS_mumu
        element.isSS_ee = false;
        
        element.RegionName = "CR_nonISR_SS_mumu";
        element.setOfChannel.clear();
        element.qFChannel.clear();
        element.setOfChannel.push_back(4);
        RegionInfo.push_back(element);
        
        element.RegionName = "CR_ISR_SS_mumu";
        element.setOfChannel.clear();
        element.qFChannel.clear();
        element.setOfChannel.push_back(10);
        RegionInfo.push_back(element);
        
        element.RegionName = "CR_SS_mumu";
        element.setOfChannel.clear();
        element.qFChannel.clear();
        element.setOfChannel.push_back(4);
        element.setOfChannel.push_back(10);
        RegionInfo.push_back(element);
        
        //setOfBG
        for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
        {
            if(RegionInfo[RegionIndex].isSS)
            {
                RegionInfo[RegionIndex].setOfBGMC.push_back("VV");
                RegionInfo[RegionIndex].setOfBGMC.push_back("Vgamma");
                
                if(RegionInfo[RegionIndex].isSS_ee) RegionInfo[RegionIndex].setOfBGData.push_back("charge flip");
                RegionInfo[RegionIndex].setOfBGData.push_back("fake lepton");
            }
            else
            {
                RegionInfo[RegionIndex].setOfBGMC.push_back("Zee");
                RegionInfo[RegionIndex].setOfBGMC.push_back("Zmumu");
                RegionInfo[RegionIndex].setOfBGMC.push_back("Ztautau");
                RegionInfo[RegionIndex].setOfBGMC.push_back("ttbar");
                RegionInfo[RegionIndex].setOfBGMC.push_back("Wt");
                RegionInfo[RegionIndex].setOfBGMC.push_back("VV");
                RegionInfo[RegionIndex].setOfBGMC.push_back("Vgamma");
            }
        }
    }
    
    
    if(doOptimize)
    {
        //Significance
        int CutVarIndex = 13;
        
        //cut value
        double jetptCut[] = {20,25,30,35,40,45,50,70,100,150,200};
        //double jetptCut[] = {25,30,35,40,45,50,70,100,150,200};
        const int cutN = sizeof(jetptCut)/sizeof(jetptCut[0]);
        
        for(int RegionIndex=4;RegionIndex<=4;RegionIndex++)
        //for(int RegionIndex=10;RegionIndex<=10;RegionIndex++)
        //for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
        {
            double uncertainty[] = {0.1,0.2,0.3};
            const int uncertaintyN = sizeof(uncertainty)/sizeof(uncertainty[0]);
            double Significance[uncertaintyN][cutN];
            //uncertainty,cut value
            
            double sumOfEvent[cutN][BGMCGroupData.size()+4][2];
            //cut value,sample,expN/error
            
            int SigID = 18;
            
            std::vector<TChain*> tree2BGMC;
            initializeTree2(tree2BGMC,RegionInfo[RegionIndex].setOfChannel,BGMCSampleID,channel);
            std::vector<TChain*> tree2Sig;
            initializeTree2(tree2Sig,RegionInfo[RegionIndex].setOfChannel,SigSampleID,channel);
            
            //calculate sumOfEvent and Significance
            for(int u=0;u<uncertaintyN;u++)
            {
                for(int q=0;q<cutN;q++)
                {
                    const int VarIndex = 0;
                    //initialize histograms
                    TString title;
                    title += Var[VarIndex].VarTitle;
                    
                    TString xaxis;
                    xaxis += Var[VarIndex].VarTitle;
                    xaxis += " ";
                    xaxis += Var[VarIndex].unit;
                    
                    //h2BGMC
                    TH1F* h2BGMC[BGMCSampleID.size()];
                    std::vector<TString> hName2BGMC;
                    for(unsigned int j=0;j<BGMCSampleID.size();j++)
                    {
                        TString NameTemp = "BG_";
                        NameTemp += TString::Itoa(j,10);
                        h2BGMC[j] = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                        hName2BGMC.push_back(NameTemp);
                    }
                    
                    //h2BGGruop
                    Group BGGroup[BGMCGroupData.size()];
                    for(unsigned int j=0;j<BGMCGroupData.size();j++)
                    {
                        BGGroup[j].info = &(BGMCGroupData[j]);
                        TString NameTemp = BGMCGroupData[j].GroupName;
                        BGGroup[j].h2 = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                    }
                    
                    //h2Sig
                    TH1F* h2Sig;
                    TString hName2Sig;
                    {
                        hName2Sig = "Sig";
                        h2Sig = new TH1F(hName2Sig.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                    }
                    
                    //Draw for BG
                    for(unsigned int j=0;j<BGMCSampleID.size();j++)
                    {
                        TString temp = Var[VarIndex].VarName;
                        temp += ">>";
                        temp += hName2BGMC[j];
                        
                        TString Cut = "weight*(";
                        Cut += Var[CutVarIndex].VarName;
                        Cut += "<=";
                        Cut += TString::Itoa(jetptCut[q],10);;
                        Cut += ")";
                        tree2BGMC[j]->Draw(temp.Data(),Cut.Data());
                    }
                    
                    //Draw for Signal
                    {
                        TString temp = Var[VarIndex].VarName;
                        temp += ">>";
                        temp += hName2Sig;
                        
                        TString Cut = "weight*(";
                        Cut += Var[CutVarIndex].VarName;
                        Cut += "<=";
                        Cut += TString::Itoa(jetptCut[q],10);;
                        Cut += ")";
                        tree2Sig[SigID]->Draw(temp.Data(),Cut.Data());
                    }
                    
                    //normalization for BG
                    for(unsigned int j=0;j<BGMCSampleID.size();j++)
                    {
                        h2BGMC[j]->Scale(BGMCXS[j]/BGMCnAOD[j] *sumDataL);
                    }
                    
                    //add background
                    sumOfEvent[q][BGMCGroupData.size()][0]=0;
                    sumOfEvent[q][BGMCGroupData.size()][1]=0;
                    for(unsigned int j=0;j<BGMCGroupData.size();j++)
                    {
                        for(unsigned int k=BGMCGroupData[j].lower;k<=BGMCGroupData[j].upper;k++)
                        {
                            BGGroup[j].h2->Add(h2BGMC[k]);
                        }
                        
                        //expected number of events for BG
                        sumOfEvent[q][j][0] = BGGroup[j].h2->IntegralAndError(0,-1,sumOfEvent[q][j][1]);
                        //cout<<BGMCGroupData[j].GroupName.Data()<<": "<<sumOfEvent[q][j][0]<<" +/- "<<sumOfEvent[q][j][1]<<endl;
                        sumOfEvent[q][BGMCGroupData.size()][0] += sumOfEvent[q][j][0];
                        sumOfEvent[q][BGMCGroupData.size()][1] += sumOfEvent[q][j][1]*sumOfEvent[q][j][1];
                        
                    }
                    //if(VarIndex==0) cout<<"Total BG: "<<sumOfEvent[q][BGMCGroupData.size()][0]<<" +/- "<<TMath::Sqrt(sumOfEvent[q][BGMCGroupData.size()][1])<<endl;
                    
                    
                    //expected number of events
                    h2Sig->Scale(SigXS[SigID]/SignAOD[SigID] *sumDataL);
                    sumOfEvent[q][BGMCGroupData.size()+2][0] = h2Sig->IntegralAndError(0,-1,sumOfEvent[q][BGMCGroupData.size()+2][1]);
                    //cout<<"Signal: "<<sumOfEvent[q][BGMCGroupData.size()+2][0]<<" +/- "<<sumOfEvent[q][BGMCGroupData.size()+2][1]<<endl;
                    
                    //Significance
                    sumOfEvent[q][BGMCGroupData.size()+3][0] = RooStats::NumberCountingUtils::BinomialExpZ(sumOfEvent[q][BGMCGroupData.size()+2][0],sumOfEvent[q][BGMCGroupData.size()][0],uncertainty[u]);
                    //cout<<"Significance: "<<sumOfEvent[q][BGMCGroupData.size()+3][0]<<endl<<endl;
                    Significance[u][q] = sumOfEvent[q][BGMCGroupData.size()+3][0];
                    //cout<<"Signal: "<<sumOfEvent[q][BGMCGroupData.size()+2][0]<<", BG: "<<sumOfEvent[q][BGMCGroupData.size()][0]<<", Significance: "<<sumOfEvent[q][BGMCGroupData.size()+3][0]<<endl<<endl;
                    
                    //delete
                    //h2BGMC
                    for(unsigned int j=0;j<BGMCSampleID.size();j++)
                    {
                        delete h2BGMC[j];
                    }
                    
                    //h2BGGruop
                    for(unsigned int j=0;j<BGMCGroupData.size();j++)
                    {
                        delete BGGroup[j].h2;
                    }
                    
                    //h2Sig
                    delete h2Sig;
                }
            }
            
            //delete tree2
            for(unsigned int i=0;i<BGMCSampleID.size();i++)
            {
                delete tree2BGMC[i];
            }
            for(unsigned int i=0;i<SigSampleID.size();i++)
            {
                delete tree2Sig[i];
            }
            
            TCanvas* c2 = new TCanvas();
            TString xaxis;
            xaxis += Var[CutVarIndex].VarTitle;
            xaxis += " ";
            xaxis += Var[CutVarIndex].unit;
            {
                //Significance graph
                double min = Significance[0][0];
                double max = Significance[0][0];
                TGraph* SignificanceGraph[uncertaintyN];
                for(int u=0;u<uncertaintyN;u++)
                {
                    for(int q=0;q<cutN;q++)
                    {
                        cout<<Var[CutVarIndex].VarName.Data()<<" cut: "<<jetptCut[q]<<", Significance: "<<Significance[u][q]<<endl;
                        if(Significance[u][q] < min) min = Significance[u][q];
                        if(Significance[u][q] > max) max = Significance[u][q];
                    }
                    cout<<endl;
                    SignificanceGraph[u] = new TGraph(cutN,jetptCut,Significance[u]);
                    SignificanceGraph[u]->SetMarkerColor(u+2);
                    SignificanceGraph[u]->SetLineColor(u+2);
                }
                
                SignificanceGraph[0]->GetXaxis()->SetTitle(xaxis.Data());
                SignificanceGraph[0]->GetYaxis()->SetTitle("Significance");
                SignificanceGraph[0]->SetMinimum(min/1.1);
                SignificanceGraph[0]->SetMaximum(max*1.1);
                
                //Legend
                TLegend* leg;
                {
                    Double_t xl1, yl1, xl2, yl2;
                    xl2=0.92;
                    yl2=0.95;
                    xl1=xl2-0.3;
                    yl1=yl2-0.2;
                    
                    leg = new TLegend(xl1,yl1,xl2,yl2);
                    leg->SetNColumns(1);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(42);
                    leg->SetBorderSize(0);
                    
                    leg->AddEntry(SignificanceGraph[0],"uncertainty = 0.1","pl");
                    leg->AddEntry(SignificanceGraph[1],"uncertainty = 0.2","pl");
                    leg->AddEntry(SignificanceGraph[2],"uncertainty = 0.3","pl");
                }
                
                c2->cd();
                SignificanceGraph[0]->Draw("AL*");
                SignificanceGraph[1]->Draw("sameL*");
                SignificanceGraph[2]->Draw("sameL*");
                leg->Draw();
                //c2->WaitPrimitive();
                
                {
                    //export graph in eps format
                    TString NameTemp = "plot/significance_";
                    NameTemp += Var[CutVarIndex].VarName;
                    NameTemp += "_";
                    NameTemp += RegionInfo[RegionIndex].RegionName;
                    NameTemp += ".eps";
                    c2->Print(NameTemp,"eps");
                }
                
                //delete
                //delete SignificanceGraph
                for(int u=0;u<uncertaintyN;u++)
                {
                    delete SignificanceGraph[u];
                }
                
                //delete legend
                delete leg;
            }
            
            {
                //number event for BG and Sig
                double min = sumOfEvent[0][BGMCGroupData.size()][0];
                double max = sumOfEvent[0][BGMCGroupData.size()][0];
                double BGNumber[cutN];
                double SigNumber[cutN];
                for(int q=0;q<cutN;q++)
                {
                    BGNumber[q] = sumOfEvent[q][BGMCGroupData.size()][0];
                    SigNumber[q] = sumOfEvent[q][BGMCGroupData.size()+2][0];
                    if(BGNumber[q] < min) min = BGNumber[q];
                    if(SigNumber[q] < min) min = SigNumber[q];
                    if(BGNumber[q] > max) max = BGNumber[q];
                    if(SigNumber[q] > max) max = SigNumber[q];
                }
                TGraph* BGGraph = new TGraph(cutN,jetptCut,BGNumber);
                BGGraph->SetMarkerColor(2);
                BGGraph->SetLineColor(2);
                
                TGraph* SigGraph = new TGraph(cutN,jetptCut,SigNumber);
                SigGraph->SetMarkerColor(3);
                SigGraph->SetLineColor(3);
                
                BGGraph->GetXaxis()->SetTitle(xaxis.Data());
                BGGraph->GetYaxis()->SetTitle("Number of events");
                BGGraph->SetMinimum(0);
                BGGraph->SetMaximum(max*1.1);
                
                //Legend
                TLegend* leg;
                {
                    Double_t xl1, yl1, xl2, yl2;
                    xl2=0.5;
                    yl2=0.95;
                    xl1=xl2-0.3;
                    yl1=yl2-0.2;
                    
                    leg = new TLegend(xl1,yl1,xl2,yl2);
                    leg->SetNColumns(1);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(42);
                    leg->SetBorderSize(0);
                    
                    leg->AddEntry(BGGraph,"Total Background","pl");
                    leg->AddEntry(SigGraph,"Signal","pl");
                }
                
                BGGraph->Draw("AL*");
                SigGraph->Draw("sameL*");
                leg->Draw();
                
                {
                    //export graph in eps format
                    TString NameTemp = "plot/Number_";
                    NameTemp += Var[CutVarIndex].VarName;
                    NameTemp += "_";
                    NameTemp += RegionInfo[RegionIndex].RegionName;
                    NameTemp += ".eps";
                    c2->Print(NameTemp,"eps");
                }
                
                //delete Graph
                delete BGGraph;
                delete SigGraph;
                
                //delete legend
                delete leg;
            }
            
            //delete canvas
            delete c2;
        }
    }
    
    
    struct SampleData
    {
        TString SampleName;
        int index;
    };
    
    std::vector<SampleData> BGVVData;
    for(unsigned int i = BGMCGroupData[5].lower;i <= BGMCGroupData[5].upper;i++)
    {
        SampleData element;
        element.SampleName = BGMCSampleID[i].Data();
        element.SampleName.Remove(0,19);
        element.SampleName.ReplaceAll("_","\\_");
        element.index = i;
        BGVVData.push_back(element);
    }
    
    if(false)
    //if(true)
    {
        double sumOfEvent[BGMCGroupData.size()+SigMassSplitting.size()+2][3];
        //sample,expN/error/significance
        
        double sumOfEventVV[BGVVData.size()][2];
        //sample,expN/error
        
        //print expected number of event
        for(unsigned int RegionIndex=4;RegionIndex<=4;RegionIndex++)
        //for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
        {
            std::vector<TChain*> tree2Data;
            initializeTree2(tree2Data,RegionInfo[RegionIndex].setOfChannel,DataSampleID,channel);
            std::vector<TChain*> tree2BGMC;
            initializeTree2(tree2BGMC,RegionInfo[RegionIndex].setOfChannel,BGMCSampleID,channel);
            std::vector<TChain*> tree2Sig;
            initializeTree2(tree2Sig,RegionInfo[RegionIndex].setOfChannel,SigSampleID,channel);
            
            cout<<RegionInfo[RegionIndex].RegionName.Data()<<endl;
            sumOfEvent[BGMCGroupData.size()][0]=0;
            sumOfEvent[BGMCGroupData.size()][1]=0;
            for(unsigned int j=0;j<BGMCGroupData.size();j++)
            {
                sumOfEvent[j][0]=0;
                sumOfEvent[j][1]=0;
                for(unsigned int k=BGMCGroupData[j].lower;k<=BGMCGroupData[j].upper;k++)
                {
                    double expN = 0;
                    double error2 = 0;
                    for(int m=0;m<tree2BGMC[k]->GetEntries();m++)
                    {
                        tree2BGMC[k]->GetEntry(m);
                        expN += weight;
                        error2 += weight*weight;
                    }
                    double cross = BGMCXS[k]/BGMCnAOD[k] *sumDataL;
                    sumOfEvent[j][0] += expN *cross;
                    sumOfEvent[j][1] += error2 *cross*cross;
                    
                    if(j==5)
                    {
                        int VVindex = k - BGMCGroupData[5].lower;
                        sumOfEventVV[VVindex][0] = expN *cross;
                        sumOfEventVV[VVindex][1] = TMath::Sqrt(error2) *cross;
                        cout<<BGVVData[VVindex].SampleName.Data()<<": "<<sumOfEventVV[VVindex][0]<<" +/- "<<sumOfEventVV[VVindex][1]<<endl;
                    }
                }
                cout<<BGMCGroupData[j].GroupName.Data()<<": "<<sumOfEvent[j][0]<<" +/- "<<TMath::Sqrt(sumOfEvent[j][1])<<endl;
                sumOfEvent[BGMCGroupData.size()][0] += sumOfEvent[j][0];
                sumOfEvent[BGMCGroupData.size()][1] += sumOfEvent[j][1];
            }
            cout<<"Total BG: "<<sumOfEvent[BGMCGroupData.size()][0]<<" +/- "<<TMath::Sqrt(sumOfEvent[BGMCGroupData.size()][1])<<endl;
            
            sumOfEvent[BGMCGroupData.size()+1][0]=0;
            for(unsigned int j=0;j<DataSampleID.size();j++)
            {
                sumOfEvent[BGMCGroupData.size()+1][0] += double(tree2Data[j]->GetEntries());
            }
            cout<<"Data: "<<sumOfEvent[BGMCGroupData.size()+1][0]<<endl<<endl;
            
            for(unsigned int i=0;i<SigMassSplitting.size();i++)
            {
                double expN = 0;
                double error2 = 0;
                for(int m=0;m<tree2Sig[SigMassSplitting[i].ID]->GetEntries();m++)
                {
                    tree2Sig[SigMassSplitting[i].ID]->GetEntry(m);
                    expN += weight;
                    error2 += weight*weight;
                }
                double cross = SigXS[SigMassSplitting[i].ID]/SignAOD[SigMassSplitting[i].ID] *sumDataL;
                sumOfEvent[BGMCGroupData.size()+i+2][0] = expN *cross;
                sumOfEvent[BGMCGroupData.size()+i+2][1] = TMath::Sqrt(error2) *cross;
                sumOfEvent[BGMCGroupData.size()+i+2][2] = RooStats::NumberCountingUtils::BinomialExpZ(sumOfEvent[BGMCGroupData.size()+i+2][0],sumOfEvent[BGMCGroupData.size()][0],0.3);
                
                cout<<"Signal ("<<SigMass1[SigMassSplitting[i].ID]<<", "<<SigMass2[SigMassSplitting[i].ID]<<"): "<<sumOfEvent[BGMCGroupData.size()+i+2][0]<<" +/- "<<sumOfEvent[BGMCGroupData.size()+i+2][1]<<", Significance: "<<sumOfEvent[BGMCGroupData.size()+i+2][2]<<endl;
            }
            cout<<endl;
            
            //delete tree2
            for(unsigned int i=0;i<DataSampleID.size();i++)
            {
                delete tree2Data[i];
            }
            for(unsigned int i=0;i<BGMCSampleID.size();i++)
            {
                delete tree2BGMC[i];
            }
            for(unsigned int i=0;i<SigSampleID.size();i++)
            {
                delete tree2Sig[i];
            }
        }
    }
    
    //plot graph
    bool optimize = 0;
    unsigned int countVariable = 31;
    for(unsigned int RegionIndex=0;RegionIndex<=0;RegionIndex++)
    //for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
    {
        std::vector<TChain*> tree2Data;
        initializeTree2(tree2Data,RegionInfo[RegionIndex].setOfChannel,DataSampleID,channel);
        
        std::vector< std::vector<TChain*> > tree2BGMC;
        std::vector< std::vector<double> > BGMCGroupXS;
        std::vector< std::vector<unsigned int> > BGMCGroupnAOD;
        std::vector<Group> BGGroup;
        {
            //For MC background
            for(unsigned int j=0;j<BGMCGroupData.size();j++)
            {
                for(unsigned int k=0;k<RegionInfo[RegionIndex].setOfBGMC.size();k++)
                {
                    if(BGMCGroupData[j].GroupName == RegionInfo[RegionIndex].setOfBGMC[k])
                    {
                        std::vector<TString> BGMCGroupSampleID;
                        std::vector<double> BGMCGroupXSElement;
                        std::vector<unsigned int> BGMCGroupnAODElement;
                        for(unsigned int m=BGMCGroupData[j].lower;m<=BGMCGroupData[j].upper;m++)
                        {
                            BGMCGroupSampleID.push_back(BGMCSampleID[m]);
                            BGMCGroupXSElement.push_back(BGMCXS[m]);
                            BGMCGroupnAODElement.push_back(BGMCnAOD[m]);
                        }
                        std::vector<TChain*> tree2BGMCElement;
                        initializeTree2(tree2BGMCElement,RegionInfo[RegionIndex].setOfChannel,BGMCGroupSampleID,channel);
                        
                        tree2BGMC.push_back(tree2BGMCElement);
                        BGMCGroupXS.push_back(BGMCGroupXSElement);
                        BGMCGroupnAOD.push_back(BGMCGroupnAODElement);
                        
                        Group BGGroupElement;
                        BGGroupElement.info = &(BGMCGroupData[j]);
                        BGGroup.push_back(BGGroupElement);
                    }
                }
                
                //Z pt reweighting
                if(dorw
                   &&
                   (
                    RegionIndex==0 ||
                    RegionIndex==1 ||
                    RegionIndex==6 ||
                    RegionIndex==7
                   )
                   &&
                   (
                    BGGroup[j].info->GroupName == "Zee" ||
                    BGGroup[j].info->GroupName == "Zmumu" ||
                    BGGroup[j].info->GroupName == "Ztautau"
                   )
                  )
                {
                    for(unsigned int k=0;k<tree2BGMC[j].size();k++)
                    {
                        TString NameTemp = "tree_rw_";
                        NameTemp += TString::Itoa(k + BGGroup[j].info->lower,10);
                        TChain* ch_rw = new TChain(NameTemp.Data());
                        
                        for(unsigned int m=0;m<RegionInfo[RegionIndex].setOfChannel.size();m++)
                        {
                            TString FileName = "skimming/rw_";
                            FileName += TString::Itoa(RegionInfo[RegionIndex].setOfChannel[m],10);
                            FileName += ".root";
                            
                            ch_rw->Add(FileName.Data());
                        }
                        tree2BGMC[j][k]->AddFriend(NameTemp.Data());
                    }
                }
                
            }
            
            //For data-driven background
            for(unsigned int j=0;j<BGDataGroupData.size();j++)
            {
                for(unsigned int k=0;k<RegionInfo[RegionIndex].setOfBGData.size();k++)
                {
                    if(BGDataGroupData[j].GroupName == RegionInfo[RegionIndex].setOfBGData[k])
                    {
                        Group BGGroupElement;
                        BGGroupElement.info = &(BGDataGroupData[j]);
                        BGGroup.push_back(BGGroupElement);
                    }
                }
            }
        }
        
        std::vector<TChain*> tree2Sig;
        initializeTree2(tree2Sig,RegionInfo[RegionIndex].setOfChannel,SigSampleID,channel);
        
        std::vector<TChain*> tree2DataOS;
        if(RegionInfo[RegionIndex].isSS_ee) initializeTree2(tree2DataOS,RegionInfo[RegionIndex].qFChannel,DataSampleID,channel);
        
        /*
        //for charge filp BG
        if(RegionIndex>=12 && RegionIndex<=14)
        {
            for(unsigned int k=BGMCGroupData[0].lower;k<=BGMCGroupData[0].upper;k++)
            {
                TString NameTemp = "tree_cfw_";
                NameTemp += TString::Itoa(k,10);
                TChain* ch_cfw = new TChain(NameTemp.Data());
                
                for(unsigned int j=0;j<RegionInfo[RegionIndex].setOfChannel.size();j++)
                {
                    TString FileName = "skimming/cfw_";
                    FileName += TString::Itoa(RegionInfo[RegionIndex].setOfChannel[j],10);
                    FileName += ".root";
                    
                    ch_cfw->Add(FileName.Data());
                }
                tree2BGMC[k]->AddFriend(NameTemp.Data());
            }
        }
        */
        
        //for(unsigned int VarIndex=5;VarIndex<=5;VarIndex++)
        for(unsigned int VarIndex=countVariable;VarIndex<=countVariable;VarIndex++)
        //for(unsigned int VarIndex=0;VarIndex<Var.size();VarIndex++)
        {
            //initialize histograms
            TString title;
            title += Var[VarIndex].VarTitle;
            title += "_";
            title += RegionInfo[RegionIndex].RegionName;
            
            TString xaxis;
            xaxis += Var[VarIndex].VarTitle;
            xaxis += " ";
            xaxis += Var[VarIndex].unit;
            
            //Fill histograms from trees
            TH1F* h2DataSum;
            TH1F* h2SigSum[SigMassSplitting.size()];
            TH1F* h2Sig[SigSampleID.size()];
            double sumOfEventVV[BGVVData.size()][2];
            //sample,expN/error

            {
                //Common Cut
                TString CommonCut = "";
                if(Var[VarIndex].VarName=="jetpt"  ||
                   Var[VarIndex].VarName=="jeteta" ||
                   Var[VarIndex].VarName=="jetphi" )
                {
                    CommonCut += " && nJet>0";
                }
                
                if(Var[VarIndex].VarName=="bjetpt"  ||
                   Var[VarIndex].VarName=="bjeteta" ||
                   Var[VarIndex].VarName=="bjetphi" )
                {
                    CommonCut += " && nBJet>0";
                }
                
                if(Var[VarIndex].VarName=="cjetpt"  ||
                   Var[VarIndex].VarName=="cjeteta" ||
                   Var[VarIndex].VarName=="cjetphi" )
                {
                    CommonCut += " && nCJet>0";
                }
                
                if(Var[VarIndex].VarName=="fjetpt"  ||
                   Var[VarIndex].VarName=="fjeteta" ||
                   Var[VarIndex].VarName=="fjetphi" )
                {
                    CommonCut += " && nFJet>0";
                }
                
                CommonCut += RegionInfo[RegionIndex].Cut;
                
                //h2DataSum
                {
                    TString NameTemp = "DataSum";
                    h2DataSum = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                    h2DataSum->GetYaxis()->SetTitle("Number of events");
                    h2DataSum->SetMarkerColor(1);
                    h2DataSum->SetMarkerStyle(20);
                    h2DataSum->SetMarkerSize(1);
                    h2DataSum->SetLineColor(1);
                }
                
                //h2Data
                TH1F* h2Data[DataSampleID.size()];
                std::vector<TString> hName2Data;
                for(unsigned int j=0;j<DataSampleID.size();j++)
                {
                    TString NameTemp = "Data_";
                    NameTemp += TString::Itoa(j,10);
                    h2Data[j] = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                    hName2Data.push_back(NameTemp);
                    
                    //Fill data
                    TString temp;
                    if(Var[VarIndex].VarName=="averageMu")
                    {
                        temp = "averageMu/1.16";
                    }
                    else
                    {
                        temp = Var[VarIndex].VarName;
                    }
                    temp += ">>";
                    temp += hName2Data[j];
                    
                    TString Cut = "(1";
                    Cut += CommonCut;
                    Cut += "&& fLwt==0";
                    
                    if(optimize)
                    {
                        Cut += " && jetpt<=";
                        Cut += TString::Itoa(35,10);
                    }
                    Cut += ")";
                    tree2Data[j]->Draw(temp.Data(),Cut.Data());
                    
                    //Add data
                    h2DataSum->Add(h2Data[j]);
                }
                
                //Background
                //For MC background
                for(unsigned int j=0;j<tree2BGMC.size();j++)
                {
                    BGGroup[j].h2 = new TH1F(BGGroup[j].info->GroupName.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                    BGGroup[j].h2->GetYaxis()->SetTitle("Number of events");
                    BGGroup[j].h2->SetLineColor(j+2);
                    BGGroup[j].h2->SetFillColor(j+2);
                    
                    for(unsigned int k=0;k<tree2BGMC[j].size();k++)
                    {
                        TH1F* hTemp = new TH1F("BGMC",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                        
                        //fill histograms from trees
                        TString temp = Var[VarIndex].VarName;
                        temp += ">>BGMC";
                        
                        TString Cut = "weight";
                        
                        
                        //Z pt reweighting
                        if(dorw
                           &&
                            (
                             RegionIndex==0 ||
                             RegionIndex==1 ||
                             RegionIndex==6 ||
                             RegionIndex==7
                            )
                           &&
                            (
                             BGGroup[j].info->GroupName == "Zee" ||
                             BGGroup[j].info->GroupName == "Zmumu" ||
                             BGGroup[j].info->GroupName == "Ztautau"
                            )
                          )
                        {
                            Cut += "*rw";
                        }
                        
                        //for charge filp BG
                        if(RegionIndex>=12 && RegionIndex<=14 &&
                           BGGroup[j].info->GroupName == "Zee")
                        {
                            //Cut += "*cfw";
                        }
                        
                        
                        Cut += "*(1";
                        Cut += CommonCut;
                        
                        if(optimize)
                        {
                            Cut += " && jetpt<=";
                            Cut += TString::Itoa(35,10);
                        }
                        Cut += ")";
                        tree2BGMC[j][k]->Draw(temp.Data(),Cut.Data());
                        
                        //normalization for BG
                        hTemp->Scale(BGMCGroupXS[j][k]/BGMCGroupnAOD[j][k] *sumDataL);
                        
                        //expN for BGVV
                        if(VarIndex==countVariable && BGGroup[j].info->GroupName == "VV")
                        {
                            sumOfEventVV[k][0] = hTemp->IntegralAndError(0,-1,sumOfEventVV[k][1]);
                            cout<<BGVVData[k].SampleName.Data()<<": "<<sumOfEventVV[k][0]<<" +/- "<<sumOfEventVV[k][1]<<endl;
                        }
                        
                        BGGroup[j].h2->Add(hTemp);
                        delete hTemp;
                    }
                }
                
                //For data-driven background
                for(unsigned int j=tree2BGMC.size();j<BGGroup.size();j++)
                {
                    TString NameTemp = BGGroup[j].info->GroupName;
                    BGGroup[j].h2 = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                    BGGroup[j].h2->GetYaxis()->SetTitle("Number of events");
                    BGGroup[j].h2->SetLineColor(j+2);
                    BGGroup[j].h2->SetFillColor(j+2);
                    
                    for(unsigned int k=0;k<DataSampleID.size();k++)
                    {
                        h2Data[k]->Scale(0);
                        
                        TString temp = Var[VarIndex].VarName;
                        temp += ">>";
                        temp += hName2Data[k];
                        
                        TString Cut = "1";
                        if(BGGroup[j].info->GroupName == "charge flip")
                        {
                            Cut += "*qFwt";
                        }
                        if(BGGroup[j].info->GroupName == "fake lepton")
                        {
                            Cut += "*fLwt";
                        }
                        
                        Cut += "*(1";
                        Cut += CommonCut;
                        if(BGGroup[j].info->GroupName == "charge flip")
                        {
                            Cut += "&& fLwt==0";
                        }
                        if(BGGroup[j].info->GroupName == "fake lepton")
                        {
                            Cut += "&& fLwt!=0";
                        }
                        
                        if(optimize)
                        {
                            Cut += " && jetpt<=";
                            Cut += TString::Itoa(35,10);
                        }
                        Cut += ")";
                        
                        if(BGGroup[j].info->GroupName == "charge flip")
                        {
                            tree2DataOS[k]->Draw(temp.Data(),Cut.Data());
                        }
                        if(BGGroup[j].info->GroupName == "fake lepton")
                        {
                            tree2Data[k]->Draw(temp.Data(),Cut.Data());
                        }
                        
                        //Add MCData
                        BGGroup[j].h2->Add(h2Data[k]);
                    }
                }
                
                //delete h2Data
                for(unsigned int j=0;j<DataSampleID.size();j++)
                {
                    delete h2Data[j];
                }
                
                //h2SigSum
                for(unsigned int j=0;j<SigMassSplitting.size();j++)
                {
                    TString NameTemp = "SigSum_";
                    NameTemp += TString::Itoa(SigMassSplitting[j].MassDiff,10);
                    h2SigSum[j] = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                    h2SigSum[j]->GetYaxis()->SetTitle("Number of events");
                    h2SigSum[j]->SetLineColor(SigMassSplitting[j].colour);
                    h2SigSum[j]->SetLineStyle(1);
                }
                
                //h2Sig
                
                for(unsigned int j=0;j<SigSampleID.size();j++)
                {
                    TString NameTemp = "Sig_";
                    NameTemp += TString::Itoa(j,10);
                    h2Sig[j] = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                    
                    //Fill Signal
                    TString temp = Var[VarIndex].VarName;
                    temp += ">>";
                    temp += NameTemp;
                    
                    TString Cut = "weight*(1";
                    Cut += CommonCut;
                    
                    if(optimize)
                    {
                        Cut += " && jetpt<=";
                        Cut += TString::Itoa(35,10);
                    }
                    Cut += ")";
                    tree2Sig[j]->Draw(temp.Data(),Cut.Data());
                }
            }
            
            //Add Signal for the same mass splitting
            const int SigScale = 10;
            for(unsigned int i=0;i<SigMassSplitting.size();i++)
            {
                unsigned int AOD = 0;
                for(unsigned int j=0;j<SigSampleID.size();j++)
                {
                    if(SigMass1[j]-SigMass2[j] == SigMassSplitting[i].MassDiff)
                    {
                        h2SigSum[i]->Add(h2Sig[j]);
                        AOD += SignAOD[j];
                    }
                }
                
                //normalization for Signal
                h2SigSum[i]->Scale(SigXS[SigMassSplitting[i].ID]/AOD *sumDataL *SigScale);
            }
            
            //Add BG
            TH1F* h2BGSum = new TH1F("BGSum",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
            std::vector<Group> vBGGroup;
            for(unsigned int j=0;j<BGGroup.size();j++)
            {
                h2BGSum->Add(BGGroup[j].h2);
                vBGGroup.push_back(BGGroup[j]);
            }
            
            //sort
            std::sort(vBGGroup.begin(),vBGGroup.end(),compare2);
            
            //stack
            THStack stack;
            for(unsigned int j=0;j<vBGGroup.size();j++)
            {
                stack.Add(vBGGroup[j].h2);
            }
            
            if(VarIndex==countVariable)
            {
                double sumOfEvent[BGGroup.size()+SigMassSplitting.size()+2][3];
                //sample,expN/error/significance
                
                //expected number of events for Data
                sumOfEvent[BGGroup.size()+1][0] = h2DataSum->Integral(0,-1);
                cout<<"Data: "<<sumOfEvent[BGGroup.size()+1][0]<<endl;
                
                //expected number of events for background
                sumOfEvent[BGGroup.size()][0]=0;
                sumOfEvent[BGGroup.size()][1]=0;
                
                for(unsigned int j=0;j<BGGroup.size();j++)
                {
                    //expected number of events for BG
                    sumOfEvent[j][0] = BGGroup[j].h2->IntegralAndError(0,-1,sumOfEvent[j][1]);
                    cout<<BGGroup[j].info->GroupName.Data()<<": "<<sumOfEvent[j][0]<<" +/- "<<sumOfEvent[j][1]<<endl;
                    sumOfEvent[BGGroup.size()][0] += sumOfEvent[j][0];
                    sumOfEvent[BGGroup.size()][1] += sumOfEvent[j][1]*sumOfEvent[j][1];
                }
                cout<<"Total BG: "<<sumOfEvent[BGGroup.size()][0]<<" +/- "<<TMath::Sqrt(sumOfEvent[BGGroup.size()][1])<<endl<<endl;
        
                //expected number of events for signal
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    //expected number of events
                    h2Sig[SigMassSplitting[i].ID]->Scale(SigXS[SigMassSplitting[i].ID]/SignAOD[SigMassSplitting[i].ID] *sumDataL);
                    sumOfEvent[BGGroup.size()+i+2][0] = h2Sig[SigMassSplitting[i].ID]->IntegralAndError(0,-1,sumOfEvent[BGGroup.size()+i+2][1]);
                    
                    //Significance
                    sumOfEvent[BGGroup.size()+i+2][2] = RooStats::NumberCountingUtils::BinomialExpZ(sumOfEvent[BGGroup.size()+i+2][0],sumOfEvent[BGGroup.size()][0],0.3);
                    
                    cout<<"Signal ("<<SigMass1[SigMassSplitting[i].ID]<<", "<<SigMass2[SigMassSplitting[i].ID]<<"): "<<sumOfEvent[BGGroup.size()+i+2][0]<<" +/- "<<sumOfEvent[BGGroup.size()+i+2][1]<<", Significance: "<<sumOfEvent[BGGroup.size()+i+2][2]<<endl;
                }
                cout<<endl;
                
                //output file
                TString PathName = "latex/data/expN/";
                PathName += RegionInfo[RegionIndex].RegionName;
                PathName += ".tex";
                
                ofstream fout;
                fout.open(PathName.Data());
                fout<<setprecision(1)<<std::fixed;
                
                for(unsigned int j=0;j<BGGroup.size();j++)
                {
                    fout<<BGGroup[j].info->LatexName.Data();
                    fout<<" & $";
                    fout<<sumOfEvent[j][0];
                    fout<<"\\pm";
                    fout<<sumOfEvent[j][1];
                    fout<<"$ & \\\\"<<endl<<"\\hline"<<endl;
                }
                
                fout<<"Total BG & $";
                fout<<sumOfEvent[BGGroup.size()][0];
                fout<<"\\pm";
                fout<<TMath::Sqrt(sumOfEvent[BGGroup.size()][1]);
                fout<<"$ & \\\\"<<endl<<"\\hline"<<endl;
                
                if(RegionInfo[RegionIndex].showData)
                {
                    fout<<"Data & $";
                    fout<<sumOfEvent[BGGroup.size()+1][0];
                    fout<<"$ & \\\\"<<endl<<"\\hline"<<endl;
                }
                
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    fout<<"Signal (";
                    fout<<SigMass1[SigMassSplitting[i].ID];
                    fout<<", ";
                    fout<<SigMass2[SigMassSplitting[i].ID];
                    fout<<") & $";
                    fout<<setprecision(1)<<std::fixed;
                    fout<<sumOfEvent[BGGroup.size()+i+2][0];
                    fout<<"\\pm";
                    fout<<sumOfEvent[BGGroup.size()+i+2][1];
                    fout<<"$ &";
                    if(RegionInfo[RegionIndex].showSignificance)
                    {
                        fout<<"$";
                        fout<<setprecision(3)<<std::fixed;
                        fout<<sumOfEvent[BGGroup.size()+i+2][2];
                        fout<<"$";
                    }
                    fout<<"\\\\"<<endl<<"\\hline"<<endl;
                }
                fout.close();
                
                //output file for VV
                PathName = "latex/data/expN/BGVV_";
                PathName += RegionInfo[RegionIndex].RegionName;
                PathName += ".tex";
                
                fout.open(PathName.Data());
                fout<<setprecision(3)<<std::fixed;
                for(unsigned int j=0;j<BGVVData.size();j++)
                {
                    fout<<BGVVData[j].SampleName.Data();
                    fout<<" & $";
                    fout<<sumOfEventVV[j][0];
                    fout<<"\\pm";
                    fout<<sumOfEventVV[j][1];
                    fout<<"$ \\\\"<<endl<<"\\hline"<<endl;
                }
                fout.close();
            }
            
            {
                //adjust the max and min value for y-axis
                double min = h2DataSum->GetBinContent(h2DataSum->GetMinimumBin());
                double max = h2DataSum->GetBinContent(h2DataSum->GetMaximumBin());
                
                if(h2BGSum->GetBinContent(h2BGSum->GetMinimumBin()) < min) min = h2BGSum->GetBinContent(h2BGSum->GetMinimumBin());
                if(h2BGSum->GetBinContent(h2BGSum->GetMaximumBin()) > max) max = h2BGSum->GetBinContent(h2BGSum->GetMaximumBin());
                
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    if(h2SigSum[i]->GetBinContent(h2SigSum[i]->GetMinimumBin()) < min) min = h2SigSum[i]->GetBinContent(h2SigSum[i]->GetMinimumBin());
                }
                
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    if(h2SigSum[i]->GetBinContent(h2SigSum[i]->GetMaximumBin()) > max) max = h2SigSum[i]->GetBinContent(h2SigSum[i]->GetMaximumBin());
                }
                
                if(Var[VarIndex].log)
                {
                    h2DataSum->SetMaximum(max*100);
                    if(min<0.1)
                    {
                        h2DataSum->SetMinimum(Var[VarIndex].ymin);
                    }
                    else
                    {
                        h2DataSum->SetMinimum(min/2);
                    }
                }
                else
                {
                    //h2DataSum->SetMinimum(Var[VarIndex].ymin);
                    //h2DataSum->SetMaximum(Var[VarIndex].ymax);
                    
                    if(min>0) h2DataSum->SetMinimum(min*0.9);
                    else h2DataSum->SetMinimum(min*1.1);
                    h2DataSum->SetMaximum(max*1.1);
                }
            }
            
            //Legend
            TLegend* leg;
            {
                Double_t xl1, yl1, xl2, yl2;
                if(RegionInfo[RegionIndex].showData)
                {
                    xl2=0.92;
                    yl2=0.95;
                    xl1=xl2-0.3;
                    yl1=yl2-0.3;
                }
                else
                {
                    xl2=0.92;
                    yl2=0.95;
                    xl1=xl2-0.3;
                    yl1=yl2-0.2;
                }
                leg = new TLegend(xl1,yl1,xl2,yl2);
                leg->SetNColumns(2);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetBorderSize(0);
                if(RegionInfo[RegionIndex].showData)
                {
                    leg->AddEntry(h2DataSum,"Data","p");
                }
                for(unsigned int j=0;j<BGGroup.size();j++)
                {
                    leg->AddEntry(BGGroup[j].h2,BGGroup[j].info->LegendName.Data(),"fl");
                }
                
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    //TString NameTemp = "C1N2#rightarrow WZ";
                    //TString NameTemp = "C1N2#rightarrow slep";
                    TString NameTemp = "";
                    NameTemp += "(";
                    NameTemp += TString::Itoa(SigMass1[SigMassSplitting[i].ID],10);
                    NameTemp += ", ";
                    NameTemp += TString::Itoa(SigMass2[SigMassSplitting[i].ID],10);
                    NameTemp += ") x";
                    NameTemp += TString::Itoa(SigScale,10);
                    leg->AddEntry(h2SigSum[i],NameTemp.Data(),"l");
                }
            }
            
            TH1F* h2Ratio = nullptr;
            TPad* pad1 = nullptr;
            TPad* pad2 = nullptr;
            if(RegionInfo[RegionIndex].showData)
            {
                //size for two pads
                const double size1 = 0.65;
                const double size2 = 0.35;
                const double scale1 = 1/size1;
                const double scale2 = 1/size2;
                
                //adjust the title
                const double x_label = h2DataSum->GetXaxis()->GetLabelSize();
                const double y_label = h2DataSum->GetYaxis()->GetLabelSize();
                const double x_title = h2DataSum->GetXaxis()->GetTitleSize();
                const double y_title = h2DataSum->GetYaxis()->GetTitleSize();
                const double y_offset = h2DataSum->GetYaxis()->GetTitleOffset();
                
                h2DataSum->GetXaxis()->SetLabelSize(x_label*scale1);
                h2DataSum->GetYaxis()->SetLabelSize(y_label*scale1);
                h2DataSum->GetXaxis()->SetTitleSize(x_title*scale1);
                h2DataSum->GetYaxis()->SetTitleSize(y_title*scale1);
                h2DataSum->GetYaxis()->SetTitleOffset(y_offset/scale1);
                
                //ratio plot
                TString NameTemp = "Ratio";
                h2Ratio = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                h2Ratio->GetXaxis()->SetTitle(xaxis.Data());
                h2Ratio->GetYaxis()->SetTitle("Data/MC");
                h2Ratio->GetYaxis()->CenterTitle();
                h2Ratio->SetMarkerSize(1.0);
                h2Ratio->SetLineColor(1);
                
                h2Ratio->Sumw2();
                h2Ratio->Add(h2DataSum);
                h2Ratio->Divide(h2BGSum);
                
                h2Ratio->GetXaxis()->SetLabelSize(x_label*scale2);
                h2Ratio->GetYaxis()->SetLabelSize(y_label*scale2);
                h2Ratio->GetYaxis()->SetRangeUser(0.5,1.5);
                h2Ratio->GetXaxis()->SetTitleSize(x_title*scale2);
                h2Ratio->GetYaxis()->SetTitleSize(y_title*scale2);
                h2Ratio->GetYaxis()->SetTitleOffset(y_offset/scale2);
                h2Ratio->GetYaxis()->SetNdivisions(8);
                
                //Two pads
                pad1 = new TPad("pad1","pad1",0,1-size1,1,1);
                pad2 = new TPad("pad2","pad2",0,0,1,size2);
                pad1->SetBottomMargin(0.1);
                pad2->SetBottomMargin(0.4);
                pad2->SetTopMargin(0.1);
                pad2->SetGridy();
            }
            else
            {
                //x-asix title for pad1
                h2DataSum->GetXaxis()->SetTitle(xaxis.Data());
                
                //pad1
                pad1 = new TPad("pad1","pad1",0,0,1,1);
            }
            pad1->SetLogy(Var[VarIndex].log);
            
            //Draw for pad1
            gStyle->SetOptStat(0);
            TCanvas* c2 = new TCanvas();
            c2->cd();
            pad1->Draw();
            pad1->cd();
            h2DataSum->Draw("axis");
            stack.Draw("histsame");
            for(unsigned int i=0;i<SigMassSplitting.size();i++)
            {
                h2SigSum[i]->Draw("histsame");
            }
            if(RegionInfo[RegionIndex].showData)
            {
                h2DataSum->Draw("Esame");
            }
            h2DataSum->Draw("sameaxis");
            leg->Draw();
            
            {
                //text
                ATLASLabel(0.3,0.88,"Internal");
                
                TLatex lt2;
                lt2.DrawLatexNDC(0.3,0.83, "#sqrt{#it{s}} = 13 TeV, 10.6 fb^{-1}");
                lt2.SetTextSize(lt2.GetTextSize());
                
                TLatex lt1;
                lt1.DrawLatexNDC(0.3,0.78,RegionInfo[RegionIndex].RegionName.Data());
                lt1.SetTextSize(lt1.GetTextSize());
            }
            
            if(RegionInfo[RegionIndex].showData)
            {
                //Draw for pad2
                c2->cd();
                pad2->Draw();
                pad2->cd();
                h2Ratio->Draw();
                
                TLine l;
                TLine* l1 = l.DrawLine(h2Ratio->GetXaxis()->GetXmin(), 1., h2Ratio->GetXaxis()->GetXmax(), 1.);
                l1->SetLineStyle(2);
                l1->SetLineWidth(2);
            }
            
            {
                //export histograms in eps format
                TString NameTemp = "plot/";
                NameTemp += Var[VarIndex].VarName;
                NameTemp += "_";
                NameTemp += RegionInfo[RegionIndex].RegionName;
                NameTemp += ".eps";
                c2->Print(NameTemp,"eps");
            }
            
            //c2->WaitPrimitive();
            
            //delete
            //h2DataSum
            delete h2DataSum;
            
            //h2 for BGGruop
            for(unsigned int j=0;j<BGGroup.size();j++)
            {
                delete BGGroup[j].h2;
            }
            
            //h2BGSum
            delete h2BGSum;
            
            //h2Sig
            for(unsigned int j=0;j<SigSampleID.size();j++)
            {
                delete h2Sig[j];
            }
            
            //h2SigSum
            for(unsigned int j=0;j<SigMassSplitting.size();j++)
            {
            	delete h2SigSum[j];
            }
            
            //legend
            delete leg;
            
            if(RegionInfo[RegionIndex].showData)
            {
                //h2Ratio
                delete h2Ratio;
                
                //pad2
                delete pad2;
            }
            
            //pad1
            delete pad1;
            
            //canvas
            delete c2;
        }
        
        //delete tree2
        for(unsigned int i=0;i<DataSampleID.size();i++)
        {
            delete tree2Data[i];
            if(RegionInfo[RegionIndex].isSS_ee) delete tree2DataOS[i];
        }
        for(unsigned int j=0;j<tree2BGMC.size();j++)
        {
            for(unsigned int k=0;k<tree2BGMC[j].size();k++)
            {
                //Z pt reweighting
                if(dorw
                   &&
                   (
                    RegionIndex==0 ||
                    RegionIndex==1 ||
                    RegionIndex==6 ||
                    RegionIndex==7
                   )
                   &&
                   (
                    BGGroup[j].info->GroupName == "Zee" ||
                    BGGroup[j].info->GroupName == "Zmumu" ||
                    BGGroup[j].info->GroupName == "Ztautau"
                   )
                  )
                {
                    TString NameTemp = "tree_rw_";
                    NameTemp += TString::Itoa(k + BGGroup[j].info->lower,10);
                    delete tree2BGMC[j][k]->GetFriend(NameTemp.Data());
                }
                
                /*
                //for charge filp BG
                if(RegionIndex>=12 && RegionIndex<=14 &&
                   BGGroup[j].info->GroupName == "Zee")
                {
                    TString NameTemp = "tree_cfw_";
                    NameTemp += TString::Itoa(i,10);
                    delete tree2BGMC[i]->GetFriend(NameTemp.Data());
                }
                */
                
                delete tree2BGMC[j][k];
            }
        }
        for(unsigned int i=0;i<SigSampleID.size();i++)
        {
            delete tree2Sig[i];
        }
    }
    
    //latex for tables
    //latex for expN
    for(unsigned int ISR=0;ISR<=6;ISR+=6)
    {
        TString PathName = "latex/data/";
        PathName += "expN_";
        if(ISR==0) PathName += "non";
        PathName += "ISR.tex";
        
        ofstream fout;
        fout.open(PathName.Data());
        
        for(unsigned int SixChannel=0;SixChannel<6;SixChannel++)
        {
            fout<<"\\begin{frame}"<<endl;
            fout<<"\\frametitle{Expected number of events (For ";
            if(ISR==0) fout<<"non";
            fout<<"ISR)}"<<endl;
            
            fout<<"For ";
            if(SixChannel%3 == 0) fout<<"ee";
            else if(SixChannel%3 == 1) fout<<"$\\mu\\mu$";
            else if(SixChannel%3 == 2) fout<<"e$\\mu$";
            fout<<" channel, ";
            if(SixChannel<=2) fout<<"opposite";
            else fout<<"same";
            fout<<" sign \\\\"<<endl;
            
            fout<<"\\vspace{5mm}"<<endl;
            fout<<"\\begin{tabular}{|c|c|c|}"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"& Number of events & Significance \\\\"<<endl;
            fout<<"\\hline"<<endl;
            
            fout<<"\\input{data/expN/"<<channel[SixChannel+ISR].Data()<<".tex}"<<endl;
            
            fout<<"\\end{tabular}"<<endl;
            fout<<"\\end{frame}"<<endl<<endl;
        }
        
        fout.close();
        
        //For VV
        PathName = "latex/data/";
        PathName += "expN_BGVV_";
        if(ISR==0) PathName += "non";
        PathName += "ISR.tex";
        fout.open(PathName.Data());
        
        for(unsigned int SixChannel=0;SixChannel<6;SixChannel++)
        {
            fout<<"\\begin{frame}"<<endl;
            fout<<"\\frametitle{Expected number of events for VV (For ";
            if(ISR==0) fout<<"non";
            fout<<"ISR)}"<<endl;
            
            fout<<"For ";
            if(SixChannel%3 == 0) fout<<"ee";
            else if(SixChannel%3 == 1) fout<<"$\\mu\\mu$";
            else if(SixChannel%3 == 2) fout<<"e$\\mu$";
            fout<<" channel, ";
            if(SixChannel<=2) fout<<"opposite";
            else fout<<"same";
            fout<<" sign \\\\"<<endl;
            
            fout<<"\\vspace{5mm}"<<endl;
            fout<<"\\begin{tabular}{|c|c|}"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"& Number of events \\\\"<<endl;
            fout<<"\\hline"<<endl;
            
            fout<<"\\input{data/expN/BGVV_"<<channel[SixChannel+ISR].Data()<<".tex}"<<endl;
            
            fout<<"\\end{tabular}"<<endl;
            fout<<"\\end{frame}"<<endl<<endl;
        }
        
        fout.close();
    }
    
    //latex for expN for CR
    {
        TString PathName = "latex/data/expN_CR.tex";
        ofstream fout;
        fout.open(PathName.Data());
        
        for(unsigned int RegionIndex=0;RegionIndex<6;RegionIndex++)
        {
            TString latexName = RegionInfo[12+RegionIndex].RegionName;
            latexName.ReplaceAll("_","\\_");
            
            fout<<"\\begin{frame}"<<endl;
            fout<<"\\frametitle{Expected number of events (For ";
            fout<<latexName.Data();
            fout<<")}"<<endl;
            
            fout<<"For ";
            fout<<latexName.Data();
            fout<<",\\\\"<<endl;
            
            fout<<"\\vspace{5mm}"<<endl;
            fout<<"\\begin{tabular}{|c|c|c|}"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"& Number of events & Significance \\\\"<<endl;
            fout<<"\\hline"<<endl;
            
            fout<<"\\input{data/expN/"<<RegionInfo[12+RegionIndex].RegionName.Data()<<".tex}"<<endl;
            
            fout<<"\\end{tabular}"<<endl;
            fout<<"\\end{frame}"<<endl<<endl;
        }
        
        fout.close();
    }
    
    //latex for plot
    //plot_nonISR.tex and plot_ISR.tex
    for(unsigned int ISR=0;ISR<=6;ISR+=6)
    {
        TString PathName = "latex/data/";
        PathName += "plot_";
        if(ISR==0) PathName += "non";
        PathName += "ISR.tex";
        
        ofstream fout;
        fout.open(PathName.Data());
        
        for(unsigned int VarIndex=0;VarIndex<Var.size();VarIndex++)
        {
            if(Var[VarIndex].VarName=="averageMu") continue;
            if(Var[VarIndex].VarName=="nVtx") continue;
            if( ISR==0 &&
                (Var[VarIndex].VarName=="bjetpt"  ||
                 Var[VarIndex].VarName=="bjeteta" ||
                 Var[VarIndex].VarName=="bjetphi" ||
                 Var[VarIndex].VarName=="cjetpt"  ||
                 Var[VarIndex].VarName=="cjeteta" ||
                 Var[VarIndex].VarName=="cjetphi" )
              ) continue;
            
            fout<<"\\begin{frame}"<<endl;
            
            fout<<"\\frametitle{"<<Var[VarIndex].VarTitle.Data()<<" (For ";
            if(ISR==0) fout<<"non";
            fout<<"ISR)}"<<endl;
            
            fout<<"\\Wider[5em]{"<<endl;
            for(unsigned int SixChannel=0;SixChannel<6;SixChannel++)
            {
                fout<<"\\includegraphics[width=0.33\\textwidth]{\\PathToPlot/"
                    <<Var[VarIndex].VarName.Data()<<"_"<<channel[SixChannel+ISR].Data()<<"}";
                if(SixChannel==2) fout<<" \\\\";
                fout<<endl;
            }
            //fout<<"\\caption{"<<Var[VarIndex].latexName.Data()<<" for ee channel (left), $\\mu\\mu$ channel (middle) and e$\\mu$ channel (right), for opposite side (top) and same sign (bottom).}"<<endl;
            fout<<"}"<<endl;
            
            fout<<"\\end{frame}"<<endl<<endl;
        }
        
        fout.close();
    }
    
    //plot_CR.tex
    {
        TString PathName = "latex/data/plot_CR.tex";
        ofstream fout;
        fout.open(PathName.Data());
        
        for(unsigned int VarIndex=0;VarIndex<Var.size();VarIndex++)
        {
            if(Var[VarIndex].VarName=="averageMu") continue;
            if(Var[VarIndex].VarName=="nVtx") continue;
            if(
               (Var[VarIndex].VarName=="bjetpt"  ||
                Var[VarIndex].VarName=="bjeteta" ||
                Var[VarIndex].VarName=="bjetphi" ||
                Var[VarIndex].VarName=="cjetpt"  ||
                Var[VarIndex].VarName=="cjeteta" ||
                Var[VarIndex].VarName=="cjetphi" )
               ) continue;
            
            fout<<"\\begin{frame}"<<endl;
            
            fout<<"\\frametitle{"<<Var[VarIndex].VarTitle.Data()<<" (For CR)}"<<endl;;
            
            fout<<"\\Wider[5em]{"<<endl;
            for(unsigned int RegionIndex=0;RegionIndex<6;RegionIndex++)
            {
                fout<<"\\includegraphics[width=0.33\\textwidth]{\\PathToPlot/"
                <<Var[VarIndex].VarName.Data()<<"_"<<RegionInfo[12+RegionIndex].RegionName.Data()<<"}";
                if(RegionIndex==2) fout<<" \\\\";
                fout<<endl;
            }
            //fout<<"\\caption{"<<Var[VarIndex].latexName.Data()<<" for ee channel (left), $\\mu\\mu$ channel (middle) and e$\\mu$ channel (right), for opposite side (top) and same sign (bottom).}"<<endl;
            fout<<"}"<<endl;
            
            fout<<"\\end{frame}"<<endl<<endl;
        }
        
        fout.close();
    }
    
    //plot_special.tex
    {
        TString PathName = "latex/data/plot_special.tex";
        ofstream fout;
        fout.open(PathName.Data());
        for(unsigned int VarIndex=0;VarIndex<Var.size();VarIndex++)
        {
            for(unsigned int ISR=0;ISR<=6;ISR+=6)
            {
                if(!
                   (Var[VarIndex].VarName=="pt1"   ||
                    Var[VarIndex].VarName=="pt2"   ||
                    Var[VarIndex].VarName=="mll"   ||
                    Var[VarIndex].VarName=="ptll"  ||
                    Var[VarIndex].VarName=="MET"   ||
                    Var[VarIndex].VarName=="mTtwo" ||
                    Var[VarIndex].VarName=="mT1"   ||
                    Var[VarIndex].VarName=="mT2"   ||
                    Var[VarIndex].VarName=="HT"    )
                   )continue;
                
                if( ISR==0 &&
                   (Var[VarIndex].VarName=="bjetpt"  ||
                    Var[VarIndex].VarName=="bjeteta" ||
                    Var[VarIndex].VarName=="bjetphi" ||
                    Var[VarIndex].VarName=="cjetpt"  ||
                    Var[VarIndex].VarName=="cjeteta" ||
                    Var[VarIndex].VarName=="cjetphi" )
                   ) continue;
                
                fout<<"\\begin{frame}"<<endl;
                
                fout<<"\\frametitle{"<<Var[VarIndex].VarTitle.Data()<<" (For ";
                if(ISR==0) fout<<"non";
                fout<<"ISR)}"<<endl;
                
                fout<<"\\Wider[5em]{"<<endl;
                for(unsigned int SixChannel=0;SixChannel<6;SixChannel++)
                {
                    fout<<"\\includegraphics[width=0.33\\textwidth]{\\PathToPlot/"
                    <<Var[VarIndex].VarName.Data()<<"_"<<channel[SixChannel+ISR].Data()<<"}";
                    if(SixChannel==2) fout<<" \\\\";
                    fout<<endl;
                }
                //fout<<"\\caption{"<<Var[VarIndex].latexName.Data()<<" for ee channel (left), $\\mu\\mu$ channel (middle) and e$\\mu$ channel (right), for opposite side (top) and same sign (bottom).}"<<endl;
                fout<<"}"<<endl;
                
                fout<<"\\end{frame}"<<endl<<endl;
            }
        }
        fout.close();
    }
}
