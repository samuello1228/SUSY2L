#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>

#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
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

#include <TGraph2D.h>

const bool dorw = 0;
const bool direct = 1;
const bool fitting = 0;
const bool simple = 1;
const bool combined = 0;

const bool docfw = 0;

const bool doOptimize = 0;
const unsigned int SigOptimizingIndex = 0;

const bool useDani = 0;

TString setupTag = "v7";

// Cutflow Attention
const bool doVVCount = 0;

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

struct ChannelData
{
    TString ChannelName;
    bool isSS;
    bool isSS_qF;
    unsigned int qFChannel;
    std::vector<TString> setOfBGMC;
    std::vector<TString> setOfBGData;
};

void initializeTree1new(TChain*& tree, TString& SampleName, TString& channel)
{
    tree = new TChain("tree");
    
    TString fileName = "skimming/skimming.";
    fileName += SampleName;
    fileName += "_";
    fileName += channel;
    fileName += ".root";
    //cout<<fileName<<endl;
    
    tree->Add(fileName.Data());
    
    tree->SetBranchAddress("ptll", &ptll, &b_ptll);
    
    tree->SetBranchAddress("pt1", &pt1, &b_pt1);
    tree->SetBranchAddress("pt2", &pt2, &b_pt2);
    tree->SetBranchAddress("eta1", &eta1, &b_eta1);
    tree->SetBranchAddress("eta2", &eta2, &b_eta2);
    
}

void initializeTree2(std::vector<TChain*>& tree2,std::vector<unsigned int>& SetOfChannel, std::vector<TString>& SampleName, std::vector<ChannelData>& ChannelInfo)
{
    for(unsigned int i=0;i<SampleName.size();i++)
    {
        TChain* treeTemp = new TChain("tree");
        for(unsigned int j=0;j<SetOfChannel.size();j++)
        {
            TString fileName = "skimming/skimming.";
            fileName += SampleName[i];
            fileName += "_";
            fileName += ChannelInfo[SetOfChannel[j]].ChannelName;
            fileName += ".root";
            
            treeTemp->Add(fileName.Data());
        }
        
        treeTemp->SetBranchAddress("weight", &weight, &b_weight);
        tree2.push_back(treeTemp);
    }
}

void getnwAOD(std::vector<double>& BGMCnwAOD,std::vector<unsigned int>& SetOfChannel, std::vector<TString>& SampleName, std::vector<ChannelData>& ChannelInfo)
{
    for(unsigned int i=0;i<SampleName.size();i++)
    {
        TString fileName = "skimming/skimming.";
        fileName += SampleName[i];
        fileName += "_";
        fileName += ChannelInfo[SetOfChannel[0]].ChannelName;
        fileName += ".root";
        
        TFile* file = new TFile(fileName.Data(),"READ");
        TH1D *h1 = (TH1D*) file->Get("hist");
        double nAOD = h1->GetBinContent(2);
        delete file;
        
        if(nAOD == 0)
        {
            cout<<"nAOD is 0. The sample: "<<SampleName[i]<<" is missing."<<endl;
            nAOD=1;
        }
        BGMCnwAOD.push_back(nAOD);
    }
}

void getnAOD(std::vector<int>& BGMCnAOD,std::vector<unsigned int>& SetOfChannel, std::vector<TString>& SampleName, std::vector<ChannelData>& ChannelInfo)
{
    for(unsigned int i=0;i<SampleName.size();i++)
    {
        TString fileName = "skimming/skimming.";
        fileName += SampleName[i];
        fileName += "_";
        fileName += ChannelInfo[SetOfChannel[0]].ChannelName;
        fileName += ".root";
        
        TFile* file = new TFile(fileName.Data(),"READ");
        TH1D *h1 = (TH1D*) file->Get("hist");
        int nAOD = h1->GetBinContent(1);
        delete file;
        
        if(nAOD == 0)
        {
            cout<<"nAOD is 0. The sample: "<<SampleName[i]<<" is missing."<<endl;
            nAOD=1;
        }
        BGMCnAOD.push_back(nAOD);
    }
}

struct MCBGInfo
{
    TString SampleName;
    unsigned int SampleID;
    double XS; //cross section in pb
};

void fillBGMCIndex(std::vector<MCBGInfo> const& BGMCSampleInfo, std::vector<unsigned int>& BGMCIndex,unsigned int lower, unsigned int upper)
{
    for(unsigned int i=lower;i<=upper;i++)
    {
        bool isFound = false;
        for(unsigned int j=0;j<BGMCSampleInfo.size();j++)
        {
            if(i == BGMCSampleInfo[j].SampleID)
            {
                BGMCIndex.push_back(j);
                isFound = true;
            }
        }
        if(!isFound) cout<<"The sample ID "<<i<<" cannot be found."<<endl;
    }
}

double GetSignificance(double nSig, double nBG, double nBGError2 = 0)
{
    //return RooStats::NumberCountingUtils::BinomialExpZ(nSig,nBG,0.3);
    
    double Error = sqrt((0.25*0.25) + nBGError2/(nBG*nBG));
    return RooStats::NumberCountingUtils::BinomialExpZ(nSig,nBG,Error);
}

struct VarData
{
    TString VarName;
    TString VarFormula;
    TString VarTitle;
    TString unit;
    TString latexName;
    int CutDirection;
    /*
    0: none
    1: x > 100 GeV
    -1: x < 100 GeV
    */
    
    int bin;
    double xmin;
    double xmax;
    
    bool log;
    double ymin;
    double ymax;
};

unsigned int findVarIndex(TString const& VarName, std::vector<VarData>& Var)
{
    for(unsigned int i=0;i<Var.size();i++)
    {
        if(VarName == Var[i].VarName)
        {
            return i;
        }
    }
    return 0;
}

struct GroupData
{
    TString GroupName;
    TString LegendName;
    TString LatexName;
    std::vector<unsigned int> BGMCIndex;
    int colour;
    int unweighted;
    double weighted;
    double error;
};

struct Group
{
    GroupData* info;
    TH1D* h2;
    TH2D* h3;
};

bool compare2(Group Group1,Group Group2)
{
    return Group1.h2->GetSumOfWeights() < Group2.h2->GetSumOfWeights();
}

bool NextBin(int& bin1, int& bin2, int const nBin, bool const OnlyHasLowerCut, bool const OnlyHasUpperCut)
{
    if(bin2 == nBin || OnlyHasLowerCut)
    {
        if(OnlyHasUpperCut) return false;
        
        if(bin1 == nBin)
        {
            bin1 = nBin+1;
        }
        else if(bin1 == nBin+1) //For OnlyHasLowerCut
        {
            return false;
        }
        else
        {
            bin1++;
        }
        
        bin2 = nBin+1;
    }
    else if(bin2 == nBin+1)
    {
        if(bin1 == nBin+1)
        {
            return false;
        }
        else
        {
            bin2 = bin1;
        }
    }
    else
    {
        bin2++;
    }
    
    return true;
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
    std::vector<ChannelData> ChannelInfo;
    {
        ChannelData element;
        //TString ISR[2] = {"nonISR","ISR"};
        TString JET[3] = {"jet1","jet23","jet4"};
        TString sign[2] = {"OS","SS"};
        TString lepton[3] = {"ee","mumu","emu"};
        //for(int i=0;i<2;i++)
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<2;j++)
            {
                for(int k=0;k<3;k++)
                {
                    element.ChannelName = "";
                    //element.ChannelName += ISR[i];
                    element.ChannelName += JET[i];
                    element.ChannelName += "_";
                    element.ChannelName += sign[j];
                    element.ChannelName += "_";
                    element.ChannelName += lepton[k];
                    
                    element.isSS = j==1;
                    element.isSS_qF = (j==1 && (k==0 || k==2));
                    if(element.isSS_qF) element.qFChannel = ChannelInfo.size() - 3;
                    
                    
                    element.setOfBGMC.clear();
                    element.setOfBGData.clear();
                    if(element.isSS)
                    {
                        ///* Cutflow Attention
                        //element.setOfBGData.push_back("fake lepton");
                        if(element.isSS_qF)
                        {
                            //element.setOfBGData.push_back("charge flip");
                        }
                        //*/
                        
                        /*
                        element.setOfBGMC.push_back("Zjets");
                        element.setOfBGMC.push_back("Wjets");
                        
                        element.setOfBGMC.push_back("ttbar");
                        element.setOfBGMC.push_back("singletop");
                        element.setOfBGMC.push_back("ttV");
                        element.setOfBGMC.push_back("multitop");
                        
                        element.setOfBGMC.push_back("VV");
                        element.setOfBGMC.push_back("Vgamma");
                        element.setOfBGMC.push_back("VVV");
                        element.setOfBGMC.push_back("Higgs");
                        */
                        
                        ///*
                        //SR
                        element.setOfBGMC.push_back("WZ");
                        element.setOfBGMC.push_back("WW");
                        element.setOfBGMC.push_back("ZZ");
                        element.setOfBGMC.push_back("ttV");
                        element.setOfBGMC.push_back("Rare");
                        element.setOfBGMC.push_back("Fakes");
                        element.setOfBGMC.push_back("charge flip");
                        //*/
                        
                        /*
                        //Dani
                        element.setOfBGMC.push_back("Fakes");
                        element.setOfBGMC.push_back("VV_DL");
                        element.setOfBGMC.push_back("ttV");
                        element.setOfBGMC.push_back("VVV");
                        element.setOfBGMC.push_back("Higgs");
                        element.setOfBGMC.push_back("multitop");
                        element.setOfBGMC.push_back("DY");
                        */
                        
                        /*
                        //VR_Fake
                        element.setOfBGMC.push_back("WZ");
                        element.setOfBGMC.push_back("WW");
                        element.setOfBGMC.push_back("ZZ");
                        element.setOfBGMC.push_back("ttV");
                        element.setOfBGMC.push_back("Rare");
                        element.setOfBGMC.push_back("Samuel_Fake");
                        element.setOfBGMC.push_back("charge flip");
                        */
                    }
                    else
                    {
                        ///* Cutflow Attention
                        element.setOfBGMC.push_back("Zjets");
                        element.setOfBGMC.push_back("Wjets");
                        element.setOfBGMC.push_back("top");
                        element.setOfBGMC.push_back("VV");
                        element.setOfBGMC.push_back("Vgamma");
                        element.setOfBGMC.push_back("VVV");
                        //*/
                    }

                    ChannelInfo.push_back(element);
                }
            }
        }
    }
    
    //For Data
    std::vector<TString> DataSampleName;
    DataSampleName.reserve(20);
    double sumDataL = 0; //cross section in pb-1
    
    {
        //read DataSample.txt
        ifstream fin;
        fin.open("DataSample.txt");
        while(!fin.eof())
        {
            TString SampleNameTemp;
            fin>>SampleNameTemp;
            if(fin.eof()) break;
            
            DataSampleName.push_back(SampleNameTemp);
            
            double DataLTemp;
            fin>>DataLTemp;
            sumDataL += DataLTemp;
            
            long nAOD;
            fin>>nAOD;
            fin>>nAOD;
        }
        fin.close();
    }
    //DataSampleName.clear(); sumDataL = 36100; //Cutflow Attention
    cout<<"Total Luminosity: "<<sumDataL<<endl;
    
    //For BGMC
    std::vector<MCBGInfo> BGMCSampleInfo;
    {
        MCBGInfo element;
        
        //read BGSample.txt
        ifstream fin;
        fin.open("BGSample.txt");
        while(!fin.eof())
        {
            TString SampleNameTemp = "mc15_13TeV.";
            TString SampleNameTemp2;
            unsigned int SampleIDTemp;
            fin>>SampleIDTemp;
            if(fin.eof()) break;
            element.SampleID = SampleIDTemp;
            
            SampleNameTemp += TString::Itoa(SampleIDTemp,10);
            
            fin>>SampleNameTemp2;
            SampleNameTemp += ".";
            SampleNameTemp += SampleNameTemp2;
            element.SampleName = SampleNameTemp;
            
            double BGMCXSTemp;
            double BGMCXSTemp2;
            fin>>BGMCXSTemp2;
            fin>>BGMCXSTemp;
            BGMCXSTemp2 *= BGMCXSTemp;
            fin>>BGMCXSTemp;
            BGMCXSTemp2 *= BGMCXSTemp;
            element.XS = BGMCXSTemp2;
            
            fin>>BGMCXSTemp;
            BGMCSampleInfo.push_back(element);
            
            long nAOD;
            fin>>nAOD;
            fin>>nAOD;
        }
        fin.close();
    }
    cout<<"Total number of BG files: "<<BGMCSampleInfo.size()<<endl;
    
    //For Signal MC
    struct PlotSigInfo
    {
        double MassDiff;
        unsigned int ID;
        TString IDName;
        int colour;
        int linestyle;
        int unweighted;
        double weighted;
        double error;
        double significance;
        double scale;
    };
    std::vector<PlotSigInfo> SigMassSplitting;
    {
        PlotSigInfo element;
        
        ///*
        //Dani
        element.MassDiff = 150;   element.ID = 20; element.colour = 1;   element.linestyle = 2; element.scale = 1;  SigMassSplitting.push_back(element);
        element.MassDiff = 150;   element.ID = 7;  element.colour = 920; element.linestyle = 5; element.scale = 1;  SigMassSplitting.push_back(element);
        element.MassDiff = 175;   element.ID = 3;  element.colour = 922; element.linestyle = 9; element.scale = 1;  SigMassSplitting.push_back(element);
        //*/
    }
    
    struct SigInfo
    {
        TString SampleName;
        double XS; //cross section in pb
        double nwAOD;
        double Mass1;
        double Mass2;
        int unweighted;
        double weighted;
        double error;
        double significance;
        double significance2; //for combined significance
    };
    
    std::vector<SigInfo> SigSampleInfo;
    {
        SigInfo element;
        
        //read SigSample.txt
        ifstream fin;
        fin.open("SigSample.txt");
        while(!fin.eof())
        {
            TString SampleNameTemp = "mc15_13TeV.";
            //for 125
            TString SampleNameTemp2;
            fin>>SampleNameTemp2;
            if(fin.eof()) break;
            
            SampleNameTemp += SampleNameTemp2;
            
            fin>>SampleNameTemp2;
            SampleNameTemp += ".";
            SampleNameTemp += SampleNameTemp2;
            
            element.SampleName = SampleNameTemp;
            
            double SigMass;
            fin>>SigMass;
            element.Mass1 = SigMass;
            fin>>SigMass;
            element.Mass2 = SigMass;
            
            double SigXSTemp2;
            fin>>SigXSTemp2;
            
            double SigXSTemp;
            fin>>SigXSTemp;
            //SigXSTemp2 *= SigXSTemp;
            
            fin>>SigXSTemp;
            SigXSTemp2 *= SigXSTemp;
            
            fin>>SigXSTemp;
            SigXSTemp2 *= SigXSTemp;

            //cout<<SigXSTemp2<<endl;
            element.XS = SigXSTemp2;
            SigSampleInfo.push_back(element);
            
            long nAOD;
            fin>>nAOD;
            fin>>nAOD;
        }
        fin.close();
    }
    
    for(unsigned int i=0;i<SigMassSplitting.size();i++)
    {
        cout<<"Mass splitting: "<<SigMassSplitting[i].MassDiff<<endl;
        cout<<"index:"<<SigMassSplitting[i].ID<<endl;
        cout<<"Name: "<<SigSampleInfo[SigMassSplitting[i].ID].SampleName.Data()<<endl;
        cout<<"MassDiff: "<<SigSampleInfo[SigMassSplitting[i].ID].Mass1 - SigSampleInfo[SigMassSplitting[i].ID].Mass2<<endl;
        cout<<"XS: "<<SigSampleInfo[SigMassSplitting[i].ID].XS<<endl<<endl;
        
        cout<<"All samples with the same mass splitting "<<SigMassSplitting[i].MassDiff<<" GeV:"<<endl;
        for(unsigned int j=0;j<SigSampleInfo.size();j++)
        {
            if(SigSampleInfo[j].Mass1 - SigSampleInfo[j].Mass2 == SigMassSplitting[i].MassDiff) cout<<"index:"<<j<<" "<<SigSampleInfo[j].SampleName.Data()<<endl;
        }
        cout<<endl;
    }
    
    for(unsigned int i=0;i<SigMassSplitting.size();i++)
    {
        TString GroupName = "(";
        GroupName += TString::Format("%.1f",SigSampleInfo[SigMassSplitting[i].ID].Mass1);
        GroupName += ",";
        GroupName += TString::Format("%.1f",SigSampleInfo[SigMassSplitting[i].ID].Mass2);
        GroupName += ")";
        SigMassSplitting[i].IDName = GroupName;
    }
    
    //Get number of events in AOD
    for(unsigned int i=0;i<SigSampleInfo.size();i++)
    {
        TString NameTemp = "skimming/skimming.";
        NameTemp += SigSampleInfo[i].SampleName;
        cout<<NameTemp<<": ";
        NameTemp += "_";
        NameTemp += ChannelInfo[0].ChannelName;
        NameTemp += ".root";
        
        TFile* file = new TFile(NameTemp.Data(),"READ");
        TH1D *h1 = (TH1D*) file->Get("hist");
        SigSampleInfo[i].nwAOD = h1->GetBinContent(2);
        cout<<SigSampleInfo[i].nwAOD<<endl;
        
        delete file;
        
        SigSampleInfo[i].significance2 = 0;
    }
    
    //Group for MC background
    std::vector<GroupData> BGMCGroupData;
    {
        GroupData element;
        
        //Z+jets
        element.GroupName = "Zjets"; element.LegendName = "Z+jets"; element.LatexName = "Z+jets";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,364100,364141);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,364198,364215);
        element.colour = 2; BGMCGroupData.push_back(element);
        
        element.GroupName = "Zee"; element.LegendName = "Z#rightarrow ee"; element.LatexName = "Z$\\rightarrow ee$";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,364114,364127);
        element.colour = 2; BGMCGroupData.push_back(element);
        
        element.GroupName = "Zmumu"; element.LegendName = "Z#rightarrow #mu#mu"; element.LatexName = "Z$\\rightarrow\\mu\\mu$";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,364100,364113);
        element.colour = 2; BGMCGroupData.push_back(element);
        
        element.GroupName = "Ztautau"; element.LegendName = "Z#rightarrow #tau#tau"; element.LatexName = "Z$\\rightarrow\\tau\\tau$";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,364128,364141);
        element.colour = 2; BGMCGroupData.push_back(element);
        
        //Drell Yan
        element.GroupName = "DY"; element.LegendName = "Drell Yan"; element.LatexName = "Drell Yan";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,364198,364215);
        element.colour = 2; BGMCGroupData.push_back(element);
        
        //W+jets
        element.GroupName = "Wjets"; element.LegendName = "W+jets"; element.LatexName = "W+jets";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,364156,364197);
        element.colour = 3; BGMCGroupData.push_back(element);
        
        //Top
        element.GroupName = "top"; element.LegendName = "top"; element.LatexName = "top";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410000,410000);
        element.colour = 4; BGMCGroupData.push_back(element);
        
        element.GroupName = "ttbar"; element.LegendName = "t#bar{t}"; element.LatexName = "$t\\bar{t}$";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410000,410000);
        element.colour = 881; BGMCGroupData.push_back(element);

        element.GroupName = "singletop"; element.LegendName = "single top"; element.LatexName = "single top";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410011,410012);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410025,410026);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410015,410016);
        element.colour = 30; BGMCGroupData.push_back(element);
        
        element.GroupName = "Wt"; element.LegendName = "Wt"; element.LatexName = "Wt";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410015,410016);
        element.colour = 30; BGMCGroupData.push_back(element);
        
        element.GroupName = "ttV"; element.LegendName = "t#bar{t}+V"; element.LatexName = "$t\\bar{t}+V$";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410218,410220);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410155,410155);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410081,410081);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,407321,407321);
        element.colour = kOrange-2; BGMCGroupData.push_back(element);
        
        element.GroupName = "multitop"; element.LegendName = "multi top"; element.LatexName = "multi top";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,304014,304014);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410080,410080);
        element.colour = 635; BGMCGroupData.push_back(element);
        
        //VV
        element.GroupName = "VV";     element.LegendName = "VV"; element.LatexName = "VV";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,361069,361073);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,361077,361077);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363356,363356);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363358,363360);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363489,363492);
        element.colour = 801; BGMCGroupData.push_back(element);
        
        element.GroupName = "VVReal"; element.LegendName = "VV"; element.LatexName = "VV";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,361069,361073);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363490,363491);
        element.colour = 801; BGMCGroupData.push_back(element);
        
        element.GroupName = "VV_qF"; element.LegendName = "VV_qF"; element.LatexName = "VV\\_qF";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,361077,361077);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363356,363356);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363358,363358);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363492,363492);
        element.colour = 801; BGMCGroupData.push_back(element);
        
        element.GroupName = "VV_DL"; element.LegendName = "VV_DL"; element.LatexName = "VV\\_DL";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,361069,361073);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,361077,361077);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363356,363356);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363358,363358);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363490,363492);
        element.colour = 801; BGMCGroupData.push_back(element);
        
        element.GroupName = "WZ"; element.LegendName = "WZ"; element.LatexName = "WZ";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,361071,361071);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363491,363491);
        element.colour = kRed+1; BGMCGroupData.push_back(element);
        
        element.GroupName = "WW"; element.LegendName = "WW"; element.LatexName = "WW";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,361069,361070);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,361077,361077);
        element.colour = kRed-9; BGMCGroupData.push_back(element);
        
        element.GroupName = "ZZ"; element.LegendName = "ZZ"; element.LatexName = "ZZ";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,361072,361073);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363490,363490);
        element.colour = kOrange-3; BGMCGroupData.push_back(element);
        
        {
            element.GroupName = "Rare"; element.LegendName = "Rare"; element.LatexName = "Rare";
            element.BGMCIndex.clear();
            
            //higgs
            fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,341177,341177);
            fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,341270,341270);
            fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,342284,342285);
            
            //VVV
            fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,407311,407315);
            
            //multi top
            fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,304014,304014);
            fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410080,410080);
            
            //Two Rare
            fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,343272,343272);
            fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410215,410215);
            element.colour = kViolet-8; BGMCGroupData.push_back(element);
        }
        
        //Vgamma
        element.GroupName = "Vgamma"; element.LegendName = "V + #gamma"; element.LatexName = "V$+\\gamma$";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,301535,301536);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,301890,301898);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,301899,301907);
        element.colour = 5; BGMCGroupData.push_back(element);
        
        element.GroupName = "Wgamma"; element.LegendName = "W + #gamma"; element.LatexName = "W$+\\gamma$";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,301890,301898);
        element.colour = 5; BGMCGroupData.push_back(element);
        
        element.GroupName = "Zgamma"; element.LegendName = "Z + #gamma"; element.LatexName = "Z$+\\gamma$";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,301535,301536);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,301899,301907);
        element.colour = 5; BGMCGroupData.push_back(element);
        
        //VVV
        element.GroupName = "VVV"; element.LegendName = "VVV"; element.LatexName = "VVV";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,407311,407315);
        element.colour = 607; BGMCGroupData.push_back(element);
        
        //Higgs
        element.GroupName = "Higgs"; element.LegendName = "Higgs"; element.LatexName = "Higgs";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,341079,341079);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,341122,341122);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,341195,341195);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,342178,342178);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,341080,341080);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,341155,341155);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,341206,341206);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,342189,342189);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,342284,342285);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,341177,341177);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,341270,341271);
        element.colour = 7; BGMCGroupData.push_back(element);
        
        //Fakes
        element.GroupName = "Fakes"; element.LegendName = "Fakes"; element.LatexName = "Fakes";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,999999,999999);
        element.colour = 3; BGMCGroupData.push_back(element);
        
        //charge flip
        element.GroupName = "charge flip"; element.LegendName = "charge flip"; element.LatexName = "charge flip";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,999998,999998);
        element.colour = kOrange-5; BGMCGroupData.push_back(element);
        
        //Samuel Fake
        element.GroupName = "Samuel_Fake"; element.LegendName = "Fakes"; element.LatexName = "Fakes";
        element.BGMCIndex.clear();
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,999997,999997);
        element.colour = kGray; BGMCGroupData.push_back(element);
        
        /*
        //charge flip
        element.GroupName = "charge flip"; element.LegendName = "charge flip"; element.LatexName = "charge flip";
        element.BGMCIndex.clear();
        //VV
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,361077,361077);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363356,363356);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363358,363358);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,363492,363492);
        //top
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410000,410000);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,410015,410016);
        //Zgamma
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,301535,301536);
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,301899,301907);
        //Z+jets
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,364100,364141);
        //Drell Yan
        fillBGMCIndex(BGMCSampleInfo,element.BGMCIndex,364198,364215);
        element.colour = 2; BGMCGroupData.push_back(element);
        */
    }
    
    //Group for data-driven background
    std::vector<GroupData> BGDataGroupData;
    {
        GroupData element;
        element.GroupName = "charge flip"; element.LegendName = element.GroupName; element.LatexName = element.GroupName;
        element.colour = 2; BGDataGroupData.push_back(element);
        
        element.GroupName = "fake lepton"; element.LegendName = element.GroupName; element.LatexName = element.GroupName;
        element.colour = 3; BGDataGroupData.push_back(element);
    }
    
    //plotting
    //Variables for plotting
    std::vector<VarData> Var;
    {
        VarData element;
        
        element.VarName = "pt1";           element.VarTitle = "p_{T}^{l1}";                     element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        //element.bin=5;          element.xmin=25;                element.xmax=300;
        element.bin=20;         element.xmin=25;                element.xmax=250;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$p_T^{l1}$";
        element.CutDirection=1;
        Var.push_back(element);
        
        element.VarName = "pt2";           element.VarTitle = "p_{T}^{l2}";                     element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=25;                element.xmax=250;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$p_T^{l2}$";
        element.CutDirection=1;
        Var.push_back(element);
        
        element.VarName = "eta1";          element.VarTitle = "#eta^{l1}";                      element.unit = "";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=-3;                element.xmax=3;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$\\eta^{l1}$";
        element.CutDirection=0;
        Var.push_back(element);
        
        element.VarName = "eta2";          element.VarTitle = "#eta^{l2}";                      element.unit = "";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=-3;                element.xmax=3;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$\\eta^{l2}$";
        element.CutDirection=0;
        Var.push_back(element);
        
        element.VarName = "phi1";          element.VarTitle = "#phi^{l1}";                      element.unit = "";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=-TMath::Pi();      element.xmax=TMath::Pi();
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$\\phi^{l1}$";
        element.CutDirection=0;
        Var.push_back(element);
        
        element.VarName = "mll";           element.VarTitle = "m_{ll}";                         element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=0;                 element.xmax=250;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$m_{ll}$";
        element.CutDirection=0;
        Var.push_back(element);
        
        element.VarName = "ptll";          element.VarTitle = "p_{T}^{ll}";                     element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=0;                 element.xmax=400;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$p_T^{ll}$";
        element.CutDirection=1;
        Var.push_back(element);
        
        element.VarName = "MET";           element.VarTitle = "E_{T}^{miss}";                   element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        //element.bin=5;          element.xmin=0;                 element.xmax=500;
        element.bin=20;         element.xmin=0;                 element.xmax=200;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$E_{\\text{T}}^{\\text{miss}}$";
        element.CutDirection=1;
        Var.push_back(element);
        
        element.VarName = "mTtwo";         element.VarTitle = "m_{T2}";                         element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=0;                 element.xmax=150;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$m_{T2}$";
        element.CutDirection=1;
        Var.push_back(element);
        
        element.VarName = "mt1";           element.VarTitle = "m_{T}^{l1}";                     element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        //element.bin=5;          element.xmin=0;                 element.xmax=300;
        element.bin=20;         element.xmin=0;                 element.xmax=250;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$m_{\\text{T}}^{l1}$";
        element.CutDirection=1;
        Var.push_back(element);
        
        element.VarName = "mt2";           element.VarTitle = "m_{T}^{l2}";                     element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=0;                 element.xmax=250;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$m_{\\text{T}}^{l2}$";
        element.CutDirection=1;
        Var.push_back(element);
        
        element.VarName = "averageMu";     element.VarTitle = "averageMu";                      element.unit = "";
        element.VarFormula = element.VarName;
        element.bin=35;         element.xmin=0;                 element.xmax=35;
        element.log=0;          element.ymin=0;                 element.ymax=0;
        element.latexName = "Average number of interactions per bunch crossing";
        element.CutDirection=0;
        Var.push_back(element);
        
        element.VarName = "jetpt";         element.VarTitle = "p_{T} of the leading jet";       element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=0;                 element.xmax=300;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$p_T$ of the leading jet";
        element.CutDirection=1;
        Var.push_back(element);
        
        element.VarName = "jeteta";        element.VarTitle = "#eta of the leading jet";        element.unit = "";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=-3;                element.xmax=3;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$\\eta$ of the leading jet";
        element.CutDirection=0;
        Var.push_back(element);
        
        element.VarName = "nJet";          element.VarTitle = "Number of jets";                 element.unit = "";
        element.VarFormula = element.VarName;
        //element.bin=6;          element.xmin=0;                 element.xmax=6;
        element.bin=15;         element.xmin=0;                 element.xmax=15;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = element.VarTitle;
        element.CutDirection=0;
        Var.push_back(element);
        
        element.VarName = "nBJet";         element.VarTitle = "Number of b-jets";               element.unit = "";
        element.VarFormula = element.VarName;
        element.bin=6;          element.xmin=0;                 element.xmax=6;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = element.VarTitle;
        element.CutDirection=-1;
        Var.push_back(element);
        
        element.VarName = "nCJet";         element.VarTitle = "Number of central jets";         element.unit = "";
        element.VarFormula = element.VarName;
        element.bin=6;          element.xmin=0;                 element.xmax=6;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = element.VarTitle;
        element.CutDirection=0;
        Var.push_back(element);
        
        element.VarName = "l12_dPhi";      element.VarTitle = "|#Delta#phi_{ll}|";              element.unit = "";
        element.VarFormula = "fabs(l12_dPhi)";
        element.bin=20;         element.xmin=0;                 element.xmax=TMath::Pi();
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$|\\Delta\\phi_{ll}|$";
        element.CutDirection=0;
        Var.push_back(element);
        
        element.VarName = "l12_MET_dPhi";  element.VarTitle = "|#Delta#phi_{ll,MET}|";          element.unit = "";
        element.VarFormula = "fabs(l12_MET_dPhi)";
        element.bin=20;         element.xmin=0;                 element.xmax=TMath::Pi();
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$|\\Delta\\phi_{ll,\\text{MET}}|$";
        element.CutDirection=0;
        Var.push_back(element);
        
        element.VarName = "jet0_MET_dPhi"; element.VarTitle = "|#Delta#phi_{jet0,MET}|";        element.unit = "";
        element.VarFormula = "fabs(jet0_MET_dPhi)";
        element.bin=20;         element.xmin=0;                 element.xmax=TMath::Pi();
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$|\\Delta\\phi_{\\text{jet0,MET}}|$";
        element.CutDirection=0;
        Var.push_back(element);
        
        element.VarName = "dEta";          element.VarTitle = "|#Delta#eta_{ll}|";              element.unit = "";
        element.VarFormula = "fabs(eta1-eta2)";
        element.bin=20;         element.xmin=0;                 element.xmax=5;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$|\\Delta\\eta_{ll}|$";
        element.CutDirection=-1;
        Var.push_back(element);
        
        element.VarName = "METRel";        element.VarTitle = "E_{T}^{miss,rel}";               element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=0;                 element.xmax=200;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$E_{\\text{T}}^{\\text{miss,rel}}$";
        element.CutDirection=1;
        Var.push_back(element);
        
        element.VarName = "meff";          element.VarTitle = "m_{eff}";                        element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        //element.bin=5;          element.xmin=200;               element.xmax=900;
        element.bin=20;         element.xmin=0;                 element.xmax=600;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$m_{\\text{eff}}$";
        element.CutDirection=1;
        Var.push_back(element);
        
        element.VarName = "mtm";           element.VarTitle = "m_{T}^{max}";                    element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=0;                 element.xmax=250;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$m_{\\text{T}}^{\\text{max}}$";
        element.CutDirection=1;
        Var.push_back(element);
        
        element.VarName = "mlj";           element.VarTitle = "m_{lj}/m_{ljj}";                 element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=0;                 element.xmax=300;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$m_{lj}$/$m_{ljj}$";
        element.CutDirection=-1;
        Var.push_back(element);
        
        element.VarName = "mjj";           element.VarTitle = "m_{jj}";                         element.unit = "[GeV]";
        element.VarFormula = element.VarName;
        element.bin=20;         element.xmin=0;                 element.xmax=300;
        element.log=1;          element.ymin=1.0/element.bin;   element.ymax=1;
        element.latexName = "$m_{jj}$";
        element.CutDirection=-1;
        Var.push_back(element);
    }
    SetAtlasStyle();
    
    struct AdditionalCutData
    {
        TString Cut;
        TString RelatedVariable;
    };
    
    struct OptimizingCutData
    {
        double lower;
        double upper;
        double min;
        double max;
        int nBin;
        TString RelatedVariable;
        bool OnlyHasLowerCut;
        bool OnlyHasUpperCut;
    };
    
    struct RegionGroupData
    {
        TString GroupName;
        unsigned int lower;
        unsigned int upper;
        bool showData;
        bool showSignificance;
    };
    
    struct RegionData
    {
        TString RegionName;
        TString RegionShortName;
        std::vector<unsigned int> setOfChannel;
        TString Cut;
        vector<AdditionalCutData> AdditionalCut;
        vector< vector< vector<OptimizingCutData> > > OptimizingCut;
    };
    
    if(dorw)
    {
        //Z pt reweighting
        int VarIndex = 0;
        for(unsigned int i=0;i<Var.size();i++)
        {
            if(Var[i].VarName == "ptll") VarIndex = i;
        }
        
        std::vector<RegionData> RegionInfo;
        {
            RegionData element;
            
            /*
            element.RegionName = "OS_ee";
            element.setOfChannel.clear();
            element.setOfChannel.push_back(0);
            element.setOfChannel.push_back(6);
            RegionInfo.push_back(element);
            
            element.RegionName = "OS_mumu";
            element.setOfChannel.clear();
            element.setOfChannel.push_back(1);
            element.setOfChannel.push_back(7);
            RegionInfo.push_back(element);
            */
            
            ///*
            element.RegionName = ChannelInfo[0].ChannelName;
            element.setOfChannel.clear();
            element.setOfChannel.push_back(0);
            RegionInfo.push_back(element);
            
            element.RegionName = ChannelInfo[1].ChannelName;
            element.setOfChannel.clear();
            element.setOfChannel.push_back(1);
            RegionInfo.push_back(element);
            
            element.RegionName = ChannelInfo[6].ChannelName;
            element.setOfChannel.clear();
            element.setOfChannel.push_back(6);
            RegionInfo.push_back(element);
            
            element.RegionName = ChannelInfo[7].ChannelName;
            element.setOfChannel.clear();
            element.setOfChannel.push_back(7);
            RegionInfo.push_back(element);
            //*/
        }

        //calculate ratio plot
        TH1D* h2Ratio_rw[RegionInfo.size()];
        for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
        {
            const unsigned int channelRepresentative = RegionInfo[RegionIndex].setOfChannel[0];
            
            std::vector<TChain*> tree2Data;
            initializeTree2(tree2Data,RegionInfo[RegionIndex].setOfChannel,DataSampleName,ChannelInfo);
            
            std::vector< std::vector<TChain*> > tree2BGMC;
            std::vector< std::vector<double> > BGMCGroupXS;
            std::vector< std::vector<double> > BGMCGroupnwAOD;
            std::vector<Group> BGGroup;
            {
                //For MC background
                for(unsigned int k=0;k<ChannelInfo[channelRepresentative].setOfBGMC.size();k++)
                {
                    for(unsigned int j=0;j<BGMCGroupData.size();j++)
                    {
                        if(BGMCGroupData[j].GroupName == ChannelInfo[channelRepresentative].setOfBGMC[k])
                        {
                            std::vector<TString> BGMCGroupSampleName;
                            std::vector<double> BGMCGroupXSElement;
                            for(unsigned int m=0;m<BGMCGroupData[j].BGMCIndex.size();m++)
                            {
                                BGMCGroupSampleName.push_back(BGMCSampleInfo[BGMCGroupData[j].BGMCIndex[m]].SampleName);
                                BGMCGroupXSElement.push_back(BGMCSampleInfo[BGMCGroupData[j].BGMCIndex[m]].XS);
                            }
                            
                            std::vector<TChain*> tree2BGMCElement;
                            initializeTree2(tree2BGMCElement,RegionInfo[RegionIndex].setOfChannel,BGMCGroupSampleName,ChannelInfo);
                            
                            std::vector<double> BGMCGroupnwAODElement;
                            getnwAOD(BGMCGroupnwAODElement,RegionInfo[RegionIndex].setOfChannel,BGMCGroupSampleName,ChannelInfo);
                            
                            tree2BGMC.push_back(tree2BGMCElement);
                            BGMCGroupXS.push_back(BGMCGroupXSElement);
                            BGMCGroupnwAOD.push_back(BGMCGroupnwAODElement);
                            
                            Group BGGroupElement;
                            BGGroupElement.info = &(BGMCGroupData[j]);
                            BGGroup.push_back(BGGroupElement);
                            
                        }
                    }
                }
            }

            TString title = Var[VarIndex].VarTitle;
            
            TString xaxis;
            xaxis += Var[VarIndex].VarTitle;
            xaxis += " ";
            xaxis += Var[VarIndex].unit;
            
            //h2DataSum
            TH1D* h2DataSum = new TH1D("DataSum",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
            
            //h2Data
            for(unsigned int j=0;j<DataSampleName.size();j++)
            {
                TH1D* hTemp = new TH1D("Data",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                
                //fill histograms from trees
                TString temp;
                temp += Var[VarIndex].VarFormula;
                temp += ">>Data";
                tree2Data[j]->Draw(temp.Data(),"fLwt==0");
                
                //add data
                h2DataSum->Add(hTemp);
                delete hTemp;
            }
            
            //For MC background
            for(unsigned int j=0;j<tree2BGMC.size();j++)
            {
                BGGroup[j].h2 = new TH1D(BGGroup[j].info->GroupName.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                for(unsigned int k=0;k<tree2BGMC[j].size();k++)
                {
                    TH1D* hTemp = new TH1D("BGMC",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                    
                    //fill histograms from trees
                    TString temp;
                    temp += Var[VarIndex].VarFormula;
                    temp += ">>BGMC";
                    tree2BGMC[j][k]->Draw(temp.Data(),"weight");
                    
                    //normalization for BG
                    hTemp->Scale(BGMCGroupXS[j][k]/BGMCGroupnwAOD[j][k] *sumDataL);
                    
                    //add BG
                    BGGroup[j].h2->Add(hTemp);
                    delete hTemp;
                }
            }
            
            //ratio plot for reweighting
            {
                TString NameTemp = "Ratio_";
                NameTemp += RegionInfo[RegionIndex].RegionName;
                h2Ratio_rw[RegionIndex] = new TH1D(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                h2Ratio_rw[RegionIndex]->GetXaxis()->SetTitle(xaxis.Data());
                h2Ratio_rw[RegionIndex]->GetYaxis()->SetTitle("(Data - NonZjets)/Zjets");
                h2Ratio_rw[RegionIndex]->SetMarkerColor(RegionIndex+1);
                h2Ratio_rw[RegionIndex]->SetMarkerSize(1);
                h2Ratio_rw[RegionIndex]->SetLineColor(RegionIndex+1);
                h2Ratio_rw[RegionIndex]->Sumw2();
            }
            h2Ratio_rw[RegionIndex]->Add(h2DataSum);
            
            //Z+jets
            TH1D* h2Zjets = new TH1D("Zjets",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
            TH1D* h2nonZjets = new TH1D("nonZjets",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);

            //Non Z+jets
            for(unsigned int j=0;j<tree2BGMC.size();j++)
            {
                if(BGGroup[j].info->GroupName == "Zee"     ||
                   BGGroup[j].info->GroupName == "Zmumu"   ||
                   BGGroup[j].info->GroupName == "Ztautau" )
                {
                    h2Zjets->Add(BGGroup[j].h2);
                }
                else
                {
                    h2nonZjets->Add(BGGroup[j].h2);
                }
            }
            h2nonZjets->Scale(-1);
            h2Ratio_rw[RegionIndex]->Add(h2nonZjets);
            h2Ratio_rw[RegionIndex]->Divide(h2Zjets);
            
            //delete
            //h2DataSum
            delete h2DataSum;
            
            //h2BGGruop
            for(unsigned int j=0;j<BGGroup.size();j++)
            {
                delete BGGroup[j].h2;
            }
            
            //Z+jets
            delete h2Zjets;
            delete h2nonZjets;
            
            //delete tree2
            for(unsigned int i=0;i<DataSampleName.size();i++)
            {
                delete tree2Data[i];
            }
            
            for(unsigned int j=0;j<tree2BGMC.size();j++)
            {
                for(unsigned int k=0;k<tree2BGMC[j].size();k++)
                {
                    delete tree2BGMC[j][k];
                }
            }
        }
        
        TF1* fun[RegionInfo.size()];
        //simple fit
        if(fitting && simple)
        {
            for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
            {
                TString NameTemp = "fun_";
                NameTemp += TString::Itoa(RegionIndex,10);
                fun[RegionIndex] = new TF1(NameTemp.Data(),"pol2",Var[VarIndex].xmin,Var[VarIndex].xmax);
                fun[RegionIndex]->SetLineColor(RegionIndex+1);
                fun[RegionIndex]->SetLineStyle(2);
                h2Ratio_rw[RegionIndex]->Fit(NameTemp.Data(),"R");
            }
        }
        
        //combined fit
        if(fitting && combined)
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
            
            for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
            {
                const unsigned int channelRepresentative = RegionInfo[RegionIndex].setOfChannel[0];
                TString NameTemp = "plot/";
                leg->AddEntry(h2Ratio_rw[RegionIndex],ChannelInfo[channelRepresentative].ChannelName.Data(),"p");
            }
        }
        
        TCanvas* c2 = new TCanvas();
        c2->cd();
        gStyle->SetOptStat(0);
        h2Ratio_rw[0]->SetMinimum(0);
        h2Ratio_rw[0]->SetMaximum(2);
        h2Ratio_rw[0]->Draw();
        for(unsigned int RegionIndex=1;RegionIndex<RegionInfo.size();RegionIndex++)
        {
            h2Ratio_rw[RegionIndex]->Draw("same");
        }
        leg->Draw();
        
        {
            //export histograms in eps format
            TString NameTemp = "plot/";
            NameTemp += Var[VarIndex].VarName;
            NameTemp += "_rw";
            NameTemp += ".eps";
            c2->Print(NameTemp,"eps");
        }
        
        //add weight in the tree
        for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
        {
            for(unsigned int i=0;i<RegionInfo[RegionIndex].setOfChannel.size();i++)
            {
                TString FileName = "skimming/rw_";
                const int ChannelIndex = RegionInfo[RegionIndex].setOfChannel[i];
                FileName += TString::Itoa(ChannelIndex,10);
                FileName += ".root";
                TFile* frw = new TFile(FileName.Data(),"RECREATE");
                
                for(unsigned int j=0;j<BGMCGroupData.size();j++)
                {
                    if(BGMCGroupData[j].GroupName == "Zee"     ||
                       BGMCGroupData[j].GroupName == "Zmumu"   ||
                       BGMCGroupData[j].GroupName == "Ztautau" )
                    {
                        for(unsigned int k=0;k<BGMCGroupData[j].BGMCIndex.size();k++)
                        {
                            TChain* tree1 = nullptr;
                            initializeTree1new(tree1,BGMCSampleInfo[BGMCGroupData[j].BGMCIndex[k]].SampleName,ChannelInfo[ChannelIndex].ChannelName);
                            
                            TString NameTemp = "tree_rw_";
                            NameTemp += TString::Itoa(BGMCGroupData[j].BGMCIndex[k],10);
                            TTree* tree2 = new TTree(NameTemp.Data(),NameTemp.Data());
                            tree2->Branch("rw",&rw,"rw/D");
                            
                            for(int m=0;m<tree1->GetEntries();m++)
                            {
                                tree1->GetEntry(m);
                                if(fitting)
                                {
                                    if(simple) rw=fun[RegionIndex]->Eval(ptll);
                                    else if(combined) rw=fun[1]->Eval(ptll);
                                }
                                if(direct)
                                {
                                    rw = h2Ratio_rw[RegionIndex]->GetBinContent(h2Ratio_rw[RegionIndex]->FindBin(ptll));
                                }
                                tree2->Fill();
                            }
                            delete tree1;
                            
                            frw->cd();
                            tree2->Write();
                        }
                    }
                }
                delete frw;
            }
        }
        
        //delete h2Ratio_rw
        for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
        {
            delete h2Ratio_rw[RegionIndex];
        }
        //delete legend
        delete leg;
        //delete canvas
        delete c2;
        
        if(fitting)
        {
            //delete functions
            for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
            {
                delete fun[RegionIndex];
            }
        }
    }
    
    if(docfw)
    //if(false)
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
            
            for(unsigned int j=0;j<BGMCGroupData.size();j++)
            {
                if(BGMCGroupData[j].GroupName == "Zee")
                {
                    for(unsigned int k=0;k<BGMCGroupData[j].BGMCIndex.size();k++)
                    {
                        TChain* tree1 = nullptr;
                        initializeTree1new(tree1,BGMCSampleInfo[BGMCGroupData[j].BGMCIndex[k]].SampleName,ChannelInfo[ChannelIndex].ChannelName);
                        
                        TString NameTemp = "tree_cfw_";
                        NameTemp += TString::Itoa(BGMCGroupData[j].BGMCIndex[k],10);
                        TTree* tree2 = new TTree(NameTemp.Data(),NameTemp.Data());
                        tree2->Branch("cfw",&cfw,"cfw/D");
                        
                        for(int m=0;m<tree1->GetEntries();m++)
                        {
                            tree1->GetEntry(m);
                            
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
                            tree2->Fill();
                        }
                        delete tree1;
                        
                        fcfw->cd();
                        tree2->Write();
                    }
                }
            }
            delete fcfw;
        }
        delete f_cf_data;
        delete f_cf_mc;
    }
    
    std::vector<RegionGroupData> RegionGroup;
    std::vector<RegionData> RegionInfo;
    {
        RegionGroupData GroupElement;
        RegionData element;
        AdditionalCutData AdditionalCutElement;
        vector< vector<OptimizingCutData> > OptimizingCutElement1;
        vector<OptimizingCutData> OptimizingCutElement2;
        OptimizingCutData OptimizingCutElement3;
        
        //CR_OS
        GroupElement.GroupName = "CR_OS";
        GroupElement.lower = RegionInfo.size();
        GroupElement.showData = true;
        GroupElement.showSignificance = false;
        element.Cut = "";
        element.AdditionalCut.clear();
        
        for(unsigned int ChannelIndex=0;ChannelIndex<ChannelInfo.size();ChannelIndex++)
        {
            if((ChannelIndex>=3 && ChannelIndex<=5) || (ChannelIndex>=9 && ChannelIndex<=11)) continue;
            
            element.RegionName = "CR_";
            element.RegionName += ChannelInfo[ChannelIndex].ChannelName;
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(ChannelIndex);
            RegionInfo.push_back(element);
        }
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
        
        //SR_SS
        GroupElement.GroupName = "SR_SS";
        GroupElement.lower = RegionInfo.size();
        GroupElement.showData = false;
        GroupElement.showSignificance = true;
        element.Cut = "";
        element.AdditionalCut.clear();
        
        for(unsigned int ChannelIndex=0;ChannelIndex<ChannelInfo.size();ChannelIndex++)
        {
            if(ChannelIndex<=2 || (ChannelIndex>=6 && ChannelIndex<=8)) continue;
            
            element.RegionName = "SR_";
            element.RegionName += ChannelInfo[ChannelIndex].ChannelName;
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(ChannelIndex);
            RegionInfo.push_back(element);
        }
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
        
        //CR_OS_1B
        GroupElement.GroupName = "CR_OS_1B";
        GroupElement.lower = RegionInfo.size();
        GroupElement.showData = true;
        GroupElement.showSignificance = false;
        element.Cut = " && nBJet == 1";
        element.AdditionalCut.clear();
        
        for(unsigned int ChannelIndex=0;ChannelIndex<ChannelInfo.size();ChannelIndex++)
        {
            if((ChannelIndex>=3 && ChannelIndex<=5) || (ChannelIndex>=9 && ChannelIndex<=11)) continue;
            
            element.RegionName = "CR_";
            element.RegionName += ChannelInfo[ChannelIndex].ChannelName;
            element.RegionName += "_1B";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(ChannelIndex);
            RegionInfo.push_back(element);
        }
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
        
        //CR_SS_1B
        GroupElement.GroupName = "CR_SS_1B";
        GroupElement.lower = RegionInfo.size();
        GroupElement.showData = true;
        GroupElement.showSignificance = false;
        element.Cut = " && nBJet == 1";
        element.AdditionalCut.clear();
        
        for(unsigned int ChannelIndex=0;ChannelIndex<ChannelInfo.size();ChannelIndex++)
        {
            if(ChannelIndex<=2 || (ChannelIndex>=6 && ChannelIndex<=8)) continue;
            
            element.RegionName = "CR_";
            element.RegionName += ChannelInfo[ChannelIndex].ChannelName;
            element.RegionName += "_1B";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(ChannelIndex);
            RegionInfo.push_back(element);
        }
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
        
        //CR_OS_2B
        GroupElement.GroupName = "CR_OS_2B";
        GroupElement.lower = RegionInfo.size();
        GroupElement.showData = true;
        GroupElement.showSignificance = false;
        element.Cut = " && nBJet == 2";
        element.AdditionalCut.clear();
        
        for(unsigned int ChannelIndex=0;ChannelIndex<ChannelInfo.size();ChannelIndex++)
        {
            if((ChannelIndex>=3 && ChannelIndex<=5) || (ChannelIndex>=9 && ChannelIndex<=11)) continue;
            
            element.RegionName = "CR_";
            element.RegionName += ChannelInfo[ChannelIndex].ChannelName;
            element.RegionName += "_2B";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(ChannelIndex);
            RegionInfo.push_back(element);
        }
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
        
        //CR_SS_2B
        GroupElement.GroupName = "CR_SS_2B";
        GroupElement.lower = RegionInfo.size();
        GroupElement.showData = true;
        GroupElement.showSignificance = false;
        element.Cut = " && nBJet == 2";
        element.AdditionalCut.clear();
        
        for(unsigned int ChannelIndex=0;ChannelIndex<ChannelInfo.size();ChannelIndex++)
        {
            if(ChannelIndex<=2 || (ChannelIndex>=6 && ChannelIndex<=8)) continue;
            
            element.RegionName = "CR_";
            element.RegionName += ChannelInfo[ChannelIndex].ChannelName;
            element.RegionName += "_2B";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(ChannelIndex);
            RegionInfo.push_back(element);
        }
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
        
        //CR_SS_mumu_low_mT2
        GroupElement.GroupName = "CR_SS_mumu_low_mT2";
        GroupElement.lower = RegionInfo.size();
        
        GroupElement.showData = true;
        GroupElement.showSignificance = false;
        element.Cut = " && mTtwo<30";
        element.AdditionalCut.clear();
        
        element.RegionName = "CR_nonISR_SS_mumu_low_mT2";
        element.setOfChannel.clear();
        element.setOfChannel.push_back(4);
        RegionInfo.push_back(element);
        
        element.RegionName = "CR_ISR_SS_mumu_low_mT2";
        element.setOfChannel.clear();
        element.setOfChannel.push_back(10);
        RegionInfo.push_back(element);
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
        
        //CR_SS_ee_Zmass
        GroupElement.GroupName = "CR_SS_ee_Zmass";
        GroupElement.lower = RegionInfo.size();
        
        GroupElement.showData = true;
        GroupElement.showSignificance = false;
        element.Cut = "";
        element.AdditionalCut.clear();
        
        AdditionalCutElement.Cut = " && fabs(mll - 91.2) < 10";
        AdditionalCutElement.RelatedVariable = "mll";
        element.AdditionalCut.push_back(AdditionalCutElement);
        
        element.RegionName = "CR_SS_ee_Zmass";
        element.setOfChannel.clear();
        element.setOfChannel.push_back(3);
        element.setOfChannel.push_back(9);
        RegionInfo.push_back(element);
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
        
        //SR_SS_0B
        GroupElement.GroupName = "SR_SS_0B";
        GroupElement.lower = RegionInfo.size();
        GroupElement.showData = false;
        GroupElement.showSignificance = true;
        element.Cut = " && nBJet == 0";
        element.AdditionalCut.clear();
        
        for(unsigned int ChannelIndex=0;ChannelIndex<ChannelInfo.size();ChannelIndex++)
        {
            if(ChannelIndex<=2 || (ChannelIndex>=6 && ChannelIndex<=8)) continue;
            
            element.RegionName = "SR_";
            element.RegionName += ChannelInfo[ChannelIndex].ChannelName;
            element.RegionName += "_0B";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(ChannelIndex);
            RegionInfo.push_back(element);
        }
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
        
        //Signal region for run 1
        GroupElement.GroupName = "SR_SS_run1";
        GroupElement.lower = RegionInfo.size();
        GroupElement.showData = false;
        GroupElement.showSignificance = true;
        
        //ee_1
        {
            element.RegionName = "SR_SS_ee_1_run1";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(3);
            element.setOfChannel.push_back(9);
            
            element.Cut = " && nCJet == 1";
            element.Cut += " && nBJet == 0";
            
            element.AdditionalCut.clear();
            
            AdditionalCutElement.Cut = " && pt1 > 30";
            AdditionalCutElement.RelatedVariable = "pt1";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && pt2 > 20";
            AdditionalCutElement.RelatedVariable = "pt2";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && fabs(mll - 91.2) > 10";
            AdditionalCutElement.RelatedVariable = "mll";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && METRel > 55";
            AdditionalCutElement.RelatedVariable = "METRel";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && meff > 200";
            AdditionalCutElement.RelatedVariable = "meff";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && mlj < 90";
            AdditionalCutElement.RelatedVariable = "mlj";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            RegionInfo.push_back(element);
        }
        
        //mumu_1
        {
            element.RegionName = "SR_SS_mumu_1_run1";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(4);
            element.setOfChannel.push_back(10);
            
            element.Cut = " && nCJet == 1";
            element.Cut += " && nBJet == 0";
            
            element.AdditionalCut.clear();
            
            AdditionalCutElement.Cut = " && pt1 > 30";
            AdditionalCutElement.RelatedVariable = "pt1";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && pt2 > 20";
            AdditionalCutElement.RelatedVariable = "pt2";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && fabs(eta1 - eta2) < 1.5";
            AdditionalCutElement.RelatedVariable = "dEta";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && meff > 200";
            AdditionalCutElement.RelatedVariable = "meff";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && mtm > 110";
            AdditionalCutElement.RelatedVariable = "mtm";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && mlj < 90";
            AdditionalCutElement.RelatedVariable = "mlj";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            RegionInfo.push_back(element);
        }
        
        //emu_1
        {
            element.RegionName = "SR_SS_emu_1_run1";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(5);
            element.setOfChannel.push_back(11);
            
            element.Cut = " && nCJet == 1";
            element.Cut += " && nBJet == 0";
            
            element.AdditionalCut.clear();
            
            AdditionalCutElement.Cut = " && pt1 > 30";
            AdditionalCutElement.RelatedVariable = "pt1";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && pt2 > 30";
            AdditionalCutElement.RelatedVariable = "pt2";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && fabs(eta1 - eta2) < 1.5";
            AdditionalCutElement.RelatedVariable = "dEta";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && meff > 200";
            AdditionalCutElement.RelatedVariable = "meff";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && mtm > 110";
            AdditionalCutElement.RelatedVariable = "mtm";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && mlj < 90";
            AdditionalCutElement.RelatedVariable = "mlj";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            RegionInfo.push_back(element);
        }
        
        //ee_2
        {
            element.RegionName = "SR_SS_ee_2_run1";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(3);
            element.setOfChannel.push_back(9);
            
            element.Cut = " && (nCJet == 2 || nCJet == 3)";
            element.Cut += " && nBJet == 0";
            
            element.AdditionalCut.clear();
            
            AdditionalCutElement.Cut = " && pt1 > 30";
            AdditionalCutElement.RelatedVariable = "pt1";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && pt2 > 20";
            AdditionalCutElement.RelatedVariable = "pt2";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && fabs(mll - 91.2) > 10";
            AdditionalCutElement.RelatedVariable = "mll";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && METRel > 30";
            AdditionalCutElement.RelatedVariable = "METRel";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && mtm > 110";
            AdditionalCutElement.RelatedVariable = "mtm";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && mlj < 120";
            AdditionalCutElement.RelatedVariable = "mlj";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            RegionInfo.push_back(element);
        }
        
        //mumu_2
        {
            element.RegionName = "SR_SS_mumu_2_run1";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(4);
            element.setOfChannel.push_back(10);
            
            element.Cut = " && (nCJet == 2 || nCJet == 3)";
            element.Cut += " && nBJet == 0";
            
            element.AdditionalCut.clear();
            
            AdditionalCutElement.Cut = " && pt1 > 30";
            AdditionalCutElement.RelatedVariable = "pt1";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && pt2 > 30";
            AdditionalCutElement.RelatedVariable = "pt2";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && fabs(eta1 - eta2) < 1.5";
            AdditionalCutElement.RelatedVariable = "dEta";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && meff > 200";
            AdditionalCutElement.RelatedVariable = "meff";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && mlj < 120";
            AdditionalCutElement.RelatedVariable = "mlj";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            RegionInfo.push_back(element);
        }
        
        //emu_2
        {
            element.RegionName = "SR_SS_emu_2_run1";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(5);
            element.setOfChannel.push_back(11);
            
            element.Cut = " && (nCJet == 2 || nCJet == 3)";
            element.Cut += " && nBJet == 0";
            
            element.AdditionalCut.clear();
            
            AdditionalCutElement.Cut = " && pt1 > 30";
            AdditionalCutElement.RelatedVariable = "pt1";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && pt2 > 30";
            AdditionalCutElement.RelatedVariable = "pt2";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && fabs(eta1 - eta2) < 1.5";
            AdditionalCutElement.RelatedVariable = "dEta";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && meff > 200";
            AdditionalCutElement.RelatedVariable = "meff";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && mtm > 110";
            AdditionalCutElement.RelatedVariable = "mtm";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            AdditionalCutElement.Cut = " && mlj < 120";
            AdditionalCutElement.RelatedVariable = "mlj";
            element.AdditionalCut.push_back(AdditionalCutElement);
            
            RegionInfo.push_back(element);
        }
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
        
        
        //Signal region for optimization
        GroupElement.GroupName = "SR_SS_opt";
        GroupElement.lower = RegionInfo.size();
        GroupElement.showData = false;
        GroupElement.showSignificance = true;
        
        
        {
            element.OptimizingCut.clear();
            OptimizingCutElement1.clear();
            
            OptimizingCutElement2.clear();
            OptimizingCutElement3.RelatedVariable = "pt1";
            OptimizingCutElement3.min = 25;
            OptimizingCutElement3.max = 250;
            OptimizingCutElement3.nBin = 45;
            OptimizingCutElement3.lower = 25;
            OptimizingCutElement3.upper = -1;
            OptimizingCutElement2.push_back(OptimizingCutElement3);
            OptimizingCutElement1.push_back(OptimizingCutElement2);
            
            OptimizingCutElement2.clear();
            OptimizingCutElement3.RelatedVariable = "pt2";
            OptimizingCutElement3.min = 25;
            OptimizingCutElement3.max = 250;
            OptimizingCutElement3.nBin = 45;
            OptimizingCutElement3.lower = 25;
            OptimizingCutElement3.upper = -1;
            OptimizingCutElement2.push_back(OptimizingCutElement3);
            OptimizingCutElement1.push_back(OptimizingCutElement2);
            
            OptimizingCutElement2.clear();
            OptimizingCutElement3.RelatedVariable = "dEta";
            OptimizingCutElement3.min = 0;
            OptimizingCutElement3.max = 5;
            OptimizingCutElement3.nBin = 10;
            OptimizingCutElement3.lower = 0;
            OptimizingCutElement3.upper = -1;
            //OptimizingCutElement3.upper = 1.5;
            OptimizingCutElement2.push_back(OptimizingCutElement3);
            OptimizingCutElement1.push_back(OptimizingCutElement2);
            
            OptimizingCutElement2.clear();
            //OptimizingCutElement3.RelatedVariable = "METRel";
            OptimizingCutElement3.RelatedVariable = "MET";
            OptimizingCutElement3.min = 0;
            OptimizingCutElement3.max = 200;
            OptimizingCutElement3.nBin = 40;
            OptimizingCutElement3.lower = 0;
            OptimizingCutElement3.upper = -1;
            OptimizingCutElement2.push_back(OptimizingCutElement3);
            OptimizingCutElement1.push_back(OptimizingCutElement2);
            
            OptimizingCutElement2.clear();
            OptimizingCutElement3.RelatedVariable = "meff";
            OptimizingCutElement3.min = 0;
            OptimizingCutElement3.max = 600;
            OptimizingCutElement3.nBin = 60;
            OptimizingCutElement3.lower = 0;
            //OptimizingCutElement3.lower = 200;
            OptimizingCutElement3.upper = -1;
            OptimizingCutElement2.push_back(OptimizingCutElement3);
            OptimizingCutElement1.push_back(OptimizingCutElement2);
            
            OptimizingCutElement2.clear();
            //OptimizingCutElement3.RelatedVariable = "mtm";
            OptimizingCutElement3.RelatedVariable = "mt1"; //Dani
            OptimizingCutElement3.min = 0;
            OptimizingCutElement3.max = 300;
            OptimizingCutElement3.nBin = 60;
            OptimizingCutElement3.lower = 0;
            OptimizingCutElement3.upper = -1;
            OptimizingCutElement2.push_back(OptimizingCutElement3);
            OptimizingCutElement1.push_back(OptimizingCutElement2);
            
            OptimizingCutElement2.clear();
            OptimizingCutElement3.RelatedVariable = "mlj";
            OptimizingCutElement3.min = 0;
            OptimizingCutElement3.max = 300;
            OptimizingCutElement3.nBin = 60;
            OptimizingCutElement3.lower = 0;
            OptimizingCutElement3.upper = -1;
            OptimizingCutElement2.push_back(OptimizingCutElement3);
            OptimizingCutElement1.push_back(OptimizingCutElement2);
            
            OptimizingCutElement2.clear();
            OptimizingCutElement3.RelatedVariable = "mTtwo";
            OptimizingCutElement3.min = 0;
            OptimizingCutElement3.max = 150;
            OptimizingCutElement3.nBin = 30;
            OptimizingCutElement3.lower = 0;
            OptimizingCutElement3.upper = -1;
            OptimizingCutElement2.push_back(OptimizingCutElement3);
            OptimizingCutElement1.push_back(OptimizingCutElement2);
            
            /*
            OptimizingCutElement2.clear();
            OptimizingCutElement3.RelatedVariable = "ptll";
            OptimizingCutElement3.min = 0;
            OptimizingCutElement3.max = 300;
            OptimizingCutElement3.nBin = 30;
            OptimizingCutElement3.lower = 0;
            OptimizingCutElement3.upper = -1;
            OptimizingCutElement2.push_back(OptimizingCutElement3);
            OptimizingCutElement1.push_back(OptimizingCutElement2);
            */
            
            for(unsigned int i=0;i<OptimizingCutElement1.size();i++)
            {
                for(unsigned int j=0;j<OptimizingCutElement1[i].size();j++)
                {
                    OptimizingCutElement1[i][j].OnlyHasLowerCut =
                    (OptimizingCutElement1[i][j].RelatedVariable == "pt1" ||
                     OptimizingCutElement1[i][j].RelatedVariable == "pt2" ||
                     OptimizingCutElement1[i][j].RelatedVariable == "ptll" ||
                     OptimizingCutElement1[i][j].RelatedVariable == "mTtwo" ||
                     OptimizingCutElement1[i][j].RelatedVariable == "METRel" ||
                     OptimizingCutElement1[i][j].RelatedVariable == "MET" ||
                     OptimizingCutElement1[i][j].RelatedVariable == "meff" ||
                     OptimizingCutElement1[i][j].RelatedVariable == "mt1" ||
                     OptimizingCutElement1[i][j].RelatedVariable == "mtm" );
                    
                    OptimizingCutElement1[i][j].OnlyHasUpperCut =
                    (OptimizingCutElement1[i][j].RelatedVariable == "dEta" ||
                     OptimizingCutElement1[i][j].RelatedVariable == "mlj" );
                }
            }
            
            element.OptimizingCut.push_back(OptimizingCutElement1);
        }
        
        ///*
        //jet1
        {
            element.RegionName = "SR_SS_jet1_opt";
            element.RegionShortName = "SRjet1";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(3);
            element.setOfChannel.push_back(4);
            element.setOfChannel.push_back(5);
            
            element.Cut = "";
            
            element.AdditionalCut.clear();
            RegionInfo.push_back(element);
        }
        
        //jet23
        {
            element.RegionName = "SR_SS_jet23_opt";
            element.RegionShortName = "SRjet23";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(9);
            element.setOfChannel.push_back(10);
            element.setOfChannel.push_back(11);
            
            element.Cut = "";
            
            element.AdditionalCut.clear();
            RegionInfo.push_back(element);
        }
        //*/
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
        
        
        //Signal region for preselection
        GroupElement.GroupName = "SR_SS_pre";
        GroupElement.lower = RegionInfo.size();
        GroupElement.showData = false;
        GroupElement.showSignificance = true;
        
        ///*
        //jet1
        {
            element.RegionName = "SR_SS_jet1_pre";
            element.RegionShortName = "SRjet1";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(3);
            element.setOfChannel.push_back(4);
            element.setOfChannel.push_back(5);
            
            element.Cut = "";
            
            element.AdditionalCut.clear();
            element.OptimizingCut.clear();
            
            RegionInfo.push_back(element);
        }
        
        //jet23
        {
            element.RegionName = "SR_SS_jet23_pre";
            element.RegionShortName = "SRjet23";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(9);
            element.setOfChannel.push_back(10);
            element.setOfChannel.push_back(11);
            
            element.Cut = "";
            
            element.AdditionalCut.clear();
            element.OptimizingCut.clear();
            
            RegionInfo.push_back(element);
        }
        //*/
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
        
        
        //Validation region for Fake
        GroupElement.GroupName = "VR_Fake";
        GroupElement.lower = RegionInfo.size();
        GroupElement.showData = true;
        GroupElement.showSignificance = false;
        
        //ee
        {
            element.RegionName = "VR_Fake_ee";
            element.RegionShortName = "ee";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(3);
            element.setOfChannel.push_back(9);
            element.setOfChannel.push_back(15);
            
            element.Cut = "";
            element.Cut += " && meff > 200";
            element.Cut += " && fabs(mll - 91.2) > 10";
            element.Cut += " && MET > 30";
            
            element.AdditionalCut.clear();
            element.OptimizingCut.clear();
            
            RegionInfo.push_back(element);
        }
        
        //mumu
        {
            element.RegionName = "VR_Fake_mumu";
            element.RegionShortName = "mumu";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(4);
            element.setOfChannel.push_back(10);
            element.setOfChannel.push_back(16);
            
            element.Cut = "";
            element.Cut += " && meff > 200";
            element.Cut += " && fabs(mll - 91.2) > 10";
            element.Cut += " && MET > 30";
            
            element.AdditionalCut.clear();
            element.OptimizingCut.clear();
            
            RegionInfo.push_back(element);
        }
        
        //emu
        {
            element.RegionName = "VR_Fake_emu";
            element.RegionShortName = "emu";
            
            element.setOfChannel.clear();
            element.setOfChannel.push_back(5);
            element.setOfChannel.push_back(11);
            element.setOfChannel.push_back(17);
            
            element.Cut = "";
            element.Cut += " && meff > 200";
            //element.Cut += " && fabs(mll - 91.2) > 10";
            element.Cut += " && MET > 30";
            
            element.AdditionalCut.clear();
            element.OptimizingCut.clear();
            
            RegionInfo.push_back(element);
        }
        
        GroupElement.upper = RegionInfo.size() -1;
        RegionGroup.push_back(GroupElement);
    }
    
    struct SampleData
    {
        TString SampleName;
        int index;
    };
    
    std::vector<SampleData> BGVVData;
    if(doVVCount)
    {
        unsigned int VVGroupIndex = 0;
        for(unsigned int i=0;i<BGMCGroupData.size();i++)
        {
            if(BGMCGroupData[i].GroupName == "VV") VVGroupIndex = i;
        }
        for(unsigned int i=0;i<BGMCGroupData[VVGroupIndex].BGMCIndex.size();i++)
        {
            SampleData element;
            element.SampleName = BGMCSampleInfo[BGMCGroupData[VVGroupIndex].BGMCIndex[i]].SampleName.Data();
            element.SampleName.Remove(0,18);
            element.index = BGMCGroupData[VVGroupIndex].BGMCIndex[i];
            BGVVData.push_back(element);
        }
    }
    
    std::vector<unsigned int> OptimizingSignal;
    /*
    for(unsigned int i=0;i<10;i++)
    {
        OptimizingSignal.push_back(i);
    }
    */
    ///*
    for(unsigned int i=0;i<SigSampleInfo.size();i++)
    {
        /*
        if(
           (SigSampleInfo[i].Mass1 == 150   && SigSampleInfo[i].Mass2 == 0     ) ||
           (SigSampleInfo[i].Mass1 == 152.5 && SigSampleInfo[i].Mass2 == 22.5  ) ||
           (SigSampleInfo[i].Mass1 == 175   && SigSampleInfo[i].Mass2 == 25    ) ||
           (SigSampleInfo[i].Mass1 == 177.5 && SigSampleInfo[i].Mass2 == 47.5  ) ||
           (SigSampleInfo[i].Mass1 == 187.5 && SigSampleInfo[i].Mass2 == 37.5  ) ||
           (SigSampleInfo[i].Mass1 == 190   && SigSampleInfo[i].Mass2 == 60    ) ||
           (SigSampleInfo[i].Mass1 == 202.5 && SigSampleInfo[i].Mass2 == 72.5  ) ||
           (SigSampleInfo[i].Mass1 == 215   && SigSampleInfo[i].Mass2 == 85    ) ||
           (SigSampleInfo[i].Mass1 == 227.5 && SigSampleInfo[i].Mass2 == 97.5  ) ||
           (SigSampleInfo[i].Mass1 == 237.5 && SigSampleInfo[i].Mass2 == 87.5  ) ||
           (SigSampleInfo[i].Mass1 == 240   && SigSampleInfo[i].Mass2 == 110   ) ||
           (SigSampleInfo[i].Mass1 == 250   && SigSampleInfo[i].Mass2 == 100   ) ||
           (SigSampleInfo[i].Mass1 == 300   && SigSampleInfo[i].Mass2 == 100   ) )
        */
        
        ///*
        //Dani
        if(
           //(SigSampleInfo[i].Mass1 == 225   && SigSampleInfo[i].Mass2 == 75   ) || //SRjet1
           (SigSampleInfo[i].Mass1 == 187.5 && SigSampleInfo[i].Mass2 == 37.5 ) || //SRjet23
           false )
        //*/
        {
            cout<<SigSampleInfo[i].SampleName.Data()<<endl;
            OptimizingSignal.push_back(i);
        }
    }
    //*/
    
    //plot graph
    unsigned int countVariable = 0;
    for(unsigned int i=0;i<Var.size();i++)
    {
        if(Var[i].VarName == "phi1") countVariable = i;
    }
    
    TFile* fout_plot = new TFile("plot/fout.root","recreate");
    
    for(unsigned int RegionGroupIndex=10;RegionGroupIndex<=10;RegionGroupIndex++)
    //for(unsigned int RegionGroupIndex=11;RegionGroupIndex<=11;RegionGroupIndex++)
    //for(unsigned int RegionGroupIndex=12;RegionGroupIndex<=12;RegionGroupIndex++)
    //for(unsigned int RegionGroupIndex=0;RegionGroupIndex<RegionGroup.size();RegionGroupIndex++)
    {
        //For SR
        std::vector<Group> SRBGGroup;
        std::vector<TH1D*> h2SRSig;
        std::vector<TH1D*> h2SRSignificance;
        if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
        {
            //channelRepresentative for SR
            unsigned int channelRepresentative = 0;
            for(unsigned int ChannelIndex=0;ChannelIndex<ChannelInfo.size();ChannelIndex++)
            {
                //if(ChannelInfo[ChannelIndex].ChannelName == "nonISR_SS_ee") channelRepresentative = ChannelIndex;
                if(ChannelInfo[ChannelIndex].ChannelName == "jet1_SS_ee") channelRepresentative = ChannelIndex;
            }
            
            //Find all BG
            //For MC background
            for(unsigned int j=0;j<ChannelInfo[channelRepresentative].setOfBGMC.size();j++)
            {
                for(unsigned int k=0;k<BGMCGroupData.size();k++)
                {
                    if(ChannelInfo[channelRepresentative].setOfBGMC[j] == BGMCGroupData[k].GroupName)
                    {
                        Group BGGroupElement;
                        BGGroupElement.info = &(BGMCGroupData[k]);
                        BGGroupElement.h2 = new TH1D(("SR_"+BGMCGroupData[k].GroupName).Data(),"",6,0,6);
                        BGGroupElement.h2->SetLineColor(BGMCGroupData[k].colour);
                        BGGroupElement.h2->SetFillColor(BGMCGroupData[k].colour);
                        SRBGGroup.push_back(BGGroupElement);
                    }
                }
            }
            
            //For data-driven background
            for(unsigned int j=0;j<ChannelInfo[channelRepresentative].setOfBGData.size();j++)
            {
                for(unsigned int k=0;k<BGDataGroupData.size();k++)
                {
                    if(ChannelInfo[channelRepresentative].setOfBGData[j] == BGDataGroupData[k].GroupName)
                    {
                        Group BGGroupElement;
                        BGGroupElement.info = &(BGDataGroupData[k]);
                        BGGroupElement.h2 = new TH1D(("SR_"+BGMCGroupData[k].GroupName).Data(),"",6,0,6);
                        BGGroupElement.h2->SetLineColor(BGMCGroupData[k].colour);
                        BGGroupElement.h2->SetFillColor(BGMCGroupData[k].colour);
                        SRBGGroup.push_back(BGGroupElement);
                    }
                }
            }
            
            //h2SRSig
            for(unsigned int j=0;j<SigMassSplitting.size();j++)
            {
                TH1D* hTemp = new TH1D(SigMassSplitting[j].IDName.Data(),"",6,0,6);
                hTemp->SetLineColor(SigMassSplitting[j].colour);
                hTemp->SetLineStyle(1);
                h2SRSig.push_back(hTemp);
                
                hTemp = new TH1D((SigMassSplitting[j].IDName+"_significance").Data(),"",6,0,6);
                hTemp->SetLineColor(SigMassSplitting[j].colour);
                hTemp->SetLineStyle(1);
                h2SRSignificance.push_back(hTemp);
            }
        }
        
        //For combined significance
        for(unsigned int i=0;i<SigSampleInfo.size();i++)
        {
            SigSampleInfo[i].significance2 = 0;
        }
        
        //for(unsigned int RegionIndex=RegionGroup[RegionGroupIndex].lower+0;RegionIndex<=RegionGroup[RegionGroupIndex].lower+0;RegionIndex++)
        for(unsigned int RegionIndex=RegionGroup[RegionGroupIndex].lower;RegionIndex<=RegionGroup[RegionGroupIndex].upper;RegionIndex++)
        {
            const unsigned int channelRepresentative = RegionInfo[RegionIndex].setOfChannel[0];
            
            std::vector<TChain*> tree2Data;
            initializeTree2(tree2Data,RegionInfo[RegionIndex].setOfChannel,DataSampleName,ChannelInfo);
            
            std::vector< std::vector<TChain*> > tree2BGMC;
            std::vector< std::vector<double> > BGMCGroupXS;
            std::vector< std::vector<int> > BGMCGroupnAOD;
            std::vector< std::vector<double> > BGMCGroupnwAOD;
            std::vector<Group> BGGroup;
            {
                //For MC background
                for(unsigned int k=0;k<ChannelInfo[channelRepresentative].setOfBGMC.size();k++)
                {
                    for(unsigned int j=0;j<BGMCGroupData.size();j++)
                    {
                        if(BGMCGroupData[j].GroupName == ChannelInfo[channelRepresentative].setOfBGMC[k])
                        {
                            std::vector<TString> BGMCGroupSampleName;
                            std::vector<double> BGMCGroupXSElement;
                            for(unsigned int m=0;m<BGMCGroupData[j].BGMCIndex.size();m++)
                            {
                                BGMCGroupSampleName.push_back(BGMCSampleInfo[BGMCGroupData[j].BGMCIndex[m]].SampleName);
                                BGMCGroupXSElement.push_back(BGMCSampleInfo[BGMCGroupData[j].BGMCIndex[m]].XS);
                            }
                            
                            std::vector<TChain*> tree2BGMCElement;
                            initializeTree2(tree2BGMCElement,RegionInfo[RegionIndex].setOfChannel,BGMCGroupSampleName,ChannelInfo);
                            
                            std::vector<int> BGMCGroupnAODElement;
                            getnAOD(BGMCGroupnAODElement,RegionInfo[RegionIndex].setOfChannel,BGMCGroupSampleName,ChannelInfo);
                            
                            std::vector<double> BGMCGroupnwAODElement;
                            getnwAOD(BGMCGroupnwAODElement,RegionInfo[RegionIndex].setOfChannel,BGMCGroupSampleName,ChannelInfo);
                            
                            tree2BGMC.push_back(tree2BGMCElement);
                            BGMCGroupXS.push_back(BGMCGroupXSElement);
                            BGMCGroupnAOD.push_back(BGMCGroupnAODElement);
                            BGMCGroupnwAOD.push_back(BGMCGroupnwAODElement);
                            
                            Group BGGroupElement;
                            BGGroupElement.info = &(BGMCGroupData[j]);
                            BGGroup.push_back(BGGroupElement);
                            
                            //Z pt reweighting
                            if(dorw
                               &&
                               (
                                RegionInfo[RegionIndex].RegionName == "CR_nonISR_OS_ee"   ||
                                RegionInfo[RegionIndex].RegionName == "CR_ISR_OS_ee"      ||
                                RegionInfo[RegionIndex].RegionName == "CR_nonISR_OS_mumu" ||
                                RegionInfo[RegionIndex].RegionName == "CR_ISR_OS_mumu"
                                )
                               &&
                               (
                                BGGroupElement.info->GroupName == "Zee" ||
                                BGGroupElement.info->GroupName == "Zmumu" ||
                                BGGroupElement.info->GroupName == "Ztautau"
                                )
                               )
                            {
                                for(unsigned int k=0;k<BGGroupElement.info->BGMCIndex.size();k++)
                                {
                                    TString NameTemp = "tree_rw_";
                                    NameTemp += TString::Itoa(BGGroupElement.info->BGMCIndex[k],10);
                                    TChain* ch_rw = new TChain(NameTemp.Data());
                                    
                                    for(unsigned int m=0;m<RegionInfo[RegionIndex].setOfChannel.size();m++)
                                    {
                                        TString FileName = "skimming/rw_";
                                        FileName += TString::Itoa(RegionInfo[RegionIndex].setOfChannel[m],10);
                                        FileName += ".root";
                                        
                                        ch_rw->Add(FileName.Data());
                                    }
                                    tree2BGMCElement[k]->AddFriend(NameTemp.Data());
                                }
                            }
                            
                            //for charge filp BG
                            if(docfw
                               &&
                               (
                                RegionInfo[RegionIndex].RegionName == "CR_nonISR_SS_ee" ||
                                RegionInfo[RegionIndex].RegionName == "CR_ISR_SS_ee"    )
                               &&
                               BGGroupElement.info->GroupName == "Zee"
                               )
                            {
                                for(unsigned int k=0;k<BGGroupElement.info->BGMCIndex.size();k++)
                                {
                                    TString NameTemp = "tree_cfw_";
                                    NameTemp += TString::Itoa(BGGroupElement.info->BGMCIndex[k],10);
                                    TChain* ch_cfw = new TChain(NameTemp.Data());
                                    
                                    for(unsigned int m=0;m<RegionInfo[RegionIndex].setOfChannel.size();m++)
                                    {
                                        TString FileName = "skimming/cfw_";
                                        FileName += TString::Itoa(RegionInfo[RegionIndex].setOfChannel[m],10);
                                        FileName += ".root";
                                        
                                        ch_cfw->Add(FileName.Data());
                                    }
                                    tree2BGMCElement[k]->AddFriend(NameTemp.Data());
                                }
                            }
                            
                        }
                    }
                }
                
                //For data-driven background
                for(unsigned int k=0;k<ChannelInfo[channelRepresentative].setOfBGData.size();k++)
                {
                    for(unsigned int j=0;j<BGDataGroupData.size();j++)
                    {
                        if(BGDataGroupData[j].GroupName == ChannelInfo[channelRepresentative].setOfBGData[k])
                        {
                            Group BGGroupElement;
                            BGGroupElement.info = &(BGDataGroupData[j]);
                            BGGroup.push_back(BGGroupElement);
                        }
                    }
                }
            }
            
            std::vector<TChain*> tree2Sig;
            {
                std::vector<TString> SigSampleName;
                for(unsigned int i=0;i<SigSampleInfo.size();i++)
                {
                    SigSampleName.push_back(SigSampleInfo[i].SampleName);
                }
                initializeTree2(tree2Sig,RegionInfo[RegionIndex].setOfChannel,SigSampleName,ChannelInfo);
            }
            
            std::vector<TChain*> tree2DataOS;
            if(ChannelInfo[channelRepresentative].isSS_qF)
            {
                std::vector<unsigned int> setOfQFChannel;
                for(unsigned int index=0;index<RegionInfo[RegionIndex].setOfChannel.size();index++)
                {
                    const unsigned int ChannelIndex  = RegionInfo[RegionIndex].setOfChannel[index];
                    setOfQFChannel.push_back(ChannelInfo[ChannelIndex].qFChannel);
                }
                initializeTree2(tree2DataOS,setOfQFChannel,DataSampleName,ChannelInfo);
            }
            
            //Significance optimization
            const bool doOptimize2 = doOptimize && (
            RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt"  );
            if(doOptimize2)
            {
                for(unsigned int SigIndex=SigOptimizingIndex;SigIndex<=SigOptimizingIndex;SigIndex++)
                //for(unsigned int SigIndex=0;SigIndex<RegionInfo[RegionIndex].OptimizingCut.size();SigIndex++)
                {
                    double nBGRecord3 = 0;
                    double nSigRecord3[OptimizingSignal.size()];
                    double significanceRecord3 = -999;
                    double significanceRecord4 = -999;
                    
                    bool isOptimizing = true;
                    bool isRemovingCut = false;
                    
                    while(true)
                    {
                        const vector< vector<OptimizingCutData> > OptimizingCutInfo = RegionInfo[RegionIndex].OptimizingCut[SigIndex];
                        
                        unsigned int VarIndexRecord2 = 0;
                        vector<double> lowerCutRecord2;
                        vector<double> upperCutRecord2;
                        double nBGRecord2 = 0;
                        double nSigRecord2[OptimizingSignal.size()];
                        double significanceRecord2 = 0;
                        
                        if(isOptimizing) significanceRecord2 = significanceRecord3;
                        if(isRemovingCut) significanceRecord2 = -999;
                        
                        bool isChanged = false;
                        
                        for(unsigned int VarIndex2=0;VarIndex2<OptimizingCutInfo.size();VarIndex2++)
                        {
                            if(OptimizingCutInfo[VarIndex2][0].RelatedVariable == "meff" ||
                               OptimizingCutInfo[VarIndex2][0].RelatedVariable == "dEta" )
                            {
                                //continue;
                            }
                            
                            const unsigned int VarNumber = OptimizingCutInfo[VarIndex2].size();
                            
                            if(isRemovingCut)
                            {
                                bool HasCut = false;
                                for(unsigned int i=0;i<VarNumber;i++)
                                {
                                    if(OptimizingCutInfo[VarIndex2][i].lower != OptimizingCutInfo[VarIndex2][i].min ||
                                       OptimizingCutInfo[VarIndex2][i].upper != -1)
                                    {
                                        HasCut = true;
                                    }
                                }
                                
                                if(!HasCut) continue;
                            }
                            
                            TH1D* BGGroupRaw1D[BGGroup.size()];
                            TH1D* h2Sig1D[OptimizingSignal.size()];
                            
                            TH2D* BGGroupRaw2D[BGGroup.size()];
                            TH2D* h2Sig2D[OptimizingSignal.size()];
                            {
                                TString CommonCut = RegionInfo[RegionIndex].Cut;
                                
                                for(unsigned int i=0;i<OptimizingCutInfo.size();i++)
                                {
                                    if(i != VarIndex2)
                                    {
                                        for(unsigned int j=0;j<OptimizingCutInfo[i].size();j++)
                                        {
                                            unsigned int VarIndex3 = findVarIndex(OptimizingCutInfo[i][j].RelatedVariable,Var);
                                            
                                            CommonCut += " && ";
                                            CommonCut += Var[VarIndex3].VarFormula;
                                            CommonCut += " >= ";
                                            CommonCut += OptimizingCutInfo[i][j].lower;
                                            
                                            if(OptimizingCutInfo[i][j].upper != -1)
                                            {
                                                CommonCut += " && ";
                                                CommonCut += Var[VarIndex3].VarFormula;
                                                CommonCut += " < ";
                                                CommonCut += OptimizingCutInfo[i][j].upper;
                                            }
                                        }
                                    }
                                }
                                cout<<CommonCut.Data()<<endl;
                                
                                //initialize histograms
                                //Fill histograms from trees
                                vector<TString> VarFormula;
                                for(unsigned int i=0;i<VarNumber;i++)
                                {
                                    for(unsigned int j=0;j<Var.size();j++)
                                    {
                                        if(OptimizingCutInfo[VarIndex2][i].RelatedVariable == Var[j].VarName)
                                        {
                                            VarFormula.push_back(Var[j].VarFormula);
                                        }
                                    }
                                }
                                
                                //Background
                                //For MC background
                                for(unsigned int j=0;j<tree2BGMC.size();j++)
                                {
                                    if(VarNumber==1)
                                    {
                                        BGGroup[j].h2 = new TH1D(BGGroup[j].info->GroupName.Data(),"",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max);
                                        BGGroupRaw1D[j] = new TH1D((BGGroup[j].info->GroupName+"_raw").Data(),"",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max);
                                    }
                                    else if(VarNumber==2)
                                    {
                                        BGGroup[j].h3 = new TH2D(BGGroup[j].info->GroupName.Data(),"",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max,OptimizingCutInfo[VarIndex2][1].nBin,OptimizingCutInfo[VarIndex2][1].min,OptimizingCutInfo[VarIndex2][1].max);
                                        BGGroupRaw2D[j] = new TH2D((BGGroup[j].info->GroupName+"_raw").Data(),"",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max,OptimizingCutInfo[VarIndex2][1].nBin,OptimizingCutInfo[VarIndex2][1].min,OptimizingCutInfo[VarIndex2][1].max);
                                    }
                                    
                                    for(unsigned int k=0;k<tree2BGMC[j].size();k++)
                                    {
                                        {
                                            TH1D* hTemp1D = nullptr;
                                            TH2D* hTemp2D = nullptr;
                                            TString temp;
                                            if(VarNumber==1)
                                            {
                                                hTemp1D = new TH1D("BGMC","",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max);
                                                temp = VarFormula[0];
                                            }
                                            else if(VarNumber==2)
                                            {
                                                hTemp2D = new TH2D("BGMC","",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max,OptimizingCutInfo[VarIndex2][1].nBin,OptimizingCutInfo[VarIndex2][1].min,OptimizingCutInfo[VarIndex2][1].max);
                                                temp = VarFormula[1];
                                                temp += ":";
                                                temp += VarFormula[0];
                                            }
                                            
                                            //fill histograms from trees
                                            temp += ">>BGMC";
                                            
                                            //Weight
                                            TString Cut = "weight";
                                            
                                            //Cut
                                            Cut += "*(1";
                                            Cut += CommonCut;
                                            Cut += ")";
                                            
                                            tree2BGMC[j][k]->Draw(temp.Data(),Cut.Data());
                                            
                                            //normalization for BG
                                            if(VarNumber==1)
                                            {
                                                hTemp1D->Scale(BGMCGroupXS[j][k]/BGMCGroupnwAOD[j][k] *sumDataL);
                                                BGGroup[j].h2->Add(hTemp1D);
                                                delete hTemp1D;
                                            }
                                            else if(VarNumber==2)
                                            {
                                                hTemp2D->Scale(BGMCGroupXS[j][k]/BGMCGroupnwAOD[j][k] *sumDataL);
                                                BGGroup[j].h3->Add(hTemp2D);
                                                delete hTemp2D;
                                            }
                                        }
                                        {
                                            TH1D* hTemp1D = nullptr;
                                            TH2D* hTemp2D = nullptr;
                                            TString temp;
                                            if(VarNumber==1)
                                            {
                                                hTemp1D = new TH1D("BGMCRaw","",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max);
                                                temp = VarFormula[0];
                                            }
                                            else if(VarNumber==2)
                                            {
                                                hTemp2D = new TH2D("BGMCRaw","",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max,OptimizingCutInfo[VarIndex2][1].nBin,OptimizingCutInfo[VarIndex2][1].min,OptimizingCutInfo[VarIndex2][1].max);
                                                temp = VarFormula[1];
                                                temp += ":";
                                                temp += VarFormula[0];
                                            }
                                            
                                            //fill histograms from trees
                                            temp += ">>BGMCRaw";
                                            
                                            //Weight
                                            TString Cut = "1";
                                            
                                            //Cut
                                            Cut += "*(1";
                                            Cut += CommonCut;
                                            Cut += ")";
                                            
                                            tree2BGMC[j][k]->Draw(temp.Data(),Cut.Data());
                                            
                                            if(VarNumber==1)
                                            {
                                                BGGroupRaw1D[j]->Add(hTemp1D);
                                                delete hTemp1D;
                                            }
                                            else if(VarNumber==2)
                                            {
                                                BGGroupRaw2D[j]->Add(hTemp2D);
                                                delete hTemp2D;
                                            }
                                        }
                                    }
                                }
                                
                                //For data-driven background
                                for(unsigned int j=tree2BGMC.size();j<BGGroup.size();j++)
                                {
                                    if(VarNumber==1)
                                    {
                                        BGGroup[j].h2 = new TH1D(BGGroup[j].info->GroupName.Data(),"",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max);
                                        BGGroupRaw1D[j] = new TH1D((BGGroup[j].info->GroupName+"_raw").Data(),"",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max);
                                    }
                                    else if(VarNumber==2)
                                    {
                                        BGGroup[j].h3 = new TH2D(BGGroup[j].info->GroupName.Data(),"",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max,OptimizingCutInfo[VarIndex2][1].nBin,OptimizingCutInfo[VarIndex2][1].min,OptimizingCutInfo[VarIndex2][1].max);
                                        BGGroupRaw2D[j] = new TH2D((BGGroup[j].info->GroupName+"_raw").Data(),"",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max,OptimizingCutInfo[VarIndex2][1].nBin,OptimizingCutInfo[VarIndex2][1].min,OptimizingCutInfo[VarIndex2][1].max);
                                    }
                                    
                                    for(unsigned int k=0;k<DataSampleName.size();k++)
                                    {
                                        {
                                            TH1D* hTemp1D = nullptr;
                                            TH2D* hTemp2D = nullptr;
                                            TString temp;
                                            if(VarNumber==1)
                                            {
                                                hTemp1D = new TH1D("BGData","",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max);
                                                temp = VarFormula[0];
                                            }
                                            else if(VarNumber==2)
                                            {
                                                hTemp2D = new TH2D("BGData","",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max,OptimizingCutInfo[VarIndex2][1].nBin,OptimizingCutInfo[VarIndex2][1].min,OptimizingCutInfo[VarIndex2][1].max);
                                                temp = VarFormula[1];
                                                temp += ":";
                                                temp += VarFormula[0];
                                            }
                                            
                                            //fill histograms from trees
                                            temp += ">>BGData";
                                            
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
                                                Cut += " && fLwt==0";
                                            }
                                            if(BGGroup[j].info->GroupName == "fake lepton")
                                            {
                                                Cut += " && fLwt!=0";
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
                                            if(VarNumber==1)
                                            {
                                                BGGroup[j].h2->Add(hTemp1D);
                                                delete hTemp1D;
                                            }
                                            else if(VarNumber==2)
                                            {
                                                BGGroup[j].h3->Add(hTemp2D);
                                                delete hTemp2D;
                                            }
                                        }
                                        {
                                            TH1D* hTemp1D = nullptr;
                                            TH2D* hTemp2D = nullptr;
                                            TString temp;
                                            if(VarNumber==1)
                                            {
                                                hTemp1D = new TH1D("BGDataRaw","",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max);
                                                temp = VarFormula[0];
                                            }
                                            else if(VarNumber==2)
                                            {
                                                hTemp2D = new TH2D("BGDataRaw","",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max,OptimizingCutInfo[VarIndex2][1].nBin,OptimizingCutInfo[VarIndex2][1].min,OptimizingCutInfo[VarIndex2][1].max);
                                                temp = VarFormula[1];
                                                temp += ":";
                                                temp += VarFormula[0];
                                            }
                                            
                                            //fill histograms from trees
                                            temp += ">>BGDataRaw";
                                            
                                            TString Cut = "1";
                                            
                                            Cut += "*(1";
                                            Cut += CommonCut;
                                            if(BGGroup[j].info->GroupName == "charge flip")
                                            {
                                                Cut += " && fLwt==0";
                                            }
                                            if(BGGroup[j].info->GroupName == "fake lepton")
                                            {
                                                Cut += " && fLwt!=0";
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
                                            if(VarNumber==1)
                                            {
                                                BGGroupRaw1D[j]->Add(hTemp1D);
                                                delete hTemp1D;
                                            }
                                            else if(VarNumber==2)
                                            {
                                                BGGroupRaw2D[j]->Add(hTemp2D);
                                                delete hTemp2D;
                                            }
                                        }
                                    }
                                }
                                
                                //h2Sig
                                for(unsigned int j=0;j<OptimizingSignal.size();j++)
                                {
                                    TString NameTemp = "Sig_";
                                    NameTemp += TString::Itoa(j,10);
                                    
                                    TString temp;
                                    if(VarNumber==1)
                                    {
                                        h2Sig1D[j] = new TH1D(NameTemp.Data(),"",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max);
                                        temp = VarFormula[0];
                                    }
                                    else if(VarNumber==2)
                                    {
                                        h2Sig2D[j] = new TH2D(NameTemp.Data(),"",OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].min,OptimizingCutInfo[VarIndex2][0].max,OptimizingCutInfo[VarIndex2][1].nBin,OptimizingCutInfo[VarIndex2][1].min,OptimizingCutInfo[VarIndex2][1].max);
                                        temp = VarFormula[1];
                                        temp += ":";
                                        temp += VarFormula[0];
                                    }
                                    
                                    //Fill Signal
                                    temp += ">>";
                                    temp += NameTemp;
                                    
                                    TString Cut = "weight*(1";
                                    Cut += CommonCut;
                                    Cut += ")";
                                    
                                    tree2Sig[OptimizingSignal[j]]->Draw(temp.Data(),Cut.Data());
                                    
                                    if(VarNumber==1)
                                    {
                                        h2Sig1D[j]->Scale(SigSampleInfo[OptimizingSignal[j]].XS /SigSampleInfo[OptimizingSignal[j]].nwAOD *sumDataL);
                                    }
                                    else if(VarNumber==2)
                                    {
                                        h2Sig2D[j]->Scale(SigSampleInfo[OptimizingSignal[j]].XS /SigSampleInfo[OptimizingSignal[j]].nwAOD *sumDataL);
                                    }
                                }
                            }
                            
                            int lowerBinRecord1[VarNumber];
                            int upperBinRecord1[VarNumber];
                            double nBGRecord1 = -1;
                            double nSigRecord1[OptimizingSignal.size()];
                            double significanceRecord1 = -999;
                            bool isPrint = false;
                            //isPrint = (OptimizingCutInfo[VarIndex2][0].RelatedVariable == "dEta");
                            
                            int bin1[VarNumber];
                            int bin2[VarNumber];
                            for(unsigned int i=0;i<VarNumber;i++)
                            {
                                bin1[i] = 1;
                                bin2[i] = OptimizingCutInfo[VarIndex2][i].nBin+1;
                            }
                            
                            while(true)
                            {
                                double nBG = 0;
                                double nBGError2 = 0;
                                double nSig[OptimizingSignal.size()];
                                double nSigError[OptimizingSignal.size()];
                                bool skip = false;
                                
                                if(isPrint)
                                {
                                    if(VarNumber==1)
                                    {
                                        cout<<"bin1[0]: "<<bin1[0]<<", bin2[0]: "<<bin2[0]<<endl;
                                    }
                                    else if(VarNumber==2)
                                    {
                                        cout<<"bin1[0]: "<<bin1[0]<<", bin2[0]: "<<bin2[0]<<", bin1[1]: "<<bin1[1]<<", bin2[1]: "<<bin2[1]<<endl;
                                    }
                                }
                                
                                if(!skip)
                                {
                                    for(unsigned int j=0;j<BGGroup.size();j++)
                                    {
                                        if(BGGroup[j].info->GroupName == "VV" ||
                                           BGGroup[j].info->GroupName == "ttV")
                                        {
                                            double nBGRaw = 0;
                                            if(VarNumber==1)
                                            {
                                                nBGRaw = BGGroupRaw1D[j]->Integral(bin1[0],bin2[0]);
                                            }
                                            else if(VarNumber==2)
                                            {
                                                nBGRaw = BGGroupRaw2D[j]->Integral(bin1[0],bin2[0],bin1[1],bin2[1]);
                                            }
                                            
                                            if(nBGRaw<=10)
                                            {
                                                if(isPrint) cout<<"VV or ttV do not have 10 raw events."<<endl;
                                                skip = true;
                                            }
                                        }
                                        
                                        double Error = 0;
                                        double n = 0;
                                        if(VarNumber==1)
                                        {
                                            n = BGGroup[j].h2->IntegralAndError(bin1[0],bin2[0],Error);
                                        }
                                        else if(VarNumber==2)
                                        {
                                            n = BGGroup[j].h3->IntegralAndError(bin1[0],bin2[0],bin1[1],bin2[1],Error);
                                        }
                                        
                                        if(isPrint && !skip) cout<<BGGroup[j].info->GroupName<<": "<<n<<" +/- "<<Error<<endl;
                                        if(n>0)
                                        {
                                            nBG += n;
                                        }
                                        else
                                        {
                                            if(BGGroup[j].info->GroupName != "multitop" &&
                                               BGGroup[j].info->GroupName != "Vgamma"   &&
                                               BGGroup[j].info->GroupName != "Wjets"    &&
                                               BGGroup[j].info->GroupName != "DY"       &&
                                               BGGroup[j].info->GroupName != "Zjets"   )
                                            {
                                                if(isPrint) cout<<"Group of BG is not positive: "<<BGGroup[j].info->GroupName<<endl;
                                                skip = true;
                                            }
                                        }
                                        nBGError2 += Error * Error;
                                    }
                                    if(isPrint && skip) cout<<endl;
                                    if(isPrint && !skip) cout<<"nBG: "<<nBG<<", nBGError: "<<sqrt(nBGError2)<<endl;
                                }
                                
                                for(unsigned int i=0;i<OptimizingSignal.size();i++)
                                {
                                    //expected number of events for signal
                                    if(VarNumber==1)
                                    {
                                        nSig[i] = h2Sig1D[i]->IntegralAndError(bin1[0],bin2[0],nSigError[i]);
                                    }
                                    else if(VarNumber==2)
                                    {
                                        nSig[i] = h2Sig2D[i]->IntegralAndError(bin1[0],bin2[0],bin1[1],bin2[1],nSigError[i]);
                                    }
                                    
                                    if(nSig[i] <= 1) skip = true;
                                    else if(nSigError[i] == 0) skip = true;
                                    else if(nSig[i]/nSigError[i]<=2) skip = true;
                                    
                                    if(isPrint && skip) cout<<"Signal "<<i<<" has low statistics."<<endl<<endl;
                                    if(isPrint && !skip) cout<<"nSig["<<i<<"]: "<<nSig[i]<<", nSigError["<<i<<"]: "<<nSigError[i]<<endl;
                                    
                                    if(skip) break;
                                }
                                
                                //Significance
                                if(!skip)
                                {
                                    double significanceTemp = 0;
                                    for(unsigned int i=0;i<OptimizingSignal.size();i++)
                                    {
                                        significanceTemp += GetSignificance(nSig[i],nBG,nBGError2);
                                    }
                                    significanceTemp /= OptimizingSignal.size();
                                    if(isPrint) cout<<"Significance: "<<significanceTemp;
                                    
                                    if(isOptimizing)
                                    {
                                        if(significanceTemp > significanceRecord1)
                                        {
                                            for(unsigned int i=0;i<VarNumber;i++)
                                            {
                                                lowerBinRecord1[i] = bin1[i];
                                                upperBinRecord1[i] = bin2[i];
                                            }
                                            
                                            nBGRecord1 = nBG;
                                            for(unsigned int i=0;i<OptimizingSignal.size();i++)
                                            {
                                                nSigRecord1[i] = nSig[i];
                                            }
                                            significanceRecord1 = significanceTemp;
                                            if(isPrint) cout<<", accepted."<<endl<<endl;
                                        }
                                        else
                                        {
                                            if(isPrint) cout<<endl<<endl;
                                        }
                                    }
                                    
                                    if(isRemovingCut)
                                    {
                                        nBGRecord1 = nBG;
                                        for(unsigned int i=0;i<OptimizingSignal.size();i++)
                                        {
                                            nSigRecord1[i] = nSig[i];
                                        }
                                        significanceRecord1 = significanceTemp;
                                        
                                        if(significanceRecord1 > significanceRecord3)
                                        {
                                            cout<<"Error1"<<endl;
                                            return;
                                        }
                                    }
                                }
                                else if(isRemovingCut)
                                {
                                    cout<<"Error2"<<endl;
                                    //return;
                                }
                                
                                if(isOptimizing)
                                {
                                    if(VarNumber==1)
                                    {
                                        bool isContinue = NextBin(bin1[0],bin2[0],OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].OnlyHasLowerCut,OptimizingCutInfo[VarIndex2][0].OnlyHasUpperCut);
                                        if(!isContinue) break;
                                    }
                                    else if(VarNumber==2)
                                    {
                                        bool isContinue = NextBin(bin1[1],bin2[1],OptimizingCutInfo[VarIndex2][1].nBin,OptimizingCutInfo[VarIndex2][1].OnlyHasLowerCut,OptimizingCutInfo[VarIndex2][1].OnlyHasUpperCut);
                                        if(!isContinue)
                                        {
                                            isContinue = NextBin(bin1[0],bin2[0],OptimizingCutInfo[VarIndex2][0].nBin,OptimizingCutInfo[VarIndex2][0].OnlyHasLowerCut,OptimizingCutInfo[VarIndex2][0].OnlyHasUpperCut);
                                            if(isContinue)
                                            {
                                                bin1[1] = 1;
                                                bin2[1] = -1;
                                            }
                                            else
                                            {
                                                break;
                                            }
                                        }
                                    }
                                }
                                
                                if(isRemovingCut)
                                {
                                    break;
                                }
                            }
                            
                            double lowerCutRecord1[VarNumber];
                            double upperCutRecord1[VarNumber];
                            if(isOptimizing)
                            {
                                for(unsigned int i=0;i<VarNumber;i++)
                                {
                                    double BinSize = (OptimizingCutInfo[VarIndex2][i].max - OptimizingCutInfo[VarIndex2][i].min)/OptimizingCutInfo[VarIndex2][i].nBin;
                                    
                                    lowerCutRecord1[i] = OptimizingCutInfo[VarIndex2][i].min + (lowerBinRecord1[i]-1) * BinSize;
                                    
                                    if(upperBinRecord1[i] != OptimizingCutInfo[VarIndex2][i].nBin+1) upperCutRecord1[i] = OptimizingCutInfo[VarIndex2][i].min + upperBinRecord1[i] * BinSize;
                                    else upperCutRecord1[i] = -1;
                                    
                                    cout<<OptimizingCutInfo[VarIndex2][i].RelatedVariable.Data()<<": ";
                                    cout<<"lowerBinRecord1: "<<lowerBinRecord1[i];
                                    cout<<", upperBinRecord1: "<<upperBinRecord1[i];
                                    cout<<", lowerCutRecord1: "<<lowerCutRecord1[i];
                                    cout<<", upperCutRecord1: "<<upperCutRecord1[i];
                                    cout<<endl;
                                }
                            }
                            if(isRemovingCut)
                            {
                                for(unsigned int i=0;i<VarNumber;i++)
                                {
                                    cout<<OptimizingCutInfo[VarIndex2][i].RelatedVariable.Data()<<", ";
                                }
                            }
                            cout<<"nBG: "<<nBGRecord1<<", nSig: "<<nSigRecord1[0]<<", significanceRecord1: "<<significanceRecord1<<endl;
                            
                            if(isOptimizing)
                            {
                                if(significanceRecord1 > significanceRecord2)
                                {
                                    bool isSameCut = true;
                                    for(unsigned int i=0;i<VarNumber;i++)
                                    {
                                        if(lowerCutRecord1[i] != OptimizingCutInfo[VarIndex2][i].lower) isSameCut = false;
                                        if(upperCutRecord1[i] != OptimizingCutInfo[VarIndex2][i].upper) isSameCut = false;
                                    }
                                    
                                    if(isSameCut)
                                    {
                                        cout<<"The cut is the same, but the significance increase becuase we are using different variable to draw."<<endl;
                                    }
                                    else
                                    {
                                        cout<<"The cut is changed. significanceRecord2: "<<significanceRecord2<<", significanceRecord1: "<<significanceRecord1<<endl;
                                        
                                        VarIndexRecord2 = VarIndex2;
                                        lowerCutRecord2.clear();
                                        upperCutRecord2.clear();
                                        for(unsigned int i=0;i<VarNumber;i++)
                                        {
                                            lowerCutRecord2.push_back(lowerCutRecord1[i]);
                                            upperCutRecord2.push_back(upperCutRecord1[i]);
                                        }
                                        nBGRecord2 = nBGRecord1;
                                        for(unsigned int i=0;i<OptimizingSignal.size();i++)
                                        {
                                            nSigRecord2[i] = nSigRecord1[i];
                                        }
                                        significanceRecord2 = significanceRecord1;
                                        
                                        isChanged = true;
                                    }
                                }
                            }
                            
                            if(isRemovingCut)
                            {
                                if(significanceRecord1 != -999 && (significanceRecord4 - significanceRecord1) < 0.05 && significanceRecord1 > significanceRecord2)
                                {
                                    cout<<"Higher significance is found. significanceRecord2: "<<significanceRecord2<<", significanceRecord1: "<<significanceRecord1<<endl;
                                    
                                    VarIndexRecord2 = VarIndex2;
                                    nBGRecord2 = nBGRecord1;
                                    for(unsigned int i=0;i<OptimizingSignal.size();i++)
                                    {
                                        nSigRecord2[i] = nSigRecord1[i];
                                    }
                                    significanceRecord2 = significanceRecord1;
                                    
                                    isChanged = true;
                                }
                            }
                            
                            cout<<endl;
                            
                            //delete
                            for(unsigned int j=0;j<BGGroup.size();j++)
                            {
                                if(VarNumber==1)
                                {
                                    delete BGGroup[j].h2;
                                    delete BGGroupRaw1D[j];
                                }
                                else if(VarNumber==2)
                                {
                                    delete BGGroup[j].h3;
                                    delete BGGroupRaw2D[j];
                                }
                            }
                            
                            for(unsigned int j=0;j<OptimizingSignal.size();j++)
                            {
                                if(VarNumber==1)
                                {
                                    delete h2Sig1D[j];
                                }
                                else if(VarNumber==2)
                                {
                                    delete h2Sig2D[j];
                                }
                            }
                        }
                        
                        if(isOptimizing)
                        {
                            if(isChanged)
                            {
                                cout<<"The following cut is applied."<<endl;
                                for(unsigned int i=0;i<OptimizingCutInfo[VarIndexRecord2].size();i++)
                                {
                                    cout<<"The "<<OptimizingCutInfo[VarIndexRecord2][i].RelatedVariable.Data()<<" cut is applied. ";
                                    cout<<"lowerCutRecord2: "<<lowerCutRecord2[i];
                                    cout<<", upperCutRecord2: "<<upperCutRecord2[i];
                                    cout<<endl;
                                }
                                cout<<"nBG: "<<nBGRecord2<<", nSig: "<<nSigRecord2[0]<<", significanceRecord2: "<<significanceRecord2<<endl<<endl;
                                
                                for(unsigned int i=0;i<OptimizingCutInfo[VarIndexRecord2].size();i++)
                                {
                                    RegionInfo[RegionIndex].OptimizingCut[SigIndex][VarIndexRecord2][i].lower = lowerCutRecord2[i];
                                    RegionInfo[RegionIndex].OptimizingCut[SigIndex][VarIndexRecord2][i].upper = upperCutRecord2[i];
                                }
                                
                                nBGRecord3 = nBGRecord2;
                                for(unsigned int i=0;i<OptimizingSignal.size();i++)
                                {
                                    nSigRecord3[i] = nSigRecord2[i];
                                }
                                significanceRecord3 = significanceRecord2;
                            }
                            else
                            {
                                cout<<"The cut is the same."<<endl;
                                cout<<"nBG: "<<nBGRecord3<<", nSig: "<<nSigRecord3[0]<<", significanceRecord3: "<<significanceRecord3<<endl<<endl;
                                
                                cout<<RegionInfo[RegionIndex].RegionName.Data()<<": "<<endl;
                                
                                for(unsigned int i=0;i<OptimizingCutInfo.size();i++)
                                {
                                    for(unsigned int j=0;j<OptimizingCutInfo[i].size();j++)
                                    {
                                        unsigned int VarIndex = findVarIndex(OptimizingCutInfo[i][j].RelatedVariable,Var);
                                        
                                        cout<<OptimizingCutInfo[i][j].lower;
                                        cout<<" <= ";
                                        cout<<Var[VarIndex].VarFormula.Data();
                                        cout<<" < ";
                                        cout<<OptimizingCutInfo[i][j].upper;
                                        cout<<endl;
                                    }
                                }
                                cout<<endl;
                                
                                cout<<"Start to remove cuts."<<endl<<endl;
                                isOptimizing = false;
                                isRemovingCut = true;
                                significanceRecord4 = significanceRecord3;
                            }
                        }
                        else if(isRemovingCut)
                        {
                            if(isChanged)
                            {
                                cout<<"The following cut is removed."<<endl;
                                for(unsigned int i=0;i<OptimizingCutInfo[VarIndexRecord2].size();i++)
                                {
                                    cout<<OptimizingCutInfo[VarIndexRecord2][i].RelatedVariable.Data()<<": ";
                                    cout<<"lower cut: "<<OptimizingCutInfo[VarIndexRecord2][i].lower;
                                    cout<<", upper cut: "<<OptimizingCutInfo[VarIndexRecord2][i].upper;
                                    cout<<endl;
                                    
                                    RegionInfo[RegionIndex].OptimizingCut[SigIndex][VarIndexRecord2][i].lower = OptimizingCutInfo[VarIndexRecord2][i].min;
                                    RegionInfo[RegionIndex].OptimizingCut[SigIndex][VarIndexRecord2][i].upper = -1;
                                }
                                cout<<"nBG: "<<nBGRecord2<<", nSig: "<<nSigRecord2[0]<<", significanceRecord2: "<<significanceRecord2<<endl<<endl;
                                
                                nBGRecord3 = nBGRecord2;
                                for(unsigned int i=0;i<OptimizingSignal.size();i++)
                                {
                                    nSigRecord3[i] = nSigRecord2[i];
                                }
                                significanceRecord3 = significanceRecord2;
                            }
                            else
                            {
                                cout<<"No any cut can be removed."<<endl;
                                cout<<"nBG: "<<nBGRecord3<<", nSig: "<<nSigRecord3[0];
                                cout<<", significanceRecord4: "<<significanceRecord4;
                                cout<<", significanceRecord3: "<<significanceRecord3<<endl<<endl;
                                
                                cout<<RegionInfo[RegionIndex].RegionName.Data()<<": "<<endl;
                                
                                for(unsigned int i=0;i<OptimizingCutInfo.size();i++)
                                {
                                    for(unsigned int j=0;j<OptimizingCutInfo[i].size();j++)
                                    {
                                        unsigned int VarIndex = findVarIndex(OptimizingCutInfo[i][j].RelatedVariable,Var);
                                        
                                        cout<<OptimizingCutInfo[i][j].lower;
                                        cout<<" <= ";
                                        cout<<Var[VarIndex].VarFormula.Data();
                                        cout<<" < ";
                                        cout<<OptimizingCutInfo[i][j].upper;
                                        cout<<endl;
                                    }
                                }
                                cout<<endl;
                                
                                TString PathName = "latex/data/optimization/cut_txt/cut_";
                                PathName += RegionInfo[RegionIndex].RegionName;
                                PathName += "_";
                                PathName += TString::Itoa(SigIndex,10);
                                PathName += ".txt";
                                
                                ofstream fout;
                                fout.open(PathName.Data());
                                
                                for(unsigned int i=0;i<OptimizingCutInfo.size();i++)
                                {
                                    for(unsigned int j=0;j<OptimizingCutInfo[i].size();j++)
                                    {
                                        unsigned int VarIndex = findVarIndex(OptimizingCutInfo[i][j].RelatedVariable,Var);
                                        
                                        fout<<Var[VarIndex].VarName;
                                        fout<<" ";
                                        fout<<OptimizingCutInfo[i][j].lower;
                                        fout<<" ";
                                        fout<<OptimizingCutInfo[i][j].upper;
                                        fout<<endl;
                                    }
                                }
                                fout.close();
                                
                                break;
                            }
                        }
                    }
                }
            }
            
            
            //for(unsigned int VarIndex=5;VarIndex<=5;VarIndex++) //mll
            //for(unsigned int VarIndex=6;VarIndex<=6;VarIndex++) //ptll
            for(unsigned int VarIndex=countVariable;VarIndex<=countVariable;VarIndex++)
            //for(unsigned int VarIndex=0;VarIndex<Var.size();VarIndex++)
            {
                //if(Var[VarIndex].VarName!="pt1") continue;
                
                if(Var[VarIndex].VarName!="pt1"    &&
                   Var[VarIndex].VarName!="pt2"    &&
                   Var[VarIndex].VarName!="mll"    &&
                   Var[VarIndex].VarName!="dEta"   &&
                   Var[VarIndex].VarName!="MET"    &&
                   Var[VarIndex].VarName!="mTtwo"  &&
                   Var[VarIndex].VarName!="meff"   &&
                   Var[VarIndex].VarName!="mt1"    &&
                   Var[VarIndex].VarName!="mlj"    &&
                   Var[VarIndex].VarName!="nJet"   &&
                   Var[VarIndex].VarName!="phi1"   )
                {
                    continue;
                }
                
                //initialize histograms
                TString title = Var[VarIndex].VarTitle;
                
                TString xaxis = Var[VarIndex].VarTitle;
                xaxis += " ";
                xaxis += Var[VarIndex].unit;
                
                //Fill histograms from trees
                TH1D* h2DataSum;
                TH1D* h2SigSum[SigMassSplitting.size()];
                TH1D* h2Sig[SigSampleInfo.size()];
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
                    
                    // mll>15 for OS_ee and OS_mumu
                    if((
                       ChannelInfo[channelRepresentative].ChannelName == "nonISR_OS_ee"   ||
                       ChannelInfo[channelRepresentative].ChannelName == "nonISR_OS_mumu" ||
                       ChannelInfo[channelRepresentative].ChannelName == "ISR_OS_ee"      ||
                       ChannelInfo[channelRepresentative].ChannelName == "ISR_OS_mumu"    )
                       && Var[VarIndex].VarName!="mll")
                       
                    {
                        CommonCut += " && mll>15";
                    }
                    
                    CommonCut += RegionInfo[RegionIndex].Cut;
                    
                    if(doOptimize2)
                    {
                        for(unsigned int i=0;i<RegionInfo[RegionIndex].OptimizingCut[SigOptimizingIndex].size();i++)
                        {
                            for(unsigned int j=0;j<RegionInfo[RegionIndex].OptimizingCut[SigOptimizingIndex][i].size();j++)
                            {
                                if(RegionInfo[RegionIndex].OptimizingCut[SigOptimizingIndex][i][j].RelatedVariable != Var[VarIndex].VarName)
                                {
                                    unsigned int VarIndex2 = findVarIndex(RegionInfo[RegionIndex].OptimizingCut[SigOptimizingIndex][i][j].RelatedVariable,Var);
                                    
                                    CommonCut += " && ";
                                    CommonCut += Var[VarIndex2].VarFormula;
                                    CommonCut += " >= ";
                                    CommonCut += RegionInfo[RegionIndex].OptimizingCut[SigOptimizingIndex][i][j].lower;
                                    
                                    if(RegionInfo[RegionIndex].OptimizingCut[SigOptimizingIndex][i][j].upper != -1)
                                    {
                                        CommonCut += " && ";
                                        CommonCut += Var[VarIndex2].VarFormula;
                                        CommonCut += " < ";
                                        CommonCut += RegionInfo[RegionIndex].OptimizingCut[SigOptimizingIndex][i][j].upper;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
                        {
                            TString PathName = "latex/data/optimization/";
                            if(useDani) PathName += "cut_txt_Dani/";
                            else PathName += "cut_txt/";
                            PathName += "cut_";
                            PathName += RegionInfo[RegionIndex].RegionName;
                            PathName += "_";
                            PathName += TString::Itoa(SigOptimizingIndex,10);
                            PathName += ".txt";
                            
                            ifstream fin;
                            fin.open(PathName.Data());
                            
                            for(unsigned int i=0;i<RegionInfo[RegionIndex].OptimizingCut[SigOptimizingIndex].size();i++)
                            {
                                for(unsigned int j=0;j<RegionInfo[RegionIndex].OptimizingCut[SigOptimizingIndex][i].size();j++)
                                {
                                    TString VarName;
                                    double cut;
                                    
                                    fin>>VarName;
                                    if(VarName != Var[VarIndex].VarName)
                                    {
                                        unsigned int VarIndex2 = findVarIndex(RegionInfo[RegionIndex].OptimizingCut[SigOptimizingIndex][i][j].RelatedVariable,Var);
                                        
                                        CommonCut += " && ";
                                        CommonCut += Var[VarIndex2].VarFormula;
                                        CommonCut += " >= ";
                                        
                                        fin>>cut;
                                        CommonCut += cut;
                                        
                                        fin>>cut;
                                        if(cut != -1)
                                        {
                                            CommonCut += " && ";
                                            CommonCut += Var[VarIndex2].VarFormula;
                                            CommonCut += " < ";
                                            CommonCut += cut;
                                        }
                                    }
                                    else
                                    {
                                        fin>>cut;
                                        fin>>cut;
                                    }
                                }
                            }
                            
                            fin.close();
                        }
                        else
                        {
                            for(unsigned int i=0;i<RegionInfo[RegionIndex].AdditionalCut.size();i++)
                            {
                                if(RegionInfo[RegionIndex].AdditionalCut[i].RelatedVariable != Var[VarIndex].VarName)
                                {
                                    CommonCut += RegionInfo[RegionIndex].AdditionalCut[i].Cut;
                                }
                            }
                        }
                    }
                    //cout<<CommonCut.Data()<<endl;
                    
                    //h2DataSum
                    {
                        h2DataSum = new TH1D("DataSum",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                        h2DataSum->GetYaxis()->SetTitle("Number of events");
                        h2DataSum->SetMarkerColor(1);
                        h2DataSum->SetMarkerStyle(20);
                        h2DataSum->SetMarkerSize(1);
                        h2DataSum->SetLineColor(1);
                    }
                    
                    //h2Data
                    TH1D* h2Data[DataSampleName.size()];
                    std::vector<TString> hName2Data;
                    for(unsigned int j=0;j<DataSampleName.size();j++)
                    {
                        TString NameTemp = "Data_";
                        NameTemp += TString::Itoa(j,10);
                        h2Data[j] = new TH1D(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                        hName2Data.push_back(NameTemp);
                        
                        //Fill data
                        TString temp;
                        if(Var[VarIndex].VarName=="averageMu")
                        {
                            temp = "averageMu/1.09";
                        }
                        else
                        {
                            temp = Var[VarIndex].VarFormula;
                        }
                        temp += ">>";
                        temp += hName2Data[j];
                        
                        TString Cut = "(1";
                        Cut += CommonCut;
                        Cut += " && fLwt==0";
                        Cut += ")";
                        tree2Data[j]->Draw(temp.Data(),Cut.Data());
                        
                        //Add data
                        h2DataSum->Add(h2Data[j]);
                    }
                    
                    //Background
                    //For MC background
                    for(unsigned int j=0;j<tree2BGMC.size();j++)
                    {
                        BGGroup[j].h2 = new TH1D(BGGroup[j].info->GroupName.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                        BGGroup[j].h2->GetYaxis()->SetTitle("Number of events");
                        BGGroup[j].h2->SetLineColor(BGGroup[j].info->colour);
                        BGGroup[j].h2->SetFillColor(BGGroup[j].info->colour);
                        
                        BGGroup[j].info->unweighted = 0;
                        for(unsigned int k=0;k<tree2BGMC[j].size();k++)
                        {
                            TH1D* hTemp = new TH1D("BGMC",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                            
                            //fill histograms from trees
                            TString temp = Var[VarIndex].VarFormula;
                            temp += ">>BGMC";
                            
                            //Weight
                            TString Cut = "weight";
                            
                            //Z pt reweighting
                            if(dorw
                               &&
                               (
                                RegionInfo[RegionIndex].RegionName == "CR_nonISR_OS_ee"   ||
                                RegionInfo[RegionIndex].RegionName == "CR_ISR_OS_ee"      ||
                                RegionInfo[RegionIndex].RegionName == "CR_nonISR_OS_mumu" ||
                                RegionInfo[RegionIndex].RegionName == "CR_ISR_OS_mumu"
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
                            if(docfw
                               &&
                               (
                                RegionInfo[RegionIndex].RegionName == "CR_nonISR_SS_ee" ||
                                RegionInfo[RegionIndex].RegionName == "CR_ISR_SS_ee"    )
                               &&
                               BGGroup[j].info->GroupName == "Zee"
                               )
                            {
                                Cut += "*cfw";
                            }
                            
                            //Cut
                            Cut += "*(1";
                            Cut += CommonCut;
                            Cut += ")";
                            int bgCount = tree2BGMC[j][k]->Draw(temp.Data(),Cut.Data());
                            //bgCount = tree2BGMC[j][k]->GetEntries(("1"+CommonCut).Data());
                            
                            BGGroup[j].info->unweighted += bgCount;
                            
                            /*
                            if(BGGroup[j].info->GroupName == "Zjets" && bgCount>0)
                            {
                                cout<<BGMCSampleInfo[BGMCGroupData[0].lower +k].SampleName.Data()<<endl;
                                cout<<BGMCGroupXS[j][k]<<" "<<BGMCGroupnwAOD[j][k]<<endl;
                                tree2BGMC[j][k]->Scan("weight",("1"+CommonCut).Data());
                            }
                            */
                            
                            //normalization for BG
                            hTemp->Scale(BGMCGroupXS[j][k]/BGMCGroupnwAOD[j][k] *sumDataL);
                            
                            //expN for BGVV
                            if(doVVCount && VarIndex==countVariable && BGGroup[j].info->GroupName == "VV"
                               && RegionGroup[RegionGroupIndex].GroupName == "SR_SS_run1")
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
                        BGGroup[j].h2 = new TH1D(BGGroup[j].info->GroupName.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                        BGGroup[j].h2->GetYaxis()->SetTitle("Number of events");
                        BGGroup[j].h2->SetLineColor(BGGroup[j].info->colour);
                        BGGroup[j].h2->SetFillColor(BGGroup[j].info->colour);
                        
                        BGGroup[j].info->unweighted = 0;
                        for(unsigned int k=0;k<DataSampleName.size();k++)
                        {
                            h2Data[k]->Scale(0);
                            
                            TString temp = Var[VarIndex].VarFormula;
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
                                Cut += " && fLwt==0";
                            }
                            if(BGGroup[j].info->GroupName == "fake lepton")
                            {
                                Cut += " && fLwt!=0";
                            }
                            Cut += ")";
                            
                            int bgCount = 0;
                            if(BGGroup[j].info->GroupName == "charge flip")
                            {
                                bgCount = tree2DataOS[k]->Draw(temp.Data(),Cut.Data());
                            }
                            if(BGGroup[j].info->GroupName == "fake lepton")
                            {
                                bgCount = tree2Data[k]->Draw(temp.Data(),Cut.Data());
                            }
                            BGGroup[j].info->unweighted += bgCount;
                            
                            //Add MCData
                            BGGroup[j].h2->Add(h2Data[k]);
                        }
                    }
                    
                    //delete h2Data
                    for(unsigned int j=0;j<DataSampleName.size();j++)
                    {
                        delete h2Data[j];
                    }
                    
                    //h2SigSum
                    for(unsigned int j=0;j<SigMassSplitting.size();j++)
                    {
                        TString NameTemp = "SigSum_";
                        NameTemp += TString::Format("%.0u",j);
                        h2SigSum[j] = new TH1D(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                        h2SigSum[j]->GetYaxis()->SetTitle("Number of events");
                        h2SigSum[j]->SetLineColor(SigMassSplitting[j].colour);
                        h2SigSum[j]->SetLineStyle(SigMassSplitting[j].linestyle);
                        h2SigSum[j]->SetLineWidth(3);
                    }
                    
                    //h2Sig
                    for(unsigned int j=0;j<SigSampleInfo.size();j++)
                    {
                        TString NameTemp = "Sig_";
                        NameTemp += TString::Itoa(j,10);
                        h2Sig[j] = new TH1D(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                        
                        //Fill Signal
                        TString temp = Var[VarIndex].VarFormula;
                        temp += ">>";
                        temp += NameTemp;
                        
                        TString Cut = "weight*(1";
                        Cut += CommonCut;
                        Cut += ")";
                        SigSampleInfo[j].unweighted = tree2Sig[j]->Draw(temp.Data(),Cut.Data());
                        //SigSampleInfo[j].unweighted = tree2Sig[j]->GetEntries(("1"+CommonCut).Data());
                        
                        for(unsigned int i=0;i<SigMassSplitting.size();i++)
                        {
                            if(j==SigMassSplitting[i].ID) SigMassSplitting[i].unweighted = SigSampleInfo[j].unweighted;
                        }
                        
                        h2Sig[j]->Scale(SigSampleInfo[j].XS *sumDataL/SigSampleInfo[j].nwAOD);
                    }
                }
                
                //Add Signal for the same mass splitting
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    /*
                    double AOD = 0;
                    for(unsigned int j=0;j<SigSampleInfo.size();j++)
                    {
                        if(SigSampleInfo[j].Mass1 - SigSampleInfo[j].Mass2 == SigMassSplitting[i].MassDiff)
                        {
                            h2SigSum[i]->Add(h2Sig[j]);
                            AOD += SigSampleInfo[j].nwAOD;
                        }
                    }
                    
                    //preliminary normalization for h2SigSum
                    h2SigSum[i]->Scale(SigSampleInfo[SigMassSplitting[i].ID].XS *sumDataL/AOD);
                    */
                    
                    ///*
                    h2SigSum[i]->Add(h2Sig[SigMassSplitting[i].ID]);
                    //*/
                }
                
                std::vector<Group> BGGroup2 = BGGroup;
                std::sort(BGGroup2.begin(),BGGroup2.end(),compare2);
                
                const bool DoSignificancePlot = !doOptimize &&
                Var[VarIndex].CutDirection!=0 && (
                RegionGroup[RegionGroupIndex].GroupName == "SR_SS_run1");
                
                unsigned int dataN = 0;
                double total_BG_weighted = 0;
                double total_BG_error2 = 0 ;
                int total_BG_unweighted = 0;
                double averageSignificance = 0;
                const double fake_sys_relative_error = TMath::Sqrt(1.12*1.12 + 0.64*0.64) /1.76;
                double fake_sys_error = 0;
                if(VarIndex==countVariable || RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
                {
                    //expected number of events for Data
                    dataN = h2DataSum->Integral(0,-1);
                    cout<<"Data: "<<dataN<<endl;
                    
                    //expected number of events for background
                    for(unsigned int j=0;j<BGGroup2.size();j++)
                    {
                        //expected number of events for BG
                        BGGroup2[j].info->weighted = BGGroup2[j].h2->IntegralAndError(0,-1,BGGroup2[j].info->error);
                        cout<<BGGroup2[j].info->GroupName.Data()<<": "<<BGGroup2[j].info->weighted<<" +/- "<<BGGroup2[j].info->error<<" ("<<BGGroup2[j].info->unweighted<<")"<<endl;
                        
                        if(BGGroup2[j].info->weighted > 0)
                        {
                            total_BG_weighted += BGGroup2[j].info->weighted;
                        }
                        total_BG_error2 += BGGroup2[j].info->error*BGGroup2[j].info->error;
                        
                        total_BG_unweighted += BGGroup2[j].info->unweighted;
                        
                        if(RegionGroup[RegionGroupIndex].GroupName == "VR_Fake" && BGGroup2[j].info->GroupName == "Samuel_Fake")
                        {
                            fake_sys_error = BGGroup2[j].info->weighted * fake_sys_relative_error;
                        }
                    }
                    cout<<"Total BG: "<<total_BG_weighted<<" +/- "<<TMath::Sqrt(total_BG_error2)<<" ("<<total_BG_unweighted<<")"<<endl<<endl;
                    
                    //expected number of events for signal
                    for(unsigned int i=0;i<SigSampleInfo.size();i++)
                    {
                        //expected number of events
                        SigSampleInfo[i].weighted = h2Sig[i]->IntegralAndError(0,-1,SigSampleInfo[i].error);
                        
                        //Significance
                        SigSampleInfo[i].significance = GetSignificance(SigSampleInfo[i].weighted,total_BG_weighted,total_BG_error2);
                        
                        for(unsigned int j=0;j<SigMassSplitting.size();j++)
                        {
                            if(i == SigMassSplitting[j].ID)
                            {
                                SigMassSplitting[j].weighted = SigSampleInfo[i].weighted;
                                SigMassSplitting[j].error = SigSampleInfo[i].error;
                                SigMassSplitting[j].significance = SigSampleInfo[i].significance;
                                
                                cout<<SigMassSplitting[j].IDName.Data()<<": "<<SigMassSplitting[j].weighted<<" +/- "<<SigMassSplitting[j].error<<" ("<<SigMassSplitting[j].unweighted<<")"<<", Significance: "<<SigMassSplitting[j].significance<<endl;
                            }
                        }
                        
                        for(unsigned int j=0;j<OptimizingSignal.size();j++)
                        {
                            if(i == OptimizingSignal[j])
                            {
                                averageSignificance += SigSampleInfo[i].significance;
                            }
                        }
                        
                        //combined significance
                        if(VarIndex==countVariable && SigSampleInfo[i].significance>0)
                        {
                            SigSampleInfo[i].significance2 += SigSampleInfo[i].significance * SigSampleInfo[i].significance;
                        }
                    }
                    cout<<endl;
                    
                    averageSignificance /= OptimizingSignal.size();
                }
                
                //Asymmetric errors for Fake VR
                TGraphAsymmErrors* g_mc = nullptr;
                TGraphAsymmErrors* g_ratio = nullptr;
                if(RegionGroup[RegionGroupIndex].GroupName == "VR_Fake")
                {
                    int nbin = Var[VarIndex].bin;
                    g_mc = new TGraphAsymmErrors(nbin);
                    g_ratio = new TGraphAsymmErrors(nbin);
                    
                    for (int j = 0; j < nbin; j++)
                    {
                        double x = BGGroup2[0].h2->GetBinCenter(j+1);
                        double xerr_down = BGGroup2[0].h2->GetBinWidth(j+1) /2;
                        double xerr_up = BGGroup2[0].h2->GetBinWidth(j+1) /2;
                        
                        double y = 0;
                        double error2 = 0;
                        double fake_yield = 0;
                        for(unsigned int k = 0; k < BGGroup2.size(); k++)
                        {
                            double yield = 0;
                            double error = 0;
                            yield = BGGroup2[k].h2->IntegralAndError(j+1,j+1,error);
                            if(yield > 0)
                            {
                                y += yield;
                            }
                            error2 += error * error;
                            
                            if(BGGroup2[k].info->GroupName == "Samuel_Fake") fake_yield = yield;
                        }
                        
                        double fake_error = fake_yield * fake_sys_relative_error;
                        error2 += fake_error * fake_error;
                        
                        double total_error = TMath::Sqrt(error2);
                        
                        double yerr_down = total_error;
                        double yerr_up = total_error;
                        
                        g_mc->SetPoint(j,x,y);
                        g_mc->SetPointError(j,xerr_down,xerr_up,yerr_down,yerr_up);
                        //cout<<xerr_down<<", "<<x<<", "<<xerr_up<<"; "<<yerr_down<<", "<<y<<", "<<yerr_up<<endl;
                        
                        g_ratio->SetPoint(j,x,1);
                        yerr_down = total_error/y;
                        yerr_up = total_error/y;
                        g_ratio->SetPointError(j,xerr_down,xerr_up,yerr_down,yerr_up);
                    }
                    
                    g_mc->SetFillColor(kBlack);
                    g_mc->SetFillStyle(3005);
                    
                    g_ratio->SetFillColor(kYellow);
                }
                
                if(VarIndex==countVariable)
                {
                    //output file
                    TString PathName = "latex/data/expN/";
                    PathName += RegionInfo[RegionIndex].RegionName;
                    if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
                    {
                        PathName += "_";
                        PathName += TString::Itoa(SigOptimizingIndex,10);
                    }
                    PathName += ".tex";
                    
                    ofstream fout;
                    fout.open(PathName.Data());
                    fout<<setprecision(6)<<std::fixed;
                    
                    for(unsigned int j=BGGroup2.size()-1;;j--)
                    {
                        fout<<BGGroup2[j].info->LatexName.Data();
                        fout<<" & $";
                        fout<<BGGroup2[j].info->weighted;
                        fout<<"\\pm";
                        fout<<BGGroup2[j].info->error;
                        
                        if(RegionGroup[RegionGroupIndex].GroupName == "VR_Fake" && BGGroup2[j].info->GroupName == "Samuel_Fake")
                        {
                            fout<<"\\pm";
                            fout<<fake_sys_error;
                        }
                        
                        fout<<"$ (";
                        fout<<BGGroup2[j].info->unweighted;
                        fout<<") & \\\\"<<endl<<"\\hline"<<endl;
                        if(j==0) break;
                    }
                    
                    fout<<"Total BG & $";
                    fout<<total_BG_weighted;
                    fout<<"\\pm";
                    fout<<TMath::Sqrt(total_BG_error2);
                    
                    if(RegionGroup[RegionGroupIndex].GroupName == "VR_Fake")
                    {
                        fout<<"\\pm";
                        fout<<fake_sys_error;
                    }
                    
                    fout<<"$ (";
                    fout<<total_BG_unweighted;
                    fout<<") & \\\\"<<endl<<"\\hline"<<endl;
                    
                    if(RegionGroup[RegionGroupIndex].showData)
                    {
                        fout<<"Data & $";
                        fout<<dataN;
                        fout<<"$ & \\\\"<<endl<<"\\hline"<<endl;
                    }
                    
                    for(unsigned int i=0;i<SigMassSplitting.size();i++)
                    {
                        fout<<SigMassSplitting[i].IDName.Data();
                        fout<<" & $";
                        fout<<setprecision(2)<<std::fixed;
                        fout<<SigMassSplitting[i].weighted;
                        fout<<"\\pm";
                        fout<<SigMassSplitting[i].error;
                        fout<<"$ (";
                        fout<<SigMassSplitting[i].unweighted;
                        fout<<") &";
                        if(RegionGroup[RegionGroupIndex].showSignificance)
                        {
                            fout<<"$";
                            fout<<setprecision(2)<<std::fixed;
                            fout<<SigMassSplitting[i].significance;
                            fout<<"$";
                        }
                        fout<<"\\\\"<<endl<<"\\hline"<<endl;
                    }
                    fout.close();
                    
                    //output file for VV
                    if(doVVCount
                       && RegionGroup[RegionGroupIndex].GroupName == "SR_SS_run1")
                    {
                        PathName = "latex/data/expN/";
                        PathName += RegionInfo[RegionIndex].RegionName;
                        PathName += "_BGVV.tex";
                        
                        fout.open(PathName.Data());
                        fout<<setprecision(3)<<std::fixed;
                        for(unsigned int j=0;j<BGVVData.size();j++)
                        {
                            TString latexName = BGVVData[j].SampleName;
                            latexName.ReplaceAll("_","\\_");

                            fout<<latexName.Data();
                            fout<<" & $";
                            fout<<sumOfEventVV[j][0];
                            fout<<"\\pm";
                            fout<<sumOfEventVV[j][1];
                            fout<<"$ \\\\"<<endl<<"\\hline"<<endl;
                        }
                        fout.close();
                    }
                    
                    //h2SR
                    if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
                    {
                        const unsigned int RegionIndex2 = RegionIndex-RegionGroup[RegionGroupIndex].lower;
                        for(unsigned int j=0;j<BGGroup.size();j++)
                        {
                            SRBGGroup[j].h2->SetBinContent(RegionIndex2+1,BGGroup[j].info->weighted);
                        }
                        for(unsigned int j=0;j<SigMassSplitting.size();j++)
                        {
                            h2SRSig[j]->SetBinContent(RegionIndex2+1,SigMassSplitting[j].weighted);
                            h2SRSignificance[j]->SetBinContent(RegionIndex2+1,SigMassSplitting[j].significance);
                        }
                    }
                    
                    //significance for all mass point
                    if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_run1" ||
                       RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt" ||
                       RegionGroup[RegionGroupIndex].GroupName == "SR_SS_pre" )
                    {
                        PathName = "latex/data/optimization/significance/";
                        PathName += RegionInfo[RegionIndex].RegionName;
                        PathName += "_";
                        PathName += TString::Itoa(SigOptimizingIndex,10);
                        PathName += ".txt";
                        fout.open(PathName.Data());
                        
                        fout<<"mC1:mN1:nSig:nSigErr:Zn"<<endl;
                        fout<<"#Tag: "<<RegionInfo[RegionIndex].RegionName.Data()<<endl;
                        for(unsigned int j=0;j<BGGroup2.size();j++)
                        {
                            //expected number of events for BG
                            fout<<setprecision(6)<<std::fixed;
                            fout<<"#"<<BGGroup2[j].info->GroupName.Data()<<": "<<BGGroup2[j].info->weighted<<" +/- "<<BGGroup2[j].info->error<<" ("<<BGGroup2[j].info->unweighted<<")"<<endl;
                        }
                        fout<<"#Total BG: "<<total_BG_weighted<<" +/- "<<TMath::Sqrt(total_BG_error2)<<" ("<<total_BG_unweighted<<")"<<endl<<endl;
                        
                        TGraph2D* g_nSig = new TGraph2D();
                        TGraph2D* g_significance = new TGraph2D();
                        
                        for(unsigned int j=0;j<SigSampleInfo.size();j++)
                        {
                            fout<<setprecision(1)<<std::fixed;
                            fout<<SigSampleInfo[j].Mass1<<" "<<SigSampleInfo[j].Mass2<<" ";
                            fout<<setprecision(3)<<std::fixed;
                            fout<<SigSampleInfo[j].weighted<<" "<<SigSampleInfo[j].error<<" ";
                            fout<<SigSampleInfo[j].significance<<endl;
                            
                            //2D plot
                            //if(SigSampleInfo[j].significance>0)
                            {
                                g_nSig->SetPoint(g_significance->GetN(), SigSampleInfo[j].Mass1, SigSampleInfo[j].Mass2, SigSampleInfo[j].weighted);
                                g_significance->SetPoint(g_significance->GetN(), SigSampleInfo[j].Mass1, SigSampleInfo[j].Mass2, SigSampleInfo[j].significance);
                            }
                        }
                        
                        cout<<"average significance: "<<averageSignificance<<endl;
                        fout<<"#average significance: "<<averageSignificance<<endl;
                        
                        fout.close();
                        
                        if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
                        {
                            gStyle->SetPalette(1);
                            
                            TH2D* h2 = g_nSig->GetHistogram();
                            h2->GetXaxis()->SetTitle("m_{#tilde{#chi}^{#pm}_{1}/#tilde{#chi}^{0}_{2}} [GeV]");
                            h2->GetYaxis()->SetTitle("m_{#tilde{#chi}^{0}_{1}} [GeV]");
                            h2->GetZaxis()->SetTitle("nSig");
                            
                            h2->GetXaxis()->SetRangeUser(0,600);
                            h2->GetYaxis()->SetRangeUser(0,120);
                            
                            TCanvas* c2 = new TCanvas();
                            c2->cd();
                            c2->SetRightMargin(0.16);
                            h2->Draw();
                            g_nSig->Draw("colz");
                            
                            {
                                TLatex lt1;
                                lt1.DrawLatexNDC(0.05,0.05,RegionInfo[RegionIndex].RegionName.Data());
                                lt1.SetTextSize(lt1.GetTextSize());
                            }
                            
                            {
                                TLatex lt1;
                                lt1.DrawLatexNDC(0.4,0.05,setupTag.Data());
                                lt1.SetTextSize(lt1.GetTextSize());
                            }
                            
                            TString NameTemp = "plot/";
                            NameTemp += "nSig_";
                            NameTemp += RegionInfo[RegionIndex].RegionName;
                            NameTemp += "_";
                            NameTemp += TString::Itoa(SigOptimizingIndex,10);
                            NameTemp += ".eps";
                            c2->Print(NameTemp,"eps");
                            
                            delete c2;
                        }
                        
                        if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
                        {
                            gStyle->SetPalette(1);
                            
                            TH2D* h2 = g_significance->GetHistogram();
                            h2->GetXaxis()->SetTitle("m_{#tilde{#chi}^{#pm}_{1}/#tilde{#chi}^{0}_{2}} [GeV]");
                            h2->GetYaxis()->SetTitle("m_{#tilde{#chi}^{0}_{1}} [GeV]");
                            h2->GetZaxis()->SetTitle("Z_{n}");
                            
                            h2->GetXaxis()->SetRangeUser(0,600);
                            h2->GetYaxis()->SetRangeUser(0,120);
                            //h2->GetZaxis()->SetRangeUser(0,4);
                            
                            TCanvas* c2 = new TCanvas();
                            c2->cd();
                            c2->SetRightMargin(0.16);
                            h2->Draw();
                            g_significance->Draw("colz");
                            
                            {
                                TLatex lt1;
                                lt1.DrawLatexNDC(0.05,0.05,RegionInfo[RegionIndex].RegionName.Data());
                                lt1.SetTextSize(lt1.GetTextSize());
                            }
                            
                            {
                                TLatex lt1;
                                lt1.DrawLatexNDC(0.4,0.05,setupTag.Data());
                                lt1.SetTextSize(lt1.GetTextSize());
                            }
                            
                            TString NameTemp = "plot/";
                            NameTemp += "significance_";
                            NameTemp += RegionInfo[RegionIndex].RegionName;
                            NameTemp += "_";
                            NameTemp += TString::Itoa(SigOptimizingIndex,10);
                            NameTemp += ".eps";
                            c2->Print(NameTemp,"eps");
                            
                            delete c2;
                        }
                        
                        delete g_nSig;
                        delete g_significance;
                    }
                    
                    //calculate scale factor for fake BG
                    if(RegionInfo[RegionIndex].RegionName == "CR_nonISR_SS_mumu_1B" ||
                       RegionInfo[RegionIndex].RegionName == "CR_ISR_SS_mumu_1B"    ||
                       RegionGroup[RegionGroupIndex].GroupName == "CR_SS_mumu_low_mT2")
                    {
                        for(unsigned int j=0;j<BGGroup.size();j++)
                        {
                            if(BGGroup[j].info->GroupName == "fake lepton")
                            {
                                double factor = ( dataN - total_BG_weighted )/ BGGroup[j].info->weighted +1;
                                cout<<RegionInfo[RegionIndex].RegionName.Data()<<": scale factor: "<<factor<<endl;
                            }
                        }
                    }
                }
                
                //significance calculation
                TH1D* hSignificance[SigMassSplitting.size()];
                if(DoSignificancePlot)
                {
                    for(unsigned int i=0;i<SigMassSplitting.size();i++)
                    {
                        TString NameTemp = "Significance_";
                        NameTemp += TString::Format("%.0u",i);
                        hSignificance[i] = new TH1D(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                        hSignificance[i]->SetLineColor(SigMassSplitting[i].colour);
                        hSignificance[i]->SetLineWidth(2);
                    }
                    
                    for(int bin=1;bin<=Var[VarIndex].bin;bin++)
                    {
                        double nBG = 0;
                        double nSig = 0;
                        for(unsigned int j=0;j<BGGroup.size();j++)
                        {
                            if(Var[VarIndex].CutDirection == 1)
                            {
                                double n = BGGroup[j].h2->Integral(bin,-1);
                                if(n<0) n=0;
                                nBG += n;
                            }
                            else if(Var[VarIndex].CutDirection == -1)
                            {
                                double n = BGGroup[j].h2->Integral(1,bin);
                                if(n<0) n=0;
                                nBG += n;
                            }
                        }
                        
                        //expected number of events for signal
                        for(unsigned int i=0;i<SigMassSplitting.size();i++)
                        {
                            //expected number of events
                            if(Var[VarIndex].CutDirection == 1)
                            {
                                nSig = h2SigSum[i]->Integral(bin,-1);
                            }
                            else if(Var[VarIndex].CutDirection == -1)
                            {
                                nSig = h2SigSum[i]->Integral(0,bin);
                            }
                            
                            //Significance
                            //if(SigMassSplitting[i].MassDiff==100) cout<<bin<<": "<<nBG<<", "<<nSig<<", "<<GetSignificance(nSig,nBG)<<endl;
                            if(nBG>0) hSignificance[i]->SetBinContent(bin,GetSignificance(nSig,nBG));
                        }
                    }
                }
                
                //final normalization for h2SigSum for plotting
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    h2SigSum[i]->Scale(SigMassSplitting[i].scale);
                }
                
                //scale Z+jets by 1.4 for CR
                if(RegionGroup[RegionGroupIndex].GroupName == "CR_OS_1B"  ||
                   RegionGroup[RegionGroupIndex].GroupName == "CR_OS_2B"  )
                {
                    for(unsigned int i=0;i<tree2BGMC.size();i++)
                    {
                        if(
                           BGGroup[i].info->GroupName == "Zee"    ||
                           BGGroup[i].info->GroupName == "Zmumu"  ||
                           BGGroup[i].info->GroupName == "Ztautau"
                           )
                        {
                            BGGroup[i].h2->Scale(1.4);
                        }
                    }
                }
                
                //apply scale factor for fake BG
                for(unsigned int j=0;j<BGGroup.size();j++)
                {
                    if(BGGroup[j].info->GroupName == "fake lepton")
                    {
                        if(RegionInfo[RegionIndex].RegionName == "CR_nonISR_SS_mumu_1B" ||
                           RegionInfo[RegionIndex].RegionName == "CR_nonISR_SS_mumu_2B" ||
                           RegionInfo[RegionIndex].RegionName == "SR_nonISR_SS_mumu_0B" )
                        {
                            BGGroup[j].h2->Scale(4.28729);
                        }
                        
                        if(RegionInfo[RegionIndex].RegionName == "CR_ISR_SS_mumu_1B" ||
                           RegionInfo[RegionIndex].RegionName == "CR_ISR_SS_mumu_2B" ||
                           RegionInfo[RegionIndex].RegionName == "SR_ISR_SS_mumu_0B" )
                        {
                            BGGroup[j].h2->Scale(4.18533);
                        }
                        
                        if(RegionInfo[RegionIndex].RegionName == "CR_nonISR_SS_mumu_low_mT2" ||
                           RegionInfo[RegionIndex].RegionName == "SR_nonISR_SS_mumu"         )
                        {
                            BGGroup[j].h2->Scale(4.43166);
                        }
                        
                        if(RegionInfo[RegionIndex].RegionName == "CR_ISR_SS_mumu_low_mT2" ||
                           RegionInfo[RegionIndex].RegionName == "SR_ISR_SS_mumu"         )
                        {
                            BGGroup[j].h2->Scale(3.85965);
                        }
                    }
                }
                
                //Add BG
                TH1D* h2BGSum = new TH1D("BGSum",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                for(unsigned int j=0;j<BGGroup.size();j++)
                {
                    h2BGSum->Add(BGGroup[j].h2);
                }
                
                //stack
                THStack stack;
                for(unsigned int j=0;j<BGGroup2.size();j++)
                {
                    stack.Add(BGGroup2[j].h2);
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
                    
                    if(RegionGroup[RegionGroupIndex].GroupName == "VR_Fake")
                    {
                        min = h2BGSum->GetBinContent(h2BGSum->GetMinimumBin());
                        max = h2BGSum->GetBinContent(h2BGSum->GetMaximumBin());
                    }
                    
                    if(Var[VarIndex].log)
                    {
                        //h2DataSum->SetMaximum(max*2);
                        h2DataSum->SetMaximum(max*100);
                        
                        //if(doOptimize) Var[VarIndex].ymin /= 100;
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
                    if(RegionGroup[RegionGroupIndex].showData)
                    {
                        xl2=0.92;
                        yl2=0.95;
                        xl1=xl2-0.1;
                        yl1=yl2-0.3;
                    }
                    else
                    {
                        xl2=0.75;
                        yl2=0.95;
                        xl1=xl2-0.3;
                        yl1=yl2-0.2;
                    }
                    leg = new TLegend(xl1,yl1,xl2,yl2);
                    leg->SetNColumns(1);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(42);
                    leg->SetBorderSize(0);
                    if(RegionGroup[RegionGroupIndex].showData)
                    {
                        leg->AddEntry(h2DataSum,"Data","p");
                    }
                    for(unsigned int j=BGGroup2.size()-1;;j--)
                    {
                        leg->AddEntry(BGGroup2[j].h2,BGGroup2[j].info->LegendName.Data(),"fl");
                        if(j==0) break;
                    }
                    
                    for(unsigned int i=0;i<SigMassSplitting.size();i++)
                    {
                        //TString NameTemp = "C1N2#rightarrow WZ";
                        //TString NameTemp = "C1N2#rightarrow slep";
                        TString NameTemp = "";
                        NameTemp += SigMassSplitting[i].IDName;
                        NameTemp += " x";
                        NameTemp += TString::Itoa(SigMassSplitting[i].scale,10);
                        leg->AddEntry(h2SigSum[i],NameTemp.Data(),"l");
                    }
                }
                
                TH1D* h2Ratio = nullptr;
                TPad* pad1 = nullptr;
                TPad* pad2 = nullptr;
                TH1D* hPad2 = nullptr;
                if(RegionGroup[RegionGroupIndex].showData || DoSignificancePlot )
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
                    
                    if(RegionGroup[RegionGroupIndex].showData)
                    {
                        //ratio plot
                        h2Ratio = new TH1D("Ratio",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                        h2Ratio->Sumw2();
                        h2Ratio->Add(h2DataSum);
                        h2Ratio->Divide(h2BGSum);
                        
                        h2Ratio->SetMarkerSize(1.0);
                        h2Ratio->SetLineColor(1);
                        h2Ratio->GetYaxis()->SetTitle("Data/MC");
                        h2Ratio->GetYaxis()->SetRangeUser(0.5,1.5);
                        h2Ratio->GetYaxis()->SetNdivisions(8);
                    }
                    else if(RegionGroup[RegionGroupIndex].showSignificance)
                    {
                        hSignificance[0]->GetYaxis()->SetTitle("Significance");
                        hSignificance[0]->GetYaxis()->SetRangeUser(-1,3);
                        hSignificance[0]->GetYaxis()->SetNdivisions(6);
                    }
                    
                    if(RegionGroup[RegionGroupIndex].showData) hPad2 = h2Ratio;
                    else if(RegionGroup[RegionGroupIndex].showSignificance) hPad2 = hSignificance[0];
                    
                    hPad2->GetXaxis()->SetLabelSize(x_label*scale2);
                    hPad2->GetYaxis()->SetLabelSize(y_label*scale2);
                    
                    hPad2->GetXaxis()->SetTitleSize(x_title*scale2);
                    hPad2->GetYaxis()->SetTitleSize(y_title*scale2);
                    hPad2->GetYaxis()->SetTitleOffset(y_offset/scale2);
                    hPad2->GetXaxis()->SetTitle(xaxis.Data());
                    hPad2->GetYaxis()->CenterTitle();
                    
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
                stack.Draw("histsame"); //Cutflow Attention
                if(RegionGroup[RegionGroupIndex].GroupName == "VR_Fake") g_mc->Draw("e2, SAME");
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    h2SigSum[i]->Draw("histsame");
                }
                if(RegionGroup[RegionGroupIndex].showData)
                {
                    h2DataSum->Draw("Esame");
                }
                h2DataSum->Draw("sameaxis");
                leg->Draw();
                
                {
                    //text
                    ATLASLabel(0.6,0.88,"Internal");
                    
                    TLatex lt2;
                    TString NameTemp = "#sqrt{#it{s}} = 13 TeV, ";
                    NameTemp += TString::Format("%.1f",sumDataL/1000);
                    NameTemp += " fb^{-1}";
                    lt2.DrawLatexNDC(0.6,0.83, NameTemp.Data());
                    lt2.SetTextSize(lt2.GetTextSize());
                    
                    TLatex lt1;
                    lt1.DrawLatexNDC(0.6,0.78,RegionInfo[RegionIndex].RegionName.Data());
                    lt1.SetTextSize(lt1.GetTextSize());
                }
                
                //Draw vertical line for optimized cut
                if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt" ||
                   RegionGroup[RegionGroupIndex].GroupName == "SR_SS_pre" )
                {
                    TString optGroupName = "SR_SS_opt";
                    unsigned int optRegionGroupIndex = 0;
                    for(unsigned int i=0;i<RegionGroup.size();i++)
                    {
                        if(RegionGroup[i].GroupName == optGroupName)
                        {
                            optRegionGroupIndex = i;
                        }
                    }
                    
                    const unsigned int SixChannel = RegionIndex - RegionGroup[RegionGroupIndex].lower;
                    
                    TString PathName = "latex/data/optimization/";
                    if(useDani) PathName += "cut_txt_Dani/";
                    else PathName += "cut_txt/";
                    PathName += "cut_";
                    PathName += RegionInfo[ RegionGroup[optRegionGroupIndex].lower+SixChannel ].RegionName;
                    PathName += "_";
                    PathName += TString::Itoa(SigOptimizingIndex,10);
                    PathName += ".txt";
                    
                    ifstream fin;
                    fin.open(PathName.Data());
                    
                    for(unsigned int i=0;i<RegionInfo[ RegionGroup[optRegionGroupIndex].lower+SixChannel ].OptimizingCut[SigOptimizingIndex].size();i++)
                    {
                        for(unsigned int j=0;j<RegionInfo[ RegionGroup[optRegionGroupIndex].lower+SixChannel ].OptimizingCut[SigOptimizingIndex][i].size();j++)
                        {
                            TString VarName;
                            double lower;
                            double upper;
                            fin>>VarName;
                            if(VarName == Var[VarIndex].VarName)
                            {
                                fin>>lower;
                                fin>>upper;
                                
                                TString NameTemp;
                                TLine l;
                                TLine* l1 = nullptr;
                                TLine* l2 = nullptr;
                                if(lower == 0)
                                {
                                    if(upper == -1)
                                    {
                                    }
                                    else
                                    {
                                        NameTemp += Var[VarIndex].VarTitle;
                                        NameTemp += " < ";
                                        NameTemp += upper;
                                        
                                        l1 = l.DrawLine(upper, h2DataSum->GetMinimum(), upper, h2DataSum->GetMaximum());
                                        l1->SetLineStyle(9);
                                        l1->SetLineWidth(2);
                                    }
                                }
                                else
                                {
                                    if(upper == -1)
                                    {
                                        NameTemp += Var[VarIndex].VarTitle;
                                        NameTemp += " >= ";
                                        NameTemp += lower;
                                        
                                        l1 = l.DrawLine(lower, h2DataSum->GetMinimum(), lower, h2DataSum->GetMaximum());
                                        l1->SetLineStyle(9);
                                        l1->SetLineWidth(2);
                                    }
                                    else
                                    {
                                        NameTemp += lower;
                                        NameTemp += " <= ";
                                        NameTemp += Var[VarIndex].VarTitle;
                                        NameTemp += " < ";
                                        NameTemp += upper;
                                        
                                        l1 = l.DrawLine(lower, h2DataSum->GetMinimum(), lower, h2DataSum->GetMaximum());
                                        l1->SetLineStyle(9);
                                        l1->SetLineWidth(2);
                                        
                                        l2 = l.DrawLine(upper, h2DataSum->GetMinimum(), upper, h2DataSum->GetMaximum());
                                        l2->SetLineStyle(9);
                                        l2->SetLineWidth(2);
                                    }
                                }
                                
                                TLatex lt2;
                                lt2.DrawLatexNDC(0.6,0.63, NameTemp.Data());
                                lt2.SetTextSize(lt2.GetTextSize());
                                
                                break;
                            }
                            else
                            {
                                fin>>lower;
                                fin>>upper;
                            }
                        }
                    }
                    
                    fin.close();
                    
                    if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
                    {
                        TString NameTemp = TString::Format("Z_{n} = %.2f",averageSignificance);
                        TLatex lt2;
                        lt2.DrawLatexNDC(0.6,0.68, NameTemp.Data());
                        lt2.SetTextSize(lt2.GetTextSize());
                    }
                }
                
                {
                    TLatex lt1;
                    lt1.DrawLatexNDC(0.6,0.73,setupTag.Data());
                    lt1.SetTextSize(lt1.GetTextSize());
                }
                
                if(RegionGroup[RegionGroupIndex].showData || DoSignificancePlot)
                {
                    //Draw for pad2
                    c2->cd();
                    pad2->Draw();
                    pad2->cd();
                    
                    if(RegionGroup[RegionGroupIndex].showData)
                    {
                        h2Ratio->Draw();
                        
                        if(RegionGroup[RegionGroupIndex].GroupName == "VR_Fake") g_ratio->Draw("e2, SAME");
                        
                        TLine l;
                        TLine* l1 = l.DrawLine(h2Ratio->GetXaxis()->GetXmin(), 1., h2Ratio->GetXaxis()->GetXmax(), 1.);
                        l1->SetLineStyle(2);
                        l1->SetLineWidth(2);
                        
                        h2Ratio->Draw("same");
                    }
                    else if(RegionGroup[RegionGroupIndex].showSignificance)
                    {
                        hSignificance[0]->Draw();
                        for(unsigned int i=1;i<SigMassSplitting.size();i++)
                        {
                            hSignificance[i]->Draw("same");
                        }
                    }
                }
                
                {
                    //export histograms in eps format
                    TString NameTemp = "plot/";
                    NameTemp += Var[VarIndex].VarName;
                    NameTemp += "_";
                    NameTemp += RegionInfo[RegionIndex].RegionName;
                    if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
                    {
                        NameTemp += "_";
                        NameTemp += TString::Itoa(SigOptimizingIndex,10);
                    }
                    
                    if(RegionGroup[RegionGroupIndex].GroupName == "VR_Fake")
                    {
                        NameTemp += ".pdf";
                        c2->Print(NameTemp,"pdf");
                    }
                    else
                    {
                        NameTemp += ".eps";
                        c2->Print(NameTemp,"eps");
                    }
                    
                    NameTemp = Var[VarIndex].VarName;
                    NameTemp += "_";
                    NameTemp += RegionInfo[RegionIndex].RegionName;
                    fout_plot->cd();
                    c2->Write(NameTemp);
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
                for(unsigned int j=0;j<SigSampleInfo.size();j++)
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
                
                if(RegionGroup[RegionGroupIndex].showData)
                {
                    //h2Ratio
                    delete h2Ratio;
                    
                    //pad2
                    delete pad2;
                }
                
                if(DoSignificancePlot)
                {
                    for(unsigned int i=0;i<SigMassSplitting.size();i++)
                    {
                        delete hSignificance[i];
                    }
                }
                
                //Asymmetric errors
                if(RegionGroup[RegionGroupIndex].GroupName == "VR_Fake")
                {
                    delete g_mc;
                    delete g_ratio;
                }
                
                //pad1
                delete pad1;
                
                //canvas
                delete c2;
            }
            
            //2D significance plot
            if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_0B")
            {
                unsigned int VarIndex[2];
                VarIndex[0] = 6; // ptll
                VarIndex[1] = 7; // MET
                
                TString title = Var[VarIndex[0]].VarTitle;
                title += ", ";
                title += Var[VarIndex[1]].VarTitle;
                
                TString axis[2];
                for(unsigned int i=0;i<=1;i++)
                {
                    axis[i] = Var[VarIndex[i]].VarTitle;
                    axis[i] += " ";
                    axis[i] += Var[VarIndex[i]].unit;
                }
                
                TH2D* h2SigSum[SigMassSplitting.size()];
                {
                    TString CommonCut = RegionInfo[RegionIndex].Cut;
                    //Background
                    //For MC background
                    for(unsigned int j=0;j<tree2BGMC.size();j++)
                    {
                        BGGroup[j].h3 = new TH2D(BGGroup[j].info->GroupName.Data(),title.Data(),Var[VarIndex[0]].bin,Var[VarIndex[0]].xmin,Var[VarIndex[0]].xmax,
                                                                                                Var[VarIndex[1]].bin,Var[VarIndex[1]].xmin,Var[VarIndex[1]].xmax);
                        for(unsigned int k=0;k<tree2BGMC[j].size();k++)
                        {
                            TH2D* hTemp = new TH2D("BGMC",title.Data(),Var[VarIndex[0]].bin,Var[VarIndex[0]].xmin,Var[VarIndex[0]].xmax,
                                                                       Var[VarIndex[1]].bin,Var[VarIndex[1]].xmin,Var[VarIndex[1]].xmax);
                            //fill histograms from trees
                            TString temp = Var[VarIndex[1]].VarFormula;
                            temp += ":";
                            temp += Var[VarIndex[0]].VarFormula;
                            temp += ">>BGMC";
                            
                            //Weight
                            TString Cut = "weight";
                            
                            //Cut
                            Cut += "*(1";
                            Cut += CommonCut;
                            Cut += ")";
                            tree2BGMC[j][k]->Draw(temp.Data(),Cut.Data());
                            
                            //normalization for BG
                            hTemp->Scale(BGMCGroupXS[j][k]/BGMCGroupnwAOD[j][k] *sumDataL);
                            
                            BGGroup[j].h3->Add(hTemp);
                            delete hTemp;
                        }
                    }
                    
                    //For data-driven background
                    for(unsigned int j=tree2BGMC.size();j<BGGroup.size();j++)
                    {
                        BGGroup[j].h3 = new TH2D(BGGroup[j].info->GroupName.Data(),title.Data(),Var[VarIndex[0]].bin,Var[VarIndex[0]].xmin,Var[VarIndex[0]].xmax,
                                                                                                Var[VarIndex[1]].bin,Var[VarIndex[1]].xmin,Var[VarIndex[1]].xmax);
                        
                        for(unsigned int k=0;k<DataSampleName.size();k++)
                        {
                            TH2D* hTemp = new TH2D("BGData",title.Data(),Var[VarIndex[0]].bin,Var[VarIndex[0]].xmin,Var[VarIndex[0]].xmax,
                                                                         Var[VarIndex[1]].bin,Var[VarIndex[1]].xmin,Var[VarIndex[1]].xmax);
                            
                            //fill histograms from trees
                            TString temp = Var[VarIndex[1]].VarFormula;
                            temp += ":";
                            temp += Var[VarIndex[0]].VarFormula;
                            temp += ">>BGData";
                            
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
                                Cut += " && fLwt==0";
                            }
                            if(BGGroup[j].info->GroupName == "fake lepton")
                            {
                                Cut += " && fLwt!=0";
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
                            BGGroup[j].h3->Add(hTemp);
                            delete hTemp;
                        }
                    }
                    
                    //For signal
                    for(unsigned int i=0;i<SigMassSplitting.size();i++)
                    {
                        TString NameTemp = "SigSum_";
                        NameTemp += TString::Format("%.0u",i);
                        h2SigSum[i] = new TH2D(NameTemp.Data(),title.Data(),Var[VarIndex[0]].bin,Var[VarIndex[0]].xmin,Var[VarIndex[0]].xmax,
                                                                            Var[VarIndex[1]].bin,Var[VarIndex[1]].xmin,Var[VarIndex[1]].xmax);
                        unsigned int AOD = 0;
                        for(unsigned int j=0;j<SigSampleInfo.size();j++)
                        {
                            if(SigSampleInfo[j].Mass1 - SigSampleInfo[j].Mass2 != SigMassSplitting[i].MassDiff) continue;
                            TH2D* hTemp = new TH2D("signal",title.Data(),Var[VarIndex[0]].bin,Var[VarIndex[0]].xmin,Var[VarIndex[0]].xmax,
                                                                         Var[VarIndex[1]].bin,Var[VarIndex[1]].xmin,Var[VarIndex[1]].xmax);
                            //fill histograms from trees
                            TString temp = Var[VarIndex[1]].VarFormula;
                            temp += ":";
                            temp += Var[VarIndex[0]].VarFormula;
                            temp += ">>signal";
                            
                            TString Cut = "weight*(1";
                            Cut += CommonCut;
                            Cut += ")";
                            tree2Sig[j]->Draw(temp.Data(),Cut.Data());
                            
                            h2SigSum[i]->Add(hTemp);
                            delete hTemp;
                            
                            AOD += SigSampleInfo[j].nwAOD;
                        }
                        
                        //normalization for h2SigSum
                        h2SigSum[i]->Scale(SigSampleInfo[SigMassSplitting[i].ID].XS *sumDataL/AOD);
                    }
                }
                
                TH2D* hSignificance[SigMassSplitting.size()];
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    TString NameTemp = "Significance_";
                    NameTemp += TString::Format("%.0u",i);
                    hSignificance[i] = new TH2D(NameTemp.Data(),title.Data(),Var[VarIndex[0]].bin,Var[VarIndex[0]].xmin,Var[VarIndex[0]].xmax,
                                                                             Var[VarIndex[1]].bin,Var[VarIndex[1]].xmin,Var[VarIndex[1]].xmax);
                    hSignificance[i]->GetXaxis()->SetTitle(axis[0].Data());
                    hSignificance[i]->GetYaxis()->SetTitle(axis[1].Data());
                }
                
                for(int bin1=1;bin1<=Var[VarIndex[0]].bin;bin1++)
                {
                    for(int bin2=1;bin2<=Var[VarIndex[1]].bin;bin2++)
                    {
                        double nBG = 0;
                        double nSig = 0;
                        for(unsigned int j=0;j<BGGroup.size();j++)
                        {
                            nBG += BGGroup[j].h3->Integral(bin1,-1,bin2,-1);
                            //nBG += BGGroup[j].h3->Integral(bin1,-1,0,-1);
                        }
                        
                        //expected number of events for signal
                        for(unsigned int i=0;i<SigMassSplitting.size();i++)
                        {
                            //expected number of events
                            nSig = h2SigSum[i]->Integral(bin1,-1,bin2,-1);
                            //nSig = h2SigSum[i]->Integral(bin1,-1,0,-1);
                            
                            //Significance
                            //if(SigMassSplitting[i].MassDiff==100 && bin2==1) cout<<bin1<<": "<<nBG<<", "<<nSig<<", "<<GetSignificance(nSig,nBG)<<endl;
                            if(nBG>0) hSignificance[i]->SetBinContent(bin1,bin2,GetSignificance(nSig,nBG));
                        }
                    }
                }
                
                for(unsigned int j=0;j<BGGroup.size();j++)
                {
                    delete BGGroup[j].h3;
                }
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    delete h2SigSum[i];
                }
                
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    TCanvas* c2 = new TCanvas();
                    c2->cd();
                    c2->SetRightMargin(0.16);
                    hSignificance[i]->Draw("colz");
                    
                    {
                        TString NameTemp = RegionInfo[RegionIndex].RegionName;
                        NameTemp += ": ";
                        NameTemp += SigMassSplitting[i].IDName;
                        
                        TLatex lt1;
                        lt1.DrawLatexNDC(0.05,0.05,NameTemp.Data());
                    }
                    
                    TString NameTemp = "plot/";
                    NameTemp += "significance_";
                    NameTemp += Var[VarIndex[0]].VarName;
                    NameTemp += "_";
                    NameTemp += Var[VarIndex[1]].VarName;
                    NameTemp += "_";
                    NameTemp += RegionInfo[RegionIndex].RegionName;
                    NameTemp += "_";
                    NameTemp += TString::Format("%.0u",i);
                    NameTemp += ".eps";
                    c2->Print(NameTemp,"eps");
                    
                    delete hSignificance[i];
                    delete c2;
                }
            }
        
            //delete tree2
            for(unsigned int i=0;i<DataSampleName.size();i++)
            {
                delete tree2Data[i];
                if(ChannelInfo[channelRepresentative].isSS_qF) delete tree2DataOS[i];
            }
            for(unsigned int j=0;j<tree2BGMC.size();j++)
            {
                for(unsigned int k=0;k<tree2BGMC[j].size();k++)
                {
                    //Z pt reweighting
                    if(dorw
                       &&
                       (
                        RegionInfo[RegionIndex].RegionName == "CR_nonISR_OS_ee"   ||
                        RegionInfo[RegionIndex].RegionName == "CR_ISR_OS_ee"      ||
                        RegionInfo[RegionIndex].RegionName == "CR_nonISR_OS_mumu" ||
                        RegionInfo[RegionIndex].RegionName == "CR_ISR_OS_mumu"
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
                        NameTemp += TString::Itoa(BGGroup[j].info->BGMCIndex[k],10);
                        delete tree2BGMC[j][k]->GetFriend(NameTemp.Data());
                    }
                    
                    //for charge filp BG
                    if(docfw
                       &&
                       (
                        RegionInfo[RegionIndex].RegionName == "CR_nonISR_SS_ee" ||
                        RegionInfo[RegionIndex].RegionName == "CR_ISR_SS_ee"    )
                       &&
                       BGGroup[j].info->GroupName == "Zee"
                       )
                    {
                        TString NameTemp = "tree_cfw_";
                        NameTemp += TString::Itoa(BGGroup[j].info->BGMCIndex[k],10);
                        delete tree2BGMC[j][k]->GetFriend(NameTemp.Data());
                    }
                    
                    delete tree2BGMC[j][k];
                }
            }
            for(unsigned int i=0;i<SigSampleInfo.size();i++)
            {
                delete tree2Sig[i];
            }
        }
        
        if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
        {
            //Legend
            TLegend* leg;
            {
                double xl2=0.92;
                double yl2=0.95;
                double xl1=xl2-0.3;
                double yl1=yl2-0.2;
                
                leg = new TLegend(xl1,yl1,xl2,yl2);
                leg->SetNColumns(2);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetBorderSize(0);
                
                for(unsigned int j=0;j<SRBGGroup.size();j++)
                {
                    leg->AddEntry(SRBGGroup[j].h2,SRBGGroup[j].info->LegendName.Data(),"fl");
                }
                
                for(unsigned int j=0;j<SigMassSplitting.size();j++)
                {
                    TString NameTemp = SigMassSplitting[j].IDName;
                    leg->AddEntry(h2SRSig[j],NameTemp.Data(),"l");
                }
            }
            
            std::vector<Group> SRBGGroup2 = SRBGGroup;
            std::sort(SRBGGroup2.begin(),SRBGGroup2.end(),compare2);
            
            //stack
            THStack stackSR;
            /*
            stackSR.Add(SRBGGroup[5].h2); //multi top
            stackSR.Add(SRBGGroup[3].h2); //single top
            stackSR.Add(SRBGGroup[0].h2); //Z+jets
            stackSR.Add(SRBGGroup[1].h2); //W+jets
            stackSR.Add(SRBGGroup[9].h2); //Higgs
            stackSR.Add(SRBGGroup[4].h2); //ttV
            stackSR.Add(SRBGGroup[8].h2); //VVV
            stackSR.Add(SRBGGroup[7].h2); //Vgamma
            stackSR.Add(SRBGGroup[2].h2); //ttbar
            stackSR.Add(SRBGGroup[6].h2); //VV
            */
            
            for(unsigned int i=0;i<SRBGGroup2.size();i++)
            {
                stackSR.Add(SRBGGroup2[i].h2);
            }
            
            for(unsigned int RegionIndex=RegionGroup[RegionGroupIndex].lower;RegionIndex<=RegionGroup[RegionGroupIndex].upper;RegionIndex++)
            {
                const unsigned int RegionIndex2 = RegionIndex - RegionGroup[RegionGroupIndex].lower;
                TString NameTemp = RegionInfo[RegionIndex].RegionShortName;
                h2SRSignificance[0]->GetXaxis()->SetBinLabel(RegionIndex2+1,NameTemp);
            }
            
            TPad* pad1 = nullptr;
            TPad* pad2 = nullptr;
            {
                //size for two pads
                const double size1 = 0.65;
                const double size2 = 0.35;
                const double scale1 = 1/size1;
                const double scale2 = 1/size2;
                
                //adjust the title
                const double x_label = h2SRSig[0]->GetXaxis()->GetLabelSize();
                const double y_label = h2SRSig[0]->GetYaxis()->GetLabelSize();
                const double x_title = h2SRSig[0]->GetXaxis()->GetTitleSize();
                const double y_title = h2SRSig[0]->GetYaxis()->GetTitleSize();
                const double y_offset = h2SRSig[0]->GetYaxis()->GetTitleOffset();
                
                h2SRSig[0]->GetXaxis()->SetLabelSize(x_label*scale1);
                h2SRSig[0]->GetYaxis()->SetLabelSize(y_label*scale1);
                h2SRSig[0]->GetXaxis()->SetTitleSize(x_title*scale1);
                h2SRSig[0]->GetYaxis()->SetTitleSize(y_title*scale1);
                h2SRSig[0]->GetYaxis()->SetTitleOffset(y_offset/scale1);
                
                h2SRSignificance[0]->GetYaxis()->SetTitle("Significance");
                h2SRSignificance[0]->GetYaxis()->SetRangeUser(-1,4);
                h2SRSignificance[0]->GetYaxis()->SetNdivisions(6);
                
                h2SRSignificance[0]->GetXaxis()->SetLabelSize(x_label*scale2);
                h2SRSignificance[0]->GetYaxis()->SetLabelSize(y_label*scale2);
                h2SRSignificance[0]->GetXaxis()->SetTitleSize(x_title*scale2);
                h2SRSignificance[0]->GetYaxis()->SetTitleSize(y_title*scale2);
                h2SRSignificance[0]->GetYaxis()->SetTitleOffset(y_offset/scale2);
                h2SRSignificance[0]->GetYaxis()->CenterTitle();
                
                //Two pads
                pad1 = new TPad("pad1","pad1",0,1-size1,1,1);
                pad2 = new TPad("pad2","pad2",0,0,1,size2);
                pad1->SetBottomMargin(0.1);
                pad2->SetBottomMargin(0.4);
                pad2->SetTopMargin(0.1);
                pad2->SetGridy();
            }
            pad1->SetLogy(1);
            
            //Draw
            gStyle->SetOptStat(0);
            TCanvas* c2 = new TCanvas();
            c2->cd();
            pad1->Draw();
            pad1->cd();
            h2SRSig[0]->SetMinimum(1e-3);
            h2SRSig[0]->SetMaximum(1e2);
            h2SRSig[0]->Draw("axis");
            stackSR.Draw("histsame");
            for(unsigned int j=0;j<SigMassSplitting.size();j++)
            {
                h2SRSig[j]->Draw("histsame");
            }
            leg->Draw();
            
            {
                TLatex lt1;
                lt1.DrawLatexNDC(0.2,0.13,setupTag.Data());
                lt1.SetTextSize(lt1.GetTextSize());
            }
            
            c2->cd();
            pad2->Draw();
            pad2->cd();
            h2SRSignificance[0]->Draw();
            for(unsigned int j=1;j<SigMassSplitting.size();j++)
            {
                h2SRSignificance[j]->Draw("histsame");
            }
            
            //export histograms in eps format
            TString NameTemp = "plot/";
            NameTemp += RegionGroup[RegionGroupIndex].GroupName;
            NameTemp += ".eps";
            c2->Print(NameTemp,"eps");
            delete c2;
            
            //delete
            for(unsigned int j=0;j<SRBGGroup.size();j++)
            {
                delete SRBGGroup[j].h2;
            }
            for(unsigned int j=0;j<SigMassSplitting.size();j++)
            {
                delete h2SRSig[j];
                delete h2SRSignificance[j];
            }
        }
        
        //combined significance plot
        if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
        {
            TGraph2D* g_significance = new TGraph2D();
            for(unsigned int j=0;j<SigSampleInfo.size();j++)
            {
                g_significance->SetPoint(g_significance->GetN(), SigSampleInfo[j].Mass1, SigSampleInfo[j].Mass2, sqrt(SigSampleInfo[j].significance2));
                //cout<<SigSampleInfo[j].Mass1<<", "<<SigSampleInfo[j].Mass2<<", "<<sqrt(SigSampleInfo[j].significance2)<<endl;
            }
            
            TH2D* h2 = g_significance->GetHistogram();
            h2->GetXaxis()->SetTitle("m_{#tilde{#chi}^{#pm}_{1}/#tilde{#chi}^{0}_{2}} [GeV]");
            h2->GetYaxis()->SetTitle("m_{#tilde{#chi}^{0}_{1}} [GeV]");
            h2->GetZaxis()->SetTitle("Z_{n}");
            
            h2->GetXaxis()->SetRangeUser(0,350);
            h2->GetZaxis()->SetRangeUser(0,4);
            
            gStyle->SetNumberContours(255);
            {
                gStyle->SetPalette(1);
                
                TCanvas* c2 = new TCanvas();
                c2->cd();
                c2->SetRightMargin(0.16);
                g_significance->Draw("colz");
                
                {
                    TList* grL = g_significance->GetContourList(1);
                    TIter next1(grL);
                    while(true)
                    {
                        TGraph* obj = (TGraph*) next1();
                        if(obj == nullptr) break;
                        obj->SetLineWidth(2);
                        obj->SetLineStyle(2);
                        obj->Draw("same");
                        obj->SetName("s0");
                    }
                }
                TGraph* xt0 = (TGraph*) gPad->GetPrimitive("s0");
                
                {
                    TList* grL = g_significance->GetContourList(2);
                    TIter next1(grL);
                    while(true)
                    {
                        TGraph* obj = (TGraph*) next1();
                        if(obj == nullptr) break;
                        obj->SetLineWidth(2);
                        obj->SetLineStyle(5);
                        obj->Draw("same");
                        obj->SetName("s1");
                    }
                }
                TGraph* xt1 = (TGraph*) gPad->GetPrimitive("s1");
                
                {
                    TList* grL = g_significance->GetContourList(3);
                    TIter next1(grL);
                    while(true)
                    {
                        TGraph* obj = (TGraph*) next1();
                        if(obj == nullptr) break;
                        obj->SetLineWidth(2);
                        obj->SetLineStyle(9);
                        obj->Draw("same");
                        obj->SetName("s2");
                    }
                }
                TGraph* xt2 = (TGraph*) gPad->GetPrimitive("s2");
                
                TLegend* leg;
                {
                    Double_t xl1=0.2, yl1=0.7, xl2=xl1+0.2, yl2=yl1+0.2;
                    leg = new TLegend(xl1,yl1,xl2,yl2);
                    leg->SetFillStyle(0);
                    leg->SetBorderSize(0);
                    
                    leg->AddEntry(xt0, "Z_{n}=1", "l");
                    leg->AddEntry(xt1, "Z_{n}=2", "l");
                    leg->AddEntry(xt2, "Z_{n}=3", "l");
                }
                leg->Draw();
                
                {
                    TLatex lt1;
                    lt1.DrawLatexNDC(0.05,0.05,"combined");
                    lt1.SetTextSize(lt1.GetTextSize());
                }
                
                {
                    TLatex lt1;
                    lt1.DrawLatexNDC(0.4,0.05,setupTag.Data());
                    lt1.SetTextSize(lt1.GetTextSize());
                }
                
                TString NameTemp = "plot/";
                NameTemp += "combine_significance_";
                NameTemp += TString::Itoa(SigOptimizingIndex,10);
                NameTemp += ".eps";
                c2->Print(NameTemp,"eps");
                
                delete leg;
                delete c2;
            }
            
            delete g_significance;
        }
    }
    
    //delete fout_plot
    delete fout_plot;
    
    //expN_CR_OS.tex and expN_SR_SS.tex
    //expN_CR_OS_1B.tex and expN_CR_SS_1B.tex
    //expN_CR_OS_2B.tex and expN_CR_SS_2B.tex
    //expN_CR_SS_mumu_low_mT2.tex
    //expN_CR_SS_ee_Zmass.tex
    //expN_SR_SS_0B.tex
    for(unsigned int RegionGroupIndex=0;RegionGroupIndex<RegionGroup.size();RegionGroupIndex++)
    {
        if(!
           (/* RegionGroup[RegionGroupIndex].GroupName == "CR_OS"    ||
            RegionGroup[RegionGroupIndex].GroupName == "SR_SS"    ||
            RegionGroup[RegionGroupIndex].GroupName == "CR_OS_1B" ||
            RegionGroup[RegionGroupIndex].GroupName == "CR_SS_1B" ||
            RegionGroup[RegionGroupIndex].GroupName == "CR_OS_2B" ||
            RegionGroup[RegionGroupIndex].GroupName == "CR_SS_2B" ||
            RegionGroup[RegionGroupIndex].GroupName == "CR_SS_mumu_low_mT2" ||
            RegionGroup[RegionGroupIndex].GroupName == "CR_SS_ee_Zmass" ||
            RegionGroup[RegionGroupIndex].GroupName == "SR_SS_0B" ||
            RegionGroup[RegionGroupIndex].GroupName == "SR_SS_run1" || */
            RegionGroup[RegionGroupIndex].GroupName == "SR_SS_pre")
           ) continue;
        
        TString PathName = "latex/data/expN_";
        PathName += RegionGroup[RegionGroupIndex].GroupName;
        PathName += ".tex";
        
        ofstream fout;
        fout.open(PathName.Data());
        
        for(unsigned int RegionIndex=RegionGroup[RegionGroupIndex].lower;RegionIndex<=RegionGroup[RegionGroupIndex].upper;RegionIndex++)
        {
            const unsigned int SixChannel = RegionIndex - RegionGroup[RegionGroupIndex].lower;
            
            if(RegionGroup[RegionGroupIndex].GroupName == "CR_OS_1B" ||
               RegionGroup[RegionGroupIndex].GroupName == "CR_OS_2B" ||
               ( (RegionGroup[RegionGroupIndex].GroupName == "CR_SS_1B" ||
                  RegionGroup[RegionGroupIndex].GroupName == "CR_SS_2B" ||
                  RegionGroup[RegionGroupIndex].GroupName == "SR_SS_0B" ) &&
               !(SixChannel == 1 || SixChannel == 4) )) continue;
            
            TString latexName = RegionInfo[RegionIndex].RegionShortName;
            latexName.ReplaceAll("_","\\_");
            
            fout<<"\\begin{frame}{";
            fout<<"For ";
            fout<<latexName.Data();
            fout<<"}"<<endl;
            
            fout<<"\\vspace{5mm}"<<endl;
            fout<<"\\begin{tabular}{|c|c|c|}"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"& Number of events & Significance \\\\"<<endl;
            fout<<"\\hline"<<endl;
            
            fout<<"\\input{data/expN/"<<RegionInfo[RegionIndex].RegionName.Data()<<".tex}"<<endl;
            
            fout<<"\\end{tabular}"<<endl;
            fout<<"\\end{frame}"<<endl<<endl;
        }
        
        fout.close();
    }
    
    //For VV
    if(doVVCount)
    {
        for(unsigned int RegionGroupIndex=0;RegionGroupIndex<RegionGroup.size();RegionGroupIndex++)
        {
            if(!
               (RegionGroup[RegionGroupIndex].GroupName == "SR_SS_run1")
               )continue;
            
            TString PathName = "latex/data/expN_";
            PathName += RegionGroup[RegionGroupIndex].GroupName;
            PathName += "_BGVV.tex";
            
            ofstream fout;
            fout.open(PathName.Data());
            
            for(unsigned int RegionIndex=RegionGroup[RegionGroupIndex].lower;RegionIndex<=RegionGroup[RegionGroupIndex].upper;RegionIndex++)
            {
                TString latexName = RegionInfo[RegionIndex].RegionName;
                latexName.ReplaceAll("_","\\_");
                
                fout<<"\\begin{frame}{Expected number of events for VV \\\\ ";
                fout<<"For ";
                fout<<latexName.Data();
                fout<<"}"<<endl;
                
                fout<<"\\vspace{5mm}"<<endl;
                fout<<"\\begin{tabular}{|c|c|}"<<endl;
                fout<<"\\hline"<<endl;
                fout<<"& Number of events \\\\"<<endl;
                fout<<"\\hline"<<endl;
                
                fout<<"\\input{data/expN/"<<RegionInfo[RegionIndex].RegionName.Data()<<"_BGVV.tex}"<<endl;
                
                fout<<"\\end{tabular}"<<endl;
                fout<<"\\end{frame}"<<endl<<endl;
            }
            
            fout.close();
        }
    }
    
    //plot_CR_SS_mumu_low_mT2.tex
    //plot_CR_SS_ee_Zmass.tex
    for(unsigned int RegionGroupIndex=0;RegionGroupIndex<RegionGroup.size();RegionGroupIndex++)
    {
        if(! //Cutflow Attention
           (//RegionGroup[RegionGroupIndex].GroupName == "CR_SS_mumu_low_mT2" ||
            RegionGroup[RegionGroupIndex].GroupName == "CR_SS_ee_Zmass"     )
           ) continue;
        
        TString PathName = "latex/data/plot_";
        PathName += RegionGroup[RegionGroupIndex].GroupName;
        PathName += ".tex";
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
            
            TString latexName = RegionGroup[RegionGroupIndex].GroupName;
            latexName.ReplaceAll("_","\\_");
            
            fout<<"\\begin{frame}{For ";
            fout<<latexName.Data();
            fout<<" \\\\ ";
            fout<<Var[VarIndex].latexName.Data();
            fout<<"}"<<endl;
            
            fout<<"\\Wider{"<<endl;
            for(unsigned int RegionIndex=0;RegionIndex<=0;RegionIndex++)
            {
                fout<<"\\includegraphics[width=0.8\\textwidth]{"
                <<Var[VarIndex].VarName.Data()<<"_"<<RegionInfo[RegionGroup[RegionGroupIndex].lower +RegionIndex].RegionName.Data()<<"}";
                fout<<endl;
            }
            fout<<"}"<<endl;
            fout<<"\\end{frame}"<<endl<<endl;
        }
        
        fout.close();
    }
    
    //plot_CR_OS.tex and plot_SR_SS.tex
    //plot_CR_OS_1B.tex and plot_CR_SS_1B.tex
    //plot_CR_OS_2B.tex and plot_CR_SS_2B.tex
    //plot_SR_SS_0B.tex
    for(unsigned int RegionGroupIndex=0;RegionGroupIndex<RegionGroup.size();RegionGroupIndex++)
    {
        if(! //Cutflow Attention
           (/*RegionGroup[RegionGroupIndex].GroupName == "CR_OS"    ||
            RegionGroup[RegionGroupIndex].GroupName == "SR_SS"    ||
            RegionGroup[RegionGroupIndex].GroupName == "CR_OS_1B" ||
            RegionGroup[RegionGroupIndex].GroupName == "CR_SS_1B" ||
            RegionGroup[RegionGroupIndex].GroupName == "CR_OS_2B" ||
            RegionGroup[RegionGroupIndex].GroupName == "CR_SS_2B" ||
            RegionGroup[RegionGroupIndex].GroupName == "SR_SS_0B" ||
            RegionGroup[RegionGroupIndex].GroupName == "SR_SS_run1" ||*/
            RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt" ||
            RegionGroup[RegionGroupIndex].GroupName == "SR_SS_pre" )
           ) continue;
        
        const unsigned int RegionN = RegionGroup[RegionGroupIndex].upper - RegionGroup[RegionGroupIndex].lower +1;
        
        TString PathName = "latex/data/plot_";
        PathName += RegionGroup[RegionGroupIndex].GroupName;
        if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
        {
            PathName += "_";
            PathName += TString::Itoa(SigOptimizingIndex,10);
        }
        PathName += ".tex";
        
        ofstream fout;
        fout.open(PathName.Data());
        
        unsigned int Nminus1 = 0;
        vector<TString> VarPlot;
        //if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
        if(true)
        {
            //N-1 plot
            VarPlot.push_back("pt1");
            VarPlot.push_back("pt2");
            VarPlot.push_back("dEta");
            VarPlot.push_back("MET");
            VarPlot.push_back("meff");
            VarPlot.push_back("mt1");
            VarPlot.push_back("mlj");
            VarPlot.push_back("mTtwo");
            Nminus1 = 8;
            
            /*
            //N plot
            VarPlot.push_back("nBJet");
            VarPlot.push_back("nCJet");
            VarPlot.push_back("eta1");
            VarPlot.push_back("eta2");
            */
            
            /*
            VarPlot.push_back("mt1");
            VarPlot.push_back("mt2");
            VarPlot.push_back("jetpt");
            VarPlot.push_back("jeteta");
            VarPlot.push_back("nJet");
            VarPlot.push_back("l12_dPhi");
            VarPlot.push_back("l12_MET_dPhi");
            VarPlot.push_back("jet0_MET_dPhi");
            VarPlot.push_back("mjj");
            VarPlot.push_back("phi1");
            */
        }
        else
        {
            for(unsigned int VarIndex=0;VarIndex<Var.size();VarIndex++)
            {
                VarPlot.push_back(Var[VarIndex].VarName);
            }
        }
        
        for(unsigned int i=0;i<VarPlot.size();i++)
        {
            unsigned int VarIndex = findVarIndex(VarPlot[i],Var);
            
            if(Var[VarIndex].VarName=="averageMu") continue;
            if(Var[VarIndex].VarName=="nVtx") continue;
            
            if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
            {
                if(i == 0) fout<<"\\subsection{``N-1'' plots}"<<endl;
                if(i == Nminus1) fout<<"\\subsection{``N'' plots}"<<endl;
            }
            
            fout<<"\\begin{frame}{";
            fout<<Var[VarIndex].latexName.Data();
            fout<<"}"<<endl;
            
            fout<<"\\Wider{"<<endl;
            for(unsigned int RegionIndex=RegionGroup[RegionGroupIndex].lower;RegionIndex<=RegionGroup[RegionGroupIndex].upper;RegionIndex++)
            {
                const unsigned int SixChannel = RegionIndex - RegionGroup[RegionGroupIndex].lower;
                
                if(SixChannel<=2 &&
                   (Var[VarIndex].VarName=="bjetpt"  ||
                    Var[VarIndex].VarName=="bjeteta" ||
                    Var[VarIndex].VarName=="bjetphi" ||
                    Var[VarIndex].VarName=="cjetpt"  ||
                    Var[VarIndex].VarName=="cjeteta" ||
                    Var[VarIndex].VarName=="cjetphi" )
                   ) continue;
                
                fout<<"\\includegraphics[width=";
                if(RegionN == 3 || RegionN == 4) fout<<"0.4";
                else fout<<"0.33";
                fout<<"\\textwidth]{";
                fout<<Var[VarIndex].VarName.Data();
                fout<<"_";
                fout<<RegionInfo[RegionIndex].RegionName.Data();
                if(RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
                {
                    fout<<"_";
                    fout<<TString::Itoa(SigOptimizingIndex,10);
                }
                fout<<"}";
                if(RegionN == 3 || RegionN == 4)
                {
                    if(SixChannel==1) fout<<" \\\\";
                }
                else
                {
                    if(SixChannel==2) fout<<" \\\\";
                }
                fout<<endl;
            }
            fout<<"}"<<endl;
            fout<<"\\end{frame}"<<endl<<endl;
        }
        
        fout.close();
    }
    
    //plot_significances.tex
    for(unsigned int RegionGroupIndex=0;RegionGroupIndex<RegionGroup.size();RegionGroupIndex++)
    {
        /*if(! //Cutflow Attention
           (RegionGroup[RegionGroupIndex].GroupName == "SR_SS_0B" )
           )*/ continue;
        
        TString PathName = "latex/data/plot_significances_";
        PathName += RegionGroup[RegionGroupIndex].GroupName;
        PathName += ".tex";
        
        ofstream fout;
        fout.open(PathName.Data());
        
        unsigned int VarIndex[2];
        VarIndex[0] = 6; // ptll
        VarIndex[1] = 7; // MET
        for(unsigned int i=0;i<SigMassSplitting.size();i++)
        {
            TString latexName = Var[VarIndex[0]].VarName;
            latexName += "_";
            latexName += Var[VarIndex[1]].VarName;
            latexName += "_";
            latexName += RegionGroup[RegionGroupIndex].GroupName;
            latexName += "_";
            latexName += TString::Format("%.0u",i);
            
            latexName.ReplaceAll("_","\\_");
            
            fout<<"\\begin{frame}{For ";
            fout<<latexName.Data();
            fout<<"}"<<endl;
            
            fout<<"\\Wider{"<<endl;
            for(unsigned int RegionIndex=RegionGroup[RegionGroupIndex].lower;RegionIndex<=RegionGroup[RegionGroupIndex].upper;RegionIndex++)
            {
                const unsigned int SixChannel = RegionIndex - RegionGroup[RegionGroupIndex].lower;
                
                fout<<"\\includegraphics[width=0.33\\textwidth]{";
                fout<<"significance_ptll_MET_";
                fout<<RegionInfo[RegionIndex].RegionName;
                fout<<"_";
                fout<<TString::Format("%.0u",i);
                fout<<"}";
                
                if(SixChannel==2) fout<<" \\\\";
                fout<<endl;
            }
            
            fout<<"}"<<endl;
            fout<<"\\end{frame}"<<endl<<endl;
        }
        fout.close();
    }
    
    //opt_*.tex
    for(unsigned int RegionGroupIndex=0;RegionGroupIndex<RegionGroup.size();RegionGroupIndex++)
    {
        if(!
           (RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
           ) continue;
        
        TString PathName = "latex/data/optimization/opt_";
        PathName += RegionGroup[RegionGroupIndex].GroupName;
        PathName += "_";
        PathName += TString::Itoa(SigOptimizingIndex,10);
        PathName += ".tex";
        
        ofstream fout;
        fout.open(PathName.Data());
        
        for(unsigned int RegionIndex=RegionGroup[RegionGroupIndex].lower;RegionIndex<=RegionGroup[RegionGroupIndex].upper;RegionIndex++)
        {
            TString latexName = RegionInfo[RegionIndex].RegionName;
            latexName.ReplaceAll("_","\\_");
            
            fout<<"\\begin{frame}{optimization}"<<endl;
            fout<<"\\tiny"<<endl;
            fout<<"\\input{data/optimization/cut_latex/cut_";
            fout<<RegionInfo[RegionIndex].RegionName.Data();
            fout<<"_";
            fout<<SigOptimizingIndex;
            fout<<".tex}"<<endl;

            fout<<"\\begin{tabular}{|c|c|c|}"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"& Number of events & Significance \\\\"<<endl;
            fout<<"\\hline"<<endl;
            
            fout<<"\\input{data/expN/";
            fout<<RegionInfo[RegionIndex].RegionName.Data();
            fout<<"_";
            fout<<SigOptimizingIndex;
            fout<<".tex}"<<endl;
            fout<<"\\end{tabular}"<<endl;
            fout<<"\\end{frame}"<<endl<<endl;
        }
        
        fout.close();
    }
    
    //cut_latex
    {
        TString optGroupName = "SR_SS_opt";
        unsigned int optRegionGroupIndex = 0;
        for(unsigned int i=0;i<RegionGroup.size();i++)
        {
            if(RegionGroup[i].GroupName == optGroupName)
            {
                optRegionGroupIndex = i;
            }
        }
        
        for(unsigned int RegionIndex=RegionGroup[optRegionGroupIndex].lower;RegionIndex<=RegionGroup[optRegionGroupIndex].upper;RegionIndex++)
        {
            ifstream fin;
            {
                TString PathName = "latex/data/optimization/";
                if(useDani) PathName += "cut_txt_Dani/";
                else PathName += "cut_txt/";
                PathName += "cut_";
                PathName += RegionInfo[RegionIndex].RegionName;
                PathName += "_";
                PathName += TString::Itoa(SigOptimizingIndex,10);
                PathName += ".txt";
                fin.open(PathName.Data());
            }
            
            ofstream fout;
            {
                TString PathName = "latex/data/optimization/cut_latex/cut_";
                PathName += RegionInfo[RegionIndex].RegionName;
                PathName += "_";
                PathName += TString::Itoa(SigOptimizingIndex,10);
                PathName += ".tex";
                fout.open(PathName.Data());
            }
            
            TString latexName = RegionInfo[RegionIndex].RegionName;
            latexName.ReplaceAll("_","\\_");
            fout<<latexName.Data();
            fout<<": \\\\"<<endl;
            
            for(unsigned int i=0;i<RegionInfo[RegionIndex].OptimizingCut[SigOptimizingIndex].size();i++)
            {
                for(unsigned int j=0;j<RegionInfo[RegionIndex].OptimizingCut[SigOptimizingIndex][i].size();j++)
                {
                    TString VarName;
                    double lower;
                    double upper;
                    fin>>VarName;
                    fin>>lower;
                    fin>>upper;
                    
                    unsigned int VarIndex = findVarIndex(VarName,Var);
                    
                    if(lower == 0)
                    {
                        if(upper == -1)
                        {
                        }
                        else
                        {
                            // x < a
                            fout<<Var[VarIndex].latexName.Data();
                            fout<<" $<";
                            fout<<upper;
                            fout<<"$";
                            fout<<" \\\\"<<endl;
                        }
                    }
                    else
                    {
                        if(upper == -1)
                        {
                            // x >= a
                            fout<<Var[VarIndex].latexName.Data();
                            fout<<" $\\geq ";
                            fout<<lower;
                            fout<<"$";
                            fout<<" \\\\"<<endl;
                        }
                        else
                        {
                            // a <= x < b
                            fout<<"$";
                            fout<<lower;
                            fout<<" \\leq$ ";
                            fout<<Var[VarIndex].latexName.Data();
                            fout<<" $<";
                            fout<<upper;
                            fout<<"$";
                            fout<<" \\\\"<<endl;
                        }
                    }
                }
            }
            fout.close();
            fin.close();
        }
    }
    
    //cut table
    {
        TString optGroupName = "SR_SS_opt";
        unsigned int optRegionGroupIndex = 0;
        for(unsigned int i=0;i<RegionGroup.size();i++)
        {
            if(RegionGroup[i].GroupName == optGroupName)
            {
                optRegionGroupIndex = i;
            }
        }
        
        ofstream fout;
        {
            TString PathName = "latex/data/optimization/cut_table_";
            PathName += RegionGroup[optRegionGroupIndex].GroupName;
            PathName += "_";
            PathName += TString::Itoa(SigOptimizingIndex,10);
            PathName += ".tex";
            fout.open(PathName.Data());
        }
        
        const unsigned int RegionN = RegionGroup[optRegionGroupIndex].upper - RegionGroup[optRegionGroupIndex].lower +1;
        std::vector<ifstream*> fin;
        for(unsigned int RegionIndex=RegionGroup[optRegionGroupIndex].lower;RegionIndex<=RegionGroup[optRegionGroupIndex].upper;RegionIndex++)
        {
            const unsigned int SixChannel = RegionIndex - RegionGroup[optRegionGroupIndex].lower;
            TString PathName = "latex/data/optimization/";
            if(useDani) PathName += "cut_txt_Dani/";
            else PathName += "cut_txt/";
            PathName += "cut_";
            PathName += RegionInfo[RegionIndex].RegionName;
            PathName += "_";
            PathName += TString::Itoa(SigOptimizingIndex,10);
            PathName += ".txt";
            
            ifstream* finElement = new ifstream(PathName.Data());
            fin.push_back(finElement);
        }
        
        if(RegionN == 6)
        {
            fout<<"\\begin{tabular}{|c|c|c|c|c|c|c|}"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"& SR$ee$-1 & SR$ee$-2 & SR$\\mu\\mu$-1 & SR$\\mu\\mu$-2 & SR$e\\mu$-1 & SR$e\\mu$-2 \\\\"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"Lepton flavours & $ee$ & $ee$ & $\\mu\\mu$ & $\\mu\\mu$ & $e\\mu$ & $e\\mu$ \\\\"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"$n_{jets}$ & 1 & 2 or 3 & 1 & 2 or 3 & 1 & 2 or 3 \\\\"<<endl;
            fout<<"\\hline"<<endl;
        }
        else if(RegionN == 2)
        {
            fout<<"\\begin{tabular}{|c|c|c|}"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"& SRjet1 & SRjet23 \\\\"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"$n_{jets}$ & 1 & 2 or 3 \\\\"<<endl;
            fout<<"\\hline"<<endl;
        }
        else if(RegionN == 1)
        {
            fout<<"\\begin{tabular}{|c|c|}"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"& SRjet123 \\\\"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"$n_{jets}$ & 1 to 3 \\\\"<<endl;
            fout<<"\\hline"<<endl;
        }
        else if(RegionN == 4)
        {
            fout<<"\\begin{tabular}{|c|c|c|c|c|}"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"& SRjet1\\_off & SRjet23\\_off & SRjet1\\_on & SRjet23\\_on \\\\"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"$n_{jets}$ & 1 & 2 or 3 & 1 & 2 or 3 \\\\"<<endl;
            fout<<"\\hline"<<endl;
        }
        
        for(unsigned int i=0;i<RegionInfo[ RegionGroup[optRegionGroupIndex].lower ].OptimizingCut[SigOptimizingIndex].size();i++)
        {
            for(unsigned int j=0;j<RegionInfo[ RegionGroup[optRegionGroupIndex].lower ].OptimizingCut[SigOptimizingIndex][i].size();j++)
            {
                unsigned int VarIndex = findVarIndex(RegionInfo[ RegionGroup[optRegionGroupIndex].lower ].OptimizingCut[SigOptimizingIndex][i][j].RelatedVariable,Var);
                
                if(Var[VarIndex].VarName == "pt1") fout<<"Leading lepton $p_T$";
                else if(Var[VarIndex].VarName == "pt2") fout<<"Sub-leading lepton $p_T$";
                else fout<<Var[VarIndex].latexName.Data();
                
                fout<<" ";
                if(Var[VarIndex].unit != "") fout<<Var[VarIndex].unit.Data()<<" ";
                
                unsigned int SixChannel = 0;
                while(true)
                {
                    TString VarName;
                    double lower;
                    double upper;
                    (*(fin[SixChannel]))>>VarName;
                    (*(fin[SixChannel]))>>lower;
                    (*(fin[SixChannel]))>>upper;
                    
                    fout<<"& ";
                    if(lower == 0)
                    {
                        if(upper == -1)
                        {
                            fout<<"- ";
                        }
                        else
                        {
                            // x < a
                            fout<<"$<";
                            fout<<upper;
                            fout<<"$ ";
                        }
                    }
                    else
                    {
                        if(upper == -1)
                        {
                            // x >= a
                            fout<<"$\\geq ";
                            fout<<lower;
                            fout<<"$ ";
                        }
                        else
                        {
                            // a <= x < b
                        }
                    }
                    
                    if(RegionN == 6)
                    {
                        if(SixChannel == 0) SixChannel = 3;
                        else if(SixChannel == 3) SixChannel = 1;
                        else if(SixChannel == 1) SixChannel = 4;
                        else if(SixChannel == 4) SixChannel = 2;
                        else if(SixChannel == 2) SixChannel = 5;
                        else if(SixChannel == 5) break;
                    }
                    else
                    {
                        if(SixChannel == RegionN -1)
                        {
                            break;
                        }
                        else
                        {
                            SixChannel++;
                        }
                    }
                }
                
                fout<<"\\\\"<<endl;
                fout<<"\\hline"<<endl;
            }
        }
        
        for(unsigned int SixChannel=0;SixChannel<RegionN;SixChannel++)
        {
            fin[SixChannel]->close();
            delete fin[SixChannel];
        }
        fout.close();
    }
    
    //nSig_SR_SS_opt_0.tex
    //significance_SR_SS_opt_0.tex
    for(unsigned int RegionGroupIndex=0;RegionGroupIndex<RegionGroup.size();RegionGroupIndex++)
    {
        if(!
           (RegionGroup[RegionGroupIndex].GroupName == "SR_SS_opt")
           ) continue;
        
        const unsigned int RegionN = RegionGroup[RegionGroupIndex].upper - RegionGroup[RegionGroupIndex].lower +1;
        
        for(unsigned int i=0;i<=1;i++)
        {
            TString PathName = "latex/data/optimization/";
            if(i==0) PathName += "nSig";
            if(i==1) PathName += "significance";
            PathName += "_";
            PathName += RegionGroup[RegionGroupIndex].GroupName;
            PathName += "_";
            PathName += TString::Itoa(SigOptimizingIndex,10);
            PathName += ".tex";
            
            ofstream fout;
            fout.open(PathName.Data());
            
            fout<<"\\begin{frame}{";
            if(i==0) fout<<"signal yield plot";
            if(i==1) fout<<"significance plot";
            fout<<"}"<<endl;
            fout<<"\\Wider{"<<endl;
            for(unsigned int RegionIndex=RegionGroup[RegionGroupIndex].lower;RegionIndex<=RegionGroup[RegionGroupIndex].upper;RegionIndex++)
            {
                const unsigned int SixChannel = RegionIndex - RegionGroup[RegionGroupIndex].lower;
                
                fout<<"\\includegraphics[width=";
                if(RegionN == 3 || RegionN == 4) fout<<"0.4";
                else fout<<"0.33";
                fout<<"\\textwidth]{";
                if(i==0) fout<<"nSig";
                if(i==1) fout<<"significance";
                fout<<"_";
                fout<<RegionInfo[RegionIndex].RegionName.Data();
                fout<<"_";
                fout<<TString::Itoa(SigOptimizingIndex,10);
                fout<<"}";
                if(RegionN == 3 || RegionN == 4)
                {
                    if(SixChannel==1) fout<<" \\\\";
                }
                else
                {
                    if(SixChannel==2) fout<<" \\\\";
                }
                fout<<endl;
            }
            fout<<"}"<<endl;
            fout<<"\\end{frame}"<<endl;
            
            fout.close();
        }
    }
    
    /*
    {
        TString PathName = "latex/data/plot_zpt.tex";
        ofstream fout;
        fout.open(PathName.Data());
        for(unsigned int VarIndex=0;VarIndex<Var.size();VarIndex++)
        {
            if(Var[VarIndex].VarName=="averageMu") continue;
            
            for(unsigned int lepton=0;lepton<=1;lepton++)
            {
                fout<<"\\begin{frame}"<<endl;
                
                fout<<"\\frametitle{"<<Var[VarIndex].latexName.Data()<<" (For ";
                if(lepton==0) fout<<"ee channel";
                if(lepton==1) fout<<"mumu channel";
                fout<<")}"<<endl;
                
                fout<<"\\Wider{"<<endl;
                for(unsigned int ISR=0;ISR<=6;ISR+=6)
                {
                    fout<<"\\includegraphics[width=0.5\\textwidth]{../plot_nozpt/"
                    <<Var[VarIndex].VarName.Data()<<"_"<<ChannelInfo[lepton+ISR].ChannelName.Data()<<"}";
                    fout<<endl;
                    
                    fout<<"\\includegraphics[width=0.5\\textwidth]{../plot_zpt/"
                    <<Var[VarIndex].VarName.Data()<<"_"<<ChannelInfo[lepton+ISR].ChannelName.Data()<<"}";
                    if(ISR==0) fout<<" \\\\";
                    fout<<endl;
                }
                fout<<"}"<<endl;
                fout<<"\\end{frame}"<<endl<<endl;
            }
        }
        fout.close();
    }
    */
}
