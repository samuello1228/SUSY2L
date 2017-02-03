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

#include <TGraph2D.h>

bool dorw = 0;
bool docfw = 0;
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

struct ChannelData
{
    TString ChannelName;
    bool isSS;
    bool isSS_ee;
    unsigned int qFChannel;
    std::vector<TString> setOfBGMC;
    std::vector<TString> setOfBGData;
};

void initializeTree1new(TChain*& tree, TString& SampleID, TString& channel)
{
    tree = new TChain("tree");
    
    TString fileName = "skimming/skimming.";
    fileName += SampleID;
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

void initializeTree2(std::vector<TChain*>& tree2,std::vector<unsigned int>& SetOfChannel, std::vector<TString>& SampleID, std::vector<ChannelData>& ChannelInfo)
{
    for(unsigned int i=0;i<SampleID.size();i++)
    {
        TChain* treeTemp = new TChain("tree");
        for(unsigned int j=0;j<SetOfChannel.size();j++)
        {
            TString fileName = "skimming/skimming.";
            fileName += SampleID[i];
            fileName += "_";
            fileName += ChannelInfo[SetOfChannel[j]].ChannelName;
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
    std::vector<ChannelData> ChannelInfo;
    {
        ChannelData element;
        TString ISR[2] = {"nonISR","ISR"};
        TString sign[2] = {"OS","SS"};
        TString lepton[3] = {"ee","mumu","emu"};
        for(int i=0;i<2;i++)
        {
            for(int j=0;j<2;j++)
            {
                for(int k=0;k<3;k++)
                {
                    element.ChannelName = "";
                    element.ChannelName += ISR[i];
                    element.ChannelName += "_";
                    element.ChannelName += sign[j];
                    element.ChannelName += "_";
                    element.ChannelName += lepton[k];
                    
                    element.isSS = j==1;
                    element.isSS_ee = (j==1 && k==0);
                    if(element.isSS_ee) element.qFChannel = ChannelInfo.size() - 3;
                    
                    
                    element.setOfBGMC.clear();
                    element.setOfBGData.clear();
                    if(element.isSS)
                    {
                        if(element.isSS_ee)
                        {
                            if(docfw)
                            {
                                element.setOfBGMC.push_back("Zee");
                            }
                            else
                            {
                                element.setOfBGData.push_back("charge flip");
                            }
                        }
                        
                        element.setOfBGData.push_back("fake lepton");
                        
                        element.setOfBGMC.push_back("VV");
                        element.setOfBGMC.push_back("Vgamma");
                    }
                    else
                    {
                        element.setOfBGMC.push_back("Zee");
                        element.setOfBGMC.push_back("Zmumu");
                        element.setOfBGMC.push_back("Ztautau");
                        element.setOfBGMC.push_back("ttbar");
                        element.setOfBGMC.push_back("VV");
                        element.setOfBGMC.push_back("Vgamma");
                    }

                    ChannelInfo.push_back(element);
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
        NameTemp += ChannelInfo[0].ChannelName;
        NameTemp += ".root";
        
        TFile* file = new TFile(NameTemp.Data(),"READ");
        TH1F *h1 = (TH1F*) file->Get("hist");
        unsigned int nAOD = h1->GetBinContent(1);
        cout<<nAOD<<endl;
        BGMCnAOD.push_back(nAOD);
        
        delete file;
    }
    
    //For Signal MC
    struct SigInfo
    {
        int MassDiff;
        unsigned int ID;
        TString IDName;
        int colour;
        int statCount;
    };
    std::vector<SigInfo> SigMassSplitting;
    {
        SigInfo element;
        
        element.MassDiff = 20;    element.ID = 2;   element.colour = 6;   SigMassSplitting.push_back(element);
        element.MassDiff = 50;    element.ID = 18;  element.colour = 7;   SigMassSplitting.push_back(element);
        element.MassDiff = 100;   element.ID = 19;  element.colour = 8;   SigMassSplitting.push_back(element);
        element.MassDiff = 200;   element.ID = 20;  element.colour = 9;   SigMassSplitting.push_back(element);
        element.MassDiff = 300;   element.ID = 21;  element.colour = 14;  SigMassSplitting.push_back(element);
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
    
    for(unsigned int i=0;i<SigMassSplitting.size();i++)
    {
        TString GroupName = "Signal(";
        GroupName += TString::Itoa(SigMass1[SigMassSplitting[i].ID],10);
        GroupName += ",";
        GroupName += TString::Itoa(SigMass2[SigMassSplitting[i].ID],10);
        GroupName += ")";
        SigMassSplitting[i].IDName = GroupName;
    }
    
    //Get number of events in AOD
    std::vector<unsigned int> SignAOD;
    for(unsigned int i=0;i<SigSampleID.size();i++)
    {
        TString NameTemp = "skimming/skimming.";
        NameTemp += SigSampleID[i];
        cout<<NameTemp<<": ";
        NameTemp += "_";
        NameTemp += ChannelInfo[0].ChannelName;
        NameTemp += ".root";
        
        TFile* file = new TFile(NameTemp.Data(),"READ");
        TH1F *h1 = (TH1F*) file->Get("hist");
        unsigned int nAOD = h1->GetBinContent(1);
        cout<<nAOD<<endl;
        SignAOD.push_back(nAOD);
        
        delete file;
    }
    
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
        
        element.GroupName = "VV"; element.LegendName = "VV"; element.LatexName = "VV";
        element.lower = 73;  element.upper = 86; BGMCGroupData.push_back(element);
        
        element.GroupName = "Vgamma"; element.LegendName = "V + #gamma"; element.LatexName = "V$+\\gamma$";
        element.lower = 87;  element.upper = 106;BGMCGroupData.push_back(element);
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
        element.bin=40;         element.xmin=0;                 element.xmax=250;
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
        
        element.VarName = "jetpt";      element.VarTitle = "pT of the leading jet";             element.unit = "[GeV]";
        element.bin=40;         element.xmin=20;                element.xmax=300;
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\text{p}_{\\text{T}}$ of the leading jet";
        //Var.push_back(element);
        
        element.VarName = "l12_MET_dPhi";element.VarTitle = "phi difference between l12 and MET";        element.unit = "";
        element.bin=40;         element.xmin=-TMath::Pi();      element.xmax=TMath::Pi();
        element.log=1;          element.ymin=1e-1;              element.ymax=1;
        element.latexName = "$\\phi$ difference between l12 and MET";
        Var.push_back(element);
    }
    SetAtlasStyle();
    
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
    
    if(dorw)
    //if(false)
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
            
            element.setOfBGMC.push_back("Zee");
            element.setOfBGMC.push_back("Zmumu");
            element.setOfBGMC.push_back("Ztautau");
            element.setOfBGMC.push_back("ttbar");
            element.setOfBGMC.push_back("VV");
            element.setOfBGMC.push_back("Vgamma");
            
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
        }

        //calculate ratio plot
        TH1F* h2Ratio_rw[2];
        for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
        {
            std::vector<TChain*> tree2Data;
            initializeTree2(tree2Data,RegionInfo[RegionIndex].setOfChannel,DataSampleID,ChannelInfo);
            
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
                            initializeTree2(tree2BGMCElement,RegionInfo[RegionIndex].setOfChannel,BGMCGroupSampleID,ChannelInfo);
                            
                            tree2BGMC.push_back(tree2BGMCElement);
                            BGMCGroupXS.push_back(BGMCGroupXSElement);
                            BGMCGroupnAOD.push_back(BGMCGroupnAODElement);
                            
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
            TH1F* h2DataSum = new TH1F("DataSum",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
            
            //h2Data
            for(unsigned int j=0;j<DataSampleID.size();j++)
            {
                TH1F* hTemp = new TH1F("Data",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                
                //fill histograms from trees
                TString temp;
                temp += Var[VarIndex].VarName;
                temp += ">>Data";
                tree2Data[j]->Draw(temp.Data(),"fLwt==0");
                
                //add data
                h2DataSum->Add(hTemp);
                delete hTemp;
            }
            
            //For MC background
            for(unsigned int j=0;j<tree2BGMC.size();j++)
            {
                BGGroup[j].h2 = new TH1F(BGGroup[j].info->GroupName.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                for(unsigned int k=0;k<tree2BGMC[j].size();k++)
                {
                    TH1F* hTemp = new TH1F("BGMC",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                    
                    //fill histograms from trees
                    TString temp;
                    temp += Var[VarIndex].VarName;
                    temp += ">>BGMC";
                    tree2BGMC[j][k]->Draw(temp.Data(),"weight");
                    
                    //normalization for BG
                    hTemp->Scale(BGMCGroupXS[j][k]/BGMCGroupnAOD[j][k] *sumDataL);
                    
                    //add BG
                    BGGroup[j].h2->Add(hTemp);
                    delete hTemp;
                }
            }
            
            //ratio plot for reweighting
            {
                TString NameTemp = "Ratio_";
                NameTemp += RegionInfo[RegionIndex].RegionName;
                h2Ratio_rw[RegionIndex] = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                h2Ratio_rw[RegionIndex]->GetXaxis()->SetTitle(xaxis.Data());
                h2Ratio_rw[RegionIndex]->GetYaxis()->SetTitle("(Data - NonZjets)/Zjets");
                h2Ratio_rw[RegionIndex]->SetMarkerColor(RegionIndex+1);
                h2Ratio_rw[RegionIndex]->SetMarkerSize(1);
                h2Ratio_rw[RegionIndex]->SetLineColor(RegionIndex+1);
                h2Ratio_rw[RegionIndex]->Sumw2();
            }
            h2Ratio_rw[RegionIndex]->Add(h2DataSum);
            
            //Z+jets
            TH1F* h2Zjets = new TH1F("Zjets",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
            TH1F* h2nonZjets = new TH1F("nonZjets",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);

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
            for(unsigned int i=0;i<DataSampleID.size();i++)
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
        
        TF1* fun[2];
        //simple fit
        if(simple)
        {
            for(int RegionIndex=0;RegionIndex<=1;RegionIndex++)
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
                        for(unsigned int k=BGMCGroupData[j].lower;k<=BGMCGroupData[j].upper;k++)
                        {
                            TChain* tree1 = nullptr;
                            initializeTree1new(tree1,BGMCSampleID[k],ChannelInfo[ChannelIndex].ChannelName);
                            
                            TString NameTemp = "tree_rw_";
                            NameTemp += TString::Itoa(k,10);
                            TTree* tree2 = new TTree(NameTemp.Data(),NameTemp.Data());
                            tree2->Branch("rw",&rw,"rw/D");
                            
                            for(int m=0;m<tree1->GetEntries();m++)
                            {
                                tree1->GetEntry(m);
                                if(simple)
                                {
                                    rw=fun[RegionIndex]->Eval(ptll);
                                }
                                if(combined)
                                {
                                    rw=fun[1]->Eval(ptll);
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
        
        //delete functions
        delete fun[0];
        delete fun[1];
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
                    for(unsigned int k=BGMCGroupData[j].lower;k<=BGMCGroupData[j].upper;k++)
                    {
                        TChain* tree1 = nullptr;
                        initializeTree1new(tree1,BGMCSampleID[k],ChannelInfo[ChannelIndex].ChannelName);
                        
                        TString NameTemp = "tree_cfw_";
                        NameTemp += TString::Itoa(k,10);
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
    
    std::vector<RegionData> RegionInfo;
    {
        RegionData element;
        
        for(unsigned int ChannelIndex=0;ChannelIndex<ChannelInfo.size();ChannelIndex++)
        {
            element.RegionName = ChannelInfo[ChannelIndex].ChannelName;
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
            
            if(element.isSS)
            {
                //element.Cut = " && ( mll<76.18 || mll>106.18 )";
                element.Cut = "";
            }
            else
            {
                element.Cut = "";
            }
            
            RegionInfo.push_back(element);
        }
        
        
        //Control Region
        element.isSS = true;
        element.showData = true;
        element.showSignificance = false;
        element.Cut = " && mll>81.18 && mll<101.18";
        
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
        
        //Signal region
        TString ISR[2] = {"nonISR","ISR"};
        TString MT[2] = {"mT_0_100","mT_100_inf"};
        TString PTLL[2] = {"ptll_0_50","ptll_50_inf"};
        TString MET[3] = {"MET_0_100","MET_100_150","MET_150_inf"};
        TString lepton[3] = {"ee","mumu","emu"};
        
        TString MTCut[2] = {"mtm<100","mtm>=100"};
        TString PTLLCut[2] = {"ptll<50","ptll>=50"};
        TString METCut[3] = {"MET<100","MET>=100 && MET<150","MET>=150"};
        
        for(int i=0;i<2;i++)
        {
            for(int j=0;j<2;j++)
            {
                for(int k=0;k<2;k++)
                {
                    if(j==1 && k==1) continue;
                    for(int l=0;l<3;l++)
                    {
                        for(int m=0;m<3;m++)
                        {
                            element.RegionName = "";
                            element.Cut = "";
                            
                            element.RegionName += ISR[i];
                            
                            element.RegionName += "_";
                            element.RegionName += MT[j];
                            element.Cut += " && ";
                            element.Cut += MTCut[j];
                            
                            if(j==0)
                            {
                                element.RegionName += "_";
                                element.RegionName += PTLL[k];
                                element.Cut += " && ";
                                element.Cut += PTLLCut[k];
                            }
                            else
                            {
                                element.RegionName += "_ptll_no_cut";
                            }
                            
                            element.RegionName += "_";
                            element.RegionName += MET[l];
                            element.Cut += " && ";
                            element.Cut += METCut[l];
                            
                            element.RegionName += "_";
                            element.RegionName += lepton[m];
                            
                            element.isSS = true;
                            element.isSS_ee = m==0;
                            
                            const unsigned int ChannelIndex = 3 +6*i +m;
                            element.setOfChannel.clear();
                            element.setOfChannel.push_back(ChannelIndex);
                            element.qFChannel.clear();
                            if(element.isSS_ee) element.qFChannel.push_back(ChannelIndex-3);
                            
                            element.showData = !element.isSS;
                            element.showSignificance = element.isSS;
                            
                            RegionInfo.push_back(element);
                        }
                    }
                }
            }
        }
        
        //setOfBG
        for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
        {
            //cout<<RegionInfo[RegionIndex].RegionName.Data()<<", Cut: "<<RegionInfo[RegionIndex].Cut.Data()<<endl;
            if(RegionInfo[RegionIndex].isSS)
            {
                if(RegionInfo[RegionIndex].isSS_ee)
                {
                    if(docfw)
                    {
                        RegionInfo[RegionIndex].setOfBGMC.push_back("Zee");
                    }
                    else
                    {
                        RegionInfo[RegionIndex].setOfBGData.push_back("charge flip");
                    }
                }
                
                RegionInfo[RegionIndex].setOfBGData.push_back("fake lepton");
                
                RegionInfo[RegionIndex].setOfBGMC.push_back("VV");
                RegionInfo[RegionIndex].setOfBGMC.push_back("Vgamma");
            }
            else
            {
                RegionInfo[RegionIndex].setOfBGMC.push_back("Zee");
                RegionInfo[RegionIndex].setOfBGMC.push_back("Zmumu");
                RegionInfo[RegionIndex].setOfBGMC.push_back("Ztautau");
                RegionInfo[RegionIndex].setOfBGMC.push_back("ttbar");
                RegionInfo[RegionIndex].setOfBGMC.push_back("VV");
                RegionInfo[RegionIndex].setOfBGMC.push_back("Vgamma");
            }
        }
    }
    
    //Significance optimization
    if(doOptimize)
    {
        int CutVarIndex = 0;
        for(unsigned int i=0;i<Var.size();i++)
        {
            if(Var[i].VarName == "jetpt") CutVarIndex = i;
        }
        
        int CountVarIndex = 0;
        for(unsigned int i=0;i<Var.size();i++)
        {
            if(Var[i].VarName == "pt1") CountVarIndex = i;
        }
        
        const int SigID = 18;
        
        //cut value
        const double jetptCut[] = {20,25,30,35,40,45,50,70,100,150,200};
        //const double jetptCut[] = {25,30,35,40,45,50,70,100,150,200};
        const int cutN = sizeof(jetptCut)/sizeof(jetptCut[0]);
        
        for(int RegionIndex=4;RegionIndex<=4;RegionIndex++)
        //for(int RegionIndex=10;RegionIndex<=10;RegionIndex++)
        //for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
        {
            std::vector<TChain*> tree2Data;
            initializeTree2(tree2Data,RegionInfo[RegionIndex].setOfChannel,DataSampleID,ChannelInfo);
            
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
                            initializeTree2(tree2BGMCElement,RegionInfo[RegionIndex].setOfChannel,BGMCGroupSampleID,ChannelInfo);
                            
                            tree2BGMC.push_back(tree2BGMCElement);
                            BGMCGroupXS.push_back(BGMCGroupXSElement);
                            BGMCGroupnAOD.push_back(BGMCGroupnAODElement);
                            
                            Group BGGroupElement;
                            BGGroupElement.info = &(BGMCGroupData[j]);
                            BGGroup.push_back(BGGroupElement);
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
            initializeTree2(tree2Sig,RegionInfo[RegionIndex].setOfChannel,SigSampleID,ChannelInfo);
            
            std::vector<TChain*> tree2DataOS;
            if(RegionInfo[RegionIndex].isSS_ee) initializeTree2(tree2DataOS,RegionInfo[RegionIndex].qFChannel,DataSampleID,ChannelInfo);
 
            const double uncertainty[] = {0.1,0.2,0.3};
            const int uncertaintyN = sizeof(uncertainty)/sizeof(uncertainty[0]);
            double Significance[uncertaintyN][cutN];
            //uncertainty,cut value
            
            double sumOfEvent[cutN][BGGroup.size()+4][2];
            //cut value,sample,expN/error
            
            //calculate sumOfEvent and Significance
            for(int u=0;u<uncertaintyN;u++)
            {
                for(int q=0;q<cutN;q++)
                {
                    //initialize histograms
                    TString title = Var[CountVarIndex].VarTitle;
                    
                    TString xaxis;
                    xaxis += Var[CountVarIndex].VarTitle;
                    xaxis += " ";
                    xaxis += Var[CountVarIndex].unit;
                    
                    TH1F* h2Sig;
                    
                    //Background
                    //For MC background
                    for(unsigned int j=0;j<tree2BGMC.size();j++)
                    {
                        BGGroup[j].h2 = new TH1F(BGGroup[j].info->GroupName.Data(),title.Data(),Var[CountVarIndex].bin,Var[CountVarIndex].xmin,Var[CountVarIndex].xmax);
                        BGGroup[j].h2->GetYaxis()->SetTitle("Number of events");
                        BGGroup[j].h2->SetLineColor(j+2);
                        BGGroup[j].h2->SetFillColor(j+2);
                        
                        for(unsigned int k=0;k<tree2BGMC[j].size();k++)
                        {
                            TH1F* hTemp = new TH1F("BGMC",title.Data(),Var[CountVarIndex].bin,Var[CountVarIndex].xmin,Var[CountVarIndex].xmax);
                            
                            //fill histograms from trees
                            TString temp = Var[CountVarIndex].VarName;
                            temp += ">>BGMC";
                            
                            TString Cut = "weight*(";
                            Cut += Var[CutVarIndex].VarName;
                            Cut += "<=";
                            Cut += TString::Itoa(jetptCut[q],10);;
                            Cut += ")";

                            tree2BGMC[j][k]->Draw(temp.Data(),Cut.Data());
                            
                            //normalization for BG
                            hTemp->Scale(BGMCGroupXS[j][k]/BGMCGroupnAOD[j][k] *sumDataL);
                            
                            BGGroup[j].h2->Add(hTemp);
                            delete hTemp;
                        }
                    }
                    
                    //For data-driven background
                    for(unsigned int j=tree2BGMC.size();j<BGGroup.size();j++)
                    {
                        BGGroup[j].h2 = new TH1F(BGGroup[j].info->GroupName.Data(),title.Data(),Var[CountVarIndex].bin,Var[CountVarIndex].xmin,Var[CountVarIndex].xmax);
                        BGGroup[j].h2->GetYaxis()->SetTitle("Number of events");
                        BGGroup[j].h2->SetLineColor(j+2);
                        BGGroup[j].h2->SetFillColor(j+2);
                        
                        for(unsigned int k=0;k<DataSampleID.size();k++)
                        {
                            TH1F* hTemp = new TH1F("BGData",title.Data(),Var[CountVarIndex].bin,Var[CountVarIndex].xmin,Var[CountVarIndex].xmax);
                            
                            //fill histograms from trees
                            TString temp = Var[CountVarIndex].VarName;
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
                            if(BGGroup[j].info->GroupName == "charge flip")
                            {
                                Cut += " && fLwt==0";
                            }
                            if(BGGroup[j].info->GroupName == "fake lepton")
                            {
                                Cut += " && fLwt!=0";
                            }
                            
                            Cut += " && ";
                            Cut += Var[CutVarIndex].VarName;
                            Cut += "<=";
                            Cut += TString::Itoa(jetptCut[q],10);;
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
                            BGGroup[j].h2->Add(hTemp);
                            delete hTemp;
                        }
                    }
                    
                    //h2Sig
                    {
                        h2Sig = new TH1F("Sig",title.Data(),Var[CountVarIndex].bin,Var[CountVarIndex].xmin,Var[CountVarIndex].xmax);

                        TString temp = Var[CountVarIndex].VarName;
                        temp += ">>Sig";
                        
                        TString Cut = "weight*(";
                        Cut += Var[CutVarIndex].VarName;
                        Cut += "<=";
                        Cut += TString::Itoa(jetptCut[q],10);;
                        Cut += ")";
                        tree2Sig[SigID]->Draw(temp.Data(),Cut.Data());
                        
                        //normalization for Signal
                        h2Sig->Scale(SigXS[SigID]/SignAOD[SigID] *sumDataL);
                    }
                    
                    //expected number of events
                    sumOfEvent[q][BGGroup.size()][0]=0;
                    sumOfEvent[q][BGGroup.size()][1]=0;
                    for(unsigned int j=0;j<BGGroup.size();j++)
                    {
                        //expected number of events for BG
                        sumOfEvent[q][j][0] = BGGroup[j].h2->IntegralAndError(0,-1,sumOfEvent[q][j][1]);
                        //cout<<BGGroup[j].info->GroupName.Data()<<": "<<sumOfEvent[q][j][0]<<" +/- "<<sumOfEvent[q][j][1]<<endl;
                        sumOfEvent[q][BGGroup.size()][0] += sumOfEvent[q][j][0];
                        sumOfEvent[q][BGGroup.size()][1] += sumOfEvent[q][j][1]*sumOfEvent[q][j][1];
                    }
                    //cout<<"Total BG: "<<sumOfEvent[q][BGGroup.size()][0]<<" +/- "<<TMath::Sqrt(sumOfEvent[q][BGGroup.size()][1])<<endl;
                    
                    //expected number of events for Signal
                    sumOfEvent[q][BGGroup.size()+2][0] = h2Sig->IntegralAndError(0,-1,sumOfEvent[q][BGGroup.size()+2][1]);
                    //cout<<"Signal: "<<sumOfEvent[q][BGGroup.size()+2][0]<<" +/- "<<sumOfEvent[q][BGGroup.size()+2][1]<<endl;
                    
                    //Significance
                    sumOfEvent[q][BGGroup.size()+3][0] = RooStats::NumberCountingUtils::BinomialExpZ(sumOfEvent[q][BGGroup.size()+2][0],sumOfEvent[q][BGGroup.size()][0],uncertainty[u]);
                    //cout<<"Significance: "<<sumOfEvent[q][BGGroup.size()+3][0]<<endl<<endl;
                    Significance[u][q] = sumOfEvent[q][BGGroup.size()+3][0];
                    //cout<<"Signal: "<<sumOfEvent[q][BGGroup.size()+2][0]<<", BG: "<<sumOfEvent[q][BGGroup.size()][0]<<", Significance: "<<sumOfEvent[q][BGGroup.size()+3][0]<<endl<<endl;
                    
                    //delete
                    //h2BGGruop
                    for(unsigned int j=0;j<BGGroup.size();j++)
                    {
                        delete BGGroup[j].h2;
                    }
                    
                    //h2Sig
                    delete h2Sig;
                }
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
                    delete tree2BGMC[j][k];
                }
            }
            for(unsigned int i=0;i<SigSampleID.size();i++)
            {
                delete tree2Sig[i];
            }
            
            //Significance graph
            TCanvas* c2 = new TCanvas();
            TString xaxis;
            xaxis += Var[CutVarIndex].VarTitle;
            xaxis += " ";
            xaxis += Var[CutVarIndex].unit;
            {
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
                if(min > 0)
                {
                    SignificanceGraph[0]->SetMinimum(min/1.1);
                }
                else
                {
                    SignificanceGraph[0]->SetMinimum(min*1.1);
                }
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
                double min = sumOfEvent[0][BGGroup.size()][0];
                double max = sumOfEvent[0][BGGroup.size()][0];
                double BGNumber[cutN];
                double SigNumber[cutN];
                for(int q=0;q<cutN;q++)
                {
                    BGNumber[q] = sumOfEvent[q][BGGroup.size()][0];
                    SigNumber[q] = sumOfEvent[q][BGGroup.size()+2][0];
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
    {
        unsigned int VVGroupIndex = 0;
        for(unsigned int i=0;i<BGMCGroupData.size();i++)
        {
            if(BGMCGroupData[i].GroupName == "VV") VVGroupIndex = i;
        }
        for(unsigned int i = BGMCGroupData[VVGroupIndex].lower;i <= BGMCGroupData[VVGroupIndex].upper;i++)
        {
            SampleData element;
            element.SampleName = BGMCSampleID[i].Data();
            element.SampleName.Remove(0,19);
            element.SampleName.ReplaceAll("_","\\_");
            element.index = i;
            BGVVData.push_back(element);
        }
    }
    
    //plot graph
    const bool optimize = 0;
    unsigned int countVariable = 0;
    for(unsigned int i=0;i<Var.size();i++)
    {
        if(Var[i].VarName == "l12_MET_dPhi") countVariable = i;
    }
    
    //expN for SR
    TString PathName_SR = "latex/data/expN/SR.txt";
    fstream fout_SR;
    fout_SR.open(PathName_SR.Data(), ios::out);
    fout_SR<<std::setw(46)<<"The Name of Signal Region";
    for(unsigned int i=0;i<BGMCGroupData.size();i++)
    {
        if(BGMCGroupData[i].GroupName == "VV") fout_SR<<std::setw(15)<<BGMCGroupData[i].GroupName.Data();
    }
    for(unsigned int i=0;i<BGMCGroupData.size();i++)
    {
        if(BGMCGroupData[i].GroupName == "Vgamma") fout_SR<<std::setw(15)<<BGMCGroupData[i].GroupName.Data();
    }
    
    for(unsigned int i=0;i<BGDataGroupData.size();i++)
    {
        if(BGDataGroupData[i].GroupName == "charge flip") fout_SR<<std::setw(15)<<BGDataGroupData[i].GroupName.Data();
    }
    for(unsigned int i=0;i<BGDataGroupData.size();i++)
    {
        if(BGDataGroupData[i].GroupName == "fake lepton") fout_SR<<std::setw(15)<<BGDataGroupData[i].GroupName.Data();
    }

    fout_SR<<std::setw(15)<<"Total BG";
    
    for(unsigned int i=0;i<SigMassSplitting.size();i++)
    {
        fout_SR<<std::setw(22)<<SigMassSplitting[i].IDName.Data();
    }
    fout_SR<<endl;
    fout_SR.close();
    
    //h2SR
    std::vector< std::vector<TH1F*> > h2SRBG;
    std::vector< std::vector<TH1F*> > h2SRSig;
    for(unsigned int i=0;i<4;i++)
    {
        std::vector<TH1F*> element;
        TH1F* hTemp = nullptr;
        TString LeptonName;
        if(i==0) LeptonName = "_ee";
        else if(i==1) LeptonName = "_mumu";
        else if(i==2) LeptonName = "_emu";
        else if(i==3) LeptonName = "_combine";
        
        //h2SRBG
        //VV
        for(unsigned int j=0;j<BGMCGroupData.size();j++)
        {
            if(BGMCGroupData[j].GroupName == "VV") hTemp = new TH1F((BGMCGroupData[j].GroupName+LeptonName).Data(),"",18,0,18);
        }
        hTemp->SetLineColor(2);
        hTemp->SetFillColor(2);
        element.push_back(hTemp);
        
        //Vgamma
        for(unsigned int j=0;j<BGMCGroupData.size();j++)
        {
            if(BGMCGroupData[j].GroupName == "Vgamma") hTemp = new TH1F((BGMCGroupData[j].GroupName+LeptonName).Data(),"",18,0,18);
        }
        hTemp->SetLineColor(4);
        hTemp->SetFillColor(4);
        element.push_back(hTemp);
        
        //charge flip
        if(i==0 || i==3)
        {
            for(unsigned int j=0;j<BGDataGroupData.size();j++)
            {
                if(BGDataGroupData[j].GroupName == "charge flip") hTemp = new TH1F((BGDataGroupData[j].GroupName+LeptonName).Data(),"",18,0,18);
            }
            hTemp->SetLineColor(3);
            hTemp->SetFillColor(3);
            element.push_back(hTemp);
        }
        
        //fake lepton
        for(unsigned int j=0;j<BGDataGroupData.size();j++)
        {
            if(BGDataGroupData[j].GroupName == "fake lepton") hTemp = new TH1F((BGDataGroupData[j].GroupName+LeptonName).Data(),"",18,0,18);
        }
        hTemp->SetLineColor(5);
        hTemp->SetFillColor(5);
        element.push_back(hTemp);
        
        h2SRBG.push_back(element);
        element.clear();
        
        //h2SRSig
        for(unsigned int j=0;j<SigMassSplitting.size();j++)
        {
            hTemp = new TH1F((SigMassSplitting[j].IDName+LeptonName).Data(),"",18,0,18);
            hTemp->SetLineColor(SigMassSplitting[j].colour);
            hTemp->SetLineStyle(1);
            element.push_back(hTemp);
        }
        h2SRSig.push_back(element);
    }
    
    for(unsigned int RegionIndex=6;RegionIndex<=6;RegionIndex++)
    //for(unsigned int RegionIndex=0;RegionIndex<=17;RegionIndex++)
    //for(unsigned int RegionIndex=18;RegionIndex<=71;RegionIndex++)
    //for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
    {
        std::vector<TChain*> tree2Data;
        initializeTree2(tree2Data,RegionInfo[RegionIndex].setOfChannel,DataSampleID,ChannelInfo);
        
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
                        initializeTree2(tree2BGMCElement,RegionInfo[RegionIndex].setOfChannel,BGMCGroupSampleID,ChannelInfo);
                        
                        tree2BGMC.push_back(tree2BGMCElement);
                        BGMCGroupXS.push_back(BGMCGroupXSElement);
                        BGMCGroupnAOD.push_back(BGMCGroupnAODElement);
                        
                        Group BGGroupElement;
                        BGGroupElement.info = &(BGMCGroupData[j]);
                        BGGroup.push_back(BGGroupElement);
                        
                        //Z pt reweighting
                        if(dorw
                           &&
                           (
                            RegionInfo[RegionIndex].RegionName == "nonISR_OS_ee"   ||
                            RegionInfo[RegionIndex].RegionName == "ISR_OS_ee"      ||
                            RegionInfo[RegionIndex].RegionName == "nonISR_OS_mumu" ||
                            RegionInfo[RegionIndex].RegionName == "ISR_OS_mumu"
                           )
                           &&
                           (
                            BGGroupElement.info->GroupName == "Zee" ||
                            BGGroupElement.info->GroupName == "Zmumu" ||
                            BGGroupElement.info->GroupName == "Ztautau"
                           )
                          )
                        {
                            for(unsigned int k=0;k<tree2BGMCElement.size();k++)
                            {
                                TString NameTemp = "tree_rw_";
                                NameTemp += TString::Itoa(k + BGGroupElement.info->lower,10);
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
                            RegionInfo[RegionIndex].RegionName == "CR_ISR_SS_ee"    ||
                            RegionInfo[RegionIndex].RegionName == "CR_SS_ee"
                           )
                           &&
                           BGGroupElement.info->GroupName == "Zee"
                          )
                        {
                            for(unsigned int k=0;k<tree2BGMCElement.size();k++)
                            {
                                TString NameTemp = "tree_cfw_";
                                NameTemp += TString::Itoa(k + BGGroupElement.info->lower,10);
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
        initializeTree2(tree2Sig,RegionInfo[RegionIndex].setOfChannel,SigSampleID,ChannelInfo);
        
        std::vector<TChain*> tree2DataOS;
        if(RegionInfo[RegionIndex].isSS_ee) initializeTree2(tree2DataOS,RegionInfo[RegionIndex].qFChannel,DataSampleID,ChannelInfo);
        
        //for(unsigned int VarIndex=4;VarIndex<=4;VarIndex++)
        //for(unsigned int VarIndex=countVariable;VarIndex<=countVariable;VarIndex++)
        for(unsigned int VarIndex=0;VarIndex<Var.size();VarIndex++)
        {
            //initialize histograms
            TString title = Var[VarIndex].VarTitle;
            
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
                //CommonCut += " && pt1>25 && pt2>20";
                //CommonCut += " && mll>60";
                
                //h2DataSum
                {
                    h2DataSum = new TH1F("DataSum",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
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
                    Cut += " && fLwt==0";
                    
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
                        
                        //Weight
                        TString Cut = "weight";
                        
                        //Z pt reweighting
                        if(dorw
                           &&
                            (
                             RegionInfo[RegionIndex].RegionName == "nonISR_OS_ee"   ||
                             RegionInfo[RegionIndex].RegionName == "ISR_OS_ee"      ||
                             RegionInfo[RegionIndex].RegionName == "nonISR_OS_mumu" ||
                             RegionInfo[RegionIndex].RegionName == "ISR_OS_mumu"
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
                            RegionInfo[RegionIndex].RegionName == "CR_ISR_SS_ee"    ||
                            RegionInfo[RegionIndex].RegionName == "CR_SS_ee"
                           )
                           &&
                           BGGroup[j].info->GroupName == "Zee"
                          )
                        {
                            Cut += "*cfw";
                        }
                        
                        //Cut
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
                    BGGroup[j].h2 = new TH1F(BGGroup[j].info->GroupName.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
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
                            Cut += " && fLwt==0";
                        }
                        if(BGGroup[j].info->GroupName == "fake lepton")
                        {
                            Cut += " && fLwt!=0";
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
                    double sigCount = tree2Sig[j]->Draw(temp.Data(),Cut.Data());
                    for(unsigned int i=0;i<SigMassSplitting.size();i++)
                    {
                        if(j==SigMassSplitting[i].ID) SigMassSplitting[i].statCount = sigCount;
                    }
                }
            }
            
            //Add Signal for the same mass splitting
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
                
                //preliminary normalization for h2SigSum
                h2SigSum[i]->Scale(sumDataL/AOD);
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
                    TH1F hTemp = TH1F("SignalTemp",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                    hTemp.Add(h2SigSum[i]);
                    hTemp.Scale(SigXS[SigMassSplitting[i].ID]);
                    sumOfEvent[BGGroup.size()+i+2][0] = hTemp.IntegralAndError(0,-1,sumOfEvent[BGGroup.size()+i+2][1]);
                    
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
                
                //output SR.txt file
                if(RegionIndex>=18)
                {
                    fout_SR.open(PathName_SR.Data(), ios::out | ios::app);
                    
                    fout_SR<<std::setw(46)<<RegionInfo[RegionIndex].RegionName;
                    
                    fout_SR<<setprecision(1)<<std::fixed;
                    for(unsigned int j=0;j<BGGroup.size();j++)
                    {
                        if(!RegionInfo[RegionIndex].isSS_ee && j==2) fout_SR<<std::setw(8)<<"0"<<std::setw(7)<<"";
                        
                        fout_SR<<std::setw(8)<<sumOfEvent[j][0];
                        fout_SR<<"+/-";
                        fout_SR<<std::setw(4)<<sumOfEvent[j][1];
                    }
                    
                    fout_SR<<std::setw(8)<<sumOfEvent[BGGroup.size()][0];
                    fout_SR<<"+/-";
                    fout_SR<<std::setw(4)<<TMath::Sqrt(sumOfEvent[BGGroup.size()][1]);
                    
                    for(unsigned int i=0;i<SigMassSplitting.size();i++)
                    {
                        fout_SR<<setprecision(1)<<std::fixed;
                        fout_SR<<std::setw(5)<<sumOfEvent[BGGroup.size()+i+2][0];
                        fout_SR<<"+/-";
                        fout_SR<<std::setw(3)<<sumOfEvent[BGGroup.size()+i+2][1];
                        
                        fout_SR<<setprecision(3)<<std::fixed;
                        fout_SR<<std::setw(7)<<sumOfEvent[BGGroup.size()+i+2][2];
                        fout_SR<<std::setw(4)<<SigMassSplitting[i].statCount;
                    }
                    fout_SR<<endl;
                    fout_SR.close();
                }
                
                //h2SR
                if(RegionIndex>=18)
                {
                    for(unsigned int j=0;j<BGGroup.size();j++)
                    {
                        h2SRBG[RegionIndex%3][j]->SetBinContent((RegionIndex-18)/3 +1,sumOfEvent[j][0]);
                    }
                    for(unsigned int j=0;j<SigMassSplitting.size();j++)
                    {
                        h2SRSig[RegionIndex%3][j]->SetBinContent((RegionIndex-18)/3 +1,sumOfEvent[BGGroup.size()+j+2][0]);
                    }
                }
                
                //significance for all mass point
                if(RegionIndex>=18)
                {
                    PathName = "latex/data/significance/";
                    PathName += RegionInfo[RegionIndex].RegionName;
                    PathName += ".txt";
                    fout.open(PathName.Data());
                    
                    TH2F* h2 = new TH2F("h2","h2;m_{C1/N2} [GeV];m_{N1} [GeV]",10,200,700,10,100,500);
                    TGraph2D* g2 = new TGraph2D();
                    int pointN = 0;
                    for(unsigned int j=0;j<SigSampleID.size();j++)
                    {
                        for(unsigned int k=0;k<SigMassSplitting.size();k++)
                        {
                            if(SigMass1[j]-SigMass2[j] == SigMassSplitting[k].MassDiff)
                            {
                                //expected number of events
                                double expN = h2SigSum[k]->Integral(0,-1);
                                expN *= SigXS[j];
                                
                                //Significance
                                double significance = RooStats::NumberCountingUtils::BinomialExpZ(expN,sumOfEvent[BGGroup.size()][0],0.3);
                                
                                fout<<setprecision(3)<<std::fixed;
                                fout<<SigMass1[j]<<" "<<SigMass2[j]<<" "<<significance<<endl;
                                
                                //2D plot
                                if(significance>0)
                                {
                                    g2->SetPoint(g2->GetN(), SigMass1[j], SigMass2[j], significance);
                                    pointN++;
                                }
                            }
                        }
                    }
                    fout.close();
                    
                    //if(pointN>=3)
                    if(RegionIndex==43)
                    {
                        gStyle->SetPalette(1);
                        gStyle->SetPadRightMargin(0.16);
                        TCanvas* c2 = new TCanvas();
                        c2->cd();
                        h2->Draw();
                        g2->Draw("colzsame");
                        
                        TLatex lt1;
                        lt1.DrawLatexNDC(0.2,0.9,RegionInfo[RegionIndex].RegionName.Data());
                        lt1.SetTextSize(lt1.GetTextSize()*0.3);
                        
                        TString NameTemp = "plot/";
                        NameTemp += "SR_";
                        NameTemp += RegionInfo[RegionIndex].RegionName;
                        NameTemp += ".eps";
                        c2->Print(NameTemp,"eps");
                        
                        delete c2;
                    }
                    
                    delete h2;
                    delete g2;
                }
            }
            
            //final normalization for h2SigSum for plotting
            const int SigScale = 10;
            for(unsigned int i=0;i<SigMassSplitting.size();i++)
            {
                h2SigSum[i]->Scale(SigXS[SigMassSplitting[i].ID] *SigScale);
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
                h2Ratio = new TH1F("Ratio",title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
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
                TString NameTemp = "#sqrt{#it{s}} = 13 TeV, ";
                NameTemp += TString::Itoa(sumDataL/1000,10);
                NameTemp += " fb^{-1}";
                lt2.DrawLatexNDC(0.3,0.83, NameTemp.Data());
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
                    RegionInfo[RegionIndex].RegionName == "nonISR_OS_ee"   ||
                    RegionInfo[RegionIndex].RegionName == "ISR_OS_ee"      ||
                    RegionInfo[RegionIndex].RegionName == "nonISR_OS_mumu" ||
                    RegionInfo[RegionIndex].RegionName == "ISR_OS_mumu"
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
                
                //for charge filp BG
                if(docfw
                   &&
                   (
                    RegionInfo[RegionIndex].RegionName == "CR_nonISR_SS_ee" ||
                    RegionInfo[RegionIndex].RegionName == "CR_ISR_SS_ee"    ||
                    RegionInfo[RegionIndex].RegionName == "CR_SS_ee"
                   )
                   &&
                   BGGroup[j].info->GroupName == "Zee"
                  )
                {
                    TString NameTemp = "tree_cfw_";
                    NameTemp += TString::Itoa(k + BGGroup[j].info->lower,10);
                    delete tree2BGMC[j][k]->GetFriend(NameTemp.Data());
                }
                
                delete tree2BGMC[j][k];
            }
        }
        for(unsigned int i=0;i<SigSampleID.size();i++)
        {
            delete tree2Sig[i];
        }
    }
    
    //h2SR
    THStack stackSR[4];
    
    //For combined SR
    for(unsigned int i=0;i<3;i++)
    {
        //h2SRBG
        h2SRBG[3][0]->Add(h2SRBG[i][0]);
        h2SRBG[3][1]->Add(h2SRBG[i][1]);
        if(i==0)
        {
            h2SRBG[3][2]->Add(h2SRBG[i][2]);
            h2SRBG[3][3]->Add(h2SRBG[i][3]);
        }
        else
        {
            h2SRBG[3][3]->Add(h2SRBG[i][2]);
        }
        
        //h2SRSig
        for(unsigned int j=0;j<SigMassSplitting.size();j++)
        {
            h2SRSig[3][j]->Add(h2SRSig[i][j]);
        }
    }
    
    for(unsigned int i=0;i<4;i++)
    {
        //Legend
        TLegend* leg;
        {
            Double_t xl1, yl1, xl2, yl2;
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
            
            for(unsigned int j=0;j<BGMCGroupData.size();j++)
            {
                if(BGMCGroupData[j].GroupName == "VV") leg->AddEntry(h2SRBG[i][0],BGMCGroupData[j].LegendName.Data(),"fl");
            }
            for(unsigned int j=0;j<BGMCGroupData.size();j++)
            {
                if(BGMCGroupData[j].GroupName == "Vgamma") leg->AddEntry(h2SRBG[i][1],BGMCGroupData[j].LegendName.Data(),"fl");
            }
            
            if(i==0 || i==3)
            {
                for(unsigned int j=0;j<BGDataGroupData.size();j++)
                {
                    if(BGDataGroupData[j].GroupName == "charge flip") leg->AddEntry(h2SRBG[i][2],BGDataGroupData[j].LegendName.Data(),"fl");
                }
                for(unsigned int j=0;j<BGDataGroupData.size();j++)
                {
                    if(BGDataGroupData[j].GroupName == "fake lepton") leg->AddEntry(h2SRBG[i][3],BGDataGroupData[j].LegendName.Data(),"fl");
                }
            }
            else
            {
                for(unsigned int j=0;j<BGDataGroupData.size();j++)
                {
                    if(BGDataGroupData[j].GroupName == "fake lepton") leg->AddEntry(h2SRBG[i][2],BGDataGroupData[j].LegendName.Data(),"fl");
                }
            }
            
            
            for(unsigned int j=0;j<SigMassSplitting.size();j++)
            {
                TString NameTemp = "";
                NameTemp += "(";
                NameTemp += TString::Itoa(SigMass1[SigMassSplitting[j].ID],10);
                NameTemp += ", ";
                NameTemp += TString::Itoa(SigMass2[SigMassSplitting[j].ID],10);
                NameTemp += ")";
                leg->AddEntry(h2SRSig[i][j],NameTemp.Data(),"l");
            }
        }
        
        //add h2SRBG
        if(i==0)
        {
            //ee
            stackSR[i].Add(h2SRBG[i][0]);
            stackSR[i].Add(h2SRBG[i][1]);
            stackSR[i].Add(h2SRBG[i][3]);
            stackSR[i].Add(h2SRBG[i][2]);
        }
        else if(i==1)
        {
            //mumu
            stackSR[i].Add(h2SRBG[i][1]);
            stackSR[i].Add(h2SRBG[i][2]);
            stackSR[i].Add(h2SRBG[i][0]);
        }
        else if(i==2)
        {
            //emu
            stackSR[i].Add(h2SRBG[i][1]);
            stackSR[i].Add(h2SRBG[i][0]);
            stackSR[i].Add(h2SRBG[i][2]);
        }
        else if(i==3)
        {
            //combine
            stackSR[i].Add(h2SRBG[i][1]);
            stackSR[i].Add(h2SRBG[i][0]);
            stackSR[i].Add(h2SRBG[i][2]);
            stackSR[i].Add(h2SRBG[i][3]);
        }
        
        //Draw
        gStyle->SetOptStat(0);
        TCanvas* c2 = new TCanvas();
        c2->cd();
        c2->SetLogy(1);
        stackSR[i].SetMinimum(1);
        stackSR[i].Draw("hist");
        for(unsigned int j=0;j<SigMassSplitting.size();j++)
        {
            h2SRSig[i][j]->Draw("histsame");
        }
        leg->Draw();
        
        //export histograms in eps format
        TString NameTemp = "plot/";
        NameTemp += "SR_";
        if(i==0) NameTemp += "ee";
        else if(i==1) NameTemp += "mumu";
        else if(i==2) NameTemp += "emu";
        else if(i==3) NameTemp += "combine";
        NameTemp += ".eps";
        c2->Print(NameTemp,"eps");
        delete c2;
    }
    
    for(unsigned int i=0;i<4;i++)
    {
        for(unsigned int j=0;j<h2SRBG[i].size();j++)
        {
            delete h2SRBG[i][j];
        }
        for(unsigned int j=0;j<h2SRSig[i].size();j++)
        {
            delete h2SRSig[i][j];
        }
    }
    
    //latex for tables
    //latex for expN
    for(unsigned int sign=0;sign<=3;sign+=3)
    {
        TString PathName = "latex/data/";
        PathName += "expN_";
        if(sign==0) PathName += "OS";
        else PathName += "SS";
        PathName += ".tex";
        
        ofstream fout;
        fout.open(PathName.Data());
        
        for(unsigned int SixChannel=0;SixChannel<9;SixChannel++)
        {
            if(SixChannel>=3 && SixChannel<=5) continue;
            
            fout<<"\\begin{frame}"<<endl;
            fout<<"\\frametitle{Expected number of events (For ";
            if(sign==0) fout<<"opposite sign";
            else fout<<"same sign";
            fout<<")}"<<endl;
            
            fout<<"For ";
            if(SixChannel == 0 || SixChannel == 6) fout<<"ee";
            else if(SixChannel == 1 || SixChannel == 7) fout<<"$\\mu\\mu$";
            else if(SixChannel == 2 || SixChannel == 8) fout<<"e$\\mu$";
            fout<<" channel, ";
            if(SixChannel<=2) fout<<"non-ISR";
            else fout<<"ISR";
            fout<<"\\\\"<<endl;
            
            fout<<"\\vspace{5mm}"<<endl;
            fout<<"\\begin{tabular}{|c|c|c|}"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"& Number of events & Significance \\\\"<<endl;
            fout<<"\\hline"<<endl;
            
            fout<<"\\input{data/expN/"<<ChannelInfo[sign+SixChannel].ChannelName.Data()<<".tex}"<<endl;
            
            fout<<"\\end{tabular}"<<endl;
            fout<<"\\end{frame}"<<endl<<endl;
        }
        
        fout.close();
        
        //For VV
        PathName = "latex/data/";
        PathName += "expN_BGVV_";
        if(sign==0) PathName += "OS";
        else PathName += "SS";
        PathName += ".tex";
        fout.open(PathName.Data());
        
        for(unsigned int SixChannel=0;SixChannel<9;SixChannel++)
        {
            if(SixChannel>=3 && SixChannel<=5) continue;
            
            fout<<"\\begin{frame}"<<endl;
            fout<<"\\frametitle{Expected number of events for VV (For ";
            if(sign==0) fout<<"opposite sign";
            else fout<<"same sign";
            fout<<")}"<<endl;
            
            fout<<"For ";
            if(SixChannel == 0 || SixChannel == 6) fout<<"ee";
            else if(SixChannel == 1 || SixChannel == 7) fout<<"$\\mu\\mu$";
            else if(SixChannel == 2 || SixChannel == 8) fout<<"e$\\mu$";
            fout<<" channel, ";
            if(SixChannel<=2) fout<<"non-ISR";
            else fout<<"ISR";
            fout<<"\\\\"<<endl;
            
            fout<<"\\vspace{5mm}"<<endl;
            fout<<"\\begin{tabular}{|c|c|}"<<endl;
            fout<<"\\hline"<<endl;
            fout<<"& Number of events \\\\"<<endl;
            fout<<"\\hline"<<endl;
            
            fout<<"\\input{data/expN/BGVV_"<<ChannelInfo[sign+SixChannel].ChannelName.Data()<<".tex}"<<endl;
            
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
    //plot_OS.tex and plot_SS.tex
    for(unsigned int sign=0;sign<=3;sign+=3)
    {
        TString PathName = "latex/data/";
        PathName += "plot_";
        if(sign==0) PathName += "OS";
        else PathName += "SS";
        PathName += ".tex";
        
        ofstream fout;
        fout.open(PathName.Data());
        
        for(unsigned int VarIndex=0;VarIndex<Var.size();VarIndex++)
        {
            if(Var[VarIndex].VarName=="averageMu") continue;
            if(Var[VarIndex].VarName=="nVtx") continue;

            fout<<"\\begin{frame}"<<endl;
            
            fout<<"\\frametitle{"<<Var[VarIndex].VarTitle.Data()<<" (For ";
            if(sign==0) fout<<"opposite sign";
            else fout<<"same sign";
            fout<<")}"<<endl;
            
            fout<<"\\Wider[5em]{"<<endl;
            for(unsigned int SixChannel=0;SixChannel<9;SixChannel++)
            {
                if(SixChannel>=3 && SixChannel<=5) continue;
                
                if(SixChannel<=2 &&
                   (Var[VarIndex].VarName=="bjetpt"  ||
                    Var[VarIndex].VarName=="bjeteta" ||
                    Var[VarIndex].VarName=="bjetphi" ||
                    Var[VarIndex].VarName=="cjetpt"  ||
                    Var[VarIndex].VarName=="cjeteta" ||
                    Var[VarIndex].VarName=="cjetphi" )
                   ) continue;
                
                fout<<"\\includegraphics[width=0.33\\textwidth]{\\PathToPlot/"
                    <<Var[VarIndex].VarName.Data()<<"_"<<ChannelInfo[sign+SixChannel].ChannelName.Data()<<"}";
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
            for(unsigned int sign=0;sign<=3;sign+=3)
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
                
                fout<<"\\begin{frame}"<<endl;
                
                fout<<"\\frametitle{"<<Var[VarIndex].VarTitle.Data()<<" (For ";
                if(sign==0) fout<<"opposite sign";
                else fout<<"same sign";
                fout<<")}"<<endl;
                
                fout<<"\\Wider[5em]{"<<endl;
                for(unsigned int SixChannel=0;SixChannel<9;SixChannel++)
                {
                    if(SixChannel>=3 && SixChannel<=5) continue;
                    
                    if(SixChannel<=2 &&
                       (Var[VarIndex].VarName=="bjetpt"  ||
                        Var[VarIndex].VarName=="bjeteta" ||
                        Var[VarIndex].VarName=="bjetphi" ||
                        Var[VarIndex].VarName=="cjetpt"  ||
                        Var[VarIndex].VarName=="cjeteta" ||
                        Var[VarIndex].VarName=="cjetphi" )
                       ) continue;
                    
                    fout<<"\\includegraphics[width=0.33\\textwidth]{\\PathToPlot/"
                    <<Var[VarIndex].VarName.Data()<<"_"<<ChannelInfo[sign+SixChannel].ChannelName.Data()<<"}";
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
