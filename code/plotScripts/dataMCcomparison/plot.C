#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TString.h>
#include "AtlasLabels.C"
#include "AtlasStyle.C"

struct ChannelData
{
    TString ChannelName;
};

void initializeTree2(std::vector<TChain*>& tree2,std::vector<unsigned int>& SetOfChannel, std::vector<TString>& SampleID, std::vector<ChannelData>& ChannelInfo)
{
    for(unsigned int i=0;i<SampleID.size();i++)
    {
        TChain* treeTemp = new TChain("tree");
        for(unsigned int j=0;j<SetOfChannel.size();j++)
        {
            TString fileName = SampleID[i];
            fileName += "_";
            fileName += ChannelInfo[SetOfChannel[j]].ChannelName;
            fileName += ".root";
            
            treeTemp->Add(fileName.Data());
        }
        
        tree2.push_back(treeTemp);
    }
}

void plot()
{
    double sumDataL = 36470; //cross section in pb-1
    
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
                    
                    ChannelInfo.push_back(element);
                }
            }
        }
    }
    
    //For Signal MC
    struct SigInfo
    {
        int MassDiff;
        unsigned int ID;
        TString IDName;
    };
    std::vector<SigInfo> SigMassSplitting;
    {
        SigInfo element;
        
        element.MassDiff = 20;    element.ID = 2;    SigMassSplitting.push_back(element);
        element.MassDiff = 50;    element.ID = 18;   SigMassSplitting.push_back(element);
        element.MassDiff = 100;   element.ID = 19;   SigMassSplitting.push_back(element);
        element.MassDiff = 200;   element.ID = 20;   SigMassSplitting.push_back(element);
        element.MassDiff = 300;   element.ID = 21;   SigMassSplitting.push_back(element);
    }
    
    std::vector<TString> SigSampleID1;
    std::vector<TString> SigSampleID2;
    std::vector<double> SigXS1; //cross section in pb
    std::vector<double> SigXS2; //cross section in pb
    std::vector<int> SigMass1;
    std::vector<int> SigMass2;
    SigSampleID1.reserve(20);
    SigSampleID2.reserve(20);
    SigXS1.reserve(20);
    SigXS2.reserve(20);
    SigMass1.reserve(20);
    SigMass2.reserve(20);
    
    {
        //read SigSample.txt
        ifstream fin;
        fin.open("SigSample_diff.txt");
        
        //old
        for(unsigned int i=1;i<=66;i++)
        {
            TString SampleIDTemp;
            fin>>SampleIDTemp;
            
            TString SampleNameTemp;
            fin>>SampleNameTemp;
            
            SampleIDTemp += ".";
            SampleIDTemp += SampleNameTemp;
            
            SigSampleID1.push_back("skimming_signal_old/skimming."+SampleIDTemp);
            SigSampleID2.push_back("skimming_signal_new/skimming."+SampleIDTemp);
            
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
            
            //cout<<SigXSTemp2<<endl;
            SigXS1.push_back(SigXSTemp2);
        }
        
        //new
        for(unsigned int i=1;i<=66;i++)
        {
            TString SampleIDTemp;
            fin>>SampleIDTemp;
            fin>>SampleIDTemp;
            
            int SigMass;
            fin>>SigMass;
            fin>>SigMass;
            
            double SigXSTemp2;
            fin>>SigXSTemp2;
            
            double SigXSTemp;
            fin>>SigXSTemp;
            SigXSTemp2 *= SigXSTemp;
            
            //cout<<SigXSTemp2<<endl;
            SigXS2.push_back(SigXSTemp2);
        }
        fin.close();
    }
    
    for(unsigned int i=0;i<SigMassSplitting.size();i++)
    {
        cout<<"Mass splitting: "<<SigMassSplitting[i].MassDiff<<endl;
        cout<<"index:"<<SigMassSplitting[i].ID<<endl;
        cout<<"Name: "<<SigSampleID1[SigMassSplitting[i].ID].Data()<<endl;
        cout<<"MassDiff: "<<SigMass1[SigMassSplitting[i].ID]-SigMass2[SigMassSplitting[i].ID]<<endl;
        cout<<"XS: "<<SigXS1[SigMassSplitting[i].ID]<<endl<<endl;
        
        cout<<"All samples with the same mass splitting "<<SigMassSplitting[i].MassDiff<<" GeV:"<<endl;
        for(unsigned int j=0;j<SigSampleID1.size();j++)
        {
            if(SigMass1[j]-SigMass2[j] == SigMassSplitting[i].MassDiff) cout<<"index:"<<j<<" "<<SigSampleID1[j]<<endl;
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
    for(unsigned int i=0;i<SigSampleID1.size();i++)
    {
        TString NameTemp = SigSampleID1[i];
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
    };
    
    std::vector<RegionData> RegionInfo;
    {
        RegionData element;
        
        for(unsigned int ChannelIndex=0;ChannelIndex<ChannelInfo.size();ChannelIndex++)
        {
            element.RegionName = ChannelInfo[ChannelIndex].ChannelName;
            element.setOfChannel.clear();
            element.setOfChannel.push_back(ChannelIndex);

            RegionInfo.push_back(element);
        }
    }
    
    //plot graph
    unsigned int countVariable = 0;
    for(unsigned int i=0;i<Var.size();i++)
    {
        if(Var[i].VarName == "l12_MET_dPhi") countVariable = i;
    }
    
    for(unsigned int RegionIndex=0;RegionIndex<RegionInfo.size();RegionIndex++)
    {
        std::vector<TChain*> tree2Sig;
        initializeTree2(tree2Sig,RegionInfo[RegionIndex].setOfChannel,SigSampleID1,ChannelInfo);
        
        //for(unsigned int VarIndex=4;VarIndex<=4;VarIndex++)
        for(unsigned int VarIndex=countVariable;VarIndex<=countVariable;VarIndex++)
        //for(unsigned int VarIndex=0;VarIndex<Var.size();VarIndex++)
        {
            //initialize histograms
            TString title = Var[VarIndex].VarTitle;
            
            TString xaxis;
            xaxis += Var[VarIndex].VarTitle;
            xaxis += " ";
            xaxis += Var[VarIndex].unit;
            
            TH1F* h2SigSum[SigMassSplitting.size()];
            TH1F* h2Sig[SigSampleID1.size()];
            
            //h2SigSum
            for(unsigned int j=0;j<SigMassSplitting.size();j++)
            {
                TString NameTemp = "SigSum_";
                NameTemp += TString::Itoa(SigMassSplitting[j].MassDiff,10);
                h2SigSum[j] = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                h2SigSum[j]->GetYaxis()->SetTitle("Number of events");
                h2SigSum[j]->SetLineColor(2);
                h2SigSum[j]->SetLineStyle(1);
            }
            
            //h2Sig
            for(unsigned int j=0;j<SigSampleID1.size();j++)
            {
                TString NameTemp = "Sig_";
                NameTemp += TString::Itoa(j,10);
                h2Sig[j] = new TH1F(NameTemp.Data(),title.Data(),Var[VarIndex].bin,Var[VarIndex].xmin,Var[VarIndex].xmax);
                
                //Fill Signal
                TString temp = Var[VarIndex].VarName;
                temp += ">>";
                temp += NameTemp;
                
                TString Cut = "weight";
                
                tree2Sig[j]->Draw(temp.Data(),Cut.Data());
            }
            
            //Add Signal for the same mass splitting
            for(unsigned int i=0;i<SigMassSplitting.size();i++)
            {
                unsigned int AOD = 0;
                for(unsigned int j=0;j<SigSampleID1.size();j++)
                {
                    if(SigMass1[j]-SigMass2[j] == SigMassSplitting[i].MassDiff)
                    {
                        h2SigSum[i]->Add(h2Sig[j]);
                        AOD += SignAOD[j];
                    }
                }
                
                //normalization for h2SigSum
                h2SigSum[i]->Scale(sumDataL/AOD *SigXS1[SigMassSplitting[i].ID]);
            }
            
            for(unsigned int j=0;j<SigSampleID1.size();j++)
            {
                delete h2Sig[j];
            }
            
            if(VarIndex==countVariable)
            {
                double sumOfEvent[SigMassSplitting.size()][2];
                //sample,expN/error/significance
                
                //expected number of events for signal
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    //expected number of events
                    sumOfEvent[i][0] = h2SigSum[i]->IntegralAndError(0,-1,sumOfEvent[i][1]);
                    
                    cout<<"Signal ("<<SigMass1[SigMassSplitting[i].ID]<<", "<<SigMass2[SigMassSplitting[i].ID]<<"): "<<sumOfEvent[i][0]<<" +/- "<<sumOfEvent[i][1]<<endl;
                }
                cout<<endl;
                
                //output file
                TString PathName = "latex/data/expN_signal/";
                PathName += RegionInfo[RegionIndex].RegionName;
                PathName += ".tex";
                
                ofstream fout;
                fout.open(PathName.Data());
                fout<<setprecision(1)<<std::fixed;
 
                for(unsigned int i=0;i<SigMassSplitting.size();i++)
                {
                    fout<<"Signal (";
                    fout<<SigMass1[SigMassSplitting[i].ID];
                    fout<<", ";
                    fout<<SigMass2[SigMassSplitting[i].ID];
                    fout<<") & $";
                    fout<<setprecision(1)<<std::fixed;
                    fout<<sumOfEvent[i][0];
                    fout<<"\\pm";
                    fout<<sumOfEvent[i][1];
                    fout<<"$ &";
                    fout<<"\\\\"<<endl<<"\\hline"<<endl;
                }
                fout.close();
            }
            
            
            for(unsigned int i=0;i<SigMassSplitting.size();i++)
            {
                //Legend
                TLegend* leg;
                {
                    Double_t xl1, yl1, xl2, yl2;
                    
                    xl2=0.92;
                    yl2=0.95;
                    xl1=xl2-0.3;
                    yl1=yl2-0.2;
                    
                    leg = new TLegend(xl1,yl1,xl2,yl2);
                    leg->SetNColumns(2);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(42);
                    leg->SetBorderSize(0);
                    
                    
                    TString NameTemp = "";
                    NameTemp += "(";
                    NameTemp += TString::Itoa(SigMass1[SigMassSplitting[i].ID],10);
                    NameTemp += ", ";
                    NameTemp += TString::Itoa(SigMass2[SigMassSplitting[i].ID],10);
                    NameTemp += ")";
                    leg->AddEntry(h2SigSum[i],NameTemp.Data(),"l");
                    
                }
                
                TH1F* h2Ratio = nullptr;
                TPad* pad1 = nullptr;
                TPad* pad2 = nullptr;
                
                //x-asix title for pad1
                h2SigSum[i]->GetXaxis()->SetTitle(xaxis.Data());
                
                //pad1
                pad1 = new TPad("pad1","pad1",0,0,1,1);
                
                pad1->SetLogy(Var[VarIndex].log);
                
                //Draw for pad1
                gStyle->SetOptStat(0);
                TCanvas* c2 = new TCanvas();
                c2->cd();
                pad1->Draw();
                pad1->cd();
                
                h2SigSum[i]->Draw("axis");
                h2SigSum[i]->Draw("histsame");
                
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
                
                {
                    //export histograms in eps format
                    TString NameTemp = "plot_signal/";
                    NameTemp += Var[VarIndex].VarName;
                    NameTemp += "_";
                    NameTemp += RegionInfo[RegionIndex].RegionName;
                    NameTemp += "_";
                    NameTemp += SigMassSplitting[i].IDName;
                    NameTemp += ".eps";
                    c2->Print(NameTemp,"eps");
                }
                
                //legend
                delete leg;
                
                //canvas
                delete c2;
            }
            
            //h2SigSum
            for(unsigned int j=0;j<SigMassSplitting.size();j++)
            {
                delete h2SigSum[j];
            }
        }
        
        for(unsigned int i=0;i<SigSampleID1.size();i++)
        {
            delete tree2Sig[i];
        }
    }

}

