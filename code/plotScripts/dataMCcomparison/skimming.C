#include <iostream>
#include <fstream>
#include <vector>

#include "obj_def.h"
#include <TInterpreter.h>

#include <TFile.h>
#include <TChain.h>
#include <TH1F.h>

Int_t ID1;
Int_t ID2;
float pt1;
float pt2;
float eta1;
float eta2;
float phi1;
float mll;
float ptll;
float MET;
float METRel;
float mTtwo;
float mt1;
float mt2;
float mtm;

float l12_dPhi;
float l12_MET_dPhi;
float jet0_MET_dPhi;
float l12_jet0_dPhi;

Int_t nJet;
float jetpt;
float jeteta;
float jetphi;
Int_t nBJet;
float bjetpt;
float bjeteta;
float bjetphi;
Int_t nCJet;
float cjetpt;
float cjeteta;
float cjetphi;
Int_t nFJet;
float fjetpt;
float fjeteta;
float fjetphi;

float weight;
float qFwt;
float fLwt;
float averageMu;

float meff;
float mlj;
float R2;
float mjj;

struct nEvent
{
    TString name;
    int n;
    double nw;
    double nAOD;
};
void skimming2(TString const& SamplePath,TString const& tag,TString const& SampleName,std::vector<nEvent>& nSS)
{
    //get the "evt2l"
    TChain *tree1 = new TChain("evt2l");
    {
        /*
        TString fileName = SamplePath;
        fileName += "user.clo.";
        fileName += tag;
        fileName += ".";
        fileName += SampleName;*/
        //fileName += "_myOutput.root/*.root*";

        ///*
        TString fileName = SamplePath;
        fileName += "all/user.*.";
        fileName += SampleName;
        fileName += ".*.myOutput.root";
        //*/
        
        /*
        TString fileName = SamplePath;
        fileName += SampleName;
        fileName += "*.root";
        //fileName += ".merge.DAOD_SUSY2.e3836_s2726_r7772_r7676_p2879.root";
        */
        
        //TString fileName = "/Users/samuel/Atlas/ntuple/test.root";
        //TString fileName = "~/test.root";
        
        cout<<fileName.Data()<<endl;
        tree1->Add(fileName.Data());
    }
    
    EVT evt;
    tree1->SetBranchAddress("evt", &evt);

    vector< L_PAR > leps;
    vector< L_PAR >* vleps = &leps;
    tree1->SetBranchAddress("leps", &vleps);

    R_PAR l12;
    tree1->SetBranchAddress("l12", &l12);

    vector< J_PAR > jets;
    vector< J_PAR >* vjets = &jets;
    tree1->SetBranchAddress("jets", &vjets);

    /*
    vector< TR_PAR > truths;
    vector< TR_PAR >* vtruths = &truths;
    tree1->SetBranchAddress("truths", &vtruths);
    */

    SIGNATURE sig;
    tree1->SetBranchAddress("sig", &sig);

    /*
    tree1->Show(0);
    tree1->GetEntry(0);
    cout << "weight = " << evt.weight << endl;
    cout << "the size of the leptons is " << leps.size() << endl;
    cout << "the lepton 0 pt= " << leps[0].pt << endl;
    cout << "l12.m = " << l12.m << endl;
    cout << "the size of the jets is " << jets.size() << endl;
    cout << "the jet 0 pt= " << jets[0].pt << endl;
    //cout << "the size of the truths is " << truths.size() << endl;
    //cout << "the truth lepton 0 pt= " << truths[0].pt << endl;
    cout << "sig.mT2 = " << sig.mT2 << endl;
    return;
    */

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
    
    //create files for storing the resulting trees and histograms
    TFile *f2[channel.size()];
    for(unsigned int j=0;j<channel.size();j++)
    {
        TString fileName = "skimming/skimming.";
        //TString fileName = "skimming_signal_old/skimming.";
        //TString fileName = "skimming_signal_new/skimming.";
        fileName += SampleName;
        fileName += "_";
        fileName += channel[j];
        fileName += ".root";
        f2[j] = new TFile(fileName.Data(),"RECREATE");
    }
    
    //trees
    TTree *tree2[channel.size()];
    for(unsigned int j=0;j<channel.size();j++)
    {
        TString treeName = "tree_";
        treeName += channel[j];
        
        f2[j]->cd();
        tree2[j] = new TTree(treeName.Data(),treeName.Data());
        tree2[j]->Branch("pt1",&pt1,"pt1/F");
        tree2[j]->Branch("pt2",&pt2,"pt2/F");
        tree2[j]->Branch("eta1",&eta1,"eta1/F");
        tree2[j]->Branch("eta2",&eta2,"eta2/F");
        tree2[j]->Branch("phi1",&phi1,"phi1/F");
        tree2[j]->Branch("mll",&mll,"mll/F");
        tree2[j]->Branch("ptll",&ptll,"ptll/F");
        tree2[j]->Branch("MET",&MET,"MET/F");
        tree2[j]->Branch("METRel",&METRel,"METRel/F");
        tree2[j]->Branch("mTtwo",&mTtwo,"mTtwo/F");
        tree2[j]->Branch("mt1",&mt1,"mt1/F");
        tree2[j]->Branch("mt2",&mt2,"mt2/F");
        tree2[j]->Branch("mtm",&mtm,"mtm/F");
        tree2[j]->Branch("meff",&meff,"meff/F");
        tree2[j]->Branch("mlj",&mlj,"mlj/F");
        tree2[j]->Branch("mjj",&mjj,"mjj/F");

        tree2[j]->Branch("l12_dPhi",&l12_dPhi,"l12_dPhi/F");
        tree2[j]->Branch("l12_MET_dPhi",&l12_MET_dPhi,"l12_MET_dPhi/F");
        tree2[j]->Branch("jet0_MET_dPhi",&jet0_MET_dPhi,"jet0_MET_dPhi/F");
        tree2[j]->Branch("l12_jet0_dPhi",&l12_jet0_dPhi,"l12_jet0_dPhi/F");

        tree2[j]->Branch("nJet",&nJet,"nJet/I");
        tree2[j]->Branch("jetpt",&jetpt,"jetpt/F");
        tree2[j]->Branch("jeteta",&jeteta,"jeteta/F");
        tree2[j]->Branch("nBJet",&nBJet,"nBJet/I");
        tree2[j]->Branch("nCJet",&nCJet,"nCJet/I");

        tree2[j]->Branch("weight",&weight,"weight/F");
        tree2[j]->Branch("qFwt",&qFwt,"qFwt/F");
        tree2[j]->Branch("fLwt",&fLwt,"fLwt/F");
        tree2[j]->Branch("averageMu",&averageMu,"averageMu/F");
    }
    
    //histograms
    TH1F* h2[channel.size()];
    for(unsigned int j=0;j<channel.size();j++)
    {
        TString hName2 = "hist_";
        hName2 += channel[j];
        
        f2[j]->cd();
        h2[j] = new TH1F(hName2.Data(), "cut flow", 20, 0, 20);
        h2[j]->GetXaxis()->SetBinLabel(1,"nAOD");
        h2[j]->GetXaxis()->SetBinLabel(2,"nwAOD");
        h2[j]->GetXaxis()->SetBinLabel(3,"ntuple");
        h2[j]->GetXaxis()->SetBinLabel(4,"trigger");
        h2[j]->GetXaxis()->SetBinLabel(5,"=2SigLep");
        h2[j]->GetXaxis()->SetBinLabel(6,"fake");
        h2[j]->GetXaxis()->SetBinLabel(7,"pt1");
        h2[j]->GetXaxis()->SetBinLabel(8,"pt2");
        h2[j]->GetXaxis()->SetBinLabel(9,channel[j].Data());
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
                h2[j]->Fill("nAOD",h1->GetBinContent(1));
                h2[j]->Fill("nwAOD",h1->GetBinContent(2));
                h2[j]->Fill("ntuple",tree1->GetEntries());
            }
            delete f1;
        }
    }
    
    nEvent element;
    element.name = SampleName;
    element.n = 0;
    element.nw = 0;
    element.nAOD = h2[0]->GetBinContent(2);
    
    //loop over all entries
    //for(int j=0;j<=10;j++)
    for(int j=0;j<tree1->GetEntries();j++)
    {
        if(j%100000==0)
        {
            cout<<"number of event: " <<j<<endl;
        }
        tree1->GetEntry(j);
        
        //trigger
        if((sig.trigCode & sig.trigMask)==0)
        {
            continue;
        }
        for(unsigned int m=0;m<channel.size();m++)
        {
            h2[m]->Fill("trigger",1);
        }
        
        fLwt = evt.fLwt;
        
        if(!evt.isMC && fLwt!=0)
        {
            for(unsigned int m=0;m<channel.size();m++)
            {
                h2[m]->Fill("fake",1);
            }
        }
        else
        {
            //exact 2 signal leptons
            if(leps.size()!=2) continue;
            if(!(leps[0].lFlag & 1<<1)) continue;
            if(!(leps[1].lFlag & 1<<1)) continue;
            
            for(unsigned int m=0;m<channel.size();m++)
            {
                h2[m]->Fill("=2SigLep",1);
            }
        }
        
        int sigIndex[2];
        sigIndex[0] = 0;
        sigIndex[1] = 1;
        ID1 = leps[sigIndex[0]].ID;
        ID2 = leps[sigIndex[1]].ID;
        pt1 = leps[sigIndex[0]].pt;
        pt2 = leps[sigIndex[1]].pt;
        eta1 = leps[sigIndex[0]].eta;
        eta2 = leps[sigIndex[1]].eta;
        //pt of leading lepton
        //if(int(abs(ID1)/1000) == 11 && !(pt1>25 && fabs(eta1)<2.47)) continue;
        //if(int(abs(ID1)/1000) == 13 && !(pt1>20 && fabs(eta1)<2.4)) continue;
        for(unsigned int m=0;m<channel.size();m++)
        {
            h2[m]->Fill("pt1",1);
        }
        
        //pt of subleading lepton
        //if(int(abs(ID2)/1000) == 11 && !(pt2>15 && fabs(eta2)<2.47)) continue;
        //if(int(abs(ID2)/1000) == 13 && !(pt2>10 && fabs(eta2)<2.4)) continue;
        for(unsigned int m=0;m<channel.size();m++)
        {
            h2[m]->Fill("pt2",1);
        }
        
        phi1 = leps[sigIndex[0]].phi;
        mll = l12.m;
        ptll = l12.pt;
        MET = sig.Met;
        METRel = sig.MetRel;
        mTtwo = sig.mT2;
        mt1 = leps[sigIndex[0]].mT;
        mt2 = leps[sigIndex[1]].mT;
        if(mt1>mt2) mtm = mt1;
        else mtm = mt2;
        meff = sig.HT + sig.Met;
        mlj = sig.mlj;
        mjj = sig.mjj;
        //R2 = MET/(MET + pt1 + pt2);
        l12_dPhi = l12.dPhi;
        l12_MET_dPhi = l12.MET_dPhi;
        if(jets.size()>=1) jet0_MET_dPhi = jets[0].MET_dPhi;
        else jet0_MET_dPhi = -999;
        l12_jet0_dPhi = l12.jet0_dPhi;
        weight = evt.weight * evt.pwt * evt.ElSF * evt.MuSF * evt.BtagSF * evt.JvtSF;
        
        qFwt = evt.qFwt;
        
        averageMu = evt.averageMu;
        
        //SF or DF
        int channelIndex = 0;
        if(int(abs(ID1)/1000) == 11)
        {
            if(int(abs(ID2)/1000) == 11)
            {
            }
            else if(int(abs(ID2)/1000) == 13)
            {
                channelIndex += 2;
            }
        }
        else if(int(abs(ID1)/1000) == 13)
        {
            if(int(abs(ID2)/1000) == 11)
            {
                channelIndex += 2;
            }
            else if(int(abs(ID2)/1000) == 13)
            {
                channelIndex += 1;
            }
        }
        
        //jets
        nJet = jets.size();
        nBJet = 0;
        nCJet = 0;
        //nFJet = 0;
        int nISR = 0;
        int leadingJetIndex = 0;
        //int leadingBJetIndex = 0;
        //int leadingCJetIndex = 0;
        //int leadingFJetIndex = 0;
        for(unsigned int k=0;k<jets.size();k++)
        {
            //signal jets
            
            //B-jets
            //if((jets[k].jFlag & 1<<5) && jets[k].pt > 20)
            if(jets[k].jFlag & 1<<5)
            {
                nBJet++;
                //if(nBJet==1) leadingBJetIndex = k;
            }
            
            if(fabs(jets[k].eta) < 2.4)
            {
                //ISR
                if(jets[k].pt > 40) nISR++;
                
                
                //Central jets
                if(!(jets[k].jFlag & 1<<5))
                {
                    //no b-tagged
                    if(jets[k].pt > 20)
                    {
                        //Central light jets
                        nCJet++;
                        //if(nCJet==1) leadingCJetIndex = k;
                    }
                }
                
            }
            else
            {
                /*
                //Forward jets
                if(jets[k].pt > 30)
                {
                    nFJet++;
                    if(nFJet==1) leadingFJetIndex = k;
                }
                */
            }
        }
        
        if(nJet>0)
        {
            jetpt = jets[leadingJetIndex].pt;
            jeteta = jets[leadingJetIndex].eta;
            //jetphi = jets[leadingJetIndex].phi;
        }
        else
        {
            jetpt = 0;
            jeteta = 0;
            //jetphi = 0;
        }
       
        /* 
        if(nBJet>0)
        {
            bjetpt = jets[leadingBJetIndex].pt;
            bjeteta = jets[leadingBJetIndex].eta;
            bjetphi = jets[leadingBJetIndex].phi;
        }
        else
        {
            bjetpt = 0;
            bjeteta = 0;
            bjetphi = 0;
        }
        
        if(nCJet>0)
        {
            cjetpt = jets[leadingCJetIndex].pt;
            cjeteta = jets[leadingCJetIndex].eta;
            cjetphi = jets[leadingCJetIndex].phi;
        }
        else
        {
            cjetpt = 0;
            cjeteta = 0;
            cjetphi = 0;
        }
        
        if(nFJet>0)
        {
            fjetpt = jets[leadingFJetIndex].pt;
            fjeteta = jets[leadingFJetIndex].eta;
            fjetphi = jets[leadingFJetIndex].phi;
        }
        else
        {
            fjetpt = 0;
            fjeteta = 0;
            fjetphi = 0;
        }
        */

        //separate the sample into channels
        if(nISR==0) {}
        else if(nISR>=1) channelIndex += 6;
        //else continue;
        
        if( (ID1>0 && ID2>0) || (ID1<0 && ID2<0) )
        {
            channelIndex += 3;
            element.n++;
            element.nw += weight;
        }
        
        h2[channelIndex]->Fill(channel[channelIndex].Data(),1);
        tree2[channelIndex]->Fill();
    }
    
    //save the resulting trees and histograms
    for(unsigned int j=0;j<channel.size();j++)
    {
        f2[j]->cd();
        tree2[j]->Write("tree");
        h2[j]->Write("hist");
    }
    
    for(int k=1;k<=8;k++)
    {
        cout<<long(h2[1]->GetBinContent(k))<<endl;
    }
    cout<<endl;
    for(unsigned int j=0;j<channel.size();j++)
    {
        cout<<int(h2[j]->GetBinContent(9))<<endl;
    }
    cout<<endl;
    
    //cout<<element.name<<" "<<element.n<<" "<<element.nw<<" "<<element.nAOD<<endl;
    nSS.push_back(element);
    //delete
    for(unsigned int j=0;j<channel.size();j++)
    {
        delete tree2[j];
    }
    
    for(unsigned int j=0;j<channel.size();j++)
    {
        delete h2[j];
    }
    delete tree1;
    for(unsigned int j=0;j<channel.size();j++)
    {
        delete f2[j];
    }
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
        
        TString SampleNameTemp;
        if(type != "Data")
        {
            SampleNameTemp2 += ".";
            fin>>SampleNameTemp;
            SampleNameTemp2 += SampleNameTemp;
        }
        SampleName.push_back(SampleNameTemp2);
        
        for(int i=1;i<=skip;i++)
        {
            double temp;
            fin>>temp;
        }
    }
    fin.close();
}

void skimming()
{
    gInterpreter->GenerateDictionary("PAR0","obj_def.h"); 
    gInterpreter->GenerateDictionary("PAR","obj_def.h"); 
    gInterpreter->GenerateDictionary("vector<L_PAR>","obj_def.h;vector"); 
    gInterpreter->GenerateDictionary("vector<J_PAR>","obj_def.h;vector"); 
    //gInterpreter->GenerateDictionary("vector<TR_PAR>","obj_def.h;vector"); 

    TString SamplePath = "/eos/atlas/user/c/clo/ntuple/";
    //TString SamplePath = "/eos/atlas/user/d/dzhang/susy_ntuples/";
    //TString SamplePath = "/srv/SUSY/ntuple/";
    //TString SamplePath = "/Users/samuel/Atlas/ntuple/";
    
    //SamplePath += "AnalysisBase-02-04-26-a73a6eda/"; TString tag = "v8.10";
    //SamplePath += "AnalysisBase-02-04-26-4dcc2f47/"; TString tag = "v8.12";
    //SamplePath += "AnalysisBase-02-04-26-da7031fc/"; TString tag = "v8.13";
    //SamplePath += "AnalysisBase-02-04-29-f86dc244/"; TString tag = "v9.0";
    //SamplePath += "AnalysisBase-02-04-29-f334c9b6/"; TString tag = "v9.1";
    //SamplePath += "AnalysisBase-02-04-30-f15e6058/"; TString tag = "v9.3";
    //SamplePath += "AnalysisBase-02-04-30-71c02737/"; TString tag = "v9.3.1";
    //SamplePath += "v19.MC/data-myOutput/"; TString tag = "";
    //SamplePath += "v19.MC.2/data-myOutput/"; TString tag = "";
    //SamplePath += "AnalysisBase-02-04-31-2cf44a2c/"; TString tag = "";
    //SamplePath += "AnalysisBase-02-04-31-35a76aa2/"; TString tag = "";
    //SamplePath += "AnalysisBase-02-04-31-ccd99030/"; TString tag = "";
    //SamplePath += "AnalysisBase-02-04-31-8bc21113/"; TString tag = "";
    SamplePath += "AnalysisBase-02-04-31-ebcb0e23/"; TString tag = "";
    
    std::vector<nEvent> nSS;
    
    //Data
    if(false)
    {
        //SamplePath += "data/";
        //tag += "b.Data";
        std::vector<TString> DataSampleName;
        DataSampleName.reserve(20);
        GetSampleName(DataSampleName,"Data",1);
        //for(unsigned int i=0;i<=0;i++)
        for(unsigned int i=0;i<DataSampleName.size();i++)
        {
            skimming2(SamplePath,tag,DataSampleName[i],nSS);
        }
    }
    
    //Background
    if(false)
    {
        //SamplePath += "bkg/";
        //tag += ".MCBkg";
        std::vector<TString> BGSampleName;
        BGSampleName.reserve(20);
        GetSampleName(BGSampleName,"BG",4);
        //for(unsigned int i=114;i<=114;i++)
        //for(unsigned int i=119;i<=120;i++)
        for(unsigned int i=0;i<BGSampleName.size();i++)
        {
            BGSampleName[i] = "mc15_13TeV." + BGSampleName[i];
            skimming2(SamplePath,tag,BGSampleName[i],nSS);
        }
    }
    
    //Signal
    //if(false)
    {
        //SamplePath += "sig/";
        //tag += ".MCSig";
        std::vector<TString> SigSampleName;
        SigSampleName.reserve(20);
        GetSampleName(SigSampleName,"Sig",5);
        for(unsigned int i=0;i<=0;i++)
        //for(unsigned int i=0;i<SigSampleName.size();i++)
        {
            SigSampleName[i] = "mc15_13TeV." + SigSampleName[i];
            skimming2(SamplePath,tag,SigSampleName[i],nSS);
        }
    }
    
    for(unsigned int i=0;i<nSS.size();i++)
    {
        cout<<nSS[i].name<<" "<<nSS[i].n<<" "<<nSS[i].nw<<" "<<nSS[i].nAOD<<endl;
    }
}
