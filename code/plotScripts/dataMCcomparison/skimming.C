#include <iostream>
#include <fstream>
#include <vector>
#include "evt2l.C"
#include "PP1_evt2l.C"
#include <TH1F.h>

Int_t ID1;
Int_t ID2;
Double_t pt1;
Double_t pt2;
Double_t eta1;
Double_t eta2;
Double_t phi1;
Double_t mll;
Double_t ptll;
Double_t MET;
Double_t mTtwo;
Double_t mt1;
Double_t mt2;
Double_t mtm;

Double_t l12_dPhi;
Double_t l12_MET_dPhi;
Double_t jets_MET_dPhi;

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

Double_t HT;
Double_t R2;

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
        fileName += tag;
        fileName += ".";
        fileName += SampleName;
        fileName += "_myOutput.root/*.root*";
        
        fileName = "/Users/samuel/Atlas/ntuple/test.root";
        
        cout<<fileName.Data()<<endl;
        tree1->Add(fileName.Data());
        if(isPP1) tree1P->Add(fileName.Data());
    }
    evt2l *evts = new evt2l(tree1);
    PP1_evt2l *evtsP = nullptr;
    if(isPP1) evtsP = new PP1_evt2l(tree1P);
    
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
        //fileName += SampleName;
        fileName += "test";
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
        tree2[j]->Branch("pt1",&pt1,"pt1/D");
        tree2[j]->Branch("pt2",&pt2,"pt2/D");
        tree2[j]->Branch("eta1",&eta1,"eta1/D");
        tree2[j]->Branch("eta2",&eta2,"eta2/D");
        tree2[j]->Branch("phi1",&phi1,"phi1/D");
        tree2[j]->Branch("mll",&mll,"mll/D");
        tree2[j]->Branch("ptll",&ptll,"ptll/D");
        tree2[j]->Branch("MET",&MET,"MET/D");
        tree2[j]->Branch("mTtwo",&mTtwo,"mTtwo/D");
        tree2[j]->Branch("mt1",&mt1,"mt1/D");
        tree2[j]->Branch("mt2",&mt2,"mt2/D");
        tree2[j]->Branch("mtm",&mtm,"mtm/D");

        tree2[j]->Branch("l12_dPhi",&l12_dPhi,"l12_dPhi/D");
        tree2[j]->Branch("l12_MET_dPhi",&l12_MET_dPhi,"l12_MET_dPhi/D");
        tree2[j]->Branch("jets_MET_dPhi",&jets_MET_dPhi,"jets_MET_dPhi/D");

        tree2[j]->Branch("nJet",&nJet,"nJet/I");
        tree2[j]->Branch("jetpt",&jetpt,"jetpt/D");
        tree2[j]->Branch("jeteta",&jeteta,"jeteta/D");
        tree2[j]->Branch("nBJet",&nBJet,"nBJet/I");
        tree2[j]->Branch("nCJet",&nCJet,"nCJet/I");

        tree2[j]->Branch("weight",&weight,"weight/D");
        tree2[j]->Branch("qFwt",&qFwt,"qFwt/D");
        tree2[j]->Branch("fLwt",&fLwt,"fLwt/D");
        tree2[j]->Branch("averageMu",&averageMu,"averageMu/D");
    }
    
    //histograms
    TH1F* h2[channel.size()];
    for(unsigned int j=0;j<channel.size();j++)
    {
        TString hName2 = "hist_";
        hName2 += channel[j];
        
        f2[j]->cd();
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
    
    nEvent element;
    element.name = SampleName;
    element.n = 0;
    element.nw = 0;
    element.nAOD = h2[0]->GetBinContent(1);
    
    //loop over all entries
    //for(int j=0;j<=10;j++)
    for(int j=0;j<tree1->GetEntries();j++)
    {
        if(j%100000==0)
        {
            cout<<"number of event: " <<j<<endl;
        }
        evts->GetEntry(j);
        if(isPP1) evtsP->GetEntry(j);
        
        //trigger
        if((evts->sig_trigCode & evts->sig_trigMask)==0)
        {
            //continue;
        }
        for(unsigned int m=0;m<channel.size();m++)
        {
            h2[m]->Fill("trigger",1);
        }
        
        if(isPP1) fLwt = evtsP->evt_fLwt;
        else fLwt = evts->evt_fLwt;
        
        if(!evts->evt_isMC && fLwt!=0)
        {
            for(unsigned int m=0;m<channel.size();m++)
            {
                h2[m]->Fill("fake",1);
            }
        }
        else
        {
            //exact 2 signal leptons
            //if(evts->leps_!=2) continue;
            if(!(evts->leps_lFlag[0] & 1<<1)) continue;
            if(!(evts->leps_lFlag[1] & 1<<1)) continue;
            
            for(unsigned int m=0;m<channel.size();m++)
            {
                h2[m]->Fill("=2SigLep",1);
            }
        }
        
        int sigIndex[2];
        sigIndex[0] = 0;
        sigIndex[1] = 1;
        ID1 = evts->leps_ID[sigIndex[0]];
        ID2 = evts->leps_ID[sigIndex[1]];
        pt1 = evts->leps_pt[sigIndex[0]];
        pt2 = evts->leps_pt[sigIndex[1]];
        eta1 = evts->leps_eta[sigIndex[0]];
        eta2 = evts->leps_eta[sigIndex[1]];
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
        
        phi1 = evts->leps_phi[sigIndex[0]];
        mll = evts->l12_m;
        ptll = evts->l12_pt;
        MET = evts->sig_Met;
        mTtwo = evts->sig_mT2;
        mt1 = evts->leps_mT[sigIndex[0]];
        mt2 = evts->leps_mT[sigIndex[1]];
        if(mt1>mt2) mtm = mt1;
        else mtm = mt2;
        //HT = evts->sig_HT;
        //R2 = MET/(MET + pt1 + pt2);
        l12_dPhi = evts->l12_dPhi;
        l12_MET_dPhi = evts->l12_MET_dPhi;
        jets_MET_dPhi = evts->jets_MET_dPhi[0];
        weight = evts->evt_weight * evts->evt_pwt * evts->evt_ElSF * evts->evt_MuSF;
        
        if(isPP1) qFwt = evtsP->evt_qFwt;
        else qFwt = evts->evt_qFwt;
        
        averageMu = evts->evt_averageMu;
        
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
        nJet = 0;
        nBJet = 0;
        nCJet = 0;
        //nFJet = 0;
        int nISR = 0;
        int leadingJetIndex = 0;
        //int leadingBJetIndex = 0;
        //int leadingCJetIndex = 0;
        //int leadingFJetIndex = 0;
        for(int k=0;k<evts->jets_;k++)
        {
            //signal jets
            if(!(evts->jets_jFlag[k] & 1<<1)) continue;
            nJet++;
            if(nJet==1) leadingJetIndex = k;
            
            //B-jets
            //if((evts->jets_jFlag[k] & 1<<5) && evts->jets_pt[k] > 20)
            if(evts->jets_jFlag[k] & 1<<5)
            {
                nBJet++;
                //if(nBJet==1) leadingBJetIndex = k;
            }
            
            if(fabs(evts->jets_eta[k]) < 2.4)
            {
                //ISR
                if(evts->jets_pt[k] > 40) nISR++;
                
                
                //Central jets
                if(!(evts->jets_jFlag[k] & 1<<5))
                {
                    //no b-tagged
                    if(evts->jets_pt[k] > 20)
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
                if(evts->jets_pt[k] > 30)
                {
                    nFJet++;
                    if(nFJet==1) leadingFJetIndex = k;
                }
                */
            }
        }
        
        if(nJet>0)
        {
            jetpt = evts->jets_pt[leadingJetIndex];
            jeteta = evts->jets_eta[leadingJetIndex];
            //jetphi = evts->jets_phi[leadingJetIndex];
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
            bjetpt = evts->jets_pt[leadingBJetIndex];
            bjeteta = evts->jets_eta[leadingBJetIndex];
            bjetphi = evts->jets_phi[leadingBJetIndex];
        }
        else
        {
            bjetpt = 0;
            bjeteta = 0;
            bjetphi = 0;
        }
        
        if(nCJet>0)
        {
            cjetpt = evts->jets_pt[leadingCJetIndex];
            cjeteta = evts->jets_eta[leadingCJetIndex];
            cjetphi = evts->jets_phi[leadingCJetIndex];
        }
        else
        {
            cjetpt = 0;
            cjeteta = 0;
            cjetphi = 0;
        }
        
        if(nFJet>0)
        {
            fjetpt = evts->jets_pt[leadingFJetIndex];
            fjeteta = evts->jets_eta[leadingFJetIndex];
            fjetphi = evts->jets_phi[leadingFJetIndex];
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
    
    //cout<<element.name<<" "<<element.n<<" "<<element.nw<<" "<<element.nAOD<<endl;
    nSS.push_back(element);
    //delete
    delete evts;
    if(isPP1) delete evtsP;
    for(unsigned int j=0;j<channel.size();j++)
    {
        delete tree2[j];
    }
    
    for(unsigned int j=0;j<channel.size();j++)
    {
        delete h2[j];
    }
    delete tree1;
    if(isPP1) delete tree1P;
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

void skimming()
{
    TString SamplePath = "/eos/atlas/user/c/clo/ntuple/";
    //TString SamplePath = "/srv/SUSY/ntuple/";
    //TString SamplePath = "/Users/samuel/Atlas/ntuple/";
    
    //SamplePath += "AnalysisBase-02-04-26-a73a6eda/"; TString tag = "v8.10";
    //SamplePath += "AnalysisBase-02-04-26-4dcc2f47/"; TString tag = "v8.12";
    //SamplePath += "AnalysisBase-02-04-26-da7031fc/"; TString tag = "v8.13";
    //SamplePath += "AnalysisBase-02-04-29-f86dc244/"; TString tag = "v9.0";
    //SamplePath += "AnalysisBase-02-04-29-f334c9b6/"; TString tag = "v9.1";
    //SamplePath += "AnalysisBase-02-04-30-f15e6058/"; TString tag = "v9.3";
    SamplePath += "AnalysisBase-02-04-30-71c02737/"; TString tag = "v9.3.1";
    
    std::vector<nEvent> nSS;
    
    //Data
    //if(true)
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
            skimming2(SamplePath,tag,DataSampleName[i],0,nSS);
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
        for(unsigned int i=68;i<=69;i++)
        //for(unsigned int i=0;i<BGSampleName.size();i++)
        {
            BGSampleName[i] = "mc15_13TeV." + BGSampleName[i];
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
            SigSampleName[i] = "mc15_13TeV." + SigSampleName[i];
            skimming2(SamplePath,tag,SigSampleName[i],false,nSS);
        }
    }
    
    for(unsigned int i=0;i<nSS.size();i++)
    {
        cout<<nSS[i].name<<" "<<nSS[i].n<<" "<<nSS[i].nw<<" "<<nSS[i].nAOD<<endl;
    }
}
