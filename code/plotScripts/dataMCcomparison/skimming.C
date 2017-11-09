#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "obj_def.h"
#include <TInterpreter.h>

#include <TFile.h>
#include <TChain.h>
#include <TH1D.h>
#include <TLorentzVector.h>

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

struct GroupData
{
    TString GroupName;
    TString LegendName;
    TString LatexName;
    unsigned int lower;
    unsigned int upper;
    int colour;
    int unweighted;
    double weighted;
    double error;
};

const bool recalculate_mlj = false;
const float mass_el = 0.000510998;
const float mass_mu = 0.105658;

const bool cutflow = false;
void initializehCutflow(vector<TH1D*>& hSRCutflow)
{
    for(unsigned int j=0;j<6;j++)
    {
        TString hName = "cutflow_";
        hName += TString::Itoa(j,10);
        TH1D* hTemp = new TH1D(hName.Data(), "cutflow", 100, 0, 100);
        hSRCutflow.push_back(hTemp);
        
        hSRCutflow[j]->GetXaxis()->SetBinLabel(1,"nAOD");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(2,"nwAOD");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(3,"trigger");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(4,"trigger,w");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(5,"=2BaseLep and =2SigLep");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(6,"=2BaseLep and =2SigLep,w");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(7,"SS");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(8,"SS,w");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(9,"flavour");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(10,"flavour,w");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(11,"jet");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(12,"jet,w");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(13,"bjet");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(14,"bjet,w");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(15,"Zmass");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(16,"Zmass,w");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(17,"pt1");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(18,"pt1,w");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(19,"pt2");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(20,"pt2,w");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(21,"dEta");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(22,"dEta,w");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(23,"meff");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(24,"meff,w");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(25,"maxmt");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(26,"maxmt,w");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(27,"mlj");
        hSRCutflow[j]->GetXaxis()->SetBinLabel(28,"mlj,w");
    }
}

void deletehCutflow(vector<TH1D*>& hSRCutflow)
{
    for(unsigned int j=0;j<6;j++)
    {
        delete hSRCutflow[j];
    }
}

void printhCutflow(vector<TH1D*>& hSRCutflow)
{
    for(unsigned int i=0;i<6;i++)
    {
        for(unsigned int j=1;j<=28;j++)
        {
            cout<<hSRCutflow[i]->GetXaxis()->GetBinLabel(j)<<": ";
            if(j%2 == 1) cout<<int(hSRCutflow[i]->GetBinContent(j))<<endl;
            if(j%2 == 0) cout<<std::fixed<<std::setprecision(6)<<hSRCutflow[i]->GetBinContent(j)<<endl;
        }
        cout<<endl;
    }
}

void skimming2(TString const& SamplePath,TString const& tag,TString const& SampleName, double XS,vector<TH1D*>& hSRCutflow)
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
        fileName += ".*.myOutput.root*";
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
        tree2[j]->Branch("nFJet",&nFJet,"nFJet/I");
        
        tree2[j]->Branch("weight",&weight,"weight/F");
        tree2[j]->Branch("qFwt",&qFwt,"qFwt/F");
        tree2[j]->Branch("fLwt",&fLwt,"fLwt/F");
        tree2[j]->Branch("averageMu",&averageMu,"averageMu/F");
    }
    
    //histograms
    TH1D* h2[channel.size()];
    for(unsigned int j=0;j<channel.size();j++)
    {
        TString hName2 = "hist_";
        hName2 += channel[j];
        
        f2[j]->cd();
        h2[j] = new TH1D(hName2.Data(), "cut flow", 20, 0, 20);
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
            TH1D *h1 = (TH1D*) f1->Get("hCutFlow");
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
    
    if(h2[0]->GetBinContent(1)==0) cout<<"The sample is missing: "<<SampleName.Data()<<endl;
    const double commonWeight = XS /h2[0]->GetBinContent(2) *(32861.6+3212.96);
    
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
            if(!cutflow) continue;
        }
        for(unsigned int m=0;m<channel.size();m++)
        {
            h2[m]->Fill("trigger",1);
        }
        
        weight = evt.weight * evt.pwt * evt.ElSF * evt.MuSF * evt.BtagSF * evt.JvtSF;
        fLwt = evt.fLwt;
        const double TotalWeight = weight * commonWeight;
        
        int sigIndex[2];
        sigIndex[0] = 0;
        sigIndex[1] = 1;
        pt1 = leps[sigIndex[0]].pt;
        pt2 = leps[sigIndex[1]].pt;
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
            if(cutflow && pt1<=25) continue;
            if(cutflow && pt2<=25) continue;
            
            for(unsigned int m=0;m<channel.size();m++)
            {
                h2[m]->Fill("=2SigLep",1);
            }
            
            if(cutflow)
            {
                for(unsigned int m=0;m<6;m++)
                {
                    hSRCutflow[m]->Fill("=2BaseLep and =2SigLep",1);
                    hSRCutflow[m]->Fill("=2BaseLep and =2SigLep,w",TotalWeight);
                }
            }
        }
        
        ID1 = leps[sigIndex[0]].ID;
        ID2 = leps[sigIndex[1]].ID;
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
        
        //OS or SS
        int channelIndex = 0;
        if( (ID1>0 && ID2>0) || (ID1<0 && ID2<0) )
        {
            channelIndex += 3;
            if(cutflow)
            {
                for(unsigned int m=0;m<6;m++)
                {
                    hSRCutflow[m]->Fill("SS",1);
                    hSRCutflow[m]->Fill("SS,w",TotalWeight);
                }
            }
        }
        else if(cutflow) continue;
        
        //SF or DF
        TLorentzVector l1;
        TLorentzVector l2;
        if(int(abs(ID1)/1000) == 11)
        {
            if(recalculate_mlj) l1.SetPtEtaPhiM(leps[sigIndex[0]].pt,leps[sigIndex[0]].eta,leps[sigIndex[0]].phi,mass_el);
            if(int(abs(ID2)/1000) == 11)
            {
                if(recalculate_mlj) l2.SetPtEtaPhiM(leps[sigIndex[1]].pt,leps[sigIndex[1]].eta,leps[sigIndex[1]].phi,mass_el);
            }
            else if(int(abs(ID2)/1000) == 13)
            {
                if(recalculate_mlj) l2.SetPtEtaPhiM(leps[sigIndex[1]].pt,leps[sigIndex[1]].eta,leps[sigIndex[1]].phi,mass_mu);
                channelIndex += 2;
            }
        }
        else if(int(abs(ID1)/1000) == 13)
        {
            if(recalculate_mlj) l1.SetPtEtaPhiM(leps[sigIndex[0]].pt,leps[sigIndex[0]].eta,leps[sigIndex[0]].phi,mass_mu);
            if(int(abs(ID2)/1000) == 11)
            {
                if(recalculate_mlj) l2.SetPtEtaPhiM(leps[sigIndex[1]].pt,leps[sigIndex[1]].eta,leps[sigIndex[1]].phi,mass_el);
                channelIndex += 2;
            }
            else if(int(abs(ID2)/1000) == 13)
            {
                if(recalculate_mlj) l2.SetPtEtaPhiM(leps[sigIndex[1]].pt,leps[sigIndex[1]].eta,leps[sigIndex[1]].phi,mass_mu);
                channelIndex += 1;
            }
        }
        
        int SRIndex = -1;
        nJet = jets.size();
        if(cutflow)
        {
            if(channelIndex==3)
            {
                hSRCutflow[0]->Fill("flavour",1);
                hSRCutflow[0]->Fill("flavour,w",TotalWeight);
                hSRCutflow[3]->Fill("flavour",1);
                hSRCutflow[3]->Fill("flavour,w",TotalWeight);
                if(nJet == 1)
                {
                    SRIndex = 0;
                }
                else if(nJet == 2 || nJet == 3)
                {
                    SRIndex = 3;
                }
                else continue;
            }
            else if(channelIndex==4)
            {
                hSRCutflow[1]->Fill("flavour",1);
                hSRCutflow[1]->Fill("flavour,w",TotalWeight);
                hSRCutflow[4]->Fill("flavour",1);
                hSRCutflow[4]->Fill("flavour,w",TotalWeight);
                if(nJet == 1)
                {
                    SRIndex = 1;
                }
                else if(nJet == 2 || nJet == 3)
                {
                    SRIndex = 4;
                }
                else continue;
            }
            else if(channelIndex==5)
            {
                hSRCutflow[2]->Fill("flavour",1);
                hSRCutflow[2]->Fill("flavour,w",TotalWeight);
                hSRCutflow[5]->Fill("flavour",1);
                hSRCutflow[5]->Fill("flavour,w",TotalWeight);
                if(nJet == 1)
                {
                    SRIndex = 2;
                }
                else if(nJet == 2 || nJet == 3)
                {
                    SRIndex = 5;
                }
                else continue;
            }
            
            hSRCutflow[SRIndex]->Fill("jet",1);
            hSRCutflow[SRIndex]->Fill("jet,w",TotalWeight);
        }
        
        //jets
        nBJet = 0;
        vector<TLorentzVector> cjet_Ls;
        nCJet = 0;
        nFJet = 0;
        int nISR = 0;
        int leadingJetIndex = 0;
        //int leadingBJetIndex = 0;
        //int leadingCJetIndex = 0;
        //int leadingFJetIndex = 0;
        for(unsigned int k=0;k<jets.size();k++)
        {
            //signal jets
            
            //B-jets
            if(jets[k].jFlag & 1<<5)
            {
                nBJet++;
                //if(nBJet==1) leadingBJetIndex = k;
            }
            
            if(fabs(jets[k].eta) < 2.8)
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
                        TLorentzVector cjet;
                        if(recalculate_mlj) cjet.SetPtEtaPhiM(jets[k].pt,jets[k].eta,jets[k].phi,jets[k].m);
                        if(recalculate_mlj) cjet_Ls.push_back(cjet);
                    }
                }
                
            }
            else
            {
                //Forward jets
                nFJet++;
                //if(nFJet==1) leadingFJetIndex = k;
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
        
        if(nBJet>0)
        {
            //bjetpt = jets[leadingBJetIndex].pt;
            //bjeteta = jets[leadingBJetIndex].eta;
            //bjetphi = jets[leadingBJetIndex].phi;
            if(cutflow) continue;
        }
        else
        {
            //bjetpt = 0;
            //bjeteta = 0;
            //bjetphi = 0;
            if(cutflow)
            {
                hSRCutflow[SRIndex]->Fill("bjet",1);
                hSRCutflow[SRIndex]->Fill("bjet,w",TotalWeight);
            }
        }
        
        /*
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
        
        phi1 = leps[sigIndex[0]].phi;
        mll = l12.m;
        ptll = l12.pt;
        MET = sig.Met;
        METRel = sig.MetRel;
        mTtwo = sig.mT2;
        //mt1 = leps[sigIndex[0]].mT;
        //mt2 = leps[sigIndex[1]].mT;
        mt1 = sqrt(2*pt1*MET*(1-cos(leps[0].MET_dPhi)));
        mt2 = sqrt(2*pt2*MET*(1-cos(leps[1].MET_dPhi)));
        if(mt1>mt2) mtm = mt1;
        else mtm = mt2;
        meff = sig.HT + sig.Met;
        mjj = sig.mjj;
        //R2 = MET/(MET + pt1 + pt2);
        l12_dPhi = l12.dPhi;
        l12_MET_dPhi = l12.MET_dPhi;
        if(jets.size()>=1) jet0_MET_dPhi = jets[0].MET_dPhi;
        else jet0_MET_dPhi = -999;
        l12_jet0_dPhi = l12.jet0_dPhi;
        
        qFwt = evt.qFwt;
        averageMu = evt.averageMu;
        
        if(recalculate_mlj)
        {
            mlj = -1;
            if(cjet_Ls.size()>=1 && cjet_Ls.size()<=3)
            {
                TLorentzVector JetSystem;
                if(cjet_Ls.size()==1)
                {
                    JetSystem = cjet_Ls[0];
                }
                else if(cjet_Ls.size()==2 || cjet_Ls.size()==3)
                {
                    JetSystem = cjet_Ls[0]+cjet_Ls[1];
                }
                
                double dR1 = JetSystem.DeltaR(l1);
                double dR2 = JetSystem.DeltaR(l2);
                TLorentzVector lmindR;
                if(dR1<dR2) lmindR = l1;
                else lmindR = l2;
                mlj = (JetSystem+lmindR).M();
                
                //if((mlj - sig.mlj)/sig.mlj >1e-6) cout<<(mlj - sig.mlj)/sig.mlj<<endl;
            }
        }
        else mlj = sig.mlj;
        
        if(cutflow)
        {
            if(fabs(mll-91.2)>10)
            {
                hSRCutflow[SRIndex]->Fill("Zmass",1);
                hSRCutflow[SRIndex]->Fill("Zmass,w",TotalWeight);
            }
            else continue;
            
            if(pt1>=65)
            {
                hSRCutflow[SRIndex]->Fill("pt1",1);
                hSRCutflow[SRIndex]->Fill("pt1,w",TotalWeight);
            }
            else continue;
            
            if(pt2>=25)
            {
                hSRCutflow[SRIndex]->Fill("pt2",1);
                hSRCutflow[SRIndex]->Fill("pt2,w",TotalWeight);
            }
            else continue;
            
            if(fabs(eta1-eta2)<1.5)
            {
                hSRCutflow[SRIndex]->Fill("dEta",1);
                hSRCutflow[SRIndex]->Fill("dEta,w",TotalWeight);
            }
            else continue;
            
            if(meff>=200)
            {
                hSRCutflow[SRIndex]->Fill("meff",1);
                hSRCutflow[SRIndex]->Fill("meff,w",TotalWeight);
            }
            else continue;
            
            if(mtm>=125)
            {
                hSRCutflow[SRIndex]->Fill("maxmt",1);
                hSRCutflow[SRIndex]->Fill("maxmt,w",TotalWeight);
            }
            else continue;
            
            if(mlj<105)
            {
                hSRCutflow[SRIndex]->Fill("mlj",1);
                hSRCutflow[SRIndex]->Fill("mlj,w",TotalWeight);
            }
            else continue;
        }
        
        //separate the sample into channels
        if(nISR==0) {}
        else if(nISR>=1) channelIndex += 6;
        //else continue;
        
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

void GetSampleName(std::vector<TString>& SampleName, std::vector<double>& XS, TString const type)
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
        
        if(type == "Data")
        {
            double lumi;
            fin>>lumi;
            XS.push_back(lumi);
            cout<<lumi<<endl;
        }
        else if(type == "BG")
        {
            double XS1;
            double XS2;
            
            fin>>XS2;
            
            fin>>XS1;
            XS2 *= XS1;
            
            fin>>XS1;
            XS2 *= XS1;
            
            XS.push_back(XS2);
            cout<<XS2<<endl;
            
            fin>>XS1;
        }
        else if(type == "Sig")
        {
            double XS1;
            double XS2;
            
            fin>>XS1;
            fin>>XS1;
            
            fin>>XS2;
            
            fin>>XS1;
            XS2 *= XS1;
            
            fin>>XS1;
            XS2 *= XS1;
            
            XS.push_back(XS2);
            cout<<XS2<<endl;
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
    //SamplePath += "AnalysisBase-02-04-31-12f0c92d/"; TString tag = "";
    
    //Group for MC background
    std::vector<GroupData> BGMCGroupData;
    {
        GroupData element;
        
        //Z+jets
        element.GroupName = "Zjets"; element.LegendName = "Z+jets"; element.LatexName = "Z+jets";
        element.lower = 0;  element.upper = 41; element.colour = 2; BGMCGroupData.push_back(element);
        
        //DY
        element.GroupName = "DY";    element.LegendName = "DY";     element.LatexName = "DY";
        element.lower = 42; element.upper = 59; element.colour = 2; BGMCGroupData.push_back(element);
        
        //W+jets
        element.GroupName = "Wjets"; element.LegendName = "W+jets"; element.LatexName = "W+jets";
        element.lower = 60;  element.upper = 101; element.colour = 3; BGMCGroupData.push_back(element);
        
        //Top
        element.GroupName = "ttbar"; element.LegendName = "t#bar{t}"; element.LatexName = "$t\\bar{t}$";
        element.lower = 102;  element.upper = 102; element.colour = 881; BGMCGroupData.push_back(element);
        
        element.GroupName = "singletop"; element.LegendName = "single top"; element.LatexName = "single top";
        element.lower = 103;  element.upper = 108; element.colour = 30; BGMCGroupData.push_back(element);
        
        element.GroupName = "ttV"; element.LegendName = "t#bar{t}+V"; element.LatexName = "$t\\bar{t}+V$";
        element.lower = 109;  element.upper = 114; element.colour = 600; BGMCGroupData.push_back(element);
        
        element.GroupName = "multitop"; element.LegendName = "multi top"; element.LatexName = "multi top";
        element.lower = 115;  element.upper = 116; element.colour = 635; BGMCGroupData.push_back(element);
        
        element.GroupName = "VV"; element.LegendName = "VV"; element.LatexName = "VV";
        element.lower = 117;  element.upper = 130; element.colour = 801; BGMCGroupData.push_back(element);
        
        element.GroupName = "Vgamma"; element.LegendName = "V + #gamma"; element.LatexName = "V$+\\gamma$";
        element.lower = 131;  element.upper = 150; element.colour = 5; BGMCGroupData.push_back(element);
        
        element.GroupName = "VVV"; element.LegendName = "VVV"; element.LatexName = "VVV";
        element.lower = 151;  element.upper = 155; element.colour = 607; BGMCGroupData.push_back(element);
        
        element.GroupName = "Higgs"; element.LegendName = "Higgs"; element.LatexName = "Higgs";
        element.lower = 156;  element.upper = 168; element.colour = 7; BGMCGroupData.push_back(element);
    }
    
    //Data
    if(false)
    {
        //SamplePath += "data/";
        //tag += "b.Data";
        std::vector<TString> DataSampleName;
        std::vector<double> XSTemp;
        DataSampleName.reserve(20);
        GetSampleName(DataSampleName,XSTemp,"Data");
        //for(unsigned int i=0;i<=0;i++)
        for(unsigned int i=0;i<DataSampleName.size();i++)
        {
            vector<TH1D*> hSRCutflow;
            initializehCutflow(hSRCutflow);
            skimming2(SamplePath,tag,DataSampleName[i],0,hSRCutflow);
            deletehCutflow(hSRCutflow);
        }
    }
    
    //Background
    if(false)
    {
        //SamplePath += "bkg/";
        //tag += ".MCBkg";
        std::vector<TString> BGSampleName;
        std::vector<double> BGXS;
        BGSampleName.reserve(20);
        GetSampleName(BGSampleName,BGXS,"BG");
        
        for(unsigned int i=0;i<BGMCGroupData.size();i++)
        {
            //if(BGMCGroupData[i].GroupName != "VV") continue;
            
            vector<TH1D*> hSRCutflow;
            initializehCutflow(hSRCutflow);
            
            //for(unsigned int j=114;j<=114;j++)
            //for(unsigned int j=119;j<=120;j++)
            //for(unsigned int j=BGMCGroupData[i].lower;j<=BGMCGroupData[i].lower;j++)
            for(unsigned int j=BGMCGroupData[i].lower;j<=BGMCGroupData[i].upper;j++)
            {
                BGSampleName[j] = "mc15_13TeV." + BGSampleName[j];
                skimming2(SamplePath,tag,BGSampleName[j],BGXS[j],hSRCutflow);
            }
            if(cutflow)
            {
                cout<<"cutflow for "<<BGMCGroupData[i].GroupName.Data()<<endl;
                printhCutflow(hSRCutflow);
            }
            deletehCutflow(hSRCutflow);
        }
    }
    
    //Signal
    //if(false)
    {
        //SamplePath += "sig/";
        //tag += ".MCSig";
        std::vector<TString> SigSampleName;
        std::vector<double> SigXS;
        SigSampleName.reserve(20);
        GetSampleName(SigSampleName,SigXS,"Sig");
        for(unsigned int i=0;i<=0;i++)
        //for(unsigned int i=0;i<SigSampleName.size();i++)
        {
            SigSampleName[i] = "mc15_13TeV." + SigSampleName[i];
            vector<TH1D*> hSRCutflow;
            initializehCutflow(hSRCutflow);
            skimming2(SamplePath,tag,SigSampleName[i],SigXS[i],hSRCutflow);
            if(cutflow) printhCutflow(hSRCutflow);
            deletehCutflow(hSRCutflow);
        }
    }
}
