#include <iostream>
#include <fstream>
#include <vector>
#include "evt2l.C"
#include <TH1F.h>

Int_t ID1;
Int_t ID2;
Double_t pt1;
Double_t pt2;
Double_t mll;
Double_t ptll;
Double_t MET;
Double_t mT2;
Int_t nJet;
Int_t nBJet;
Double_t weight;
Double_t averageMu;
Int_t nVtx;

bool PASS;


void skimming2(TString SamplePath,TString SampleName)
{
    //create a file for storing the resulting trees and histograms
    TString fileName="skimming/skimming.";
    fileName+=SampleName;
    fileName+=".root";
    TFile *f2 = new TFile(fileName.Data(),"RECREATE");
    
    //6 channels
    TString *channel[8];
    channel[1]= new TString("ee_OS");
    channel[2]= new TString("mumu_OS");
    channel[3]= new TString("emu_OS");
    channel[4]= new TString("ee_SS");
    channel[5]= new TString("mumu_SS");
    channel[6]= new TString("emu_SS");
    channel[7]= new TString("combine");

    //trees
    TTree *tree2[7];
    TString treeName[7];
    for(int j=1;j<=6;j++)
    {
        treeName[j]="skimming_";
        treeName[j]+= *channel[j];

        tree2[j] = new TTree(treeName[j].Data(),treeName[j].Data());
        tree2[j]->Branch("ID1",&ID1,"ID1/I");
        tree2[j]->Branch("ID2",&ID2,"ID2/I");
        tree2[j]->Branch("pt1",&pt1,"pt1/D");
        tree2[j]->Branch("pt2",&pt2,"pt2/D");
        tree2[j]->Branch("mll",&mll,"mll/D");
        tree2[j]->Branch("ptll",&ptll,"ptll/D");
        tree2[j]->Branch("MET",&MET,"MET/D");
        tree2[j]->Branch("mT2",&mT2,"mT2/D");
        tree2[j]->Branch("nJet",&nJet,"nJet/I");
        tree2[j]->Branch("nBJet",&nBJet,"nBJet/I");
        tree2[j]->Branch("weight",&weight,"weight/D");
        tree2[j]->Branch("averageMu",&averageMu,"averageMu/D");
        tree2[j]->Branch("nVtx",&nVtx,"nVtx/I");
    }
    
    //histograms
    TH1F* h2[8];
    TString hName2[8];
    for(int j=1;j<=7;j++)
    {
        hName2[j]="h_";
        hName2[j]+= *channel[j];
        
        h2[j] = new TH1F(hName2[j].Data(), "cut flow", 20, 0, 20);
        h2[j]->GetXaxis()->SetBinLabel(1,"AOD");
        h2[j]->GetXaxis()->SetBinLabel(2,"All");
        h2[j]->GetXaxis()->SetBinLabel(3,"NLepton");
        h2[j]->GetXaxis()->SetBinLabel(4,"trigger");
        h2[j]->GetXaxis()->SetBinLabel(5,"GRL");
        h2[j]->GetXaxis()->SetBinLabel(6,"Bad muon");
        h2[j]->GetXaxis()->SetBinLabel(7,"Comsic muon");
        h2[j]->GetXaxis()->SetBinLabel(8,"Bad jet");
        h2[j]->GetXaxis()->SetBinLabel(9,"PV");
        h2[j]->GetXaxis()->SetBinLabel(10,"2Lepton");
        h2[j]->GetXaxis()->SetBinLabel(11,"pt1");
        h2[j]->GetXaxis()->SetBinLabel(12,"pt2");
        h2[j]->GetXaxis()->SetBinLabel(13,"mll_60");
        h2[j]->GetXaxis()->SetBinLabel(14,channel[j]->Data());
    }
    
    //get the "evt2l" and "hCutFlow"
    TFile *f1 = TFile::Open(SamplePath.Data());
    TH1F *h1 = (TH1F*) f1->Get("hCutFlow");
    TTree *tree1 = (TTree*) f1->Get("evt2l");
    evt2l *evts = new evt2l(tree1);
    
    for(int j=1;j<=7;j++)
    {
        if(h1->GetBinContent(3)==0)
        {
            h2[j]->Fill("AOD",h1->GetBinContent(1));
            h2[j]->Fill("All",h1->GetBinContent(1));
            h2[j]->Fill("NLepton",h1->GetBinContent(2));
        }
        else
        {
            h2[j]->Fill("AOD",h1->GetBinContent(1));
            h2[j]->Fill("All",h1->GetBinContent(2));
            h2[j]->Fill("NLepton",h1->GetBinContent(3));
        }
    }
    
    //loop over all entries
    //for(int j=0;j<=10;j++)
    for(int j=0;j<tree1->GetEntries();j++)
    {
        
        if(j%100000==0)
        {
            cout<<"number of event: " <<j<<endl;
        }
        evts->GetEntry(j);
        
        //trigger
        if(evts->sig_trigCode==0)
        {
            continue;
        }
        for(int m=1;m<=7;m++)
        {
            h2[m]->Fill("trigger",1);
        }
        
        //GRL
        if(!(evts->evt_cuts==1))
        {
            continue;
        }
        for(int m=1;m<=7;m++)
        {
            h2[m]->Fill("GRL",1);
        }
     
        //remove events with bad muon
        PASS = false;
        for(int k=0;k<evts->leps_;k++)
        {
            if(abs(int(evts->leps_ID[k]/1000)) == 13 && (evts->leps_lFlag[k] & 1<<3))
            {
                PASS = true;
            }
        }
        if(PASS)
        {
            continue;
        }
        for(int m=1;m<=7;m++)
        {
            h2[m]->Fill("Bad muon",1);
        }
        
        //remove events with comsic muon
        PASS = false;
        for(int k=0;k<evts->leps_;k++)
        {
            if(abs(int(evts->leps_ID[k]/1000)) == 13 && (evts->leps_lFlag[k] & 1<<4))
            {
                PASS = true;
            }
        }
        if(PASS)
        {
            continue;
        }
        for(int m=1;m<=7;m++)
        {
            h2[m]->Fill("Comsic muon",1);
        }
        
        //remove events with bad jet
        PASS = false;
        for(int k=0;k<evts->jets_;k++)
        {
            if(evts->jets_jFlag[k] & 1<<3)
            {
                PASS = true;
            }
        }
        if(PASS)
        {
            continue;
        }
        for(int m=1;m<=7;m++)
        {
            h2[m]->Fill("Bad jet",1);
        }
        
        //PV
        if(!(evts->sig_nPV==1))
        {
            continue;
        }
        for(int m=1;m<=7;m++)
        {
            h2[m]->Fill("PV",1);
        }
        
        //only two leptons
        if(!(evts->leps_==2))
        {
            continue;
        }
        for(int m=1;m<=7;m++)
        {
            h2[m]->Fill("2Lepton",1);
        }
        
        //pt of leading lepton > 30 GeV
        if(!(evts->leps_pt[0]>30))
        {
            continue;
        }
        for(int m=1;m<=7;m++)
        {
            h2[m]->Fill("pt1",1);
        }
        
        //pt of subleading lepton > 30 GeV
        if(!(evts->leps_pt[1]>30))
        {
            continue;
        }
        for(int m=1;m<=7;m++)
        {
            h2[m]->Fill("pt2",1);
        }
        //mll > 60 GeV
        if(!(evts->l12_m>60))
        {
            continue;
        }
        for(int m=1;m<=7;m++)
        {
            h2[m]->Fill("mll_60",1);
        }
        
        
        ID1 = evts->leps_ID[0];
        ID2 = evts->leps_ID[1];
        pt1 = evts->leps_pt[0];
        pt2 = evts->leps_pt[1];
        mll = evts->l12_m;
        ptll = evts->l12_pt;
        MET = evts->sig_MetRel;
        mT2 = evts->sig_mT2;
        weight = evts->evt_weight * evts->evt_ElSF * evts->evt_MuSF;
        //weight = evts->evt_weight;
        averageMu = evts->evt_averageMu;
        nVtx = evts->sig_nVtx;
        
        //jets
        nJet = 0;
        nBJet = 0;
        for(int k=0;k<evts->jets_;k++)
        {
            if((evts->jets_jFlag[k] & 1<<1))
            {
                nJet++;
                if((evts->jets_jFlag[k] & 1<<5))
                {
                    nBJet++;
                }
            }
        }
        
        //separate the sample into 6 channels
        if( (evts->leps_ID[0]>0 && evts->leps_ID[1]<0) ||
           (evts->leps_ID[0]<0 && evts->leps_ID[1]>0) )
        {
            if(int(abs(evts->leps_ID[0])/1000) == 11)
            {
                
                if(int(abs(evts->leps_ID[1])/1000) == 11)
                {
                    h2[1]->Fill(channel[1]->Data(),1);
                    h2[7]->Fill(channel[7]->Data(),1);
                    tree2[1]->Fill();
                }
                else if(int(abs(evts->leps_ID[1])/1000) == 13)
                {
                    h2[3]->Fill(channel[3]->Data(),1);
                    h2[7]->Fill(channel[7]->Data(),1);
                    tree2[3]->Fill();
                }
            }
            else if(int(abs(evts->leps_ID[0])/1000) == 13)
            {
                if(int(abs(evts->leps_ID[1])/1000) == 11)
                {
                    h2[3]->Fill(channel[3]->Data(),1);
                    h2[7]->Fill(channel[7]->Data(),1);
                    tree2[3]->Fill();
                }
                else if(int(abs(evts->leps_ID[1])/1000) == 13)
                {
                    h2[2]->Fill(channel[2]->Data(),1);
                    h2[7]->Fill(channel[7]->Data(),1);
                    tree2[2]->Fill();
                }
            }
        }
        else
        {
            if(int(abs(evts->leps_ID[0])/1000) == 11)
            {
                
                if(int(abs(evts->leps_ID[1])/1000) == 11)
                {
                    h2[4]->Fill(channel[4]->Data(),1);
                    h2[7]->Fill(channel[7]->Data(),1);
                    tree2[4]->Fill();
                }
                else if(int(abs(evts->leps_ID[1])/1000) == 13)
                {
                    h2[6]->Fill(channel[6]->Data(),1);
                    h2[7]->Fill(channel[7]->Data(),1);
                    tree2[6]->Fill();
                }
            }
            else if(int(abs(evts->leps_ID[0])/1000) == 13)
            {
                if(int(abs(evts->leps_ID[1])/1000) == 11)
                {
                    h2[6]->Fill(channel[6]->Data(),1);
                    h2[7]->Fill(channel[7]->Data(),1);
                    tree2[6]->Fill();
                }
                else if(int(abs(evts->leps_ID[1])/1000) == 13)
                {
                    h2[5]->Fill(channel[5]->Data(),1);
                    h2[7]->Fill(channel[7]->Data(),1);
                    tree2[5]->Fill();
                }
            }
        }
    }
    
    //save the resulting trees and histograms
    f2->ReOpen("UPDATE");
    for(int j=1;j<=6;j++)
    {
        tree2[j]->Write();
        cout<<tree2[j]->GetEntries()<<endl;
    }
    cout<<endl;
    for(int j=1;j<=7;j++)
    {
        h2[j]->Write();
    }
    for(int k=1;k<=13;k++)
    {
        cout<<int(h2[1]->GetBinContent(k))<<endl;
    }
    cout<<endl;
    for(int j=1;j<=7;j++)
    {
        cout<<int(h2[j]->GetBinContent(14))<<endl;
    }
}

void GetSampleName(std::vector<TString*> &SampleName, TString type, int skip)
{
    ifstream fin;
    TString fileName = type;
    fileName+="Sample.txt";
    fin.open(fileName.Data());
    
    TString SampleNameTemp;
    TString SampleNameTemp2;
    TString *pSampleNameTemp;
    double temp=0;
    while(!fin.eof())
    {
        fin>>SampleNameTemp2;
        if(fin.eof()) break;
        SampleNameTemp2 += ".";
        fin>>SampleNameTemp;
        SampleNameTemp2 += SampleNameTemp;
        pSampleNameTemp = new TString(SampleNameTemp2);
        SampleName.push_back(pSampleNameTemp);
        
        for(int i=1;i<=skip;i++)
        {
            fin>>temp;
        }
    }
    fin.close();
}

void skimming()
{
    TString SamplePath = "/afs/cern.ch/work/c/clo/public/ntuple/AnalysisSUSY-02-03-33a-399128/v5.18.s.";
    TString SamplePath2;
    
    //Data
    std::vector<TString*> DataSampleName;
    DataSampleName.reserve(20);
    GetSampleName(DataSampleName,"Data",1);
    for(unsigned int i=0;i<DataSampleName.size();i++)
    {
        SamplePath2=SamplePath;
        SamplePath2+= *DataSampleName[i];
        SamplePath2+=".root";
        cout<<SamplePath2.Data()<<endl;
        skimming2(SamplePath2,*DataSampleName[i]);
    }
    
    //Background
    std::vector<TString*> BGSampleName;
    BGSampleName.reserve(20);
    GetSampleName(BGSampleName,"BG",4);
    //for(unsigned int i=0;i<=14;i++)
    for(unsigned int i=0;i<BGSampleName.size();i++)
    {
        SamplePath2=SamplePath;
        SamplePath2+= *BGSampleName[i];
        SamplePath2+=".root";
        cout<<SamplePath2.Data()<<endl;
        skimming2(SamplePath2,*BGSampleName[i]);
    }
    
    //Signal
    std::vector<TString*> SigSampleName;
    SigSampleName.reserve(20);
    GetSampleName(SigSampleName,"Sig",5);
    for(unsigned int i=0;i<SigSampleName.size();i++)
    {
        SamplePath2=SamplePath;
        SamplePath2+= *SigSampleName[i];
        SamplePath2+=".root";
        cout<<SamplePath2.Data()<<endl;
        skimming2(SamplePath2,*SigSampleName[i]);
    }
}
