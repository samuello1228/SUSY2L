// C++
#include <iostream>
#include <fstream>
using namespace std;
// ROOT
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TH2D.h"
// My packages
#include "../multiLepSearch/susyEvts.h"

int MCId;
Long64_t evn;

int NlepBL;
int nSigLep;
int isTruthLep1;
int isTruthLep2;
int nJets20;
int nBJets20;
Float_t Pt_l;
Float_t Pt_subl;
Float_t DeltaEtaLep;
Float_t met;
Float_t mt;
Float_t meff;
Float_t mljj_comb;
Float_t MT2;

//weight
Float_t mcweight;
Float_t puweight;
Float_t totweight;
Float_t lumiScaling;

Float_t wmu_nom;
Float_t wel_nom;
Float_t wjet_nom;
Float_t wpu_nom_bkg;
Float_t wtrig_nom;
Float_t wchflip_nom;
Float_t wttV_nom;

void checkError(float x, float y, TString name, float limit)
{
    float error = fabs((y-x)/x);
    if(error > limit)
    {
        cout<<"Event number: "<<evn<<": "<<name.Data()<<" are different: "<<x<<", "<<y<<", "<<error<<endl;
    }
}

struct Sample
{
    TString tag;
    int ID;
    double XS;
    vector<int> filter_first;
    vector<int> filter_last;
};

void checkAllVariables(susyEvts* evt, Sample& sample, bool isBG, bool isSig)
{
        //wight
        if(isBG || isSig)
        {
            checkError(mcweight, evt->evt.weight, "mc gen weight", 2e-7);
            checkError(puweight, evt->evt.pwt, "pileup weight", 2e-7);
            if(isBG) checkError(puweight, wpu_nom_bkg, "pileup weight", 2e-7);
            checkError(wel_nom * wchflip_nom, evt->evt.ElSF, "electron weight", 5e-7);
            checkError(wmu_nom, evt->evt.MuSF, "muon weight", 2e-7);
            checkError(wjet_nom, evt->evt.BtagSF * evt->evt.JvtSF, "jet weight", 2e-7);
            checkError(wtrig_nom, evt->evt.trigSF, "trigger weight", 2e-7);
            double ttV_SF = 1;
            if(sample.ID == 410155 && evt->sig.nBJet>=3) ttV_SF = 1.34;
            else if(sample.ID >= 410218 && sample.ID <= 410220 && evt->sig.nBJet>=3) ttV_SF = 1.37;
            checkError(wttV_nom, ttV_SF, "ttV weight", 2e-7);
            checkError(totweight, evt->evt.weight * evt->evt.pwt * evt->evt.trigSF * evt->evt.ElSF * evt->evt.MuSF * evt->evt.BtagSF * evt->evt.JvtSF * ttV_SF, "Total weight", 6e-7);
        }

        //number of leptons
        if(NlepBL != evt->sig.nBaseLep) cout<<"nBaseLep are different."<<endl;
        if(NlepBL != int(evt->leps.size())) cout<<"nBaseLep are different."<<endl;
        if(nSigLep != evt->sig.nSigLep) cout<<"nSigLep are different."<<endl;

        //lepTruth
        if(isBG && isTruthLep1 != evt->leps[0].lepTruth) cout<<"lepTruth1 are different."<<endl;
        if(isBG && isTruthLep2 != evt->leps[1].lepTruth) cout<<"lepTruth2 are different."<<endl;

        //number of jets
        if(nJets20 != evt->sig.nJet) cout<<"nJet are different."<<endl;
        if(nBJets20 != evt->sig.nBJet) cout<<"nBJet are different."<<endl;

        //variables
        checkError(Pt_l, evt->leps[0].pt * 1000, "pt1", 3e-6);
        checkError(Pt_subl, evt->leps[1].pt * 1000, "pt2", 4e-6);
        checkError(DeltaEtaLep, fabs(evt->leps[0].eta - evt->leps[1].eta), "dEta", 2e-7);
        checkError(met, evt->sig.Met * 1000, "MET", 7e-6);
        checkError(mt, evt->leps[0].mT * 1000, "mT", 3e-5);

        //meff
        float pt_sum = 0;
        for(int k=2;k<nSigLep;k++) pt_sum += evt->leps[k].pt;
        checkError(meff, (evt->sig.HT + evt->sig.Met + pt_sum) * 1000, "meff", 2e-6);

        //mlj
        if(nSigLep == 2 && nJets20 > 0) checkError(mljj_comb, evt->sig.mlj * 1000, "mlj", 2e-6);

        //MT2
        if(MT2>10000) checkError(MT2, evt->sig.mT2 * 1000, "MT2", 5e-6);
}

void FillHist_Dani(TH1D* hCutflow_Dani, bool isBG)
{
    double weight = totweight*lumiScaling;
    hCutflow_Dani->Fill("Z veto",1);
    hCutflow_Dani->Fill("Z veto,w",weight);

    if(NlepBL !=2) return;
    if(nSigLep !=2) return;
    hCutflow_Dani->Fill("=2BaseLep && =2SigLep",1);
    hCutflow_Dani->Fill("=2BaseLep && =2SigLep,w",weight);

    if(Pt_l < 25000) return;
    hCutflow_Dani->Fill("pt1",1);
    hCutflow_Dani->Fill("pt1,w",weight);

    if(Pt_subl < 25000) return;
    hCutflow_Dani->Fill("pt2",1);
    hCutflow_Dani->Fill("pt2,w",weight);

    if(isBG && isTruthLep1 != 1) return;
    hCutflow_Dani->Fill("isTruthLep1",1);
    hCutflow_Dani->Fill("isTruthLep1,w",weight);

    if(isBG && isTruthLep2 != 1) return;
    hCutflow_Dani->Fill("isTruthLep2",1);
    hCutflow_Dani->Fill("isTruthLep2,w",weight);

    if(nBJets20 !=0) return;
    hCutflow_Dani->Fill("nBJets20",1);
    hCutflow_Dani->Fill("nBJets20,w",weight);

    if(nJets20 ==1)
    {
        hCutflow_Dani->Fill("SR1:nJet",1);
        hCutflow_Dani->Fill("SR1:nJet,w",weight);

        if(fabs(DeltaEtaLep)>1.5) return;
        hCutflow_Dani->Fill("SR1:DeltaEtaLep",1);
        hCutflow_Dani->Fill("SR1:DeltaEtaLep,w",weight);

        if(met<100000) return;
        hCutflow_Dani->Fill("SR1:met",1);
        hCutflow_Dani->Fill("SR1:met,w",weight);

        if(mt<140000) return;
        hCutflow_Dani->Fill("SR1:mt",1);
        hCutflow_Dani->Fill("SR1:mt,w",weight);

        if(meff<260000) return;
        hCutflow_Dani->Fill("SR1:meff",1);
        hCutflow_Dani->Fill("SR1:meff,w",weight);

        if(mljj_comb>=180000) return;
        hCutflow_Dani->Fill("SR1:mljj",1);
        hCutflow_Dani->Fill("SR1:mljj,w",weight);

        if(MT2<80000) return;
        hCutflow_Dani->Fill("SR1:MT2",1);
        hCutflow_Dani->Fill("SR1:MT2,w",weight);

        //cout<<"Event number in SRjet1 in Dani: "<<evn<<endl;
    }
    else if(nJets20 ==2 || nJets20 ==3)
    {
        hCutflow_Dani->Fill("SR2:nJet",1);
        hCutflow_Dani->Fill("SR2:nJet,w",weight);

        hCutflow_Dani->Fill("SR2:DeltaEtaLep",1);
        hCutflow_Dani->Fill("SR2:DeltaEtaLep,w",weight);

        if(met<100000) return;
        hCutflow_Dani->Fill("SR2:met",1);
        hCutflow_Dani->Fill("SR2:met,w",weight);

        if(mt<120000) return;
        hCutflow_Dani->Fill("SR2:mt",1);
        hCutflow_Dani->Fill("SR2:mt,w",weight);

        if(meff<240000) return;
        hCutflow_Dani->Fill("SR2:meff",1);
        hCutflow_Dani->Fill("SR2:meff,w",weight);

        if(mljj_comb>=130000) return;
        hCutflow_Dani->Fill("SR2:mljj",1);
        hCutflow_Dani->Fill("SR2:mljj,w",weight);

        if(MT2<70000) return;
        hCutflow_Dani->Fill("SR2:MT2",1);
        hCutflow_Dani->Fill("SR2:MT2,w",weight);

        //cout<<"Event number in SRjet2 in Dani: "<<evn<<endl;
    }
}

void FillHist_Samuel(TH1D* hCutflow_Samuel, susyEvts* evt, const double& commonWeight, double& ttV_SF, bool isBG, bool isSig)
{
    double weight = 1;
    if(isBG || isSig) weight = evt->evt.weight * evt->evt.pwt * evt->evt.trigSF * evt->evt.ElSF * evt->evt.MuSF * evt->evt.BtagSF * evt->evt.JvtSF * ttV_SF * commonWeight;
    hCutflow_Samuel->Fill("Z veto",1);
    hCutflow_Samuel->Fill("Z veto,w",weight);

    if(evt->sig.nBaseLep !=2) return;
    if(evt->sig.nSigLep !=2) return;
    hCutflow_Samuel->Fill("=2BaseLep && =2SigLep",1);
    hCutflow_Samuel->Fill("=2BaseLep && =2SigLep,w",weight);

    if(evt->leps[0].pt < 25) return;
    hCutflow_Samuel->Fill("pt1",1);
    hCutflow_Samuel->Fill("pt1,w",weight);

    if(evt->leps[1].pt < 25) return;
    hCutflow_Samuel->Fill("pt2",1);
    hCutflow_Samuel->Fill("pt2,w",weight);

    if(isBG && evt->leps[0].lepTruth != 1) return;
    hCutflow_Samuel->Fill("isTruthLep1",1);
    hCutflow_Samuel->Fill("isTruthLep1,w",weight);

    if(isBG && evt->leps[1].lepTruth != 1) return;
    hCutflow_Samuel->Fill("isTruthLep2",1);
    hCutflow_Samuel->Fill("isTruthLep2,w",weight);

    if(evt->sig.nBJet !=0) return;
    hCutflow_Samuel->Fill("nBJets20",1);
    hCutflow_Samuel->Fill("nBJets20,w",weight);

    if(evt->sig.nJet ==1)
    {
        hCutflow_Samuel->Fill("SR1:nJet",1);
        hCutflow_Samuel->Fill("SR1:nJet,w",weight);

        if(fabs(evt->leps[0].eta - evt->leps[1].eta) > 1.5) return;
        hCutflow_Samuel->Fill("SR1:DeltaEtaLep",1);
        hCutflow_Samuel->Fill("SR1:DeltaEtaLep,w",weight);

        if(evt->sig.Met<100) return;
        hCutflow_Samuel->Fill("SR1:met",1);
        hCutflow_Samuel->Fill("SR1:met,w",weight);

        if(evt->leps[0].mT<140) return;
        hCutflow_Samuel->Fill("SR1:mt",1);
        hCutflow_Samuel->Fill("SR1:mt,w",weight);

        if(evt->sig.HT + evt->sig.Met < 260) return;
        hCutflow_Samuel->Fill("SR1:meff",1);
        hCutflow_Samuel->Fill("SR1:meff,w",weight);

        if(evt->sig.mlj>=180) return;
        hCutflow_Samuel->Fill("SR1:mljj",1);
        hCutflow_Samuel->Fill("SR1:mljj,w",weight);

        if(evt->sig.mT2<80) return;
        hCutflow_Samuel->Fill("SR1:MT2",1);
        hCutflow_Samuel->Fill("SR1:MT2,w",weight);

        //cout<<"Event number in SRjet1 in Samuel: "<<evt->evt.event<<endl;
    }
    else if(evt->sig.nJet ==2 || evt->sig.nJet ==3)
    {
        hCutflow_Samuel->Fill("SR2:nJet",1);
        hCutflow_Samuel->Fill("SR2:nJet,w",weight);

        hCutflow_Samuel->Fill("SR2:DeltaEtaLep",1);
        hCutflow_Samuel->Fill("SR2:DeltaEtaLep,w",weight);

        if(evt->sig.Met<100) return;
        hCutflow_Samuel->Fill("SR2:met",1);
        hCutflow_Samuel->Fill("SR2:met,w",weight);

        if(evt->leps[0].mT<120) return;
        hCutflow_Samuel->Fill("SR2:mt",1);
        hCutflow_Samuel->Fill("SR2:mt,w",weight);

        if(evt->sig.HT + evt->sig.Met < 240) return;
        hCutflow_Samuel->Fill("SR2:meff",1);
        hCutflow_Samuel->Fill("SR2:meff,w",weight);

        if(evt->sig.mlj>=130) return;
        hCutflow_Samuel->Fill("SR2:mljj",1);
        hCutflow_Samuel->Fill("SR2:mljj,w",weight);

        if(evt->sig.mT2<70) return;
        hCutflow_Samuel->Fill("SR2:MT2",1);
        hCutflow_Samuel->Fill("SR2:MT2,w",weight);

        //cout<<"Event number in SRjet2 in Samuel: "<<evt->evt.event<<endl;
    }
}

void RunSample(TChain* tree1, Sample& sample, TH1D* hCutflow_Dani, TH1D* hCutflow_Samuel, bool isBG, bool isSig)
{
    bool isMC = (isBG || isSig);

    struct evn_info
    {
        long evn;
        long Dani_index;
        long Samuel_index;
        int nJet; //For 361073
    };

    double lumiScaling_Dani = -1;
    vector<evn_info> selected_Dani;
    for (long i = 0; i < tree1->GetEntries(); i++)
    {
        if(i%100000==0) cout<<"number of event: " <<i<<endl;
        tree1->GetEntry(i);

        if(isMC && MCId != sample.ID) continue;

        //cutflow for Dani ntuple
        FillHist_Dani(hCutflow_Dani, isBG);

        //Fill selected_Dani
        evn_info evn_temp;
        evn_temp.evn = evn;
        evn_temp.Samuel_index = -1;
        evn_temp.Dani_index = i;
        evn_temp.nJet = nJets20;
        selected_Dani.push_back(evn_temp);

        if(lumiScaling_Dani == -1) lumiScaling_Dani = lumiScaling;
    }

    //Our ntuple
    //TString path2 = "/afs/cern.ch/user/c/clo/AnalysisBase-02-04-31/SUSY2L/code/run/t1/data-myOutput/test.root";
    
    ///*
    //TString path2 = "/eos/user/c/clo/ntuple/AnalysisBase-02-04-31-6ecc6eb7/user.clo.v13.5.";
    TString path2 = "/eos/user/c/clo/ntuple/AnalysisBase-02-04-31-f05df733/user.clo.v13.6.";
    path2 += sample.tag;
    path2 += "_myOutput.root/user.clo.";
    if(isMC)
    {
        path2 += "mc15_13TeV.";
        path2 += TString::Format("%i",sample.ID);
    }
    else
    {
        path2 += "data*";
    }
    path2 += ".*.myOutput.root*";
    //*/
    cout<<path2<<endl;

    TChain* tree2 = new TChain("evt2l");
    tree2->Add(path2);
    susyEvts* evt = new susyEvts(tree2);
    cout<<"Our ntuple has events: "<<tree2->GetEntries()<<endl;

    double nwAOD = 0;
    TObjArray* fileList = tree2->GetListOfFiles();
    for(int i=0;i<fileList->GetEntries();i++)
    {
        cout<<fileList->At(i)->GetTitle()<<endl;
        TFile *f1 = TFile::Open(fileList->At(i)->GetTitle());
        TH1D *h1 = (TH1D*) f1->Get("hCutFlow");
        //for(unsigned int j=1;j<80;j++) cout<<j<<": "<<h1->GetXaxis()->GetBinLabel(j)<<": "<<h1->GetBinContent(j)<<endl;;
        nwAOD += h1->GetBinContent(2);

        hCutflow_Samuel->Fill("AOD", h1->GetBinContent(1));
        hCutflow_Samuel->Fill("SUSY2", h1->GetBinContent(4));
        hCutflow_Samuel->Fill("GRL", h1->GetBinContent(8));
        hCutflow_Samuel->Fill("Trigger", h1->GetBinContent(9));
        hCutflow_Samuel->Fill("LAr+Tile+SCT+CoreFlag", h1->GetBinContent(10));
        hCutflow_Samuel->Fill("PV", h1->GetBinContent(11));
        hCutflow_Samuel->Fill("cosMuon", h1->GetBinContent(12));
        hCutflow_Samuel->Fill("BadMuon", h1->GetBinContent(13));
        hCutflow_Samuel->Fill("BadJet", h1->GetBinContent(14));
        hCutflow_Samuel->Fill(">=2BaseLep", h1->GetBinContent(70));
        hCutflow_Samuel->Fill(">=2SigLep", h1->GetBinContent(72));
        hCutflow_Samuel->Fill("SS leptons", h1->GetBinContent(74));
        hCutflow_Samuel->Fill(">=1 passOR jet", h1->GetBinContent(75));
        hCutflow_Samuel->Fill(">=1 signal jet", h1->GetBinContent(76));
        
        delete f1;
    }

    bool DoFilter = true;
    //DoFilter = false;
    DoFilter = DoFilter && (sample.filter_first.size()!=0);
    const double lumi = 36100;
    double commonWeight = 1;
    if(isMC) commonWeight = sample.XS *lumi /nwAOD;
    if(isMC && DoFilter)
    {
        commonWeight = lumiScaling_Dani;
        cout<<"Dani lumiScaling: "<<std::setprecision(20)<<lumiScaling_Dani<<endl;
        cout<<std::setprecision(6);
    }
    vector<evn_info> selected_Samuel;
    for (long i = 0; i < tree2->GetEntries(); i++)
    {
        if(i%100000==0) cout<<"number of event: " <<i<<endl;
        tree2->GetEntry(i);

        if(evt->sig.nSigLep < 2) continue;
        if(!evt->sig.isSS) continue;
        if(!evt->sig.JetCut) continue;
        if(evt->sig.isZ) continue;

        if(DoFilter)
        {
            bool DoContinue = false;
            for (unsigned int j = 0; j < sample.filter_first.size(); j++)
            {
                if(i >= sample.filter_first[j] && i <= sample.filter_last[j]) DoContinue = true;
            }
            if(DoContinue) continue;
        }

        //cutflow for My ntuple
        double ttV_SF = 1;
        if(sample.ID == 410155 && evt->sig.nBJet>=3) ttV_SF = 1.34;
        else if(sample.ID >= 410218 && sample.ID <= 410220 && evt->sig.nBJet>=3) ttV_SF = 1.37;
        FillHist_Samuel(hCutflow_Samuel, evt, commonWeight, ttV_SF, isBG, isSig);

        //Fill selected_Samuel
        evn_info evn_temp;
        evn_temp.evn = long(evt->evt.event);
        evn_temp.Samuel_index = i;
        evn_temp.Dani_index = -1;
        evn_temp.nJet = evt->sig.nJet;
        selected_Samuel.push_back(evn_temp);
    }

    //check all variables
    vector<evn_info> missing_Dani;
    vector<evn_info> missing_Samuel;

    //use Samuel event to find Dani event
    unsigned int j = 0;
    for (unsigned int i = 0; i < 0; i++)
    //for (unsigned int i = 0; i < selected_Samuel.size(); i++)
    {
        //check event number
        if(selected_Dani[j].evn != selected_Samuel[i].evn )
        {
            //search for event number
            bool isFound = false;
            if(sample.ID == 361073 && selected_Samuel[i].evn == 1588)
            {
                //Known issue: Two events have the same event number: 1588
                cout<<"Searching event number (special case for known issue): "<<selected_Samuel[i].evn<<endl;
                for (unsigned int k = 0; k < selected_Dani.size(); k++)
                {
                    if(selected_Samuel[i].evn == selected_Dani[k].evn && selected_Samuel[i].nJet == selected_Dani[k].nJet)
                    {
                        isFound = true;
                        j=k;
                        cout<<"Event number: "<<selected_Samuel[i].evn<<" is found. Dani index is "<<selected_Dani[j].Dani_index<<endl;
                        selected_Dani[j].Samuel_index = selected_Samuel[i].Samuel_index;
                        break;
                    }
                }
            }
            else
            {
                //cout<<"Searching event number: "<<selected_Samuel[i].evn<<endl;
                for (unsigned int k = 0; k < selected_Dani.size(); k++)
                {
                    if(selected_Samuel[i].evn == selected_Dani[k].evn)
                    {
                        isFound = true;
                        j=k;
                        //cout<<"Event number: "<<selected_Samuel[i].evn<<" is found. Dani index is "<<selected_Dani[j].Dani_index<<endl;
                        selected_Dani[j].Samuel_index = selected_Samuel[i].Samuel_index;
                        break;
                    }
                }
            }

            if(!isFound)
            {
                //cout<<"Event number: "<<selected_Samuel[i].evn<<" cannot be found. It will be further searched later."<<endl;
                missing_Dani.push_back(selected_Samuel[i]);
                continue;
            }
        }
        else
        {
            selected_Dani[j].Samuel_index = selected_Samuel[i].Samuel_index;
        }

        tree1->GetEntry(selected_Dani[j].Dani_index);
        tree2->GetEntry(selected_Samuel[i].Samuel_index);
        checkAllVariables(evt, sample, isBG, isSig);

        /*
        FillHist_Dani(hCutflow_Dani, isBG);
        FillHist_Samuel(hCutflow_Samuel, evt, commonWeight, ttV_SF, isBG, isSig);
        */

        j++;
    }
    cout<<endl;

    //Only for the method that use Samuel event to find Dani event
    cout<<"The missing event in Samuel ntuple: "<<endl;
    for (unsigned int k = 0; k < 0; k++)
    //for (unsigned int k = 0; k < selected_Dani.size(); k++)
    {
        if(selected_Dani[k].Samuel_index == -1)
        {
            cout<<k<<", "<<selected_Dani[k].Dani_index<<", "<<selected_Dani[k].evn<<endl;
            missing_Samuel.push_back(selected_Dani[k]);
        }
    }
    cout<<endl;

    //use Dani event to find Samuel event
    unsigned int i = 0;
    //for (unsigned int j = 0; j < 0; j++)
    for (unsigned int j = 0; j < selected_Dani.size(); j++)
    {
        //check event number
        if(selected_Dani[j].evn != selected_Samuel[i].evn )
        {
            //search for event number
            bool isFound = false;
            if(sample.ID == 361073 && selected_Dani[j].evn == 1588)
            {
                //Known issue: Two events have the same event number: 1588
                cout<<"Searching event number (special case for known issue): "<<selected_Dani[j].evn<<endl;
                for (unsigned int k = 0; k < selected_Samuel.size(); k++)
                {
                    if(selected_Dani[j].evn == selected_Samuel[k].evn && selected_Dani[j].nJet == selected_Samuel[k].nJet)
                    {
                        isFound = true;
                        i=k;
                        cout<<"Event number: "<<selected_Dani[j].evn<<" is found. Samuel index is "<<selected_Samuel[i].Samuel_index<<endl;
                        selected_Samuel[i].Dani_index = selected_Dani[j].Dani_index;
                        break;
                    }
                }
            }
            else
            {
                //cout<<"Searching event number: "<<evn<<endl;
                for (unsigned int k = 0; k < selected_Samuel.size(); k++)
                {
                    if(selected_Dani[j].evn == selected_Samuel[k].evn)
                    {
                        isFound = true;
                        i=k;
                        //cout<<"Event number: "<<selected_Dani[j].evn<<" is found. Samuel index is "<<selected_Samuel[i].Samuel_index<<endl;
                        selected_Samuel[i].Dani_index = selected_Dani[j].Dani_index;
                        break;
                    }
                }
            }

            if(!isFound)
            {
                //cout<<"Event number: "<<selected_Dani[j].evn<<" cannot be found. It will be further searched later."<<endl;
                missing_Samuel.push_back(selected_Dani[j]);
                continue;
            }
        }
        else
        {
            selected_Samuel[i].Dani_index = selected_Dani[j].Dani_index;
        }

        tree1->GetEntry(selected_Dani[j].Dani_index);
        tree2->GetEntry(selected_Samuel[i].Samuel_index);
        checkAllVariables(evt, sample, isBG, isSig);

        /*
        FillHist_Dani(hCutflow_Dani, isBG);
        FillHist_Samuel(hCutflow_Samuel, evt, commonWeight, ttV_SF, isBG, isSig);
        */

        i++;
    }
    cout<<endl;

    //Only for the method that use Dani event to find Samuel event
    cout<<"The missing event in Dani ntuple: "<<endl;
    int first = -1;
    int last = -1;
    int TotalMissing = 0;
    //for (unsigned int k = 0; k < 0; k++)
    for (unsigned int k = 0; k < selected_Samuel.size(); k++)
    {
        if(selected_Samuel[k].Dani_index == -1)
        {
            TotalMissing++;
            if(first == -1)
            {
                first = k;
                last = k;
                cout<<k<<", "<<selected_Samuel[k].Samuel_index<<", "<<selected_Samuel[k].evn<<endl;
                missing_Dani.push_back(selected_Samuel[k]);
            }
            else
            {
                if(k - last ==1)
                {
                    last = k;
                }
                else
                {
                    if(first < last)
                    {
                        cout<<"continuous events: "<<first<<" - "<<last<<" (Total: "<<(last - first + 1)<<")"<<endl;
                        cout<<"index: "<<selected_Samuel[first].Samuel_index<<" - "<<selected_Samuel[last].Samuel_index<<endl;
                    }
                    first = k;
                    last = k;
                    cout<<k<<", "<<selected_Samuel[k].Samuel_index<<", "<<selected_Samuel[k].evn<<endl;
                    missing_Dani.push_back(selected_Samuel[k]);
                }
            }
        }
    }
    if(first < last)
    {
        cout<<"continuous events: "<<first<<" - "<<last<<" (Total: "<<(last - first + 1)<<")"<<endl;
        cout<<"index: "<<selected_Samuel[first].Samuel_index<<" - "<<selected_Samuel[last].Samuel_index<<endl;
    }
    cout<<endl;

    //search for all event numbers in Samuel ntuple
    cout<<"Further searching all event number: "<<endl;
    cout<<"Total missing event in Samuel ntuple (all): "<<missing_Samuel.size()<<endl;
    if(missing_Samuel.size() != 0)
    {
        vector<evn_info> missing_Samuel_inclusive;
        vector<evn_info> missing_Samuel_exclusive;
        for (unsigned int k = 0; k < tree2->GetEntries(); k++)
        {
            tree2->GetEntry(k);
            
            for (unsigned int m = 0; m < missing_Samuel.size(); m++)
            {
                if(missing_Samuel[m].evn == long(evt->evt.event))
                {
                    cout<<"Event number: "<<missing_Samuel[m].evn<<" is found, but not in Samuel expected selection. Samuel index is "<<k<<endl;
                    missing_Samuel[m].Samuel_index = k;
                    missing_Samuel_inclusive.push_back(missing_Samuel[m]);
 
                    tree1->GetEntry(missing_Samuel[m].Dani_index);
                    checkAllVariables(evt, sample, isBG, isSig);
                    break;
                }
            }
        }
        cout<<endl;

        cout<<"The missing event in Samuel ntuple (exclusive): "<<endl;
        for (unsigned int k = 0; k < missing_Samuel.size(); k++)
        {
            if(missing_Samuel[k].Samuel_index == -1)
            {
                cout<<missing_Samuel[k].Dani_index<<", "<<missing_Samuel[k].evn<<endl;
                missing_Samuel_exclusive.push_back(missing_Samuel[k]);
            }
        }
        cout<<"Total missing event in Samuel ntuple (exclusive): "<<missing_Samuel_exclusive.size()<<endl;
        cout<<endl;
 
        cout<<"Total missing event in Samuel ntuple (inclusive): "<<missing_Samuel_inclusive.size()<<endl;
        for (unsigned int k = 0; k < missing_Samuel_inclusive.size(); k++)
        {
            tree2->GetEntry(missing_Samuel_inclusive[k].Samuel_index);
            cout<<missing_Samuel_inclusive[k].evn<<": ";
            if(evt->sig.isZ) cout<<"isZ = 1 (Samuel), isZ = 0 (Dani), ";
            cout<<"nBaseLep = "<<evt->sig.nBaseLep;
            cout<<": ";
            for(unsigned int m=0;m < evt->leps.size();m++)
            {
                cout<<int(evt->leps[m].ID/1000)<<" ";
            }
            cout<<endl;
        }
        cout<<endl;
    }

    cout<<"Total missing event in Dani ntuple: "<<TotalMissing<<endl;
    for (unsigned int k = 0; k < missing_Dani.size(); k++)
    {
        tree2->GetEntry(missing_Dani[k].Samuel_index);
        cout<<missing_Dani[k].evn<<": ";
        if(!evt->sig.isZ) cout<<"isZ = 0 (Samuel), isZ = 1 (Dani), ";
        cout<<"nBaseLep = "<<evt->sig.nBaseLep;
        cout<<": ";
        for(unsigned int m=0;m < evt->leps.size();m++)
        {
            cout<<int(evt->leps[m].ID/1000)<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    delete evt;
    delete tree2;
}

struct Process
{
    TString path_Dani;
    TString CHAIN_NAME;
    vector<Sample> samples;
};

void RunProcess(Process& process, vector<TString>& cutflowList)
{
    bool isMC = false;
    bool isBG = false;
    bool isSig = false;
    if(process.CHAIN_NAME == "data_nom") {}
    else if(process.CHAIN_NAME.Contains("C1N2")) 
    {
        isMC = true;
        isSig = true;
    }
    else 
    {
        isMC = true;
        isBG = true;
    }

    //read Dani ntuple
    TString path1 = "/eos/user/c/clo/ntuple/Dani_tree/Trees/IncludingLepTruth/";
    path1 += process.path_Dani;
    TChain* tree1 = new TChain(process.CHAIN_NAME);
    tree1->Add(path1.Data());
    cout << "Dani ntuple has events: " <<tree1->GetEntries()<<endl;

    tree1->SetBranchAddress("mcweight", &mcweight);
    tree1->SetBranchAddress("puweight", &puweight);
    tree1->SetBranchAddress("totweight", &totweight);
    tree1->SetBranchAddress("lumiScaling", &lumiScaling);

    tree1->SetBranchAddress("wmu_nom", &wmu_nom);
    tree1->SetBranchAddress("wel_nom", &wel_nom);
    tree1->SetBranchAddress("wjet_nom", &wjet_nom);
    tree1->SetBranchAddress("wpu_nom_bkg", &wpu_nom_bkg);
    tree1->SetBranchAddress("wtrig_nom", &wtrig_nom);
    tree1->SetBranchAddress("wchflip_nom", &wchflip_nom);
    tree1->SetBranchAddress("wttV_nom", &wttV_nom);

    tree1->SetBranchAddress("MCId", &MCId);
    tree1->SetBranchAddress("evn", &evn);

    tree1->SetBranchAddress("NlepBL", &NlepBL);
    tree1->SetBranchAddress("nSigLep", &nSigLep);
    if(isBG) tree1->SetBranchAddress("isTruthLep1", &isTruthLep1);
    if(isBG) tree1->SetBranchAddress("isTruthLep2", &isTruthLep2);
    tree1->SetBranchAddress("nJets20", &nJets20);
    tree1->SetBranchAddress("nBJets20", &nBJets20);
    tree1->SetBranchAddress("Pt_l", &Pt_l);
    tree1->SetBranchAddress("Pt_subl", &Pt_subl);
    tree1->SetBranchAddress("DeltaEtaLep", &DeltaEtaLep);
    tree1->SetBranchAddress("met", &met);
    tree1->SetBranchAddress("mt", &mt);
    tree1->SetBranchAddress("meff", &meff);
    tree1->SetBranchAddress("mljj_comb", &mljj_comb);
    tree1->SetBranchAddress("MT2", &MT2);
    //tree1->SetBranchAddress("", &);

    TH1D* hCutflow_Dani = new TH1D("cutflow_Dani", "cutflow_Dani", cutflowList.size(), 0, cutflowList.size());
    TH1D* hCutflow_Samuel = new TH1D("cutflow_Samuel", "cutflow_Samuel", cutflowList.size(), 0, cutflowList.size());
    for(unsigned int i=0;i<cutflowList.size();i++)    
    {
        hCutflow_Dani->GetXaxis()->SetBinLabel(i+1,cutflowList[i].Data());
        hCutflow_Samuel->GetXaxis()->SetBinLabel(i+1,cutflowList[i].Data());
    }

    for(unsigned int i=0;i<process.samples.size();i++)    
    {
        RunSample(tree1, process.samples[i], hCutflow_Dani, hCutflow_Samuel, isBG, isSig);
    }

    cout<<"Cutflow: Dani, Samuel"<<endl;
    for(unsigned int j=1;j<=cutflowList.size();j++)
    {
        cout<<hCutflow_Samuel->GetXaxis()->GetBinLabel(j)<<": ";
        if(j<=14) cout<<long(hCutflow_Samuel->GetBinContent(j))<<endl;
        else if(j%2 == 1)
        {
            cout<<int(hCutflow_Dani->GetBinContent(j))<<", "<<int(hCutflow_Samuel->GetBinContent(j))<<endl;
            if(int(hCutflow_Dani->GetBinContent(j)) != int(hCutflow_Samuel->GetBinContent(j))) cout<<"Error: unweighed yield is different."<<endl;
        }
        else if(isMC && j%2 == 0)
        {
            cout<<std::fixed<<std::setprecision(10)<<hCutflow_Dani->GetBinContent(j)<<", "<<hCutflow_Samuel->GetBinContent(j)<<endl;
            checkError(hCutflow_Dani->GetBinContent(j), hCutflow_Samuel->GetBinContent(j), "yield", 1e-6);
        }
        else
        {
            cout<<endl;
        }
    }

    delete hCutflow_Dani;
    delete hCutflow_Samuel;
    delete tree1;
}

void add_filter_sig(vector<Process>& processes, int sampleID, int first, int last)
{
    for(unsigned int i=0;i<processes.size();i++)    
    {
        if(!processes[i].CHAIN_NAME.Contains("C1N2") || processes[i].samples[0].ID != sampleID) continue;

        processes[i].samples[0].filter_first.push_back(first);
        processes[i].samples[0].filter_last.push_back(last);
    }
}

int main()
{
    std::vector<TString> cutflowList;
    cutflowList.push_back("AOD");
    cutflowList.push_back("SUSY2");
    cutflowList.push_back("GRL");
    cutflowList.push_back("Trigger");
    cutflowList.push_back("LAr+Tile+SCT+CoreFlag");
    cutflowList.push_back("PV");
    cutflowList.push_back("cosMuon");
    cutflowList.push_back("BadMuon");
    cutflowList.push_back("BadJet");
    cutflowList.push_back(">=2BaseLep");
    cutflowList.push_back(">=2SigLep");
    cutflowList.push_back("SS leptons");
    cutflowList.push_back(">=1 passOR jet");
    cutflowList.push_back(">=1 signal jet");
    cutflowList.push_back("Z veto");
    cutflowList.push_back("Z veto,w");

    cutflowList.push_back("=2BaseLep && =2SigLep");
    cutflowList.push_back("=2BaseLep && =2SigLep,w");
    cutflowList.push_back("pt1");
    cutflowList.push_back("pt1,w");
    cutflowList.push_back("pt2");
    cutflowList.push_back("pt2,w");
    cutflowList.push_back("isTruthLep1");
    cutflowList.push_back("isTruthLep1,w");
    cutflowList.push_back("isTruthLep2");
    cutflowList.push_back("isTruthLep2,w");
    cutflowList.push_back("nBJets20");
    cutflowList.push_back("nBJets20,w");

    cutflowList.push_back("SR1:nJet");
    cutflowList.push_back("SR1:nJet,w");
    cutflowList.push_back("SR1:DeltaEtaLep");
    cutflowList.push_back("SR1:DeltaEtaLep,w");
    cutflowList.push_back("SR1:met");
    cutflowList.push_back("SR1:met,w");
    cutflowList.push_back("SR1:mt");
    cutflowList.push_back("SR1:mt,w");
    cutflowList.push_back("SR1:meff");
    cutflowList.push_back("SR1:meff,w");
    cutflowList.push_back("SR1:mljj");
    cutflowList.push_back("SR1:mljj,w");
    cutflowList.push_back("SR1:MT2");
    cutflowList.push_back("SR1:MT2,w");

    cutflowList.push_back("SR2:nJet");
    cutflowList.push_back("SR2:nJet,w");
    cutflowList.push_back("SR2:DeltaEtaLep");
    cutflowList.push_back("SR2:DeltaEtaLep,w");
    cutflowList.push_back("SR2:met");
    cutflowList.push_back("SR2:met,w");
    cutflowList.push_back("SR2:mt");
    cutflowList.push_back("SR2:mt,w");
    cutflowList.push_back("SR2:meff");
    cutflowList.push_back("SR2:meff,w");
    cutflowList.push_back("SR2:mljj");
    cutflowList.push_back("SR2:mljj,w");
    cutflowList.push_back("SR2:MT2");
    cutflowList.push_back("SR2:MT2,w");

    vector<Process> processes;
    Process process;
    Sample sample;
    sample.filter_first.clear();
    sample.filter_last.clear();

    //WZ
    process.path_Dani = "background.36100.updatedVVVMcTruth_incltHtWZ.root"; process.CHAIN_NAME = "WZ_nom";
    process.samples.clear();
    sample.tag = "VV_CT10"; sample.ID = 361071; sample.XS = 0.042287*0.91; process.samples.push_back(sample);
    sample.tag = "VV_221";  sample.ID = 363491; sample.XS = 4.5877;        process.samples.push_back(sample);
    processes.push_back(process);

    //WW
    process.path_Dani = "background.36100.updatedVVVMcTruth_incltHtWZ.root"; process.CHAIN_NAME = "WW_nom";
    process.samples.clear();
    sample.tag = "VV_CT10"; sample.ID = 361069; sample.XS = 0.025765*0.91; process.samples.push_back(sample);
    sample.tag = "VV_CT10"; sample.ID = 361070; sample.XS = 0.043375*0.91; process.samples.push_back(sample);
    sample.tag = "VV_CT10"; sample.ID = 361077; sample.XS = 0.85492 *0.91; process.samples.push_back(sample);
    processes.push_back(process);

    //ZZ
    process.path_Dani = "background.36100.updatedVVVMcTruth_incltHtWZ.root"; process.CHAIN_NAME = "ZZ_nom";
    process.samples.clear();
    sample.tag = "VV_CT10"; sample.ID = 361072; sample.XS = 0.031496*0.91; process.samples.push_back(sample);
    sample.tag = "VV_CT10"; sample.ID = 361073; sample.XS = 0.02095 *0.91; process.samples.push_back(sample);

    ///*
    sample.tag = "VV_221";  sample.ID = 363490; sample.XS = 1.2557;
    sample.filter_first.push_back(1529811); sample.filter_last.push_back(1967630);
    sample.filter_first.push_back(3939243); sample.filter_last.push_back(4601130);
    process.samples.push_back(sample);
    sample.filter_first.clear();
    sample.filter_last.clear();
    //*/

    processes.push_back(process);

    //ttV
    process.path_Dani = "background.36100.updatedVVVMcTruth_incltHtWZ.root"; process.CHAIN_NAME = "ttV_nom";
    process.samples.clear();
    sample.tag = "ttV"; sample.ID = 407321; sample.XS = 0.000266 *1.34;   process.samples.push_back(sample);
    sample.tag = "ttV"; sample.ID = 410081; sample.XS = 0.0080975*1.2231; process.samples.push_back(sample);
    sample.tag = "ttV"; sample.ID = 410155; sample.XS = 0.54830  *1.10;   process.samples.push_back(sample);
    sample.tag = "ttV"; sample.ID = 410218; sample.XS = 0.036888 *1.12;   process.samples.push_back(sample);

    ///*
    sample.tag = "ttV"; sample.ID = 410219; sample.XS = 0.036895 *1.12;
    sample.filter_first.push_back(7); sample.filter_last.push_back(538779);
    process.samples.push_back(sample);
    sample.filter_first.clear();
    sample.filter_last.clear();
    //*/

    sample.tag = "ttV"; sample.ID = 410220; sample.XS = 0.036599 *1.12;   process.samples.push_back(sample);
    processes.push_back(process);

    //Rare
    process.path_Dani = "background.36100.updatedVVVMcTruth_incltHtWZ.root"; process.CHAIN_NAME = "Rare_nom";
    process.samples.clear();
    sample.tag = "higgs_HerwigppEvtGen"; sample.ID = 341177; sample.XS = 0.5085 * 0.10554;        process.samples.push_back(sample);
    sample.tag = "higgs_HerwigppEvtGen"; sample.ID = 341270; sample.XS = 0.5085 * 0.43929;        process.samples.push_back(sample);
    sample.tag = "higgs_Pythia8EvtGen";  sample.ID = 342284; sample.XS = 1.1021 * 1.25215497686;  process.samples.push_back(sample);
    sample.tag = "higgs_Pythia8EvtGen";  sample.ID = 342285; sample.XS = 0.60072 * 1.44759621787; process.samples.push_back(sample);
    sample.tag = "VVV";  sample.ID = 407311; sample.XS = 0.00010235; process.samples.push_back(sample);
    sample.tag = "VVV";  sample.ID = 407312; sample.XS = 0.00056766; process.samples.push_back(sample);
    sample.tag = "VVV";  sample.ID = 407313; sample.XS = 0.0043684;  process.samples.push_back(sample);
    sample.tag = "VVV";  sample.ID = 407314; sample.XS = 0.015846;   process.samples.push_back(sample);
    sample.tag = "VVV";  sample.ID = 407315; sample.XS = 0.0058366;  process.samples.push_back(sample);
    sample.tag = "multitop_fast";  sample.ID = 304014; sample.XS = 0.00164;            process.samples.push_back(sample);
    sample.tag = "multitop";       sample.ID = 410080; sample.XS = 0.0091622 * 1.0042; process.samples.push_back(sample);
    sample.tag = "Rare";       sample.ID = 343272; sample.XS = 0.7317;   process.samples.push_back(sample);
    sample.tag = "Rare";       sample.ID = 410215; sample.XS = 0.015558; process.samples.push_back(sample);
    processes.push_back(process);

    //Data
    process.path_Dani = "data.36100.root"; process.CHAIN_NAME = "data_nom";
    process.samples.clear();
    sample.tag = "data"; sample.ID = 0; sample.XS = 1;
    sample.filter_first.push_back(16092378); sample.filter_last.push_back(24160509);
    sample.filter_first.push_back(45674365); sample.filter_last.push_back(45674365); //evn: 636665874
    sample.filter_first.push_back(45855936); sample.filter_last.push_back(45855936); //evn: 1993066247
    sample.filter_first.push_back(46054256); sample.filter_last.push_back(46054676);
    sample.filter_first.push_back(46075865); sample.filter_last.push_back(46075865); //evn: 576063822
    process.samples.push_back(sample);
    sample.filter_first.clear();
    sample.filter_last.clear();
    processes.push_back(process);

    //signal
    process.path_Dani = "../C1N2_Wh_hall_36100.root";
    sample.tag = "MCSig";
    {
        ifstream fin;
        fin.open("../plotScripts/dataMCcomparison/SigSample.txt");
        while(!fin.eof())
        {
            fin>>sample.ID;
            //cout<<sample.ID<<endl;
            if(fin.eof()) break;

            TString SampleNameTemp;
            fin>>SampleNameTemp;
            SampleNameTemp.ReplaceAll("MGPy8EG_A14N23LO_","");
            SampleNameTemp.ReplaceAll("_2L7","");
            //cout<<SampleNameTemp<<endl;
            process.CHAIN_NAME = SampleNameTemp + "_nom";

            double skip;
            fin>>skip;
            fin>>skip;

            double XS1;
            double XS2;
            fin>>XS2;

            fin>>XS1;
            //XS2 *= XS1;

            fin>>XS1;
            XS2 *= XS1;

            fin>>XS1;
            XS2 *= XS1;

            sample.XS = XS2;
            //cout<<sample.XS<<endl;

            fin>>skip;
            fin>>skip;

            process.samples.clear();
            process.samples.push_back(sample);
            processes.push_back(process);
        }
        fin.close();
    }

    //use C1N2_Wh_hall_ext36100.root
    for(unsigned int i=0;i<processes.size();i++)    
    {
        if(!processes[i].CHAIN_NAME.Contains("C1N2")) continue;
        if( (processes[i].samples[0].ID >= 393832 && processes[i].samples[0].ID <= 393840) ||
            (processes[i].samples[0].ID >= 393842 && processes[i].samples[0].ID <= 393845) ||
            (processes[i].samples[0].ID >= 393847 && processes[i].samples[0].ID <= 393849) )
        {
            processes[i].path_Dani = "../C1N2_Wh_hall_ext36100.root";
        }
    }

    //add filter
    add_filter_sig(processes, 393841,  6254, 18856);

    for(unsigned int i=0;i<processes.size();i++)    
    {
        //if(processes[i].CHAIN_NAME != "WZ_nom") continue;
        //if(processes[i].CHAIN_NAME != "WW_nom") continue;
        //if(processes[i].CHAIN_NAME != "ZZ_nom") continue;
        //if(processes[i].CHAIN_NAME != "ttV_nom") continue;
        //if(processes[i].CHAIN_NAME != "Rare_nom") continue;
        //if(processes[i].CHAIN_NAME != "data_nom") continue;

        //if(processes[i].CHAIN_NAME != "C1N2_Wh_hall_150p0_0p0_nom") continue;
        if(!processes[i].CHAIN_NAME.Contains("C1N2")) continue;

        RunProcess(processes[i], cutflowList);
    }
}

