#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>

float Add(float x1,float x2,float x3)
{
    return sqrt(x1*x1+x2*x2+x3*x3);
}

void GetFakeHist()
{
    TFile* file = new TFile("RealFakeLepEff_WhSS_Aug2017.root","RECREATE");
    
    //fake eff electron
    {
        float bin_pt[] = {20,25,30,40};
        int nBin_pt = 3;
        float bin_eta[] = {0,1.35,2.47};
        int nBin_eta = 2;
        TH2F* h1 = new TH2F("h1","h1",nBin_pt,bin_pt,nBin_eta,bin_eta);
        
        h1->SetBinContent(1,1,0.0887); h1->SetBinError(1,1, Add(0.0152,0.0116,0.0058));
        h1->SetBinContent(2,1,0.1012); h1->SetBinError(2,1, Add(0.0232,0.0012,0.0111));
        h1->SetBinContent(3,1,0.1482); h1->SetBinError(3,1, Add(0.0256,0.0120,0.0155));
        h1->SetBinContent(4,1,0.2169); h1->SetBinError(4,1, Add(0.0329,0.0146,0.0344));
        
        h1->SetBinContent(1,2,0.1112); h1->SetBinError(1,2, Add(0.0223,0.0035,0.0045));
        h1->SetBinContent(2,2,0.1094); h1->SetBinError(2,2, Add(0.0263,0.0070,0.0052));
        h1->SetBinContent(3,2,0.0980); h1->SetBinError(3,2, Add(0.0212,0.0120,0.0053));
        h1->SetBinContent(4,2,0.0890); h1->SetBinError(4,2, Add(0.0156,0.0010,0.0100));
        
        file->cd();
        h1->Write("FakeeEff");
        delete h1;
    }
    
    //fake eff muon
    {
        float bin_pt[] = {20,25,30,40};
        int nBin_pt = 3;
        float bin_eta[] = {0,2.4};
        int nBin_eta = 1;
        TH2F* h1 = new TH2F("h1","h1",nBin_pt,bin_pt,nBin_eta,bin_eta);
        
        h1->SetBinContent(1,1,0.1388); h1->SetBinError(1,1, Add(0.0167,0.0082,0.0083));
        h1->SetBinContent(2,1,0.1698); h1->SetBinError(2,1, Add(0.0278,0.0045,0.0157));
        h1->SetBinContent(3,1,0.1204); h1->SetBinError(3,1, Add(0.0299,0.0025,0.0354));
        h1->SetBinContent(4,1,0.2252); h1->SetBinError(4,1, Add(0.0517,0.0086,0.0939));
        
        file->cd();
        h1->Write("FakeuEff");
        delete h1;
    }
    
    TFile* fPeter = new TFile("RLE_Histos.root","READ");
    //real eff electron
    {
        float bin_eta[] = {0,0.8,1.37,1.52,2.0,2.47};
        int nBin_eta = 5;
        
        TH1F* hPeter1 = (TH1F*) fPeter->Get("hEtaSigElStd_1");
        
        int nBin_pt = hPeter1->GetXaxis()->GetXbins()->GetSize() -1;
        float bin_pt[nBin_pt +1];
        for(int j=0;j<=nBin_pt;j++)
        {
            cout<<hPeter1->GetBinLowEdge(j+1)<<endl;
            bin_pt[j] = hPeter1->GetBinLowEdge(j+1);
        }
        cout<<endl;
        
        TH2F* h1 = new TH2F("h1","h1",nBin_pt,bin_pt,nBin_eta,bin_eta);
        for(int i=1;i<=nBin_eta;i++)
        {
            TString HistName = "hEtaSigElStd_";
            HistName += TString::Itoa(i,10);
            cout<<HistName.Data()<<endl;
            TH1F* hPeter = (TH1F*) fPeter->Get(HistName.Data());
            
            for(int j=1;j<=nBin_pt+1;j++)
            {
                cout<<hPeter->GetBinContent(j)<<endl;
                h1->SetBinContent(j,i,hPeter->GetBinContent(j));
                h1->SetBinError(j,i,hPeter->GetBinError(j));
            }
            cout<<endl;
        }
        file->cd();
        h1->Write("RealeEff");
        delete h1;
    }
    
    //real eff muon
    {
        float bin_eta[] = {0,0.6,1.2,1.8,2.5};
        int nBin_eta = 4;
        
        TH1F* hPeter1 = (TH1F*) fPeter->Get("hEtaSigMuStd_1");
        
        int nBin_pt = hPeter1->GetXaxis()->GetXbins()->GetSize() -1;
        float bin_pt[nBin_pt +1];
        for(int j=0;j<=nBin_pt;j++)
        {
            cout<<hPeter1->GetBinLowEdge(j+1)<<endl;
            bin_pt[j] = hPeter1->GetBinLowEdge(j+1);
        }
        cout<<endl;
        
        TH2F* h1 = new TH2F("h1","h1",nBin_pt,bin_pt,nBin_eta,bin_eta);
        for(int i=1;i<=nBin_eta;i++)
        {
            TString HistName = "hEtaSigMuStd_";
            HistName += TString::Itoa(i,10);
            cout<<HistName.Data()<<endl;
            TH1F* hPeter = (TH1F*) fPeter->Get(HistName.Data());
            
            for(int j=1;j<=nBin_pt+1;j++)
            {
                cout<<hPeter->GetBinContent(j)<<endl;
                h1->SetBinContent(j,i,hPeter->GetBinContent(j));
                h1->SetBinError(j,i,hPeter->GetBinError(j));
            }
            cout<<endl;
        }
        file->cd();
        h1->Write("RealuEff");
        delete h1;
    }
    
    delete fPeter;
    delete file;
}
