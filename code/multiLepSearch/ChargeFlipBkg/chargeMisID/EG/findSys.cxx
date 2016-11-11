#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>
#include <sstream>

#include "TFile.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TString.h"
#include "TList.h"
#include "TColor.h"
#include "TStyle.h"

TString nomFileName = "80.000000_100.000000_20.000000_20.000000_DATA.root";
TString nomHistoName = "80.0_100.0_20.0_20.0_DATA_misid";
TString outFileName = "rates_wSys.root";

bool drawEtaPt=true;
bool drawPtEta=true;

Color_t nomColor = kRed;
Color_t colors[] = {kGreen+3, kBlue, kCyan, kViolet, kGray+1, kCyan+4};

TFile* 	outFile;
TH2*	hVar;
TH2*	hStat;
TH2*	hBias;
TH2*  hAllSys;
TH2*  hSysAndStat;

double calcTotErr(double stat, double sys)
{
	return sqrt(stat*stat + sys*sys);
}

void findSys()
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(kLightTerrain);

	outFile = new TFile(outFileName, "recreate");

	//================================================
	// STATISTICAL UNCERTAINTY
	//================================================
	{
		TFile *nomFile = TFile::Open(nomFileName, "read");
		if(!nomFile)
		{
			cout << "Could not open nominal file " << nomFileName << endl;
			return;
		}

		hStat = (TH2*) nomFile->Get(nomHistoName);
		if(!hStat)
		{
			cout << "Could not find nominal histogram " << nomHistoName << endl
						<< " in " << nomFileName << endl;
			return;
		}
		hStat->SetName("hFlipProb_data");
		hStat->SetDirectory(outFile);
		nomFile->Close();
	}

	//================================================
	// SYSTEMATIC UNCERTAINTY
	//================================================
	std::vector<TH2*> histos;
	std::vector<TString> histLabels;
	TH2 *hMCTruth, *hMCLH;

	// Get all data histograms
	{
		// Get list of filenames 
		TString dirName = gSystem->pwd();
		TSystemDirectory sysDir(dirName, dirName);
		TList *fList = sysDir.GetListOfFiles();
		std::vector<TString> fNames;
		if(fList)
		{
			TSystemFile *file;
			TString fName;
			TIter next(fList);
			while ((file=(TSystemFile*)next())) {
				fName = file->GetName();
				if (!file->IsDirectory() && fName.EndsWith("DATA.root") && fName != nomFileName) 
				{
					cout << fName.Data() << endl;
					fNames.push_back(fName);
				}
			}
		}

		// Get histogram from files
		for(auto fName : fNames)
		{
			TFile *f = TFile::Open(fName);
			TIter next(f->GetListOfKeys());
			TKey *key = 0;
			while((key = (TKey*)next()))
			{
				// Get histogram
				TClass *cl = gROOT->GetClass(key->GetClassName());
				if(!cl->InheritsFrom("TH2")) continue;
				TH2* h = (TH2*) key->ReadObj();
				h->SetDirectory(0);
				histos.push_back(h);

				// Get label
				TString histName = h->GetName();
				stringstream sStream;
				sStream << histName;
				int ml, mr, sb;		
				string tmp;
				getline(sStream, tmp, '_'); ml = stoi(tmp);
				getline(sStream, tmp, '_'); mr = stoi(tmp);
				getline(sStream, tmp, '_'); sb = stoi(tmp);
				char tempLabel[100];
				sprintf(tempLabel, "%d<m<%d,SB=%d", ml, mr, sb);
				TString histLabel(tempLabel);
				histLabels.push_back(histLabel);
				cout << histLabel << endl;
			}
			f->Close();
		}
	}

	// Get MC histograms
	{
		TFile* fMCTruth = TFile::Open("MCTruth_80_100.root");
		hMCTruth = (TH2*) fMCTruth->Get("hMCTruthRate");
		hMCTruth->SetName("hFlipProb_MCtruth");

		TFile* fMCLH = TFile::Open("80.000000_100.000000_0.000000_0.000000_MC.root");
		hMCLH = (TH2*) fMCLH->Get("80.0_100.0_0.0_0.0_MC_misid"); 
		hMCLH->SetName("hFlipProb_MCLH");
		hMCTruth->SetDirectory(outFile); 
		hMCLH->SetDirectory(outFile);
	}
	hVar = (TH2*) hStat->Clone("hFlipProb_var");
	hVar->SetDirectory(outFile); 
	hBias = (TH2*) hStat->Clone("hFlipProb_bias");
	hBias->SetDirectory(outFile);
	hAllSys = (TH2*) hStat->Clone("hFlipProb_AllSys");
	hAllSys->SetDirectory(outFile);
	hSysAndStat = (TH2*) hStat->Clone("hFlipProb_SysAndStat");
	hSysAndStat->SetDirectory(outFile);

	// cout << "Finding sys error"  << endl;
	int nEtaBins(hStat->GetNbinsX()), nPtBins(hStat->GetNbinsY());
	for(int etaBin=1; etaBin<=nEtaBins; etaBin++){
	for(int ptBin =1;  ptBin<=nPtBins ;  ptBin++){
		// cout << ptBin + etaBin*nPtBins << endl;
		double nom = hVar->GetBinContent(etaBin, ptBin);

		// Find maximum difference from nominal measurement
		std::vector<double> var;
		for(auto h : histos) var.push_back(fabs(nom - h->GetBinContent(etaBin,ptBin)));
		double sysErr = var[std::distance(var.begin(), std::max_element(var.begin(), var.end()))];
		hVar->SetBinError(etaBin, ptBin, sysErr);

		// Find method bias
		hBias->SetBinError(etaBin, ptBin, fabs(hMCLH->GetBinContent(etaBin, ptBin)-hMCTruth->GetBinContent(etaBin, ptBin)));

		// Stat+sys err. See calcTotErr()
		double totalSys = calcTotErr(hVar->GetBinError(etaBin, ptBin), hBias->GetBinError(etaBin, ptBin));
		hAllSys->SetBinError(etaBin, ptBin, totalSys);

		double totalErr = totalSys + hStat->GetBinError(etaBin, ptBin);
		hSysAndStat->SetBinError(etaBin, ptBin, totalErr);
	}}	

	// {
	// 	TCanvas c;
	// 	histos[1]->Draw("colz");
	// 	c.Print("hStat.pdf");
	// }

	//================================================
	// PLOTTING
	//================================================
	const double *PtEdges  = hStat->GetYaxis()->GetXbins()->GetArray();
	const double *EtaEdges = hStat->GetXaxis()->GetXbins()->GetArray();

	TCanvas cAll("cAll", "cAll", 600, 600);
	TCanvas cErr("cErr", "cErr", 600, 600);
	if(drawEtaPt)
	{
		TCanvas *cEtaPt[nPtBins];
		TLegend *lEtaPt[nPtBins];

		TCanvas *cEtaPtErr[nPtBins];
		TLegend *lEtaPtErr[nPtBins];

		TString sPdfFilename = "EtaPt.pdf";
		TString sPdfFilenameErr = "EtaPtErr.pdf";
		cAll.Print(sPdfFilename+"[");
		cErr.Print(sPdfFilenameErr+"[");

		for(int ptBin=0; ptBin<nPtBins; ptBin++)
		{
			// ============== All variations =========
			cEtaPt[ptBin] = (TCanvas*) cAll.Clone(("pt"+to_string(PtEdges[ptBin])).c_str());
			cEtaPt[ptBin]->SetLogy();

			lEtaPt[ptBin] = new TLegend(0.1, 0.5, 0.45, 0.9);
			char sLegendHeader[100];
			sprintf(sLegendHeader, "%d < p_{T} < %d", (int) PtEdges[ptBin], (int) PtEdges[ptBin+1]);
			lEtaPt[ptBin]->SetHeader(sLegendHeader);

			TH1* hNom = 0;
			char hName[100];
			sprintf(hName, "%s_%d_pt_%d", hStat->GetName(), (int) PtEdges[ptBin], (int) PtEdges[ptBin+1]);
			hNom = hStat->ProjectionX(hName, ptBin+1, ptBin+1, "e");
			hNom->SetLineColor(nomColor); hNom->SetMarkerColor(nomColor);
			hNom->Draw("e1");

			lEtaPt[ptBin]->AddEntry(hNom, "#bf{80<m<100,SB=20}", "lpf");

			for(int i=0; i<histos.size(); i++)
			{
				char hName[200];
				sprintf(hName, "%s_%d_pt_%d", histos[i]->GetName(), (int) PtEdges[ptBin], (int) PtEdges[ptBin+1]);
				TH1* h = histos[i]->ProjectionX(hName, ptBin+1, ptBin+1, "e");
				h->SetLineColor(colors[i]); h->SetMarkerColor(colors[i]);
				h->Draw("e1 same");
				lEtaPt[ptBin]->AddEntry(h, histLabels[i], "lpf");
			}

			lEtaPt[ptBin]->Draw();

			cEtaPt[ptBin]->Print(sPdfFilename, cEtaPt[ptBin]->GetName());

			// ============= Sys err and stat err ============= 
			cEtaPtErr[ptBin] = (TCanvas*) cAll.Clone(("eterr"+to_string(PtEdges[ptBin])).c_str());
			cEtaPtErr[ptBin]->SetLogy(); 
			lEtaPtErr[ptBin] = new TLegend(0.55, 0.1, 0.9, 0.3);
			lEtaPtErr[ptBin]->SetHeader(sLegendHeader);

			sprintf(hName, "%s_%.2f_eta_%.2f", hVar->GetName(), PtEdges[ptBin], PtEdges[ptBin+1]);
			TH1* hVarProj = hVar->ProjectionX(hName, ptBin+1, ptBin+1, "e");
			hVarProj->SetLineColor(nomColor); hVarProj->SetMarkerColor(nomColor);

			sprintf(hName, "%s_%.2f_eta_%.2f", hBias->GetName(), PtEdges[ptBin], PtEdges[ptBin+1]);
			TH1* hBiasProj = hBias->ProjectionX(hName, ptBin+1, ptBin+1, "e");
			hBiasProj->SetLineColor(colors[1]); hBiasProj->SetMarkerColor(colors[1]);

			sprintf(hName, "%s_%.2f_eta_%.2f", hAllSys->GetName(), PtEdges[ptBin], PtEdges[ptBin+1]);
			TH1* hAllSysProj = hAllSys->ProjectionX(hName, ptBin+1, ptBin+1, "e");
			hAllSysProj->SetLineColor(colors[3]); hAllSysProj->SetMarkerColor(colors[3]);

			sprintf(hName, "%s_%.2f_eta_%.2f", hSysAndStat->GetName(), PtEdges[ptBin], PtEdges[ptBin+1]);
			TH1* hSysAndStatProj = hSysAndStat->ProjectionX(hName, ptBin+1, ptBin+1, "e");
			hSysAndStatProj->SetLineColor(colors[2]); hSysAndStatProj->SetMarkerColor(colors[2]);

			hSysAndStatProj->Draw("e1 "); lEtaPtErr[ptBin]->AddEntry(hSysAndStatProj, "Stat+sys", "lpf");
			hAllSysProj->Draw("e1 same"); lEtaPtErr[ptBin]->AddEntry(hAllSysProj, "#oplus systematics", "lpf");
			// hNom->Draw("e1 same"); lEtaPtErr[ptBin]->AddEntry(hNom, "Statistical", "lpf");
			hVarProj->Draw("e1 same"); lEtaPtErr[ptBin]->AddEntry(hVarProj, "Variation", "lpf");
			hBiasProj->Draw("e1 same"); lEtaPtErr[ptBin]->AddEntry(hBiasProj, "Method bias", "lpf");

			lEtaPtErr[ptBin]->Draw();
			cEtaPtErr[ptBin]->Print(sPdfFilenameErr, cEtaPtErr[ptBin]->GetName());
		}
		cAll.Print(sPdfFilename+"]");
		cErr.Print(sPdfFilenameErr+"]");
	}

	if(drawPtEta)
	{
		TCanvas *cPtEta[nEtaBins];
		TLegend *lPtEta[nEtaBins];

		TCanvas *cPtEtaErr[nEtaBins];
		TLegend *lPtEtaErr[nEtaBins];

		TString sPdfFilename = "PtEta.pdf";
		TString sPdfFilenameErr = "PtEtaErr.pdf";
		cAll.Print(sPdfFilename+"[");
		cErr.Print(sPdfFilenameErr+"[");

		for(int etaBin=0; etaBin<nEtaBins; etaBin++)
		{
			// ================= Variations =================
			cPtEta[etaBin] = (TCanvas*) cAll.Clone(("eta"+to_string(EtaEdges[etaBin])).c_str());
			cPtEta[etaBin]->SetLogy();
			cPtEta[etaBin]->SetLogx();

			lPtEta[etaBin] = new TLegend(0.55, 0.5, 0.9, 0.1);
			char sLegendHeader[100];
			sprintf(sLegendHeader, "%.2f < |#eta| < %.2f", EtaEdges[etaBin], EtaEdges[etaBin+1]);
			lPtEta[etaBin]->SetHeader(sLegendHeader);

			TH1* hNom = 0;
			char hName[100];
			sprintf(hName, "%s_%.2f_pt_%.2f", hStat->GetName(), EtaEdges[etaBin], EtaEdges[etaBin+1]);
			hNom = hStat->ProjectionY(hName, etaBin+1, etaBin+1, "e");
			hNom->SetLineColor(nomColor); hNom->SetMarkerColor(nomColor);
			hNom->Draw("e1");

			lPtEta[etaBin]->AddEntry(hNom, "#bf{80<m<100,SB=20}", "lpf");

			for(int i=0; i<histos.size(); i++)
			{
				char hName[200];
				sprintf(hName, "%s_%.2f_pt_%.2f", histos[i]->GetName(), EtaEdges[etaBin], EtaEdges[etaBin+1]);
				TH1* h = histos[i]->ProjectionY(hName, etaBin+1, etaBin+1, "e");
				h->SetLineColor(colors[i]); h->SetMarkerColor(colors[i]);
				h->Draw("e1 same");
				lPtEta[etaBin]->AddEntry(h, histLabels[i], "lpf");
			}

			lPtEta[etaBin]->Draw();

			cPtEta[etaBin]->Print(sPdfFilename, cPtEta[etaBin]->GetName());

			// ================= Sys and stat error =================
			cPtEtaErr[etaBin] = (TCanvas*) cAll.Clone(("eterr"+to_string(EtaEdges[etaBin])).c_str());
			cPtEtaErr[etaBin]->SetLogy(); cPtEtaErr[etaBin]->SetLogx();
			lPtEtaErr[etaBin] = new TLegend(0.55, 0.1, 0.9, 0.3);
			lPtEtaErr[etaBin]->SetHeader(sLegendHeader);

			sprintf(hName, "%s_%.2f_eta_%.2f", hVar->GetName(), EtaEdges[etaBin], EtaEdges[etaBin+1]);
			TH1* hVarProj = hVar->ProjectionY(hName, etaBin+1, etaBin+1, "e");
			hVarProj->SetLineColor(nomColor); hVarProj->SetMarkerColor(nomColor);

			sprintf(hName, "%s_%.2f_eta_%.2f", hBias->GetName(), EtaEdges[etaBin], EtaEdges[etaBin+1]);
			TH1* hBiasProj = hBias->ProjectionY(hName, etaBin+1, etaBin+1, "e");
			hBiasProj->SetLineColor(colors[1]); hBiasProj->SetMarkerColor(colors[1]);

			sprintf(hName, "%s_%.2f_eta_%.2f", hAllSys->GetName(), EtaEdges[etaBin], EtaEdges[etaBin+1]);
			TH1* hAllSysProj = hAllSys->ProjectionY(hName, etaBin+1, etaBin+1, "e");
			hAllSysProj->SetLineColor(colors[3]); hAllSysProj->SetMarkerColor(colors[3]);

			sprintf(hName, "%s_%.2f_eta_%.2f", hSysAndStat->GetName(), EtaEdges[etaBin], EtaEdges[etaBin+1]);
			TH1* hSysAndStatProj = hSysAndStat->ProjectionY(hName, etaBin+1, etaBin+1, "e");
			hSysAndStatProj->SetLineColor(colors[2]); hSysAndStatProj->SetMarkerColor(colors[2]);

			hSysAndStatProj->Draw("e1 "); lPtEtaErr[etaBin]->AddEntry(hSysAndStatProj, "Stat+sys", "lpf");
			hAllSysProj->Draw("e1 same"); lPtEtaErr[etaBin]->AddEntry(hAllSysProj, "#oplus systematics", "lpf");
			// hNom->Draw("e1 same"); lPtEtaErr[etaBin]->AddEntry(hNom, "Statistical", "lpf");
			hVarProj->Draw("e1 same"); lPtEtaErr[etaBin]->AddEntry(hVarProj, "Variation", "lpf");
			hBiasProj->Draw("e1 same"); lPtEtaErr[etaBin]->AddEntry(hBiasProj, "Method bias", "lpf");

			lPtEtaErr[etaBin]->Draw();
			cPtEtaErr[etaBin]->Print(sPdfFilenameErr);
		}
		cAll.Print(sPdfFilename+"]");
		cErr.Print(sPdfFilenameErr+"]");
	}
	//================================================
	// CLEANUP
	//================================================

	outFile->Write();

}
