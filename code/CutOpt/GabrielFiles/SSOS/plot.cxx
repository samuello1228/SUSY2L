#include "evt2l.C"
#include <TH1.h>
#include <TChain.h>

inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
   if ( (x != n) && (x % (n/100+1) != 0) ) return;

   float ratio    =  x/(float)n;
   unsigned int c =  ratio * w;

   std::cout << setw(3) << (int)(ratio*100) << "% [";
   for (unsigned int x=0; x<c; x++) std::cout << "=";
   for (unsigned int x=c; x<w; x++) std::cout << " ";
   std::cout << "]\r" << flush;
}


struct hists{
	TString name;
	TH1* OS;
	TH1* SS;
	TH1* comp;
	TH1* OSNorm;
	TH1* SSNorm;
	TH1* compNorm;
};

hists hMT2, hL12pt, hMET, hHT, hMTl1, hMTl2, hLLdPhi, hMll, hJetMet_dPhi, hMET_JetPt_R, hl1Pt_JetPt_R;
std::vector<hists> v_hists;
TFile* outFile;
TDirectory* currentDir;
TString outDir;

bool initializeHists(){
	hMT2.name = "mT2";
	hMT2.OS = new TH1D("hmT2_OS", ";mT2[GeV];Events/5 GeV", 20, 0, 100);
	hMT2.SS = new TH1D("hmT2_SS", ";mT2[GeV];Events/5 GeV", 20, 0, 100);
	v_hists.push_back(hMT2);

	hL12pt.name = "L12pt";
	hL12pt.OS = new TH1D("hL12pt_OS", ";l12 pt [GeV];Events/5 GeV", 20, 0, 100);
	hL12pt.SS = new TH1D("hL12pt_SS", ";l12 pt [GeV];Events/5 GeV", 20, 0, 100);
	v_hists.push_back(hL12pt);

	hMET.name = "MET";
	hMET.OS = new TH1D("hMET_OS", ";MET [GeV];Events/5 GeV", 60, 0, 300);
	hMET.SS = new TH1D("hMET_SS", ";MET [GeV];Events/5 GeV", 60, 0, 300);
	v_hists.push_back(hMET);

	hHT.name = "HT";
	hHT.OS = new TH1D("hHT_OS", ";HT [GeV];Events/20 GeV", 50, 0, 1000);
	hHT.SS = new TH1D("hHT_SS", ";HT [GeV];Events/20 GeV", 50, 0, 1000);
	v_hists.push_back(hHT);

	hMTl1.name = "MTl1";
	hMTl1.OS = new TH1D("hMTl1_OS", ";m_{T}^{l1} [GeV];Events/10 GeV", 30, 0, 300);
	hMTl1.SS = new TH1D("hMTl1_SS", ";m_{T}^{l1} [GeV];Events/10 GeV", 30, 0, 300);
	v_hists.push_back(hMTl1);

	hMTl2.name = "MTl2";
	hMTl2.OS = new TH1D("hMTl2_OS", ";m_{T}^{l2} [GeV];Events/10 GeV", 30, 0, 300);
	hMTl2.SS = new TH1D("hMTl2_SS", ";m_{T}^{l2} [GeV];Events/10 GeV", 30, 0, 300);
	v_hists.push_back(hMTl2);

	hLLdPhi.name = "LLdPhi";
	hLLdPhi.OS = new TH1D("hLLdPhi_OS", ";#Delta #phi_{ll};Events/0.2512", 25, -3.14, 3.14 );
	hLLdPhi.SS = new TH1D("hLLdPhi_SS", ";#Delta #phi_{ll};Events/0.2512", 25, -3.14, 3.14 );
	v_hists.push_back(hLLdPhi);

	hMll.name = "Mll";
	hMll.OS = new TH1D("hMll_OS", ";m_{ll} [GeV];Events/5 GeV", 30, 0, 150);
	hMll.SS = new TH1D("hMll_SS", ";m_{ll} [GeV];Events/5 GeV", 30, 0, 150);
	v_hists.push_back(hMll);

	hJetMet_dPhi.name = "JetMet_dPhi";
	hJetMet_dPhi.OS = new TH1D("hJetMet_dPhi_OS", ";#Delta #phi_{ll};Events/0.2512", 25, -3.14, 3.14);
	hJetMet_dPhi.SS = new TH1D("hJetMet_dPhi_SS", ";#Delta #phi_{ll};Events/0.2512", 25, -3.14, 3.14);
	v_hists.push_back(hJetMet_dPhi);

	hMET_JetPt_R.name = "MET_JetPt_R";
	hMET_JetPt_R.OS = new TH1D("hMET_JetPt_R_OS", ";MET/p_{T}(leading jet);Events/0.25", 20, 0, 5);
	hMET_JetPt_R.SS = new TH1D("hMET_JetPt_R_SS", ";MET/p_{T}(leading jet);Events/0.25", 20, 0, 5);
	v_hists.push_back(hMET_JetPt_R);

	hl1Pt_JetPt_R.name = "l1Pt_JetPt_R";
	hl1Pt_JetPt_R.OS = new TH1D("hl1Pt_JetPt_R_OS", ";p_{T}^{l1}/p_T(leading jet);Events/0.25", 12, 0, 3);
	hl1Pt_JetPt_R.SS = new TH1D("hl1Pt_JetPt_R_SS", ";p_{T}^{l1}/p_T(leading jet);Events/0.25", 12, 0, 3);
	v_hists.push_back(hl1Pt_JetPt_R);

	return true;
}
bool destructHists(){
	for(auto h : v_hists){
		h.OS->Write(); delete h.OS; h.OS = 0;
		h.SS->Write(); delete h.SS; h.SS = 0;
		h.comp->Write(); delete h.comp; h.comp = 0;
		h.OSNorm->Write(); delete h.OSNorm; h.OSNorm = 0;
		h.SSNorm->Write(); delete h.SSNorm; h.SSNorm = 0;
		h.compNorm->Write(); delete h.compNorm; h.compNorm = 0;
	}
	return true;
}

TChain* loadData(TString fileList){
   //return a TChain linked to the data files
   TChain* tc = new TChain("evt2l");

   if (fileList){
      std::ifstream inFiles(fileList);
      if (inFiles.is_open()){
         for( std::string line; getline( inFiles, line ); ){
         	if(tc->Add(line.c_str())) std::cout << "Added " << line << std::endl;
         }
         inFiles.close();
      }
   }
   return tc;
}
evt2l* evt2lTree;

bool hasISR()
{
	for(int i=0; i<evt2l::kMaxjets; i++)
	{
		if(evt2lTree->jets_pt[i]>20 && evt2lTree->jets_eta[i]<2.4) return true;
	}
	return false;
}

bool selectChannel(int channel)
{
	// Trigger
	if (evt2lTree->sig_trigCode==0) return false;

	// Two leptons
	if (evt2lTree->leps_!=2) return false;

	// Signal leptons
	if (!(evt2lTree->leps_lFlag[0] & 2)/2 || !(evt2lTree->leps_lFlag[1] & 2)/2) return false;

	// Channels:
	//   noISR:  0=ee,  1=emu,  2=mumu,  3=combFlav
	//     ISR: 10=ee, 11=emu, 12=mumu, 13=combFlav
	// combISR: 20=ee, 21=emu, 22=mumu, 23=combFlav
	// useISR = True if int(ch/10)==1 else False

	if(channel/10 == 0 &&  hasISR()) return false;
	if(channel/10 == 1 && !hasISR()) return false;

	if(channel%10==0 && 
		!(int(abs(evt2lTree->leps_ID[0])/1000) == 11 && int(abs(evt2lTree->leps_ID[1])/1000) == 11)) return false;
	if(channel%10==1 &&
		(int(abs(evt2lTree->leps_ID[0])/1000) == int(abs(evt2lTree->leps_ID[1])/1000))) return false;
	if(channel%10==2 &&
		(int(abs(evt2lTree->leps_ID[0])/1000) == 13 && int(abs(evt2lTree->leps_ID[1])/1000) == 13)) return false;

	// Zmass veto for same flavor
	if((int(abs(evt2lTree->leps_ID[0])/1000) == int(abs(evt2lTree->leps_ID[1])/1000)) 
		&& fabs(evt2lTree->l12_m - 91.1876)<=10) return false;

	return true;
}

void draw()
{
	// Hide info box of plot
	gStyle->SetOptStat(0);

	// Split pad
	TCanvas* c = new TCanvas(1);
	c->Divide(1,2);
	c->cd(1)->SetPad(0, 0.25, 1, 1);
	c->cd(2)->SetPad(0, 0, 1, 0.25);

	for(hists h : v_hists){

		///// NORMALIZED //////
		c->cd(1);
		h.OS->SetTitle("Normalized OS and SS "+h.name+" distributions");
		h.OS->SetLineColor(kRed);
		h.SS->SetLineColor(kBlue);		

		TLegend l(0.75, 0.75, 0.9, 0.9);
		l.AddEntry(h.OS, "OS", "l");
		l.AddEntry(h.SS, "SS", "l");
		
		h.compNorm = (TH1*) h.OS->Clone(h.name+"_norm_ratio");
		h.compNorm->Divide(h.OS->DrawNormalized(), h.SS->DrawNormalized("same"));
		l.Draw();
		h.compNorm->SetTitle(";;OS/SS");
		c->cd(2)->SetGridy();
		h.compNorm->GetXaxis()->SetLabelSize(0);
		h.compNorm->GetYaxis()->SetTitleSize(0.1);
		h.compNorm->GetYaxis()->SetTitleOffset(0.4);
		h.compNorm->GetYaxis()->SetLabelSize(0.08);
		h.compNorm->SetMarkerStyle(20);
		h.compNorm->SetMarkerColor(kBlack);
		h.compNorm->SetLineColor(kBlack);
		h.compNorm->Draw("p");

		// c->Print(outDir+"/"+h.name+"_norm.pdf", "Title:"+h.name);

		///// REGULAR //////

		// Comparison
		c->cd(2)->SetGridy();

		h.comp = (TH1*) h.OS->Clone(h.name+"_ratio");
		h.comp->Divide(h.SS);
		h.comp->SetTitle(";;OS/SS");

		h.comp->GetXaxis()->SetLabelSize(0);
		h.comp->GetYaxis()->SetTitleSize(0.1);
		h.comp->GetYaxis()->SetTitleOffset(0.4);
		h.comp->GetYaxis()->SetLabelSize(0.08);
		h.comp->SetMarkerStyle(20);
		h.comp->SetMarkerColor(kBlack);
		h.comp->Draw("p");

		c->cd(1);
		h.OS->SetLineColor(kRed);
		h.SS->SetLineColor(kBlue);

		h.OS->GetYaxis()->SetRangeUser(0, max(h.SS->GetMaximum(),h.OS->GetMaximum())*1.1);
		h.OS->SetTitle("OS and SS "+h.name+" distributions");
		h.OS->Draw();
		h.SS->Draw("same");

		l.Draw();

		// c->Print(outDir+"/"+h.name+".pdf", "Title:"+h.name);
		h.OS->SetTitle("");
	}
	// currentDir->Write();
}

int plot1(TString dm, int channel)
{
	outDir="dm" + dm + "_Channel" + to_string(channel);
	cout << outDir << endl;
	currentDir = outFile->mkdir(outDir);
	currentDir->cd();
	gSystem->mkdir(outDir); 
	if(!initializeHists()) return -1;

	TString fileList = "/home/ggallard/Documents/SUSY2L/code/CutOpt/GabrielFiles/sigFiles_"+dm+".txt";
	evt2lTree = new evt2l(loadData(fileList));

	long long nEntries = evt2lTree->fChain->GetEntries();

	for (long long i =0; i<nEntries; i++){
		loadbar(i, nEntries);
		evt2lTree->GetEntry(i);
		// if(!selectChannel(channel)) continue;

		if((evt2lTree->leps_ID[0]>0) == (evt2lTree->leps_ID[1]>0))
		{
			hMT2.SS->Fill(evt2lTree->sig_mT2);
			hL12pt.SS->Fill(evt2lTree->l12_pt);
			hMET.SS->Fill(evt2lTree->sig_MetRel);
			hHT.SS->Fill(evt2lTree->sig_HT);
			hMTl1.SS->Fill(evt2lTree->leps_mT[0]);
			hMTl2.SS->Fill(evt2lTree->leps_mT[1]);
			hLLdPhi.SS->Fill(evt2lTree->l12_dPhi);
			hMll.SS->Fill(evt2lTree->l12_m);

			if(hasISR()){
				hJetMet_dPhi.SS->Fill(evt2lTree->jets_MET_dPhi[0]);
				hMET_JetPt_R.SS->Fill(evt2lTree->sig_MetRel/evt2lTree->jets_pt[0]);
				hl1Pt_JetPt_R.SS->Fill(evt2lTree->leps_pt[0]/evt2lTree->jets_pt[0]);
			}
		}
		else
		{
			hMT2.OS->Fill(evt2lTree->sig_mT2);
			hL12pt.OS->Fill(evt2lTree->l12_pt);
			hMET.OS->Fill(evt2lTree->sig_MetRel);
			hHT.OS->Fill(evt2lTree->sig_HT);
			hMTl1.OS->Fill(evt2lTree->leps_mT[0]);
			hMTl2.OS->Fill(evt2lTree->leps_mT[1]);
			hLLdPhi.OS->Fill(evt2lTree->l12_dPhi);
			hMll.OS->Fill(evt2lTree->l12_m);

			if(hasISR()){
				hJetMet_dPhi.OS->Fill(evt2lTree->jets_MET_dPhi[0]);
				hMET_JetPt_R.OS->Fill(evt2lTree->sig_MetRel/evt2lTree->jets_pt[0]);
				hl1Pt_JetPt_R.OS->Fill(evt2lTree->leps_pt[0]/evt2lTree->jets_pt[0]);
			}
		}
	}
	// draw();
	destructHists();
	return 0;
}

int plot(){
	outFile = new TFile("SSOS.root", "recreate");

	TString dm[] = {"all", "20", "50", "100", "+"};
	int channel[] = {0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23};

	// for(int i : dm) for(int j : channel) plot1(i, j);
	plot1("20",23);
	return 0;
}


