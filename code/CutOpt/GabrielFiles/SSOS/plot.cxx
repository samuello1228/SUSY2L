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

hists hMT2, hL12pt, hMET, hHT, hMTl1, hMTl2, hLLdPhi, hMll, hJetMet_dPhi, hMET_JetPt_R, hl1Pt_JetPt_R, hParents, hParents1, hParents2;
std::vector<hists> v_hists;
TFile* outFile;
TDirectory* currentDir;
TString outDir;
double nOS, nSS;

bool initializeHists(){
	hMT2.name = "mT2";
	hMT2.OS = new TH1D("hmT2_OS", ";mT2[GeV];Events/5 GeV", 20, 0, 100);
	hMT2.SS = new TH1D("hmT2_SS", ";mT2[GeV];Events/5 GeV", 20, 0, 100);
	v_hists.push_back(hMT2);

	hL12pt.name = "L12pt";
	hL12pt.OS = new TH1D("hL12pt_OS", ";l12 pt [GeV];Events/5 GeV", 40, 0, 200);
	hL12pt.SS = new TH1D("hL12pt_SS", ";l12 pt [GeV];Events/5 GeV", 40, 0, 200);
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
	hMll.OS = new TH1D("hMll_OS", ";m_{ll} [GeV];Events/5 GeV", 40, 0, 200);
	hMll.SS = new TH1D("hMll_SS", ";m_{ll} [GeV];Events/5 GeV", 40, 0, 200);
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
	hl1Pt_JetPt_R.OS = new TH1D("hl1Pt_JetPt_R_OS", ";p_{T}^{l1}/p_T(leading jet);Events/0.25", 20, 0, 5);
	hl1Pt_JetPt_R.SS = new TH1D("hl1Pt_JetPt_R_SS", ";p_{T}^{l1}/p_T(leading jet);Events/0.25", 20, 0, 5);
	v_hists.push_back(hl1Pt_JetPt_R);

	hParents.name = "parents";
	hParents.OS = new TH1D("hParents_OS", "Parent PDGID", 32, -2, 30);
	hParents.SS = new TH1D("hParents_SS", "Parent PDGID", 32, -2, 30);
	v_hists.push_back(hParents);

	hParents1.name = "parents1";
	hParents1.OS = new TH1D("hParents1_OS", "Lep 1 Parent PDGID", 32, -2, 30);
	hParents1.SS = new TH1D("hParents1_SS", "Lep 1 Parent PDGID", 32, -2, 30);
	v_hists.push_back(hParents1);

	hParents2.name = "parents2";
	hParents2.OS = new TH1D("hParents2_OS", "Lep 2 Parent PDGID", 32, -2, 30);
	hParents2.SS = new TH1D("hParents2_SS", "Lep 2 Parent PDGID", 32, -2, 30);
	v_hists.push_back(hParents2);

	return true;
}
bool destructHists(){
	for(auto h : v_hists){
		if(h.OS) {h.OS->Write(); delete h.OS; h.OS = 0;}
		if(h.SS) {h.SS->Write(); delete h.SS; h.SS = 0;}
		if(h.comp) {h.comp->Write(); delete h.comp; h.comp = 0;}
		if(h.OSNorm) {h.OSNorm->Write(); delete h.OSNorm; h.OSNorm = 0;}
		if(h.SSNorm) {h.SSNorm->Write(); delete h.SSNorm; h.SSNorm = 0;}
		if(h.compNorm) {h.compNorm->Write(); delete h.compNorm; h.compNorm = 0;}
		v_hists.clear();
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
         	tc->Add(line.c_str());
         	// if(tc->Add(line.c_str())) std::cout << "Added " << line << std::endl;
         }
         inFiles.close();
      }
   }
   return tc;
}
evt2l* evt2lTree;

bool hasISR()
{
	for(int i=0; i<evt2lTree->jets_; i++)
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
	if (!(evt2lTree->leps_lFlag[0] & 2) || !(evt2lTree->leps_lFlag[1] & 2)) return false;

	// Channels:
	//   noISR:  0=ee,  1=emu,  2=mumu,  3=combFlav
	//     ISR: 10=ee, 11=emu, 12=mumu, 13=combFlav
	// combISR: 20=ee, 21=emu, 22=mumu, 23=combFlav
	// useISR = True if int(ch/10)==1 else False

	if(int(channel/10) == 0 &&  hasISR()) return false;
	if(int(channel/10) == 1 && !hasISR()) return false;

	if(int(channel%10)==0 && 
		!(int(abs(evt2lTree->leps_ID[0])/1000) == 11 && int(abs(evt2lTree->leps_ID[1])/1000) == 11)) return false;
	if(int(channel%10)==1 &&
		(int(abs(evt2lTree->leps_ID[0])/1000) == int(abs(evt2lTree->leps_ID[1])/1000))) return false;
	if(int(channel%10)==2 &&
		!(int(abs(evt2lTree->leps_ID[0])/1000) == 13 && int(abs(evt2lTree->leps_ID[1])/1000) == 13)) return false;

	// Zmass veto for same flavor
	if((int(abs(evt2lTree->leps_ID[0])/1000) == int(abs(evt2lTree->leps_ID[1])/1000)) 
		&& fabs(evt2lTree->l12_m - 91.1876)<=10) return false;

	return true;
}

void draw()
{
	// Hide info box of plot
	gStyle->SetOptStat(0);

	for(hists h : v_hists){
		TCanvas* c = new TCanvas(1);
		c->Divide(1,2);
		c->cd(1)->SetPad(0, 0.25, 1, 1);
		c->cd(2)->SetPad(0, 0, 1, 0.25);

		if(h.OS->GetSumOfWeights()<=0 || h.SS->GetSumOfWeights()<=0){
			c->Print(outDir+"/"+h.name+"_norm.pdf", "Title:"+h.name);
			c->Print(outDir+"/"+h.name+".pdf", "Title:"+h.name);
			continue;
		}

		///// NORMALIZED //////
		c->cd(1);
		h.OS->SetLineColor(kRed);
		h.SS->SetLineColor(kBlue);

		double yMax = max(h.OS->GetMaximum()/h.OS->GetSumOfWeights(), h.SS->GetMaximum()/h.SS->GetSumOfWeights());
		char axisTitles[150];
		sprintf(axisTitles, ";%s;%s", h.OS->GetXaxis()->GetTitle(), h.SS->GetYaxis()->GetTitle());
		TH2F* frame = new TH2F("frame", axisTitles,
			500, h.OS->GetXaxis()->GetXmin(), h.OS->GetXaxis()->GetXmax(),
			500, 0, yMax*1.1);
		frame->SetTitle("Normalized OS and SS "+h.name+" distributions");
		frame->Draw();

		TLegend l(0.75, 0.75, 0.9, 0.9);
		l.AddEntry(h.OS, "OS", "l");
		l.AddEntry(h.SS, "SS", "l");
		
		h.compNorm = (TH1*) h.OS->Clone(h.name+"_norm_ratio");
		h.compNorm->Divide(h.OS->DrawNormalized("same"), h.SS->DrawNormalized("same"));
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
		h.compNorm->GetYaxis()->SetRangeUser(0.5, 2);
		h.compNorm->Draw("p");

		TLine line;
		TLine *l1 = line.DrawLine(h.compNorm->GetXaxis()->GetXmin(), 1., h.compNorm->GetXaxis()->GetXmax(), 1.);
		l1->SetLineStyle(1);
		l1->SetLineWidth(2);

		c->Print(outDir+"/"+h.name+"_norm.pdf", "Title:"+h.name);
		delete frame; frame=0;

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
		h.comp->GetYaxis()->SetRangeUser(0.5, 2);
		h.comp->Draw("p");

		l1 = line.DrawLine(h.compNorm->GetXaxis()->GetXmin(), 1., h.compNorm->GetXaxis()->GetXmax(), 1.);
		l1->SetLineStyle(1);
		l1->SetLineWidth(2);

		c->cd(1);
		h.OS->SetLineColor(kRed);
		h.SS->SetLineColor(kBlue);

		yMax = max(h.OS->GetMaximum(), h.SS->GetMaximum());
		sprintf(axisTitles, ";%s;%s", h.OS->GetXaxis()->GetTitle(), h.SS->GetYaxis()->GetTitle());
		frame = new TH2F("frame", axisTitles,
			500, h.OS->GetXaxis()->GetXmin(), h.OS->GetXaxis()->GetXmax(),
			500, 0, yMax*1.1);
		frame->SetTitle("OS and SS "+h.name+" distributions"); 
		frame->Draw();
		h.OS->Draw("same");
		h.SS->Draw("same");

		l.Draw();

		c->Print(outDir+"/"+h.name+".pdf", "Title:"+h.name);
		delete frame; frame = 0;
		delete c;
		h.OS->SetTitle("");
	}
	currentDir->Write();
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
		if(!selectChannel(channel)) continue;

		double weight = 1;
		// weight *= evt2lTree->evt_weight;
		// weight *= evt2lTree->evt_pwt;
		// weight *= evt2lTree->evt_ElSF;
		// weight *= evt2lTree->evt_MuSF;
		// weight *= evt2lTree->evt_BtagSF;

		int pdgID;
		if((evt2lTree->leps_ID[0]>0) == (evt2lTree->leps_ID[1]>0))
		{
			hMT2.SS->Fill(evt2lTree->sig_mT2, weight);
			hL12pt.SS->Fill(evt2lTree->l12_pt, weight);
			hMET.SS->Fill(evt2lTree->sig_MetRel, weight);
			hHT.SS->Fill(evt2lTree->sig_HT, weight);
			hMTl1.SS->Fill(evt2lTree->leps_mT[0], weight);
			hMTl2.SS->Fill(evt2lTree->leps_mT[1], weight);
			hLLdPhi.SS->Fill(evt2lTree->l12_dPhi, weight);
			hMll.SS->Fill(evt2lTree->l12_m, weight);

			if(int(channel%10)!=0){
				hJetMet_dPhi.SS->Fill(evt2lTree->jets_MET_dPhi[0], weight);
				hMET_JetPt_R.SS->Fill(evt2lTree->sig_MetRel/evt2lTree->jets_pt[0], weight);
				hl1Pt_JetPt_R.SS->Fill(evt2lTree->leps_pt[0]/evt2lTree->jets_pt[0], weight);
			}

			if(evt2lTree->leps_truthI[0]<0) pdgID = -1;
			else if (evt2lTree->truths_motherI[evt2lTree->leps_truthI[0]]<0) pdgID = -2;
			else pdgID = evt2lTree->truths_pdgId[evt2lTree->truths_motherI[evt2lTree->leps_truthI[0]]];
			hParents.SS->Fill(pdgID, weight);
			hParents1.SS->Fill(pdgID, weight);

			if(evt2lTree->leps_truthI[1]<0) pdgID = -1;
			else if (evt2lTree->truths_motherI[evt2lTree->leps_truthI[1]]<0) pdgID = -2;
			else pdgID = evt2lTree->truths_pdgId[evt2lTree->truths_motherI[evt2lTree->leps_truthI[1]]];
			hParents.SS->Fill(pdgID, weight);
			hParents2.SS->Fill(pdgID, weight);
		}
		else
		{
			hMT2.OS->Fill(evt2lTree->sig_mT2, weight);
			hL12pt.OS->Fill(evt2lTree->l12_pt, weight);
			hMET.OS->Fill(evt2lTree->sig_MetRel, weight);
			hHT.OS->Fill(evt2lTree->sig_HT, weight);
			hMTl1.OS->Fill(evt2lTree->leps_mT[0], weight);
			hMTl2.OS->Fill(evt2lTree->leps_mT[1], weight);
			hLLdPhi.OS->Fill(evt2lTree->l12_dPhi, weight);
			hMll.OS->Fill(evt2lTree->l12_m, weight);

			if(int(channel%10)!=0){
				hJetMet_dPhi.OS->Fill(evt2lTree->jets_MET_dPhi[0], weight);
				hMET_JetPt_R.OS->Fill(evt2lTree->sig_MetRel/evt2lTree->jets_pt[0], weight);
				hl1Pt_JetPt_R.OS->Fill(evt2lTree->leps_pt[0]/evt2lTree->jets_pt[0], weight);
			}

			if(evt2lTree->leps_truthI[0]<0) pdgID = -1;
			else if (evt2lTree->truths_motherI[evt2lTree->leps_truthI[0]]<0) pdgID = -2;
			else pdgID = evt2lTree->truths_pdgId[evt2lTree->truths_motherI[evt2lTree->leps_truthI[0]]];
			hParents.OS->Fill(pdgID, weight);
			hParents1.OS->Fill(pdgID, weight);

			if(evt2lTree->leps_truthI[1]<0) pdgID = -1;
			else if (evt2lTree->truths_motherI[evt2lTree->leps_truthI[1]]<0) pdgID = -2;
			else pdgID = evt2lTree->truths_pdgId[evt2lTree->truths_motherI[evt2lTree->leps_truthI[1]]];
			hParents.OS->Fill(pdgID, weight);
			hParents2.OS->Fill(pdgID, weight);
		}
	}
	draw();
	nOS=hMT2.OS->GetSumOfWeights();
	nSS=hMT2.SS->GetSumOfWeights();
	destructHists();
	return 0;
}

int plot(){
	outFile = new TFile("SSOS.root", "recreate");
	ofstream fOut("plots.tex");

	TString dm[] = {"all", "20", "50", "100", "+"};
	int channel[] = {0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23};
	int channel1[] = {0, 1, 2, 3};
	int channel2[] = {3, 13, 23};

	for(TString i : dm) {
		fOut << "\\section{$\\Delta m="<< i.Data() << "$ GeV}\n";
			for(int j : channel){ 
			plot1(i, j);
			char temp[100];
			sprintf(temp, "\\drawPlots{%s}{%d}{%d}{%d}\n", i.Data(), j, (int) nOS, (int) nSS);
			fOut << temp;
		} 
		fOut << endl;
	}
	fOut.close();

	return 0;
}


