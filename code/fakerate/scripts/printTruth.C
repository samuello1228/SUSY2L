#include "../support/obj_def.h"
#include "../support/evt2l.C"

const bool DEBUG = false;

void printTruth()
{
	TChain* tc = new TChain("evt2l");
	tc->Add("/afs/cern.ch/user/g/ggallard/private/fakeLep/mkNTuples/run/v11.5.beta3.ZeeP.dRZ/data-myOutput/test.root");
	evt2l* t = new evt2l(tc);

	int nEntries = t->fChain->GetEntries();
	for(int i=0; i<nEntries; i++)
	{
		t->GetEntry(i);
		if(t->leps_!=2 || t->evt_flag==0) continue;
		if(abs(t->leps_ID[0]/1000) !=11) continue;
		if(abs(t->leps_ID[1]/1000) !=11) continue;

		cout << "Event #" << i << "  ###########" << endl;

		cout << "    Lep0 ID: " << setw(10) << t->leps_ID[0] << "  Lep0 truthI: " << t->leps_truthI[0] << '\n'
			<< "    Lep1 ID: " << setw(10) << t->leps_ID[1] << "  Lep1 truthI: " << t->leps_truthI[1] << '\n';

		cout << "  Truth record:  -------- \n"
			<< "    PdgId  : ";

		for(int j=0; j<t->truths_; j++)
			cout << setw(10) << t->truths_pdgId[j] ;
		cout << '\n'
			<< "    MotherI: ";
		for(int j=0; j<t->truths_; j++)
			cout << setw(10) << t->truths_motherI[j] ;
		cout << '\n'
			<< "    Barcode: ";
		for(int j=0; j<t->truths_; j++)
			cout << setw(10) << t->truths_barcode[j];
		cout << "\n\n";
	}
}
