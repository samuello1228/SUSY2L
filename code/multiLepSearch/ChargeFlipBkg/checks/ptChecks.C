/* ptChecks.C
 * Gabriel Gallardo
 * 20 July 2016
 */

#include "../common/common.C"
using namespace std;

TString inFileList="../common/inFileList-MCdR.txt";

int ptChecks(){
	TChain *ch = COMMON::loadData(inFileList);
	evt2l evt2lTree(ch);
	ofstream fout("ptPrint1.txt");

	for(int i=0; i<2000; i++){
		ch->GetEntry(i);
		fout << evt2lTree.leps_ << "\t";
		for(int j=0; j<kMaxleps; j++){
			fout << evt2lTree.leps_pt[j] << "\t";
		} 
		fout << endl;
	}

	fout.close();
	return 0;
}