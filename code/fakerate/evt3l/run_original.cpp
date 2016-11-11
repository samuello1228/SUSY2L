
///////////////////////////////
// Main Program
///////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>  //atoi
#include <unistd.h>  //getopt
#include "analysis.h"
#include "TChain.h"

int main(int argc, char** argv) {

  TString filelist("file.list");

  TChain ch("evt1l1j");

  std::cout << "===========================================" << "\n"
            << "           Reading Files                   " << "\n"
            << "===========================================" << std::endl;

  TString input_roottuple;
  analysis *sel = new analysis;
 
  ifstream fin;
  fin.open(filelist.Data());
  Int_t sum=0;

  while(fin>>input_roottuple){
        
    TFile *newROOTFile= new TFile(input_roottuple);
    TTree *SelEventROOTTree = (TTree *)newROOTFile->Get("evt1l1j");
    
    if(!SelEventROOTTree){
      std::cout << "Can not open ROOT file: " << input_roottuple << std::endl;
    } else {
      std::cout << "File: " << input_roottuple << ", Entries = " << SelEventROOTTree->GetEntries() << std::endl;
      ch.Add(input_roottuple);
      sum++;
    }
    
    delete SelEventROOTTree;
    delete newROOTFile;
  }
  
  std::cout << "===========================================" << "\n"
            << "       " << sum << " file(s) is(are) merged" << "      " << "\n"
            << "===========================================" << std::endl;
  
  fin.close();
  
  std::cout << "===========================================" << "\n"
            << "           Analyzing  Files                   " << "\n"
            << "===========================================" << std::endl; 
  ch.Process(sel);
}




