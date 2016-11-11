
///////////////////////////////
// Main Program
///////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>  //atoi
#include <unistd.h>  //getopt
#include "analysis.h"
#include <TChain.h>
#include <TROOT.h>
#include <TH1F.h>

using namespace std;

int main(int argc, char** argv) {

  TString filelist("file.list");
//  TChain ch("evt1l1j");

  std::cout << "===========================================" << "\n"
            << "           Reading Files                   " << "\n"
            << "===========================================" << std::endl;

  TString input_file;
//  analysis *sel = new analysis;
 
  ifstream fin;
  fin.open(filelist.Data());

  while(fin>>input_file){
   
    TChain ch("evt3l");
    analysis *sel = new analysis;
    // Open Inputfile
    TString input_roottuple ="input/";
    input_roottuple +=input_file;
    TFile *newROOTFile= new TFile(input_roottuple);
    TTree *SelEventROOTTree = (TTree *)newROOTFile->Get("evt3l");


    if(!SelEventROOTTree){
      std::cout << "Can not open ROOT file: " << input_roottuple << std::endl;
    } else {
      std::cout << "File: " << input_roottuple << ", Entries = " << SelEventROOTTree->GetEntries() << std::endl;
      ch.Add(input_roottuple);

    }

    TString filename = "output/skim_";
    filename += input_file;
    TFile *_file = new TFile(filename, "recreate");
    TH1F *hCutFlow = (TH1F *)((TTree*)newROOTFile->Get("evt3l"))->GetUserInfo()->At(0);
    hCutFlow->SetDirectory(_file);
    hCutFlow->Write();
    
    delete SelEventROOTTree;
    delete newROOTFile; 
    std::cout << "===========================================" << "\n"
              << "           Analyzing  Files                   " << "\n"
              << "===========================================" << std::endl; 
    
    ch.Process(sel);

    //_file->Write();
    _file->Close();
    delete sel;
    delete _file;
    
    //TString oldname = "output/output.root";
    //rename(oldname,newname);
    //TFile *newFile= new TFile(newname,"update");
    //hCutFlow->SetDirectory(newname);
    //hCutFlow->Write();
    //newFile->Close();
    //delete newFile;

  }
  
  fin.close();
}




