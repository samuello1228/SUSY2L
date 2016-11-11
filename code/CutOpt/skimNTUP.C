

{

  gSystem.AddDynamicPath( "../FiltreeUtil" );
  gROOT.Macro( "../FiltreeUtil/rootlogon.C" );
  
  //---------------------------Main-------------------------------
  // The input data
  TChain t("evt2l");
  t->Add( "~/sandbox/v5.0.s.361513.MadGraphPythia8EvtGen_A14NNPDF23LO_Ztautau_Np3.root");
  
  //the output files
  TFile outFile("test.root","RECREATE");
  
  // The nodes in the filter chain
  fsSrc = FiltreeStream("fsSrc");
  
  
  //set output verbose level
  fsSrc.printLv = 4;
  
  //connect the src node to the data
  fsSrc.ConnectTo(&t);
  fsSrc.SetOutput("evt2l_skim", &outFile);
  
  //All TTree Draw cmd can be used, esp.see the '$' keyword in TTree::Draw documentation
  fsSrc.AddCut("trigger" , "trigCode!=0");
  fsSrc.AddCut("GRL"     , "evt.cuts==1");
  fsSrc.AddCut("pt1A"     , "leps.pt[0]"); // unit: GeV
  fsSrc.AddCut("pt2A"     , "leps.pt[1]"); // unit: GeV
  //fsSrc.AddCut("2Lepton1", "Length$(leps)");
  fsSrc.AddCut("2Lepton" , "Length$(leps)==2");
  fsSrc.AddCut("mll_60"  , "l12.m>60");
  fsSrc.AddCut("pt1"     , "leps.pt[0]>30"); // unit: GeV
  fsSrc.AddCut("pt2"     , "leps.pt[1]>30");
  
  //fsSrc.AddSpectator("ID0", "int(leps.ID[0]/1000)")
  //fsSrc.AddSpectator("ID1", "int(leps.ID[1]/1000)")
  
  fsSrc.SetWriteBranch("*", false);
  fsSrc.SetWriteBranch("evt", true); //donno why have to keep this else crash
  fsSrc.SetWriteBranch("l12", true);
  fsSrc.SetWriteBranch("sig", true);
  fsSrc.SetWriteBranch("leps.pt", false);
  fsSrc.SetWriteBranch("leps.ID", true);
  
  // start the filter process chain
  fsSrc.StartDataFlow(10);
  //fsSrc.StartDataFlow();
  
  outFile.Write();
  outFile.Close();
}
