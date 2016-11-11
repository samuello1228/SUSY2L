import ROOT
import math
from ROOT import TCanvas, THStack, TH1F, gRandom, TPad, TLegend, TLatex
import sys, getopt
import os

from ReadMCInfo import MC_info
from ReadOutputList import output_info
from ReadPlotConfig import plot_config_info 
from ReadChannels import channel_info
from ReadDataInfo import data_info
from ReadVarInfo import var_info
#tree_name = "skimming_ee_OC"
#branch_name = "mll"
newfile = ROOT.TFile("graph/graph.root","update")
#calculate data luminosity
Ldata = 0.
for data_name in data_info:
    Ldata = Ldata + data_info[data_name]["Ldata"]
#test_info = ["Mt"]
for tree_name in channel_info:
  for branch_name in var_info:
### open root file
    c=ROOT.TCanvas(tree_name+"_"+branch_name,tree_name+"_"+branch_name)
    if(var_info[branch_name]["log"]):
        c.SetLogy()
    hist_dict = {}

    #combine same processes and scaled
    for output_file in output_info:
      file = ROOT.TFile(output_file[:-1])
      tree = file.Get(tree_name)
      hist = ROOT.TH1F("hist",tree_name+"_"+var_info[branch_name]["name"],var_info[branch_name]["nbin"],var_info[branch_name]["min"],var_info[branch_name]["max"])
      MCID = output_info[output_file]["MCID"]

      color = plot_config_info[MCID]["color"]
      name = plot_config_info[MCID]["name"]
      stack = plot_config_info[MCID]["stack"]
      process = plot_config_info[MCID]["process"]
      counter = 0.
      #calculate the scale
      if(process == "data"):
          scale = 1.
          tree.Draw(branch_name+">>hist")
      else:
          xsec = MC_info[MCID]["xsec"]
          eff = MC_info[MCID]["eff"]
          kfac = MC_info[MCID]["kfac"]
          hcutflow = file.Get("hCutFlow")
          Nmc = hcutflow.GetBinContent(2)
          scale = Ldata*xsec*eff*kfac/Nmc
          tree.Draw(branch_name+">>hist","lep_SF")
          print "MCID Scale Xsec Eff Kfac Nmc Ldata"
          print MCID,scale,xsec,eff,kfac,Nmc,Ldata


      hist.SetDirectory(0)
      hist.Scale(scale)
      total = hist.Integral(0,-1)
      print tree_name, process, total
      
      
      temp = {"name":name,"color":color,"stack":stack,"hist":hist}
      newprocess = True
      for process_prev in hist_dict:
        if(process_prev == process):
          newprocess = False
          hist_dict[process_prev]["hist"].Add(hist)
      if(newprocess):
        hist_dict[process]=temp

    #plot graph & legend
    stack_graph = ROOT.THStack("stack_graph",tree_name+"_"+var_info[branch_name]["name"])
    xl1=0.75
    yl1=0.75
    xl2=xl1+0.2
    yl2=yl1+0.15
    leg = ROOT.TLegend(xl1,yl1,xl2,yl2)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)

    max = 0

    #data plot
    for process in hist_dict:
        if(process == "data"):
            data=hist_dict["data"]["hist"]
            color = hist_dict["data"]["color"]
            data.SetLineColor(color)
            data.SetMarkerStyle(8)
            leg.AddEntry(data,"data","p")
            max=data.GetMaximum()

    #background plot
    for process in hist_dict:
      name = hist_dict[process]["name"]
      stack = hist_dict[process]["stack"]
      color = hist_dict[process]["color"]
      hist = hist_dict[process]["hist"]
      if(stack == "0"):
        hist.SetLineColor(color)
        hist.SetFillColor(color)
        stack_graph.Add(hist)
        leg.AddEntry(hist,name,"f")



    #find max and min
    if(max<stack_graph.GetMaximum()):
        max = stack_graph.GetMaximum()
    max=max*1.1
    ROOT.gStyle.SetOptStat(0)
    stack_graph.SetMaximum(max)
    stack_graph.Draw("hist")
    data.Draw("pesame")
    leg.Draw()
    unit=""
    if(var_info[branch_name]["unit"]!="/"):
        unit=" ["+var_info[branch_name]["unit"]+"]"
    stack_graph.GetXaxis().SetTitle(var_info[branch_name]["name"]+unit)
    stack_graph.GetYaxis().SetTitle("number of events")
   # title = ROOT.TLatex()
   # title.DrawLatexNDC(0.4,0.8,tree_name[9:]+"_"+branch_name)


    c.SaveAs("graph/"+tree_name+"_"+branch_name+".png")

    #substract background from data
    d=ROOT.TCanvas(tree_name+"_"+branch_name+"_sub",tree_name+"_"+branch_name+"_sub")
    for process in hist_dict:
        hist = hist_dict[process]["hist"]
        stack = hist_dict[process]["stack"]
        if(stack == "0"):
            data.Add(hist,-1)
    data.Draw("pe")
    data.SetNameTitle(tree_name+"_"+branch_name+"_sub", tree_name+"_"+branch_name+"_sub")
    newfile.cd()
    data.Write()

    #d.SaveAs("graph/"+tree_name+"_"+branch_name+"_sub.png")

    #plot MC histogram
    #MC_hist = ROOT.TH1F(tree_name+"_"+var_info[branch_name]["name"],tree_name+"_"+var_info[branch_name]["name"],var_info[branch_name]["nbin"],var_info[branch_name]["min"],var_info[branch_name]["max"])
    #for process in hist_dict:
    #    hist = hist_dict[process]["hist"]
    #    stack = hist_dict[process]["stack"]
    #    if(stack == "0"):
    #        MC_hist.Add(hist,1)
    #MC_hist.Draw("pe")
    #newfile.cd()
    #MC_hist.Write()



#newfile.Write("All")
newfile.Close()

