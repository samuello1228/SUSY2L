struct GroupData
{
    TString GroupName;
    TString LegendName;
    TString LatexName;
    unsigned int lower;
    unsigned int upper;
    int colour;
    int unweighted;
    double weighted;
    double error;
};

std::vector<GroupData> getBKG(){
        std::vector<GroupData> BGMCGroupData;
        GroupData element;

        element.GroupName = "Zjets"; element.LegendName = "Z+jets"; element.LatexName = "Z+jets";
        element.lower = 0;  element.upper = 59; element.colour = 2; BGMCGroupData.push_back(element);
        
        //W+jets
        element.GroupName = "Wjets"; element.LegendName = "W+jets"; element.LatexName = "W+jets";
        element.lower = 60;  element.upper = 101; element.colour = 3; BGMCGroupData.push_back(element);
        
        //Top
        element.GroupName = "top"; element.LegendName = "top"; element.LatexName = "top";
        element.lower = 102;  element.upper = 116; element.colour = 4; BGMCGroupData.push_back(element);
        
        element.GroupName = "ttbar"; element.LegendName = "t#bar{t}"; element.LatexName = "$t\\bar{t}$";
        element.lower = 102;  element.upper = 102; element.colour = 881; BGMCGroupData.push_back(element);

        element.GroupName = "singletop"; element.LegendName = "single top"; element.LatexName = "single top";
        element.lower = 103;  element.upper = 108; element.colour = 30; BGMCGroupData.push_back(element);
        //element.lower = 103;  element.upper = 106; element.colour = 30; BGMCGroupData.push_back(element);
        
        element.GroupName = "ttV"; element.LegendName = "t#bar{t}+V"; element.LatexName = "$t\\bar{t}+V$";
        element.lower = 109;  element.upper = 114; element.colour = 600; BGMCGroupData.push_back(element);
        
        element.GroupName = "multitop"; element.LegendName = "multi top"; element.LatexName = "multi top";
        element.lower = 115;  element.upper = 116; element.colour = 635; BGMCGroupData.push_back(element);
        
        element.GroupName = "VV"; element.LegendName = "VV"; element.LatexName = "VV";
        element.lower = 117;  element.upper = 130; element.colour = 801; BGMCGroupData.push_back(element);
        
        element.GroupName = "Vgamma"; element.LegendName = "V + #gamma"; element.LatexName = "V$+\\gamma$";
        element.lower = 131;  element.upper = 150; element.colour = 5; BGMCGroupData.push_back(element);

        element.GroupName = "VVV"; element.LegendName = "VVV"; element.LatexName = "VVV";
        element.lower = 151;  element.upper = 155; element.colour = 607; BGMCGroupData.push_back(element);

        element.GroupName = "Higgs"; element.LegendName = "Higgs"; element.LatexName = "Higgs";
        element.lower = 156;  element.upper = 168; element.colour = 7; BGMCGroupData.push_back(element);

     return BGMCGroupData;
    }
    
