#!/bin/bash
tag=v13.6

#My Data PRW
#dataPRW=GoodRunsLists/data16_13TeV/20170720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-009.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root

#My MC PRW
#mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Sep19.root,multiLepSearch/prw_MC/merged_prw_mc15c_signal_Sep13.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root

#Peter Data PRW
dataPRW=multiLepSearch/prw_Data/ilumicalc_histograms_None_276262-284484.root,multiLepSearch/prw_Data/ilumicalc_histograms_None_297730-311481_OflLumi-13TeV-009.root

#Peter MC PRW
mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Jun15.root,multiLepSearch/prw_MC/merged_prw_mc15c_signal_Mar1.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root

#SUSYTools Config File
SUSYconf=multiLepSearch/sel_conf/WHSS_0910.conf

#MC type
isMC=1

study=ss
#study=fakes_Peter
#study=fakes

if [ "$study" = "fakes_Peter" ]
then
    dataPRW=multiLepSearch/prw_Data/ilumicalc_histograms_None_297730-311481_OflLumi-13TeV-009.root
fi

truth=MCTC
#truth=dR

#For background MC
#mcPRW=dev/SUSYTools/merged_prw_mc15c_latest.root
#k=${tag}.MCBG
isMC=1

#file=../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt k=${tag}.Zjets #ShowerType = SHERPA
#file=../multiLepSearch/script/MCBGWjetsSherpa_sample_list.txt k=${tag}.Wjets #ShowerType = SHERPA
#file=../multiLepSearch/script/MCBGDYSherpa_sample_list.txt k=${tag}.DY #ShowerType = SHERPA
#file=../multiLepSearch/script/MCBGsingletop_sample_list.txt k=${tag}.singletop #ShowerType = PYTHIAEVTGEN
#file=../multiLepSearch/script/MCBGttV_sample_list.txt k=${tag}.ttV #ShowerType = PYTHIA8EVTGEN
#file=../multiLepSearch/script/MCBGmultitop_sample_list.txt k=${tag}.multitop #ShowerType = PYTHIA8EVTGEN
#file=../multiLepSearch/script/MCBGttbar_sample_list.txt k=${tag}.ttbar #ShowerType = PYTHIAEVTGEN
#file=../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt k=${tag}.Vgamma #ShowerType = SHERPA_CT10
file=../multiLepSearch/script/MCBGVVSherpa_CT10_sample_list.txt k=${tag}.VV_CT10 #ShowerType = SHERPA_CT10
#file=../multiLepSearch/script/MCBGVVSherpa_221_sample_list.txt k=${tag}.VV_221 #ShowerType = SHERPA
#file=../multiLepSearch/script/MCBGVVVSherpa_sample_list.txt k=${tag}.VVV #ShowerType = SHERPA
#file=../multiLepSearch/script/MCBGhiggs_Pythia8EvtGen_sample_list.txt k=${tag}.higgs_Pythia8EvtGen #ShowerType = Pythia8EvtGen
#file=../multiLepSearch/script/MCBGhiggs_HerwigppEvtGen_sample_list.txt k=${tag}.higgs_HerwigppEvtGen #ShowerType = HerwigppEvtGen
#file=../multiLepSearch/script/MCBGRare_sample_list.txt k=${tag}.Rare #ShowerType = PYTHIA8EVTGEN
#file=../multiLepSearch/script/MCBGZeeSherpaSelected.txt k=${tag}.Zee_NOchfSF #ShowerType = SHERPA

#For fast simulation
#isMC=2
#file=../multiLepSearch/script/MCBGmultitop_fast_sample_list.txt k=${tag}.multitop_fast #ShowerType = PYTHIA8EVTGEN

#For signal MC
#isMC=2
#file=../multiLepSearch/script/MCSig_sample_list.txt #ShowerType = PYTHIA8EVTGEN
#k=${tag}.MCSig

../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${k} -o ${k} -a ${isMC} --study ${study} --mcMatch ${truth} --doSys 0 --conf ${SUSYconf} --fast


# file="mc15_13TeV.*.MGPy8EG_A14N23LO_C1N2_Wh_hall_*_2L7.merge.DAOD_SUSY2.e6153_a766_a821_r7676_p2949"
# ../multiLepSearch/util/run_ss_selection.py --driver grid --inputDS ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${k} -o ${k} -a ${isMC} --study ss --mcMatch TruthLink --doSys 0  --conf ${SUSYconf} --fast

#For running on local files

#file=/afs/cern.ch/work/c/clo/sample/mc15_13TeV.364104.Sherpa_221_NNPDF30NNLO_Zmumu_MAXHTPTV70_140_CFilterBVeto.merge.DAOD_SUSY2.e5271_s2726_r7772_r7676_p2949/DAOD_SUSY2.11084041._000011.pool.root.1
#eventNumber=47

#file=/afs/cern.ch/work/c/clo/sample/mc15_13TeV/DAOD_SUSY2.11084160._000250.pool.root.1

#file=/afs/cern.ch/work/c/clo/sample/mc15_13TeV.361071.Sherpa_CT10_lllvjj_EW6.merge.DAOD_SUSY2.e3836_s2726_r7772_r7676_p2949/DAOD_SUSY2.11085088._000011.pool.root.1
#eventNumber=160509

#file=/afs/cern.ch/work/c/clo/sample/mc15_13TeV.410009.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_dil.merge.DAOD_SUSY2.e4511_s2608_s2183_r7725_r7676_p2949/DAOD_SUSY2.10524294._000001.pool.root.1
#eventNumber=5799815

#file=/afs/cern.ch/work/c/clo/sample/mc15_13TeV.363491.Sherpa_221_NNPDF30NNLO_lllv.merge.DAOD_SUSY2.e5332_s2726_r7772_r7676_p2949/DAOD_SUSY2.10506279._000001.pool.root.1
#eventNumber=263

#file=/eos/atlas/user/d/dzhang/Samples/Wh/mc15_13TeV.993820.MGPy8EG_A14N13LO_C1N2_Wh_2L_175_0.merge.DAOD_SUSY2.e5678_a766_a821_r7676_p2949_p2972_pUM999999/merge.DAOD_SUSY2.993820._00001.pool.root.1

#file=/eos/atlas/user/d/dzhang/Samples/Wh/mc15_13TeV.993821.MGPy8EG_A14N13LO_C1N2_Wh_2L_165_35.merge.DAOD_SUSY2.e5678_a766_a821_r7676_p2949_p2972_pUM999999/merge.DAOD_SUSY2.993821._00001.pool.root.1

#file=/eos/atlas/user/d/dzhang/Samples/Wh/mc15_13TeV.993822.MGPy8EG_A14N13LO_C1N2_Wh_2L_400_0.merge.DAOD_SUSY2.e5678_a766_a821_r7676_p2949_p2972_pUM999999/merge.DAOD_SUSY2.993822._00001.pool.root.1
#eventNumber=10566

#file=/afs/cern.ch/work/c/clo/sample/mc15_13TeV.393885.MGPy8EG_A14N23LO_C1N2_Wh_hall_400p0_0p0_2L7.merge.DAOD_SUSY2.e6153_a766_a821_r7676_p2949/DAOD_SUSY2.11792577._000001.pool.root.1
#eventNumber=263
#isMC=2

#file=/afs/cern.ch/user/c/clo/work/sample/mc15_13TeV.364162.Sherpa_221_NNPDF30NNLO_Wmunu_MAXHTPTV140_280_CVetoBVeto.merge.DAOD_SUSY2.e5340_s2726_r7772_r7676_p2949/DAOD_SUSY2.11134137._000011.pool.root.1

#mcPRW=dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root

#../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} -o t1 -w -a ${isMC} --study ${study} --mcMatch ${truth} --doSys 0 --conf ${SUSYconf} # --printEvent ${eventNumber} --cutflow | tee cutflow.log

#for(int i=1;i<60;i++) cout<<hCutFlow->GetXaxis()->GetBinLabel(i)<<": "<<hCutFlow->GetBinContent(i)<<endl;
#grep -e '^[0-9]\+$' cutflow.log |tee trigger.txt
