# ROOTCONFIG=`root-config --cflags --glibs`
# OPTFLAGS=-O3 -Wall -Wextra -std=c++11

# data: ConvertNTuples inFileList-data.txt
# 	./ConvertNTuples /afs/cern.ch/work/c/clo/public/ntuple/AnalysisSUSY-02-03-45-405704/ inFileList-data.txt /afs/cern.ch/work/g/ggallard/NTuples/converted/

# mc: ConvertNTuples inFileList-MC.txt
# 	./ConvertNTuples /afs/cern.ch/work/g/ggallard/NTuples/user.ggallard.v6.3.6.dR.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_myOutput.root/ inFileList-MC.txt /afs/cern.ch/work/g/ggallard/NTuples/convertedMC/

# ConvertNTuples: ConvertNTuples.o susyEvts.o obj_dict.o
# 	g++ $(OPTFLAGS)  $^ -o $@ $(ROOTCONFIG)

# ConvertNTuples.o: ConvertNTuples.cxx 
# 	g++ $(OPTFLAGS) -c $< $(ROOTCONFIG)

# susyEvts.o: susyEvts.cxx susyEvts.h
# 	g++ $(OPTFLAGS) -c $< $(ROOTCONFIG)

# obj_dict.o

# clean: 
# 	rm *.o ConvertNTuples

# .PHONY: clean

ROOTFLAG=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)
PKGDIR=/afs/cern.ch/user/g/ggallard/private/multiLepSearch
CXX=gcc
# OPTFLAG=-g -pt
OPTFLAGS=-O3 -std=c++11
CXXFLAGS=$(ROOTFLAG) $(OPTFLAGS) -I$(PKGDIR)
LD=g++ $(CXXFLAGS)

PROG=ConvertNTuples


data: $(PROG) ../common/inFileList-data.txt
	./$(PROG) ../common/inFileList-data.txt /afs/cern.ch/work/g/ggallard/NTuples/convertedData/

zee: $(PROG) ../common/inFileList-Zee.txt
	./$(PROG) ../common/inFileList-Zee.txt /afs/cern.ch/work/g/ggallard/NTuples/convertedZee/

zeeyl: $(PROG) ../common/inFileList-ZeeYL.txt
	./$(PROG) ../common/inFileList-ZeeYL.txt /afs/cern.ch/work/g/ggallard/NTuples/convertedZeeYL/


datayl: $(PROG) ../common/inFileList-DataYL.txt
	./$(PROG) ../common/inFileList-DataYL.txt /afs/cern.ch/work/g/ggallard/NTuples/convertedDataYL/


zeetest: $(PROG) ../common/inFileList-ZeeTest.txt
	./$(PROG) ../common/inFileList-ZeeTest.txt /afs/cern.ch/work/g/ggallard/NTuples/convertedZeeTest/

zeepowheg: $(PROG) ../common/inFileList-ZeePowheg.txt
	./$(PROG) ../common/inFileList-ZeePowheg.txt /afs/cern.ch/work/g/ggallard/NTuples/convertedZeePowheg/

$(PROG): ConvertNTuples.o obj_dict.o susyEvts.o
	$(LD) $(ROOTLIBS) $^ -o $@

ConvertNTuples.o: ConvertNTuples.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

susyEvts.o: susyEvts.cxx susyEvts.h obj_dict.o
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj_dict.o: obj_dict.C
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj_dict.C: obj_def.h Linkdef.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

clean:
	-rm *.o
	-rm *dict*
	-rm $(PROG)

.PHONY: clean


