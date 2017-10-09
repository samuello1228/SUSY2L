ROOTFLAG=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)
# PKGDIR=/home/dzhang/work/bsmSearch/ewSUSY/analysis/r20_evtSelect/multiLepSearch
PKGDIR=/net/s3_data_home/dzhang/work/bsmSearch/ewSUSY/analysis/v20_7b/SUSY2L/code/multiLepSearch
CXX=gcc
# OPTFLAG=-g -pt
OPTFLAGS=-O3
CXXFLAGS=$(ROOTFLAG) $(OPTFLAGS) -I$(PKGDIR)
LD=g++ $(CXXFLAGS)


PROG=run_over_ntuple.exe



$(PROG): test_run_over_ntuple.o obj_dict.o susyEvts.o
	$(LD) $(ROOTLIBS) $^ -o $@
# test_run_over_ntuple.o: test_run_over_ntuple.C
# 	$(CXX) $(CXXFLAGS) -c $< -o $@
susyEvts.o: $(PKGDIR)/Root/susyEvts.cxx $(PKGDIR)/multiLepSearch/susyEvts.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
# obj_dict.o: obj_dict.C
# 	$(CXX) $(CXXFLAGS) -c $< -o $@

run_%.exe: %.o obj_dict.o susyEvts.o
	$(LD) $(ROOTLIBS) $^ -o $@
%.o: %.C
	$(CXX) $(CXXFLAGS) -c $< -o $@



obj_dict.C: $(PKGDIR)/multiLepSearch/obj_def.h Linkdef.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^
clean:
	-rm *.o
	-rm *dict*
	-rm $(PROG)
