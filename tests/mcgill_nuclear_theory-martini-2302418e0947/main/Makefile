-include ../pythia8201/Makefile.inc
include Makefile.inc

INCLUDES   = -I./src -I../pythia8201/include -I/homes/zark/chanwook/lhapdf/include -I/homes/zark/chanwook/FastJet/include
PYTHIALIBDIR = pythia8201/lib/
OBJECTS	= main_ALICE_pi-_Optical.o main_CMS_dNdA.o main_singlejets.o #main_ALICE_charged.o main.o main2.o main3.o main4.o main5.o main6.o main8.o main9.o main_ALICE_charged_alpha25.o dijets_pp_celljet.o  main_Maxime_photons.o main_Maxime_photons_FIC.o main_jet_profile.o main_HQ_testing.o main_HQ_JPsiYields.o main_HQ_JPsiYieldsvsTau.o 
DIRS = src
OBJLIBS	= ./lib/libsrc.a
LIBS	= -lgsl -lgslcblas -L./lib -L/usr/local/lib -L/homes/zark/chanwook/lhapdf/lib -lLHAPDF -lm -L/homes/zark/chanwook/FastJet/lib -lfastjet
# includes gnu scientific library, pythia libs and local libs, as well as lLHAPDF
EXECUTABLE = main_ALICE_pi-_Optical main_CMS_dNdA main_singlejets #main_ALICE_charged main main2 main3 main4 main5 main6 main8 main9 main_ALICE_charged_alpha25 dijets_pp_celljet main_Maxime_photons main_Maxime_photons_FIC main_jet_profile main_HQ_testing main_HQ_JPsiYields main_HQ_JPsiYieldsvsTau
CXXFLAGS = -O3

all: main_ALICE_pi-_Optical main_CMS_dNdA main_singlejets #main_ALICE_charged main main2 main3 main4 main5 main6 main8 main9 main_ALICE_charged_alpha25 dijets_pp_celljet main_Maxime_photons main_Maxime_photons_FIC main_jet_profile main_HQ_testing main_HQ_JPsiYields main_HQ_JPsiYieldsvsTau
$(EXECUTABLE): $(OBJECTS) $(OBJLIBS) 
	$(CXX) $(CXXFLAGS) $@.o $(LIBS) -o $@ -lsrc -L../$(PYTHIALIBDIR) -lpythia8 -lpythia8lhapdf5 -ldl


main_ALICE_pi-_Optical.o: $(OBJLIBS) main_ALICE_pi-_Optical.cpp
	 $(CXX) -c $(INCLUDES) main_ALICE_pi-_Optical.cpp 

main_CMS_dNdA.o: $(OBJLIBS) main_CMS_dNdA.cpp
	 $(CXX) -c $(INCLUDES) main_CMS_dNdA.cpp 

main_singlejets.o: $(OBJLIBS) main_singlejets.cpp
	 $(CXX) -c $(INCLUDES) main_singlejets.cpp 

./lib/libsrc.a : force_look
#	cd src; $(MAKE) $(MFLAGS)

.PHONY: clean # explicitly declare 'clean' to be PHONY
clean :
	-$(RM) -f $(EXECUTABLE) $(OBJECTS) #$(OBJLIBS)
#	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true
