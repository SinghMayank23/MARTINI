SHELL = /bin/sh

-include ../pythia/config.mk
include Makefile.inc

INCLUDES   = -I./src -I../pythia/include -I../lhapdf/include -I../../fastjet-install/include
PYTHIALIBDIR = pythia/lib/archive
OBJECTS	= main.o main2.o main3.o main4.o main5.o main6.o main_ALICE_charged.o main_ALICE_charged_alpha25.o main8.o main9.o dijets_pp_celljet.o main_CMS_dNdA.o main_Maxime_photons.o main_Maxime_photons_FIC.o main_jet_profile.o main_JPsi_v2.o main_HQ_JPsiYields.o main_HQ_JPsiYieldsvsTau.o
DIRS = src
OBJLIBS	= ./lib/libsrc.a
LIBS	= -lgsl -lgslcblas -L./lib -L/usr/local/lib -L../lhapdf/lib -lLHAPDF -Wl,-rpath,/homes/jorgen/clinty/junkdir/fastjet-install/lib -lm -L/homes/jorgen/clinty/junkdir/fastjet-install/lib -lfastjet
# includes gnu scientific library, pythia libs and local libs, as well as lLHAPDF
EXECUTABLE = main main2 main3 main4 main5 main6 main_ALICE_charged main_ALICE_charged_alpha25 main8 main9 dijets_pp_celljet main_CMS_dNdA main_Maxime_photons main_Maxime_photons_FIC main_jet_profile main_JPsi_v2 main_HQ_JPsiYields main_HQ_JPsiYieldsvsTau
CXXFLAGS = -O3

all: main main2 main3 main4 main5 main6 main_ALICE_charged main_ALICE_charged_alpha25 main8 main9 dijets_pp_celljet main_CMS_dNdA main_Maxime_photons main_Maxime_photons_FIC main_jet_profile main_JPsi_v2 main_HQ_JPsiYields main_HQ_JPsiYieldsvsTau

$(EXECUTABLE): $(OBJECTS) $(OBJLIBS) 
	$(CXX) $(CXXFLAGS) $@.o $(LIBS) -o $@ -lsrc -L../$(PYTHIALIBDIR) -lpythia8 -llhapdfdummy

main.o: $(OBJLIBS) main.cpp
	$(CXX) -c $(INCLUDES) main.cpp 

main2.o: $(OBJLIBS) main2.cpp
	 $(CXX) -c $(INCLUDES) main2.cpp 

main3.o: $(OBJLIBS) main3.cpp
	 $(CXX) -c $(INCLUDES) main3.cpp 

main4.o: $(OBJLIBS) main4.cpp
	 $(CXX) -c $(INCLUDES) main4.cpp 

main5.o: $(OBJLIBS) main5.cpp
	 $(CXX) -c $(INCLUDES) main5.cpp 

main6.o: $(OBJLIBS) main6.cpp
	 $(CXX) -c $(INCLUDES) main6.cpp 

main_ALICE_charged.o: $(OBJLIBS) main_ALICE_charged.cpp
	 $(CXX) -c $(INCLUDES) main_ALICE_charged.cpp 

main_Maxime_photons.o: $(OBJLIBS) main_Maxime_photons.cpp
	 $(CXX) -c $(INCLUDES) main_Maxime_photons.cpp 

main_Maxime_photons_FIC.o: $(OBJLIBS) main_Maxime_photons_FIC.cpp
	 $(CXX) -c $(INCLUDES) main_Maxime_photons_FIC.cpp 

main_jet_profile.o: $(OBJLIBS) main_jet_profile.cpp
	 $(CXX) -c $(INCLUDES) main_jet_profile.cpp 

main_ALICE_charged_alpha25.o: $(OBJLIBS) main_ALICE_charged_alpha25.cpp
	 $(CXX) -c $(INCLUDES) main_ALICE_charged_alpha25.cpp 

main8.o: $(OBJLIBS) main8.cpp
	 $(CXX) -c $(INCLUDES) main8.cpp 

main9.o: $(OBJLIBS) main9.cpp
	 $(CXX) -c $(INCLUDES) main9.cpp 

dijets_pp_celljet.o: $(OBJLIBS) dijets_pp_celljet.cpp
	 $(CXX) -c $(INCLUDES) dijets_pp_celljet.cpp 

main_CMS_dNdA.o: $(OBJLIBS) main_CMS_dNdA.cpp
	 $(CXX) -c $(INCLUDES) main_CMS_dNdA.cpp 

main_JPsi_v2.o: $(OBJLIBS) main_JPsi_v2.cpp
	 $(CXX) -c $(INCLUDES) main_JPsi_v2.cpp 

main_HQ_JPsiYields.o: $(OBJLIBS) main_HQ_JPsiYields.cpp
	 $(CXX) -c $(INCLUDES) main_HQ_JPsiYields.cpp 

main_HQ_JPsiYieldsvsTau.o: $(OBJLIBS) main_HQ_JPsiYieldsvsTau.cpp
	 $(CXX) -c $(INCLUDES) main_HQ_JPsiYieldsvsTau.cpp 

./lib/libsrc.a : force_look
	cd src; $(MAKE) $(MFLAGS)

.PHONY: clean # explicitly declare 'clean' to be PHONY
clean :
	-$(RM) -f $(EXECUTABLE) $(OBJECTS) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true
