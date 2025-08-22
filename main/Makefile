# This Makefile is a part of the MARTINI event generator.
# Copyright (C) 2009 Bjoern Schenke.
# MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.


-include Makefile.inc

INCLUDES   = -I. \
	     -I/usr/include/gsl \
	     -I$(lhapdfpath)/include \
	     -I$(pythiapath)/include \
	     -I../src \
             $(shell gsl-config --cflags) 

OBJLIBS = $(LOCAL_LIB)/libsrc.a
LIBS	= -L$(LOCAL_LIB) -lsrc -L$(pythiapath)/lib -lpythia8 -lpythia8lhapdf5 -L$(lhapdfpath)/lib -lLHAPDF -L/usr/local/lib -lgsl -lgslcblas -lm -ldl -Wl,-rpath $(pythiapath) -Wl,-rpath $(lhapdfpath)/lib -Wl,-rpath $(pythiapath)/lib -R$(LOCAL_LIB) -R$(pythiapath)/lib -R$(lhapdfpath)/lib -R/usr/local/lib
#LIBS	= -L$(LOCAL_LIB) -lsrc -L$(pythiapath)/lib -lpythia8 -lpythia8lhapdf5 -L$(lhapdfpath)/lib -lLHAPDF -L/usr/local/lib -lgsl -lgslcblas -lm -ldl -Wl,-rpath $(pythiapath)/lib -Wl,-rpath $(lhapdfpath)/lib -Wl,-rpath $(pythiapath)/lib
#LIBS    = -L$(LOCAL_LIB) -lsrc -L$(pythiapath)/lib -lpythia8 -lpythia8lhapdf5 -L$(lhapdfpath)/lib -lLHAPDF -L/usr/local/lib -lgsl -lgslcblas -lm -ldl -Wl,-rpath $(pythiapath)/lib/lib -Wl,-rpath $(lhapdfpath)/lib -Wl,-rpath $(pythiapath)/lib
#-Wl,-rpath,/homes/jorgen/clinty/junkdir/fastjet-install/lib -lm -L/homes/jorgen/clinty/junkdir/fastjet-install/lib -lfastjet

#ROOT_FLAGS = `root-config --cflags --glibs`

#SOURCES = main_HQYieldsvsTau.cpp main_HQYields.cpp main_v2.cpp main_v2_wo_evolution.cpp

#EXECUTABLE = main_HQYieldsvsTau main_HQYields main_v2 main_v2_wo_evolution

SOURCES = main_v2.cpp main_v2_wo_evolution.cpp main_v2_initial_distributions.cpp

EXECUTABLE = main_v2 main_v2_wo_evolution main_v2_initial_distributions

OBJECTS = $(SOURCES:.cpp=.o)

all: $(EXECUTABLE)

# main_HQYieldsvsTau main_HQYields

$(EXECUTABLE): $(OBJECTS) $(OBJLIBS)
	$(CXX) $(CFLAGS) $(INCLUDES) $@.o $(LIBS) -o $@
#dummy

#main_HQYieldsvsTau.o: main_HQYieldsvsTau.cpp
#	$(CXX) $(CFLAGS) -c $(INCLUDES) -o $@ $<

#main_HQYields.o: main_HQYields.cpp
#	$(CXX) $(CFLAGS) -c $(INCLUDES) -o $@ $<

main_v2.o: main_v2.cpp
	$(CXX) $(CFLAGS) -c $(INCLUDES) -o $@ $<

main_v2_wo_evolution.o: main_v2_wo_evolution.cpp
	$(CXX) $(CFLAGS) -c $(INCLUDES) -o $@ $<

main_v2_initial_distributions.o: main_v2_initial_distributions.cpp
	$(CXX) $(CFLAGS) -c $(INCLUDES) -o $@ $<

# general rule to make a .o out of a .cpp file
.cpp.o:
	$(CXX) $(CFLAGS) -Wall -c $(INCLUDES) -o $@ $< 


.PHONY: clean
clean :
	rm -f $(OBJECTS)
#	rm -fr $(LOCAL_LIB)

force_look :
	true
