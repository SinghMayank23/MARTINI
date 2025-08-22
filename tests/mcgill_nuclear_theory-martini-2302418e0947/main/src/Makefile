# This Makefile is a part of the MARTINI event generator.
# Copyright (C) 2009 Bjoern Schenke.
# MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.

-include ../../pythia8201/Makefile.inc
include ../Makefile.inc

INCLUDES   = -I. -I../../pythia8201/include -I/homes/zark/chanwook/lhapdf/include
OBJLIBS = ../lib/libsrc.a
CFLAGS	= $(PRJCFLAGS) 
SOURCES = Glauber.cpp HydroSetup.cpp Random.cpp Import.cpp ImportLRates.cpp Rates.cpp Elastic.cpp Information.cpp Setup.cpp Testing.cpp  MARTINI.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(OBJLIBS)

../lib/libsrc.a : $(OBJECTS)
	$(AR) $(ARFLAGS) ../lib/libsrc.a $(OBJECTS)
	$(RANLIB) ../lib/libsrc.a

MARTINI.o: MARTINI.cpp
	$(CXX) -O3 -c -I../../pythia8201/include -I/homes/zark/chanwook/lhapdf/include -I. -o $@ $<

# general rule to make a .o out of a .cpp file
.cpp.o:
	$(CXX) -O3 -Wno-write-strings -c $(INCLUDES) -o $@ $<


.PHONY: clean
clean :
	-$(RM) -f $(OBJECTS)

force_look :
	true
