# This Makefile is a part of the MARTINI event generator.
# Copyright (C) 2009 Bjoern Schenke.
# MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.


-include Makefile.inc

INCLUDES   = -I. \
             -I$(pythiapath)/include \
             $(shell gsl-config --cflags) \
             -I$(lhapdfpath)/include 

OBJLIBS = $(LOCAL_LIB)/libsrc.a

#ROOT_FLAGS = `root-config --cflags --glibs`

SOURCES = Glauber.cpp HydroSetup.cpp Random.cpp Import.cpp ImportLRates.cpp \
          Rates.cpp Elastic.cpp Information.cpp Setup.cpp Testing.cpp \
          MARTINI.cpp binarycollision_info.cpp JetPhotonConversion.cpp

OBJECTS = $(SOURCES:.cpp=.o)

all: $(OBJLIBS)

$(LOCAL_LIB)/libsrc.a : $(OBJECTS)
	rm -fr $(LOCAL_LIB)
	mkdir $(LOCAL_LIB)
	ar rc $(LOCAL_LIB)/libsrc.a $(OBJECTS)
	ranlib $(LOCAL_LIB)/libsrc.a

#MARTINI.o: MARTINI.cpp
#	$(CXX) $(CFLAGS) $(ROOT_FLAGS) -c $(INCLUDES) -o $@ $<

MARTINI.o: MARTINI.cpp
	$(CXX) $(CFLAGS) -c $(INCLUDES) -o $@ $<

# general rule to make a .o out of a .cpp file
.cpp.o:
	$(CXX) $(CFLAGS) -Wall -c $(INCLUDES) -o $@ $< 


.PHONY: clean
clean :
	rm -f $(OBJECTS)
	rm -fr $(LOCAL_LIB)

force_look :
	true
