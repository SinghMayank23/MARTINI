# This Makefile is a part of the MARTINI event generator.
# Copyright (C) 2009 Bjoern Schenke.
# MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.

-include ../Makefile.inc
fastjetpath = /home/chanwook/fastjet3

INCLUDES   = -I. \
			 -I../$(pythiapath)/include \
			 $(shell gsl-config --cflags) \
             -I$(fastjetpath)/include

OBJLIBS = ../$(LOCAL_LIB)/libbackground.a

SOURCES = Background.cpp 

OBJECTS = $(SOURCES:.cpp=.o)

all: $(OBJLIBS)

../$(LOCAL_LIB)/libbackground.a : $(OBJECTS)
	ar rc ../$(LOCAL_LIB)/libbackground.a $(OBJECTS)
	ranlib ../$(LOCAL_LIB)/libbackground.a

Bacground.o: Background.cpp
	$(CXX) $(CFLAGS) -c $(INCLUDES) -o $@ $<

# general rule to make a .o out of a .cpp file
.cpp.o:
	$(CXX) $(CFLAGS) -Wall -c $(INCLUDES) -o $@ $<


.PHONY: clean
clean :
	rm -f $(OBJECTS)
	rm -f ../$(LOCAL_LIB)/libbackground.a

force_look :
	true
