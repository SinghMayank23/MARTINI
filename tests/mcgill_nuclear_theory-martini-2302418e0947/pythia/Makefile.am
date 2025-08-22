SUBDIRS = 
bin_PROGRAMS = main

main_MAIN = main.cpp
main_LDFLAGS = $(LDFLAGS)
main_SOURCES = $(main_MAIN) Random.h Random.cpp Import.h Import.cpp
main_LDADD =  
# Location of directories
main_TMPDIR=tmp
main_TOPDIR=$(shell \pwd)
main_INCDIR=include
main_SRCDIR=src
main_LIBDIR=lib
main_LIBDIRARCH=lib/archive
main_BINDIR=bin


# note that in _LDADD the order of the libraries is important. Order is such that the lowest level (tools) modules are last in the list.