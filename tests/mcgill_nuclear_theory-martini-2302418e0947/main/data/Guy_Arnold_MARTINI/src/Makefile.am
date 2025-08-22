SUBDIRS = tools main
bin_PROGRAMS = eloss

eloss_MAIN = elossMain.cxx
eloss_LDFLAGS = $(LDFLAGS)
eloss_SOURCES = $(eloss_MAIN)
eloss_LDADD = main/libmain.a tools/libtools.a 

# note that in _LDADD the order of the libraries is important. Order is such that the lowest level (tools) modules are last in the list.