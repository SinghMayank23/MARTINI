# PYTHIA configuration file.
# Generated on Fri May 27 15:03:39 EDT 2022 with the user supplied options:
# --prefix=/global/homes/s/singhm/MARTINI/pythia8209/lib
# --arch=Linux
# --with-lhapdf5=/global/homes/s/singhm/MARTINI/lhapdf-5.9.1/lib/

# Install directory prefixes.
PREFIX_BIN=/global/homes/s/singhm/MARTINI/pythia8209/lib/bin
PREFIX_INCLUDE=/global/homes/s/singhm/MARTINI/pythia8209/lib/include
PREFIX_LIB=/global/homes/s/singhm/MARTINI/pythia8209/lib/lib
PREFIX_SHARE=/global/homes/s/singhm/MARTINI/pythia8209/lib/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=false
CXX=g++
CXX_COMMON=-O2 -ansi -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname
LIB_SUFFIX=.so

# EVTGEN configuration.
EVTGEN_USE=false
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

# FASTJET3 configuration.
FASTJET3_USE=false
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

# HEPMC2 configuration.
HEPMC2_USE=false
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

# HEPMC3 configuration.
HEPMC3_USE=false
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

# LHAPDF5 configuration.
LHAPDF5_USE=true
LHAPDF5_BIN=/global/homes/s/singhm/MARTINI/lhapdf-5.9.1/lib//bin/
LHAPDF5_INCLUDE=/global/homes/s/singhm/MARTINI/lhapdf-5.9.1/lib//include
LHAPDF5_LIB=/global/homes/s/singhm/MARTINI/lhapdf-5.9.1/lib//lib
LHAPDF5_PLUGIN=LHAPDF5.h

# LHAPDF6 configuration.
LHAPDF6_USE=false
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=
LHAPDF6_PLUGIN=LHAPDF5.h

# POWHEG configuration.
POWHEG_USE=false
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

# PROMC configuration.
PROMC_USE=false
PROMC_BIN=
PROMC_INCLUDE=
PROMC_LIB=

# ROOT configuration.
ROOT_USE=false
ROOT_BIN=
ROOT_INCLUDE=
ROOT_LIB=

# GZIP configuration.
GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

# BOOST configuration.
BOOST_USE=false
BOOST_BIN=
BOOST_INCLUDE=
BOOST_LIB=
