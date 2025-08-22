# PYTHIA configuration file.
# Generated on Fri Dec 5 15:28:20 EST 2014 with the user supplied options:
# --with-lhapdf5-lib=../../lhapdf/lib
# --with-lhapdf5-include=../../lhapdf/include
# --with-fastjet3-lib=../../FastJet/lib
# --with-fastjet3-include=../../FastJet/include

# Install directory prefixes.
PREFIX_BIN=/usr/bin
PREFIX_INCLUDE=/usr/include
PREFIX_LIB=/usr/lib
PREFIX_SHARE=/usr/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=false
CXX=g++
CXX_COMMON=-O2 -ansi -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname
LIB_SUFFIX=.so

# FASTJET3 configuration.
FASTJET3_USE=true
FASTJET3_BIN=
FASTJET3_INCLUDE=../../FastJet/include
FASTJET3_LIB=../../FastJet/lib

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
LHAPDF5_BIN=
LHAPDF5_INCLUDE=../../lhapdf/include
LHAPDF5_LIB=../../lhapdf/lib

# LHAPDF6 configuration.
LHAPDF6_USE=false
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=

# ROOT configuration.
ROOT_USE=false
ROOT_BIN=
ROOT_INCLUDE=
ROOT_LIB=

# BOOST configuration.
BOOST_USE=false
BOOST_BIN=
BOOST_INCLUDE=
BOOST_LIB=

# GZIP configuration.
GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=
