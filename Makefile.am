### Makefile.am
### Process this file with automake to produce a Makefile
###
### Encore
### Nick Davis (c) 2011-2012
### nick-davis@utulsa.edu
### http://insilico.utulsa.edu

ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = foreign
AM_CXXFLAGS = @OS_CXX@ $(BOOST_CPPFLAGS)

# Random Jungle flags
AM_CXXFLAGS += -DHAVE__BOOL  -D__NOPLUGIN__

lib_LTLIBRARIES=libencore.la
libencore_la_SOURCES= regain.cpp snprank.cpp

libencore_la_HEADERS= regain.h snprank.h

libencore_la_LDFLAGS  = -release 1.0.0

if WITH_OPENMP
AM_CXXFLAGS += -fopenmp
endif

if FORCE_32BIT
AM_CXXFLAGS += -m32
endif

if FORCE_DYNAMIC

else
AM_CXXFLAGS += -static
endif

if SYS_MAC
AM_CXXFLAGS += -D_FILE_OFFSET_BITS=64 -Dfopen64=fopen
endif

if SYS_WIN
AM_CXXFLAGS += -static -static-libgcc -static-libstdc++
endif

# dir to install headers for encore library
libencore_ladir=$(includedir)/encore

# build encore binary, linking to library built above
bin_PROGRAMS = encore
encore_SOURCES =  encore.cpp
encore_LDADD= -lencore $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS)
