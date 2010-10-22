######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here


CXX:=g++
F77:=gfortran

CXXFLAGS := -O3  \
   -funroll-loops -ffast-math \
  $(INCLUDEPATH) -fomit-frame-pointer


DEBUG:= -Wall -DNO_RANGE_CHECKING -DNDEBUG    -DDEBUG_WRITE
LDFLAGS:= 

HDF_LIBS:=-L/opt/hdf5-1.8.5/lib -lhdf5
HDF_INCLUDE:=-I/opt/hdf5-1.8.5/include

#XML_LIBS:=`xml2-config --libs`
#XML_INCLUDE:=`xml2-config --cflags`
XML_LIBS:=$(shell xml2-config --libs)
XML_INCLUDE:=$(shell xml2-config --cflags)

######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM  $(INCLUDEPATH)

######################################################################
# example of changing options based upon processor type:
# if compiling on 386, don't use 486 extensions

#ifeq ("i386",$(shell uname -m'))
#  CXXFLAGS+=-mno-486
#endif
#ifeq ("i486",$(shell uname -m))
#  CXXFLAGS+=-m486
#endif
#ifeq ("i586",$(shell uname -m))
#  CXXFLAGS+=-m586  # if only gcc supported this...
#endif
