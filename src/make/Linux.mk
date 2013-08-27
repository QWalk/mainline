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

XML_LIBS:=$(shell xml2-config --libs)
XML_INCLUDE:=$(shell xml2-config --cflags)

######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM  $(INCLUDEPATH)

