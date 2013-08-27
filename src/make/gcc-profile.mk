######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here

CXX := g++
CXXFLAGS:= -O2 $(INCLUDEPATH)
#CXXFLAGS+= -DUSE_LAPACK -DUSE_BLAS
#LAPACK_LIBS := -L/home/apps/ATLAS/lib/Linux_ATHLONSSE1 -llapack -lf77blas -lcblas -latlas -static
#LAPACK_INCLUDE := -I/home/apps/ATLAS/include
#BLAS_LIBS := $(LAPACK_LIBS)
#BLAS_INCLUDE := $(LAPACK_INCLUDE)
DEBUG:=-Wall -DNDEBUG -DNO_RANGE_CHECKING -DDEBUG_WRITE  -pg  -g

######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM  $(INCLUDEPATH)

