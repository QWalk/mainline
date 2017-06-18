######################################################################
# Compiler definitions for Linux systems
# Setup should work for Ubuntu systems with ATLAS and OpenMPI installed
CXX_SERIAL:=g++

CXXFLAGS := -O2  
CXXFLAGS += -DUSE_BLAS -DUSE_LAPACK  ${INCLUDEPATH}

#LAPACK_LIBS := -L/usr/lib/atlas-base -llapack -lf77blas -lcblas -latlas -lf77blas 
LAPACK_LIBS := -L/usr/lib/atlas-base -llapack -lblas
LAPACK_INCLUDE := -I/usr/include/atlas/

BLAS_LIBS :=  $(LAPACK_LIBS) 
BLAS_INCLUDE := $(LAPACK_INCLUDE) 

DEBUG:= -Wall -DNO_RANGE_CHECKING -DNDEBUG   -D__USE_GNU -DDEBUG_WRITE
LDFLAGS:= 

######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM  $(INCLUDEPATH)

