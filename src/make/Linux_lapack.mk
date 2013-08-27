######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here

CXX:=g++
CXX_SERIAL:=g++

CXXFLAGS := -O3  \
   -funroll-loops -ffast-math \
   -fomit-frame-pointer
CXXFLAGS += -DUSE_BLAS -DUSE_LAPACK  ${INCLUDEPATH}

LAPACK_LIBS := -L/usr/lib/atlas -llapack -lf77blas -lcblas -latlas -lblas -lg2c -lm -lgfortran
LAPACK_INCLUDE := -I/usr/include/atlas
BLAS_LIBS :=  $(LAPACK_LIBS) 
BLAS_INCLUDE := $(LAPACK_INCLUDE) 

DEBUG:= -Wall -DNO_RANGE_CHECKING -DNDEBUG   -D__USE_GNU -DDEBUG_WRITE
LDFLAGS:= -static

######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM  $(INCLUDEPATH)

