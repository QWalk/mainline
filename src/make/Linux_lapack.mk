######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here

CXX:=mpicxx
CXX_SERIAL:=mpicxx
include $(HOME)/mylibs/make.inc
CXX:=mpicxx -I/usr/include/x86_64-linux-gnu/c++/4.8/ -I$(HOME)/mylibs/include

CXXFLAGS := -O3  \
   -funroll-loops -ffast-math \
   -fomit-frame-pointer
CXXFLAGS += -DUSE_MPI -DUSE_BLAS -DUSE_LAPACK  ${INCLUDEPATH}

LAPACK_LIBS := -L/usr/lib/atlas-base -llapack_atlas -lf77blas -lcblas -latlas
LAPACK_LIBS += -L/usr/lib/atlas-base/atlas/ -llapack -lblas
LAPACK_INCLUDE := -I/usr/include/atlas
BLAS_LIBS :=  $(LAPACK_LIBS) 
#BLAS_LIBS +=  -L${HOME}/mylibs -cblas
BLAS_INCLUDE := $(LAPACK_INCLUDE) 

DEBUG:= -Wall -DNO_RANGE_CHECKING -DNDEBUG   -D__USE_GNU -DDEBUG_WRITE
#LDFLAGS:= -static

######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM  $(INCLUDEPATH)

