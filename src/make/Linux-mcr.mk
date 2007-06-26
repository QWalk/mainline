######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
#CXX:=g++-3.3
CXX:=mpiicpc

CXXFLAGS := -O3 -axW -ip -unroll -fno-alias $(INCLUDEPATH)
CXXFLAGS+= -DUSE_LAPACK -DUSE_MPI -DUSE_BLAS

BLAS_LIBS := -L/usr/local/intel/mkl/lib/32 -lmkl  -lguide -lpthread
BLAS_INCLUDE :=  -I/usr/local/intel/mkl/include

LAPACK_LIBS := $(BLAS_LIBS) -lmkl_ia32 -lmkl_lapack64 -lmkl 
LAPACK_INCLUDE := $(BLAS_INCLUDE)


#CXXFLAGS := -O2 -I$(INCLUDEPATH) -DUSE_MPI

#CXX := mpiCC
#CXXFLAGS:= -O3 -fomit-frame-pointer -funroll-loops -ffast-math -I$(INCLUDEPATH) -DUSE_MPI
#CXXFLAGS:= -O2 --no_exceptions -I$(INCLUDEPATH)
#DEBUG := -Wall  -DRANGE_CHECKING 
DEBUG:=  -DNO_RANGE_CHECKING   -DNDEBUG 
#DEBUG:=-Wall -DNDEBUG -DNO_RANGE_CHECKING -pg
#LDFLAGS := -L/home/apps/mpich-1.2.5/lib -lpmpich++ -lmpich


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
