######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here

CXX:=mpiCC
CXXFLAGS := -O3 \
   -funroll-loops -ffast-math \
  $(INCLUDEPATH) -fomit-frame-pointer

CXXFLAGS+= -DUSE_LAPACK -DUSE_BLAS -DUSE_MPI

BLAS_LIBS := -L/usr/local/intel/mkl/lib/64 -lmkl
BLAS_INCLUDE :=  -I/usr/local/intel/mkl/include

LAPACK_LIBS := $(BLAS_LIBS) -lmkl_lapack  -lmkl_lapack64
LAPACK_INCLUDE := $(BLAS_INCLUDE)


#CXXFLAGS:= -O2 -I$(INCLUDEPATH)
#DEBUG := -Wall  -DRANGE_CHECKING -DDEBUG_WRITE 
#DEBUG:= -Wall -DNO_RANGE_CHECKING   -DDEBUG_WRITE
DEBUG:= -Wall -DNO_RANGE_CHECKING -DNDEBUG   -D__USE_GNU
#DEBUG:=-Wall -DNDEBUG -DNO_RANGE_CHECKING -DDEBUG_WRITE -D__USE_GNU -pg 
#LDFLAGS := -static

######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM -I $(INCLUDEPATH)

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
