######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
CXX:=mpiCC

MPIDIRECTORY:=/opt/mpich/ethernet/gcc/

CXXFLAGS:= -O3 -fomit-frame-pointer -funroll-loops -ffast-math 
CXXFLAGS += -DUSE_MPI  ${INCLUDEPATH}

#LAPACK_LIBS := -L/usr/lib -lgslcblas
#LAPACK_INCLUDE := ""

#BLAS_LIBS := $(LAPACK_LIBS)
#BLAS_INCLUDE := $(LAPACK_INCLUDE)


#DEBUG := -Wall  -DRANGE_CHECKING -DDEBUG_WRITE
DEBUG:= -Wall -DNO_RANGE_CHECKING   -DNDEBUG 
#DEBUG:=-Wall -DNDEBUG -DNO_RANGE_CHECKING -pg
LDFLAGS := "" 


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
