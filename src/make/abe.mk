######################################################################
# Compiler definitions for abe.ncsa.uiuc.edu
# Intel 64 Linux cluster, Intel compiler, MVAPICH2
######################################################################

CXX:=mpicxx
CXXFLAGS:= -O2 -ip -mp -xT -DMPICH_IGNORE_CXX_SEEK 
CXXFLAGS += -DUSE_MPI -DUSE_BLAS ${INCLUDEPATH}

#To set up the MKL environment call:  soft add +intel-mkl
#LAPACK_LIBS := -L${MKL_HOME}/lib/em64t -lmkl_lapack -lmkl -lguide -lpthread
#LAPACK_INCLUDE := ""
BLAS_LIBS := -L${MKL_HOME}/lib/em64t -lmkl -lguide -lpthread
#BLAS_INCLUDE := ""

DEBUG:= -DNO_RANGE_CHECKING -DNDEBUG 
#LDFLAGS := "" 

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
