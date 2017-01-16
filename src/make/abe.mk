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

