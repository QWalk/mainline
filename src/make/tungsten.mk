######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
CXX:=cmpic++

CXXFLAGS:= -gnu -O3 -fomit-frame-pointer -funroll-loops -ffast-math 
#CXXFLAGS:= -O3 -axW -unroll -fno-alias
CXXFLAGS += -DUSE_MPI -DUSE_BLAS ${INCLUDEPATH}

#LAPACK_LIBS := -L/usr/lib -lgslcblas
#LAPACK_INCLUDE := ""

BLAS_LIBS := -L${MKL_HOME}/lib/32 -lmkl -lguide -lpthread
#BLAS_INCLUDE := ""


DEBUG:= -DNO_RANGE_CHECKING   -DNDEBUG 
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
