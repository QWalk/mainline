 ######################################################################

CXX:=g++

CXXFLAGS := -O3   \
   -funroll-loops -ffast-math  \
  $(INCLUDEPATH) # -fomit-frame-pointer
 #not omitting the frame pointer so we get instrumentation

CXXFLAGS+= -DUSE_LAPACK -DUSE_BLAS 
# uncomment here if you've installed ATLAS through Fink, otherwise use the accelerate framework
#LAPACK_LIBS := -L/sw/lib/ -llapack -lf77blas -lcblas -latlas
#LAPACK_INCLUDE := -I/sw/include/

BLAS_LIBS := $(LAPACK_LIBS)
BLAS_INCLUDE := $(LAPACK_INCLUDE)

DEBUG:= -Wall -DNO_RANGE_CHECKING -DNDEBUG    -DDEBUG_WRITE
#DEBUG := -Wall -DRANGE_CHECKING  -DDEBUG_WRITE
LDFLAGS:=  -framework Accelerate

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
