 ######################################################################

CXX:=clang++

CXXFLAGS := -O3  -g  \
  $(INCLUDEPATH) 
#   -funroll-loops -ffast-math  \
 #not omitting the frame pointer so we get instrumentation

CXXFLAGS+= -DUSE_LAPACK -DUSE_BLAS 
# uncomment here if you've installed ATLAS through Fink, otherwise use the accelerate framework
#LAPACK_LIBS := -L/sw/lib/ -llapack -lf77blas -lcblas -latlas
#LAPACK_INCLUDE := -I/sw/include/

BLAS_LIBS := $(LAPACK_LIBS)
BLAS_INCLUDE := $(LAPACK_INCLUDE)

DEBUG:= -Wall -DNO_RANGE_CHECKING -DNDEBUG   # -DDEBUG_WRITE -DSUPERDEBUG
#DEBUG := -Wall -DRANGE_CHECKING  -DDEBUG_WRITE
LDFLAGS:=  -framework Accelerate

######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM  $(INCLUDEPATH)

