######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here

CXX:=mpicxx 
CXX_SERIAL:=g++

CXXFLAGS := -std=gnu++98 -O3 -ffast-math -fomit-frame-pointer \
	    -malign-double -fstrict-aliasing -funroll-loops 
CXXFLAGS += -DUSE_BLAS -DUSE_LAPACK -DUSE_MPI -DUSE_RESTRICT  ${INCLUDEPATH} 

LAPACK_LIBS := -L/usr/lib/atlas-base -llapack -lf77blas -lcblas -latlas -lblas 
LAPACK_INCLUDE := -I/usr/include/atlas

EINSPLINE_LIBS := #-L/usr/local/lib -leinspline
EINSPLINE_INCLUDE := #-I/usr/local/include/einspline

BLAS_LIBS :=  $(LAPACK_LIBS) 
BLAS_INCLUDE := $(LAPACK_INCLUDE) 

BLAS_LIBS += $(EINSPLINE_LIBS)
BLAS_INCLUDE += $(EINSPLINE_INCLUDE)

DEBUG:= -Wall -DNO_RANGE_CHECKING -DNDEBUG   -D__USE_GNU -DDEBUG_WRITE
LDFLAGS:= 

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
