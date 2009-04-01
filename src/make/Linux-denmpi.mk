######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
CXX:=g++

MPIDIRECTORY:=/opt/mpich/gnu_unshared

CXXFLAGS:= -O3 -fomit-frame-pointer -funroll-loops -ffast-math 

CXXFLAGS += -DUSE_MPI -DUSE_STDARG -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 \
            -DHAVE_STDARG_H=1 -DUSE_STDARG=1 -DMALLOC_RET_VOID=1 \
            -I$(MPIDIRECTORY)/include $(INCLUDEPATH) 

#CXXFLAGS+= -DUSE_LAPACK -DUSE_BLAS
#LAPACK_LIBS := -L/usr/lib -llapack -lf77blas -lcblas -latlas 
#BLAS_LIBS := $(LAPACK_LIBS)

CXXFLAGS += -DUSE_EINSPLINE -malign-double
EINSPLINE_LIBS := -L/home/apps/einspline-0.8.9/lib -leinspline
EINSPLINE_INCLUDE :=-I/home/apps/einspline-0.8.9/include/einspline

BLAS_LIBS += $(EINSPLINE_LIBS)
BLAS_INCLUDE += $(EINSPLINE_INCLUDE)

#DEBUG := -Wall  -DRANGE_CHECKING -DDEBUG_WRITE
DEBUG:= -Wall -DNO_RANGE_CHECKING   -DNDEBUG -DDEBUG_WRITE
#DEBUG:=-Wall -DNDEBUG -DNO_RANGE_CHECKING -pg
LDFLAGS := -L$(MPIDIRECTORY)/lib -lmpich  -static

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
