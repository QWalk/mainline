######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here

CXX:=g++-3.3
CXXFLAGS := -O3 \
   -funroll-loops -ffast-math \
  -I$(INCLUDEPATH) -fomit-frame-pointer

CXXFLAGS+= -DUSE_LAPACK -DUSE_BLAS

LAPACK_LIBS := -L/home/apps/ATLAS/lib/Linux_ATHLONSSE1 -llapack -lf77blas -lcblas -latlas -static
LAPACK_INCLUDE := -I/home/apps/ATLAS/include

BLAS_LIBS := $(LAPACK_LIBS)
BLAS_INCLUDE := $(LAPACK_INCLUDE)

#CXXFLAGS:= -O2 -I$(INCLUDEPATH)
#DEBUG := -Wall  -DRANGE_CHECKING -DDEBUG_WRITE 
#DEBUG:= -Wall -DNO_RANGE_CHECKING   -DDEBUG_WRITE
DEBUG:= -Wall -DNO_RANGE_CHECKING -DNDEBUG   -D__USE_GNU
#DEBUG:=-Wall -DNDEBUG -DNO_RANGE_CHECKING -DDEBUG_WRITE -D__USE_GNU -pg 
#LDFLAGS:=-L/usr/pgi/linux86/lib -lm -lstd
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
