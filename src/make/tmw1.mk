######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
#module add openmpi acml
CXX:=mpiCC
#CXX:=/opt/mpich/gnu/bin/mpiCC


CXXFLAGS:= -std=gnu++98 -O3 -ffast-math -fomit-frame-pointer \
            -fstrict-aliasing -funroll-loops  
#CXXFLAGS += -DUSE_MPI -DUSE_BLAS -DUSE_LAPACK -DUSE_EINSPLINE -DUSE_RESTRICT ${INCLUDEPATH}
CXXFLAGS += -DUSE_MPI -DUSE_LAPACK -DUSE_RESTRICT ${INCLUDEPATH}

LAPACK_LIBS := -L/opt/acml/4.3.0/gfortran64/lib -lacml -lgfortran
LAPACK_INCLUDE := -I/opt/acml/4.3.0/gfortran64/include

BLAS_LIBS := $(LAPACK_LIBS)
BLAS_INCLUDE := $(LAPACK_INCLUDE)

#EINSPLINE_LIBS := -L/opt/share/einspline/0.9.1/lib -leinspline
#EINSPLINE_INCLUDE :=-I/opt/share/einspline/0.9.1/include/einspline


BLAS_LIBS += $(EINSPLINE_LIBS)
BLAS_INCLUDE += $(EINSPLINE_INCLUDE)


#DEBUG := -Wall  -DRANGE_CHECKING -DDEBUG_WRITE
DEBUG:= -Wall -DNO_RANGE_CHECKING   -DNDEBUG 
#DEBUG:=-Wall -DNDEBUG -DNO_RANGE_CHECKING -pg
LDFLAGS := -static


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
