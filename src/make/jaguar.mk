######################################################################
# Compiler definitions for jaguar.ccs.ornl.gov (compute nodes) 
#  all compiler specific information should be declared here

#for gnu with blas and lapack issue this command: 
# module swap PrgEnv-pgi PrgEnv-gnu; module add atlas ; module remove xt-libsci ; module add xt-libsci

CXX:=CC
CXXFLAGS:= -std=gnu++98 -O3 -ffast-math -fomit-frame-pointer \
            -fstrict-aliasing -funroll-loops -march=barcelona 
CXXFLAGS += -DUSE_MPI -DUSE_BLAS -DUSE_LAPACK -DUSE_EINSPLINE -DUSE_RESTRICT ${INCLUDEPATH}

BLAS_LIBS := -L$(ATLASDIR) -lf77blas -lcblas -latlas 
BLAS_INCLUDE := -I$(ATLASINCLUDE)

EINSPLINE_LIBS := -L/ccs/proj/mat001/jaguar/lib -leinspline
EINSPLINE_INCLUDE :=-I/ccs/proj/mat001/jaguar/include/einspline

#EINSPLINE_LIBS := -L/ccs/home/bajdich/codes/lib/lib -leinspline  
#EINSPLINE_INCLUDE :=-I/ccs/home/bajdich/codes/lib/include/einspline

BLAS_LIBS += $(LAPACK_LIBS)
BLAS_INCLUDE += $(LAPACK_LIBS)

BLAS_LIBS += $(EINSPLINE_LIBS)
BLAS_INCLUDE +=$(EINSPLINE_INCLUDE)


DEBUG:= -DNO_RANGE_CHECKING   -DNDEBUG 
######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=CC -MM  $(INCLUDEPATH)

