######################################################################
# Compiler definitions for ranger.tacc.utexas.edu (compute nodes) 
#  all compiler specific information should be declared here

# NOTES:
# 1/ Intel compiler (short testing indicates that it provides faster code than
#    PGI; besides, no luck so far with compiling einspline with PGI on ranger)
# 2/ TACC documentation is not quite clear on what is the best BLAS/LAPACK
#    implementation (MKL or ACML or ...), will try MKL first
# 3/ einspline library is enabled, user is expected to have it instaled under
#    $(HOME)/einspline ; since RPATH (runtime search path) is used, it is NOT
#    neccessary to include $HOME/einspline/lib into LD_LIBRARY_PATH

# STEPS BEFORE COMPILATION:
#  module unload pgi mvapich; module load intel mvapich mkl 

# COMILATION OF EINSPLINE
#  ./configure --prefix=$HOME/einspline CC=icc CXX=icpc F77=ifort CFLAGS="-O2 -xW" CXXFLAGS="-O2 -xW" FFLAGS="-O2 -xW"

CXX:=mpicxx
CXXFLAGS := -O2 -xW 
CXXFLAGS += -DUSE_MPI
CXXFLAGS += -DUSE_BLAS -DUSE_LAPACK
#CXXFLAGS += -DUSE_EINSPLINE -DUSE_RESTRICT
CXXFLAGS += $(INCLUDEPATH)

#BLAS_LIBS :=
#BLAS_INCLUDE :=

BLAS_LIBS := -Wl,-rpath,$(TACC_MKL_LIB) -L$(TACC_MKL_LIB) -lmkl_em64t -lguide
BLAS_INCLUDE := -I$(TACC_MKL_INC)

#EINSPLINE_LIBS := -Wl,-rpath,$(HOME)/einspline/lib -L$(HOME)/einspline/lib -leinspline
#EINSPLINE_INCLUDE :=-I$(HOME)/einspline/include/einspline

DEBUG:= -DNO_RANGE_CHECKING   -DNDEBUG 
######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=icc -MM  $(INCLUDEPATH)

