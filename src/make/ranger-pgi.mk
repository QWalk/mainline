######################################################################
# Compiler definitions for ranger.tacc.utexas.edu (compute nodes) 
#  all compiler specific information should be declared here

# NOTES:
# 1/ default compiler is PGI, perhaps they have a reason so I will use it;
# 2/ TACC documentation is not quite clear on what is the best BLAS/LAPACK
#    implementation (MKL or ACML or ...), will try MKL first
# 3/ einspline library is not enabled (yet), someone has to compile it first
#    (some problems encountered in linking stage) 

# STEPS BEFORE COMPILATION: module load mkl

CXX:=mpicxx
CXXFLAGS := -tp barcelona-64 -fast 
CXXFLAGS += -DUSE_MPI -DUSE_BLAS -DUSE_LAPACK
#CXXFLAGS += -DUSE_EINSPLINE -DUSE_RESTRICT
CXXFLAGS += $(INCLUDEPATH)

BLAS_LIBS := -Wl,-rpath,$(TACC_MKL_LIB) -L$(TACC_MKL_LIB) -lmkl_em64t -lguide
BLAS_INCLUDE := -I$(TACC_MKL_INC)

#EINSPLINE_LIBS := -L/ccs/proj/mat001/jaguar/lib -leinspline
#EINSPLINE_INCLUDE :=-I/ccs/proj/mat001/jaguar/include/einspline

DEBUG:= -DNO_RANGE_CHECKING   -DNDEBUG 
######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=pgCC -MM  $(INCLUDEPATH)

