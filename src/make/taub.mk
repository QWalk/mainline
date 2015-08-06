######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
# Run the following to load the correct libraries:
# module load gcc/4.7.1
# module load openmpi/1.4-gcc
# module load intel/11.1


CXX:= mpicxx
CXX_SERIAL:=g++

CXXFLAGS := -O2
CXXFLAGS += -DUSE_MPI -DUSE_BLAS -DUSE_LAPACK  ${INCLUDEPATH}

LAPACK_LIBS :=   -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread  -lm
LAPACK_INCLUDE :=  -I$(MKLROOT)/include

BLAS_LIBS :=  $(LAPACK_LIBS)
BLAS_INCLUDE := $(LAPACK_INCLUDE)

DEBUG:= -Wall -DNO_RANGE_CHECKING -DNDEBUG   -D__USE_GNU -DDEBUG_WRITE
LDFLAGS:=

######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM  $(INCLUDEPATH)


