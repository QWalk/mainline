######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
CXX:=openmpicxx
CXXFLAGS:= -O3 -fomit-frame-pointer -funroll-loops -ffast-math 
CXXFLAGS += -DUSE_MPI $(INCLUDEPATH)

#DEBUG := -Wall  -DRANGE_CHECKING -DDEBUG_WRITE
DEBUG:= -Wall -DNO_RANGE_CHECKING   -DNDEBUG 


######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM  $(INCLUDEPATH)
