######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
#CXX:=g++-3.3

#CXX:=icpc
CXX:=/home/apps/intel/compiler80/bin/icc
CXXFLAGS := -O3 \
   -axW  -ip  -DUSE_MPI \
    -DUSE_STDARG -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 \
            -DHAVE_STDARG_H=1 -DUSE_STDARG=1 -DMALLOC_RET_VOID=1 -I/home/apps/mpich-1.2.4/include \
  -I$(INCLUDEPATH) 

#DEBUG := -Wall  -DRANGE_CHECKING -DDEBUG_WRITE 
DEBUG:= -DNO_RANGE_CHECKING -DNDEBUG 
#DEBUG:=-Wall -DNDEBUG -DNO_RANGE_CHECKING -DDEBUG_WRITE -pg 
#LDFLAGS := -lpthread -i_dynamic
LDFLAGS := -static -L/home/apps/mpich-1.2.5/lib -lmpich
######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM -I $(INCLUDEPATH)

