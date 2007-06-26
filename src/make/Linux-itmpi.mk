######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
CXX:=g++
CXXFLAGS:= -O2 -funroll-loops -ffast-math -I$(INCLUDEPATH) -DUSE_MPI \
           -I/usr/local/vmi/mpich/include 
#CXX:=pgCC
#CXXFLAGS:= -O2 --no_exceptions -I$(INCLUDEPATH)
#DEBUG := -Wall  -DRANGE_CHECKING 
DEBUG:= -Wall -DNO_RANGE_CHECKING   -DNDEBUG
#DEBUG:=-Wall -DNDEBUG -DNO_RANGE_CHECKING -pg
#LDFLAGS:=-L/usr/pgi/linux86/lib -lm -lstd
LDFLAGS:=-L/usr/local/vmi/mpich/lib/gcc -lmpich -lvmi -ldl -lpthread


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
