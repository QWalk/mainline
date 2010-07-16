######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
CXX:=mpicxx


CXXFLAGS:= -O3 -fomit-frame-pointer -funroll-loops -ffast-math 
#CXXFLAGS := -O2  -fomit-frame-pointer -funroll-loops -ffast-math
CXXFLAGS += -DUSE_MPI $(INCLUDEPATH)

#DEBUG := -Wall  -DRANGE_CHECKING -DDEBUG_WRITE
DEBUG:= -Wall -DNO_RANGE_CHECKING   -DNDEBUG 
#DEBUG:=-Wall -DNDEBUG -DNO_RANGE_CHECKING -pg
#LDFLAGS := -static


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
