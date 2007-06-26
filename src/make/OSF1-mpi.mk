######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
CXX:=cxx
CXXFLAGS:= -DUSE_MPI -ieee -underflow_to_zero  -fast -I$(INCLUDEPATH)
DEBUG:= -DNO_RANGE_CHECKING -DNDEBUG # -DDEBUG_WRITE
LDFLAGS:= -lmpi  -lm
#LDFLAGS:= -lm

######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=c++ -MM -I $(INCLUDEPATH)

