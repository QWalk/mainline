######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
CXX:=cxx
CXXFLAGS:= -ieee -underflow_to_zero  -fast -I$(INCLUDEPATH)
DEBUG:= -DNO_RANGE_CHECKING -DNDEBUG # -DDEBUG_WRITE
LDFLAGS:= -lm

######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=c++ -MM -I $(INCLUDEPATH)

