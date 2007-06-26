######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
CXX:=g++
CXXFLAGS:= -O2   -I$(INCLUDEPATH)
DEBUG:=-Wall -DNDEBUG -DNO_RANGE_CHECKING -pg


######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM -I $(INCLUDEPATH)


