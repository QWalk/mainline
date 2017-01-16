######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
CXX:=g++
CXXFLAGS:= -O2 $(INCLUDEPATH)
DEBUG := -Wall  -DRANGE_CHECKING -DDEBUG_WRITE -g 
######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM  $(INCLUDEPATH)

