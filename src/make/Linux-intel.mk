CXX:=ecpc
CXXFLAGS:= -O2 -Ob2 -mcpu=itanium -I$(INCLUDEPATH)
#DEBUG:=  -DRANGE_CHECKING #-DNDEBUG
DEBUG:= -DNO_RANGE_CHECKING -DNDEBUG
LDFLAGS:= ""
######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM -I $(INCLUDEPATH)

