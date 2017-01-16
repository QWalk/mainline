CXX:=xlC
CXXFLAGS:= -qrtti -q64 -qarch=pwr4 $(INCLUDEPATH) -O3 -qstrict
DEBUG:= -DNDEBUG
LDFLAGS:= 

######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM  $(INCLUDEPATH)

