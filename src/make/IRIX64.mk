CXX:=KCC
CXXFLAGS:= -I$(INCLUDEPATH) -O3
LDFLAGS:= -lm

#We can also use SGI's compiler, although it's more or less inferior to 
#either KCC or gcc.  
#CXX:=CC
#CXXFLAGS:= -LANG:std -I$(INCLUDEPATH) -O2
#LDFLAGS := -lm
