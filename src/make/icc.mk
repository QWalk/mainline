######################################################################
# Compiler definitions for Linux systems
#  all compiler specific information should be declared here
#also try -O3 -xW -tpp7 -ipo -ipo_obj
CXX:=icpc
CXXFLAGS := -O3 -fno-alias -xW -vec_report3  $(INCLUDEPATH)

DEBUG:= -DNO_RANGE_CHECKING -DNDEBUG -DDEBUG_WRITE 


######################################################################
# This is the invokation to generate dependencies
DEPENDMAKER:=g++ -MM -I $(INCLUDEPATH)

