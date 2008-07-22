####################################################################
# Set platform and build specific information
#   any of these can be overridden by environment variables of the
#   same name, or by command line options "make target PLATFORM=foo"
PLATFORM = $(shell uname -s)
DVLPATH = .

####################################################################
# Derived paths, all relative to the base path.  These shouldn't
# change without restructuring the make environment

MAKEPATH:=$(DVLPATH)/make
OBJPATH:=$(DVLPATH)/obj/$(PLATFORM)
SRCPATH:=$(DVLPATH)
PROJECT:=$(DVLPATH)/qwalk-$(PLATFORM)
DEPENDFILE:=$(OBJPATH)/depend

include make.include
include basis_function/make.include
include macopt/make.include
include method/make.include
include orbitals/make.include
include program_utils/make.include
include properties/make.include
include system/make.include
include wavefunction/make.include

#remove the directory prefix and .cpp extension, and add the .o extension
OBJS:= $(addprefix $(OBJPATH)/,$(addsuffix .o,$(notdir $(basename $(SOURCES)))))
NODEPENDGOALS = depend clean totalclean cleanhtml html all

####################################################################
# Special rules
all: $(PROJECT) gosling  
.PHONY: clean totalclean depend  gosling pw2lcao
#.SILENT:
####################################################################
# Included files for dependency and platform specific defines
include $(MAKEPATH)/$(PLATFORM).mk
include $(DEPENDFILE)

####################################################################
# Project and build rules

$(PROJECT):$(OBJS)
	@echo ________________________________________________________________
	@echo Linking $@
	$(CXX) $(DEBUG) $(CXXFLAGS) -o $(PROJECT) $(OBJS) $(LDFLAGS) $(BLAS_LIBS)

gosling:
	$(MAKE) -f Makefile.gosling

pw2lcao:
	$(MAKE) -f Makefile.pw2lcao

clean:
	@echo ________________________________________________________________
	@echo Cleaning object files...
	rm -f $(OBJS) core $(PROJECT)
	make -f Makefile.gosling clean

include make.rules
