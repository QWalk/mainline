PLATFORM = $(shell uname -s)
DVLPATH = $(shell pwd)
##########################using the make files of qwalk########################
MAKEPATH:=$(DVLPATH)/../make
OBJPATH:=$(DVLPATH)/obj/$(PLATFORM)

ALLOBJS:=$(OBJPATH)/converter.o $(OBJPATH)/Pseudo_writer.o $(OBJPATH)/wf_writer.o $(OBJPATH)/basis_writer.o
include $(MAKEPATH)/$(PLATFORM).mk
GAMESS2QMCOBJS:=$(OBJPATH)/gamess2qmc.o $(ALLOBJS)
GAMESSCI2QMCOBJS:=$(OBJPATH)/gamessci2qmc.o $(ALLOBJS) 
CRYSTAL2QMCOBJS:=$(OBJPATH)/crystal2qmc.o $(ALLOBJS)
SIESTA2QMCOBJS:= $(OBJPATH)/siesta2qmc.o $(ALLOBJS)
ABINIT2QMCOBJS := $(OBJPATH)/abinit2qmc.o $(ALLOBJS)
G032QMCOBJS := $(OBJPATH)/g032qmcmain.o $(OBJPATH)/g03tools.o $(ALLOBJS)
HEG2QMCOBJS := $(OBJPATH)/heg2qmc.o $(OBJPATH)/indexx.o
HEG2QMCCOBJS := $(OBJPATH)/heg2qmcc.o $(OBJPATH)/indexx.o $(OBJPATH)/setk01.o
NWCHEM2QMCOBJS := $(OBJPATH)/nwchem2qmc.o $(ALLOBJS)
SQD2QMCOBJS:=$(OBJPATH)/sqd2qmc.o $(ALLOBJS)

SUFFIX := -$(PLATFORM) 

all: nwchem2qmc$(SUFFIX) gamess2qmc$(SUFFIX) gamessci2qmc$(SUFFIX) crystal2qmc$(SUFFIX) siesta2qmc$(SUFFIX) g032qmc$(SUFFIX) heg2qmc$(SUFFIX) abinit2qmc$(SUFFIX) 
.PHONY: clean 

gamess2qmc$(SUFFIX): $(GAMESS2QMCOBJS)
	@echo ________________________________________________________________
	@echo Linking $@
	$(CXX) $(DEBUG) $(CXXFLAGS) -o gamess2qmc$(SUFFIX) $(LDFLAGS) $(BLAS_LIBS) $(GAMESS2QMCOBJS)

gamessci2qmc$(SUFFIX): $(GAMESSCI2QMCOBJS)
	@echo ________________________________________________________________
	@echo Linking $@
	$(CXX) $(DEBUG) $(CXXFLAGS) -o gamessci2qmc$(SUFFIX) $(LDFLAGS) $(BLAS_LIBS) $(GAMESSCI2QMCOBJS)

crystal2qmc$(SUFFIX): $(CRYSTAL2QMCOBJS)
	@echo ________________________________________________________________
	@echo Linking $@
	$(CXX) $(DEBUG) $(CXXFLAGS) -o crystal2qmc$(SUFFIX) $(LDFLAGS) $(BLAS_LIBS) $(CRYSTAL2QMCOBJS)

crystal2qmc_all$(SUFFIX): $(CRYSTAL2QMCALLOBJS)
	@echo ________________________________________________________________
	@echo Linking $@
	$(CXX) $(DEBUG) $(CXXFLAGS) -o crystal2qmc_all$(SUFFIX) $(LDFLAGS) $(BLAS_LIBS) $(CRYSTAL2QMCALLOBJS)

nwchem2qmc$(SUFFIX): $(NWCHEM2QMCOBJS)
	@echo ________________________________________________________________
	@echo Linking $@
	$(CXX) $(DEBUG) $(CXXFLAGS) -o nwchem2qmc$(SUFFIX) $(LDFLAGS) $(BLAS_LIBS) $(NWCHEM2QMCOBJS)

siesta2qmc$(SUFFIX): $(SIESTA2QMCOBJS)
	@echo ________________________________________________________________
	@echo Linking $@
	$(CXX) $(DEBUG) $(CXXFLAGS) -o siesta2qmc$(SUFFIX) $(LDFLAGS) $(BLAS_LIBS) $(SIESTA2QMCOBJS)

abinit2qmc$(SUFFIX) : $(ABINIT2QMCOBJS)
	@echo ________________________________________________________________
	@echo Linking $@
	$(CXX) $(DEBUG)  $(CXXFLAGS) -o abinit2qmc$(SUFFIX) $(LDFLAGS) $(BLAS_LIBS) $(ABINIT2QMCOBJS)

abinit2qmc2$(SUFFIX) : $(ABINIT2QMC2OBJS)
	@echo ________________________________________________________________  
	@echo Linking $@
	$(CXX) $(DEBUG) -I$(INCLUDEPATH) $(CXXFLAGS) -o abinit2qmc2$(SUFFIX) $(LDFLAGS) $(BLAS_LIBS) $(ABINIT2QMC2OBJS)

g032qmc$(SUFFIX) : $(G032QMCOBJS)
	@echo ________________________________________________________________
	@echo Linking $@
	$(CXX) $(DEBUG) $(CXXFLAGS) -o g032qmc$(SUFFIX) $(LDFLAGS) $(BLAS_LIBS) $(G032QMCOBJS)

heg2qmc$(SUFFIX) : $(HEG2QMCOBJS)
	@echo ________________________________________________________________
	@echo Linking $@
	$(CXX) $(DEBUG) $(CXXFLAGS) -o heg2qmc$(SUFFIX) $(LDFLAGS) $(BLAS_LIBS) $(HEG2QMCOBJS)

heg2qmcc$(SUFFIX) : $(HEG2QMCCOBJS)
	@echo ________________________________________________________________
	@echo Linking $@
	$(CXX) $(DEBUG) $(CXXFLAGS) -o heg2qmcc$(SUFFIX) $(LDFLAGS) $(BLAS_LIBS) $(HEG2QMCCOBJS) -lgfortran #-lg2c

sqd2qmc$(SUFFIX) : $(SQD2QMCOBJS)
	@echo ________________________________________________________________
	@echo Linking $@
	$(CXX) $(DEBUG) $(CXXFLAGS) -o sqd2qmc$(SUFFIX) $(LDFLAGS) $(BLAS_LIBS) $(SQD2QMCOBJS) $(XML_LIBS) $(HDF_LIBS) 

clean: 
	@echo ________________________________________________________________
	@echo Cleaning object files... 
	rm -f $(GAMESS2QMCOBJS) $(NWCHEM2QMCOBJS) $(GAMESSCI2QMCOBJS) $(CRYSTAL2QMCOBJS) $(JEEP2QMCOBJS) $(SIESTA2QMCOBJS) $(G032QMCOBJS) $(HEG2QMCOBJS) $(HEG2QMCCOBJS) $(ABINIT2QMCOBJS)$(SQD2QMCOBJS) siesta2qmc$(SUFFIX) g032qmc$(SUFFIX) gamess2qmc$(SUFFIX) gamessci2qmc$(SUFFIX) jeep2qmc$(SUFFIX) crystal2qmc$(SUFFIX) heg2qmc$(SUFFIX) heg2qmcc$(SUFFIX) abinit2qmc$(SUFFIX) sqd2qmc$(SUFFIX) nwchem2qmc$(SUFFIX) 

$(OBJPATH)/%.o:$(DVLPATH)/%.cpp
	@echo ________________________________________________________________
	@echo Compiling $<
	[ -e $(DVLPATH)/obj ] || mkdir $(DVLPATH)/obj
	[ -e $(OBJPATH) ] || mkdir $(OBJPATH)
	$(CXX) -c $(CXXFLAGS) $(BLAS_INCLUDE) $(CPPFLAGS) -o $@ $< $(XML_INCLUDE) $(HDF_INCLUDE)
        

$(OBJPATH)/%.o:$(DVLPATH)/%.f
	@echo ________________________________________________________________
	@echo Compiling $<
	[ -e $(DVLPATH)/obj ] || mkdir $(DVLPATH)/obj
	[ -e $(OBJPATH) ] || mkdir $(OBJPATH)
	$(F77) -c -o $@ $<
