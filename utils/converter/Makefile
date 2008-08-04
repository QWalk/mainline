CXXFLAGS= -g 

ALLOBJS:=converter.o Pseudo_writer.o wf_writer.o basis_writer.o

GAMESS2QMCOBJS:=gamess2qmc.o $(ALLOBJS)
GAMESSCI2QMCOBJS:=gamessci2qmc.o $(ALLOBJS) 
CRYSTAL2QMCOBJS:=crystal2qmc.o $(ALLOBJS)
SIESTA2QMCOBJS:= siesta2qmc.o $(ALLOBJS)
JEEP2QMCOBJS :=  jeep2qmc.o JeepBasis.o $(ALLOBJS)
G032QMCOBJS := g032qmcmain.o g03tools.o $(ALLOBJS)
HEG2QMCOBJS := heg2qmc.o indexx.o
HEG2QMCCOBJS := heg2qmcc.o indexx.o setk01.o

all: gamess2qmc gamessci2qmc crystal2qmc siesta2qmc g032qmc heg2qmc
.PHONY: clean 

gamess2qmc: $(GAMESS2QMCOBJS)
	$(CXX) -o gamess2qmc $(GAMESS2QMCOBJS)

gamessci2qmc: $(GAMESSCI2QMCOBJS)
	$(CXX) -o gamessci2qmc $(GAMESSCI2QMCOBJS)

crystal2qmc: $(CRYSTAL2QMCOBJS)
	$(CXX) -o crystal2qmc $(CRYSTAL2QMCOBJS)

siesta2qmc: $(SIESTA2QMCOBJS)
	$(CXX) -o siesta2qmc $(SIESTA2QMCOBJS)

jeep2qmc : $(JEEP2QMCOBJS)
	$(CXX) -o jeep2qmc $(JEEP2QMCOBJS)

g032qmc: $(G032QMCOBJS)
	$(CXX) -o g032qmc $(G032QMCOBJS)

heg2qmc : $(HEG2QMCOBJS)
	$(CXX) -o heg2qmc $(HEG2QMCOBJS)

heg2qmcc : $(HEG2QMCCOBJS)
	$(CXX) -o heg2qmcc $(HEG2QMCCOBJS) -lg2c

clean: 
	rm -f $(GAMESS2QMCOBJS) $(GAMESSCI2QMCOBJS) $(CRYSTAL2QMCOBJS) $(JEEP2QMCOBJS) $(SIESTA2QMCOBJS) $(G032QMCOBJS) $(HEG2QMCOBJS) $(HEG2QMCCOBJS) siesta2qmc g032qmc gamess2qmc gamessci2qmc jeep2qmc crystal2qmc heg2qmc heg2qmcc

