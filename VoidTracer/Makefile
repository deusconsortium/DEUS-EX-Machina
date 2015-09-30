CC=g++
CFLAGS=-c -O3 -g -fopenmp -std=c++0x
LDFLAGS=-fopenmp
SOURCES= FOFReaderLib/FOFFiles/FOFCube.cpp \
    FOFReaderLib/FOFFiles/FOFFile.cpp \
    FOFReaderLib/FOFFiles/FOFHalo.cpp \
    FOFReaderLib/FOFFiles/FOFMasst.cpp \
    FOFReaderLib/FOFFiles/FOFMultiCube.cpp \
    FOFReaderLib/FOFFiles/FOFParticles.cpp \
    FOFReaderLib/FOFFiles/FOFStrct.cpp \
    FOFReaderLib/FOFFiles/FOFExtrema.cpp \
    FOFReaderLib/FOFFiles/FOFCubeGrav.cpp \
    FOFReaderLib/FOFFiles/FOFCubeGravLevel.cpp \
    FOFReaderLib/DEUSSimulation/DEUSHalos.cpp \
    FOFReaderLib/DEUSSimulation/DEUSHalo.cpp \
    FOFReaderLib/DEUSSimulation/DEUSCubes.cpp \
    FOFReaderLib/DEUSSimulation/DEUSGrav.cpp \
    FOFReaderLib/DEUSSimulation/DEUSArea.cpp \
    FOFReaderLib/DEUSSimulation/DEUSSimulationSingleton.cpp \
    Simu.cpp \
    Tools.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: $(SOURCES) main
	
main: main.o $(OBJECTS) 
	$(CC) $(LDFLAGS) $@.o $(OBJECTS) -o ProfileTracer

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm -rf *o FOFReaderLib/*.o FOFReaderLib/FOFFiles/*.o FOFReaderLib/DEUSSimulation/*.o *.o ProfileTracer

