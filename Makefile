CC=g++
CFLAGS=-c -Wall -g
LDFLAGS= -g
SOURCES= \
akaykinsv.cpp\
akimovada.cpp\
artamonovav.cpp\
bugreevaam.cpp\
venediktovayap.cpp\
vecherskiymp.cpp\
denisovrv.cpp\
zinkinakv.cpp\
kaderovro.cpp\
kochetkovpa.cpp\
makarovaayu.cpp\
melkonyanma.cpp\
melyakinev.cpp\
negryame.cpp\
nikishkinev.cpp\
nuyanzinma.cpp\
pomelovaas.cpp\
prokopenkoas.cpp\
prokopenkods.cpp\
rodkinav.cpp\
ryabikinks.cpp\
timovkinayu.cpp\
turaevdv.cpp\
fedinda.cpp\
zhalninrv.cpp\
lab.cpp\
main.cpp


OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=vvm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o vvm

        