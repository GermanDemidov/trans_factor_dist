CC=g++
CFLAGS=-c -Wall -std=c++11
LDFLAGS=
SOURCES=main.cpp parser.cpp transcription_factor.cpp cell.cpp transcription_factors_in_cell.cpp auxuilary.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=tf_distrib

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o