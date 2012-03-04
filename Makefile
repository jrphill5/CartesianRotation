CXX=g++
CXXFLAGS=-g
LIBS=-lm
OBJECTS=main.o
PROJECT=CartesianRotation

all: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LIBS) $(OBJECTS) -o $(PROJECT)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

clean:
	rm *.o $(PROJECT)
