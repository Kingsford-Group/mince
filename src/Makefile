CC = g++
CXXFLAGS = -g -fpermissive -fopenmp $(shell pkg-config --cflags jellyfish-2.0) -I/opt/local/stow/boost-1.55.0-gcc/include -I/opt/local/include -I../include -Wall -O3 -std=c++11 -Wno-unused-local-typedefs
LDFLAGS  = $(shell pkg-config --libs jellyfish-2.0) -L/opt/local/stow/boost-1.55.0-gcc/lib -lboost_program_options -lboost_system -L/opt/local/lib -ltbb -ltbbmalloc -lz

all: mince
LibraryFormat.o: LibraryFormat.cpp ../include/LibraryFormat.hpp
	$(CC) $(CXXFLAGS) LibraryFormat.cpp -c 
naive: NaiveGraph.cpp
	$(CC) $(CXXFLAGS) NaiveGraph.cpp $(LDFLAGS) -o naive
mince: Mince.cpp LibraryFormat.o
	$(CC) $(CXXFLAGS) Mince.cpp LibraryFormat.o $(LDFLAGS) -o mince
canonicalize: Canonicalize.cpp
	$(CC) $(CXXFLAGS) Canonicalize.cpp $(LDFLAGS) -o canonicalize
mergepairs: MergePairs.cpp LibraryFormat.o
	$(CC) $(CXXFLAGS) LibraryFormat.o MergePairs.cpp $(LDFLAGS) -o mergepairs


clean:
	rm -rf *o naive mince
