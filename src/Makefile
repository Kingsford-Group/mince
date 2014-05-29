CC = g++
CXXFLAGS = -g -fpermissive -fopenmp $(shell pkg-config --cflags jellyfish-2.0) -I/opt/local/stow/boost-1.55.0-gcc/include -I/opt/local/include -I../include -Wall -O3 -std=c++11
LDFLAGS  = $(shell pkg-config --libs jellyfish-2.0) -L/home/robp/libdivsufsort-2.0.1/build/lib/ -L/opt/local/stow/boost-1.55.0-gcc/lib -lboost_program_options -lboost_system -L/opt/local/lib -ltbb -ltbbmalloc -lz

all: mince
naive: NaiveGraph.cpp
	$(CC) $(CXXFLAGS) NaiveGraph.cpp $(LDFLAGS) -o naive
mince: Mince.cpp
	$(CC) $(CXXFLAGS) Mince.cpp $(LDFLAGS) -o mince
clean:
	rm -rf *o naive mince
