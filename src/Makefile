CC = g++
CXXFLAGS = -g -fopenmp $(shell pkg-config --cflags jellyfish-2.0) -I/opt/local/stow/boost-1.55.0-gcc/include -I/opt/local/include -Wall -O3 -std=c++11
LDFLAGS  = $(shell pkg-config --libs jellyfish-2.0) -L/opt/local/stow/boost-1.55.0-gcc/lib  -lboost_system -L/opt/local/lib -ltbb -ltbbmalloc

all: bucket
naive: NaiveGraph.cpp
	$(CC) $(CXXFLAGS) NaiveGraph.cpp $(LDFLAGS) -o naive
bucket: BucketByMinimizer.cpp
	$(CC) $(CXXFLAGS) BucketByMinimizer.cpp $(LDFLAGS) -o bucket
clean:
	rm -rf *o naive bucket
