
# This affects the memory usage of the program
# we use 1 byte for every 4 bp in kmers. Ideally
# this parameter should be a multiple of 4.
# Actual maximum kmer size is 1 less.
MAX_KMER_SIZE = 32

CC = g++
CXX = g++
INCLUDES = -I.
CXXFLAGS = -c -Wall -Wno-reorder $(INCLUDES) -DMAX_KMER_SIZE=$(MAX_KMER_SIZE) -fPIC -fopenmp
LDFLAGS =
LDLIBS  = -lm -lz -lgomp

all: CXXFLAGS += -O3
all: target


debug: CXXFLAGS += -g -O0
debug: LDFLAGS += -g
debug: target

profile: CXXFLAGS += -p -g -O2
profile: LDFLAGS += -p -g
profile: clean
profile: target

target: BFCounter

OBJECTS =  CountBF.o DumpBF.o Kmer.o KmerIntPair.o KmerIterator.o hash.o fastq.o 

testread: testread.o $(OBJECTS)
	$(CC) $(INCLUDES) $(LDFLAGS) $(LDLIBS) $(OBJECTS) testread.o -o testread

BFCounter: BFCounter.o $(OBJECTS)
	$(CC) $(INCLUDES) $(LDFLAGS) $(LDLIBS) $(OBJECTS) BFCounter.o -o BFCounter -lstdc++

Naive: Naive.o $(OBJECTS)
	$(CC) $(INCLUDES) $(LDFLAGS) $(LDLIBS) $(OBJECTS) Naive.o -o Naive


BFCounter.o: BFCounter.cpp
CountBF.o: CountBF.cpp
DumpBF.o: DumpBF.cpp
KmerIntPair.o: KmerIntPair.cpp
KmerIterator.o: KmerIterator.hpp KmerIterator.cpp
fastq.o: fastq.hpp fastq.cpp 
kmer.o: kmer.hpp kmer.cpp
hash.o: hash.hpp hash.cpp	

clean:
	rm -rf *.o
	rm -rf BFCounter
