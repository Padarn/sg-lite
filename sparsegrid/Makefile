CC      = g++ -std=c++11 -Wall -g -ggdb
CFLAGS  = -I$(CURDIR)/../
LDFLAGS = 

all: sparsegrid.o filereading.o

sparsegrid.o: sparsegrid.cpp sparsegrid.hpp
	$(CC) -c $(CFLAGS) sparsegrid.cpp

.PHONY: clean

clean:
	rm *.o