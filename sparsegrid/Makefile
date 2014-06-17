CC      = g++ -std=c++11 -Wall -g -ggdb -fPIC -shared
CFLAGS  = -I$(CURDIR)/../
LDFLAGS = 

all: sparsegrid.o filereading.o

sparsegrid.o: sparsegrid.cpp sparsegrid.hpp
	$(CC) -c $(CFLAGS) sparsegrid.cpp

filereading.o: filereading.cpp filereading.hpp
	$(CC) -c $(CFLAGS) filereading.cpp

.PHONY: clean

clean:
	rm *.o