CC      = g++ -std=c++11 -Wall -g -ggdb
CFLAGS  = -I$(CURDIR)
LDFLAGS = 
OBJECTS = testmain.o sparsegrid/sparsegrid.o

main: $(OBJECTS)
	$(CC) $(OBJECTS) -o test 

testmain.o: testmain.cpp sparsegrid/filereading.hpp
	$(CC) -c testmain.cpp $(CFLAGS)

sparsegrid/sparsegrid.o: sparsegrid/sparsegrid.cpp sparsegrid/sparsegrid.hpp
	cd sparsegrid; make sparsegrid.o;

.PHONY: clean

clean:
	rm *.o
	rm test
