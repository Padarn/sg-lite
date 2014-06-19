CC      = g++ -std=c++11 -Wall -g -ggdb -fPIC
CFLAGS  = -I$(CURDIR)
LDFLAGS = 
OBJECTS = testmain.o sparsegrid/sparsegrid.o sparsegrid/filereading.o

main: $(OBJECTS)
	$(CC) $(OBJECTS) -o test 

testmain.o: testmain.cpp sparsegrid/filereading.hpp 
	$(CC) -c testmain.cpp $(CFLAGS)

sparsegrid/sparsegrid.o: sparsegrid/sparsegrid.cpp sparsegrid/sparsegrid.hpp
	cd sparsegrid; make sparsegrid.o;

sparsegrid/filereading.o: sparsegrid/filereading.cpp sparsegrid/filereading.hpp
	cd sparsegrid; make filereading.o;

.PHONY: clean

clean:
	rm *.o
	rm test
