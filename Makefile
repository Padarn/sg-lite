CC      = 
CFLAGS  = 
LDFLAGS = 
OBJECTS = 

all:  
	cd sparsegrid; make
	cd tests; make

debug:
	cd sparsegrid; make debug
	cd tests; make

release:
	cd sparsegrid; make release
	cd tests; make

clean:
	cd sparsegrid; make clean
	cd swig; make clean
	cd tests; make clean

swig:
	cd swig; make