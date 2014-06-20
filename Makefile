CC      = 
CFLAGS  = 
LDFLAGS = 
OBJECTS = 

all:  
	cd sparsegrid; make
	cd swig; make

debug:
	cd sparsegrid; make debug
	cd swig; make

release:
	cd sparsegrid; make release
	cd swig; make

clean:
	cd sparsegrid; make clean
	cd swig; make clean
