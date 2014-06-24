CC      = g++ -std=c++11 -Wall 
platform=$(shell uname)
DEBUG   =  -g -ggdb
RELEASE = -O2
INCLUDE = -I$(CURDIR)/../ -I$(CURDIR)
CFLAGS  = -fPIC
SRCS = $(wildcard *.cpp)
HEADERS = $(patsubst %.cpp,%.hpp,$(SRCS))
OBJS = $(patsubst %.cpp,%.o,$(SRCS))
EXE = $(patsubst %.cpp,%_run,$(SRCS))
LIBRARY = -lsg-lite
LIBINCLUDE = -L../sparsegrid

ifeq  ($(platform), Darwin)
CC+= -stdlib=libstdc++ 
LDFLAGS+= -install_name ${CURDIR}"/"${TARGETLIB}
endif

all: $(EXE)

%_run: %.cpp
	echo $@
	$(CC) $(FLAGS) $(CFLAGS) $(INCLUDE) $(LIBINCLUDE) $(LIBRARY) $< -o $@

.PHONY: clean

clean:
	rm *.o
	rm *.so
	rm *_run