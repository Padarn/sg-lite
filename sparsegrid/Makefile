CC      = g++ -std=c++11 -stdlib=libstdc++ -Wall 
DEBUG   =  -g -ggdb
RELEASE = -O2
INCLUDE = -I$(CURDIR)/../ -I$(CURDIR)
CFLAGS  = -fPIC
LDFLAGS =  -shared
SRCS = $(wildcard *.cpp)
HEADERS = $(patsubst %.cpp,%.hpp,$(SRCS))
OBJS = $(patsubst %.cpp,%.o,$(SRCS))
TARGETLIB = libsg-lite.so

all: ${TARGETLIB}

debug: 
	export FLAGS=$(DEBUG); make all 

release: 
	export FLAGS=$(RELEASE); make all

$(TARGETLIB): $(OBJS)
	$(CC) -o $@ $(LDFLAGS) $^

%.o: %.cpp %.hpp
	$(CC) $(FLAGS) $(CFLAGS) $(INCLUDE) -c $< >$@

.PHONY: clean

clean:
	rm *.o
	rm *.so