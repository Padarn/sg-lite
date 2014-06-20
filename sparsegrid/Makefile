CC      = g++ -std=c++11 -Wall 
platform=$(shell uname)
DEBUG   =  -g -ggdb
RELEASE = -O2
INCLUDE = -I$(CURDIR)/../ -I$(CURDIR)
CFLAGS  = -fPIC
LDFLAGS =  -shared
SRCS = $(wildcard *.cpp)
HEADERS = $(patsubst %.cpp,%.hpp,$(SRCS))
OBJS = $(patsubst %.cpp,%.o,$(SRCS))
TARGETLIB = libsg-lite.so

ifeq  ($(platform), Darwin)
CC+= -stdlib=libstdc++ 
LDFLAGS+= -install_name ${CURDIR}"/"${TARGETLIB}
endif

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