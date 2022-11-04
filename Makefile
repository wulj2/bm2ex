CC=g++
CPPFLAGS=-std=c++11 -g -Wall -O3
LINKLIB=-lhts -lz

all: bm2ex

bm2ex:
	$(CC) $(CPPFLAGS) $(LINKLIB) -o bm2ex bm2ex.cpp

.PHONY:all clean cleanlocal

cleanlocal:
	rm -rf bm2ex bm2ex.dSYM

clean:cleanlocal
