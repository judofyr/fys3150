CC = g++
PYTHON = python2
CFLAGS = -O3 -std=c++11 -Wall -pedantic
LDFLAGS = -larmadillo
DEPS = src/jacobi.hh

.PHONY: clean test all

all: build/schrodinger

build:
	mkdir -p build

build/schrodinger: src/schrodinger.cc $(DEPS) | build
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

build/tests: test/tests.cc $(DEPS) | build
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

test: build/tests
	./build/tests

clean:
	rm -rf build
