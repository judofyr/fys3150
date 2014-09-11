CC = g++
PYTHON = python2
CFLAGS = -O3 -std=c++11 -Wall -pedantic
LDFLAGS = -larmadillo
DEPS = src/jacobi.hh

.PHONY: clean test

build:
	mkdir -p build

build/tests: test/tests.cc $(DEPS) | build
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

test: build/tests
	./build/tests

clean:
	rm -rf build
