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

build/wra.txt: build/schrodinger | build
	$< 280 5 vec 0.01 > build/wra.txt
	$< 280 5 vec    1 > build/wrb.txt
	$< 280 5 vec    5 > build/wrc.txt

plot: tools/plot.py build/wra.txt
	$(PYTHON) $<

test: build/tests
	./build/tests

clean:
	rm -rf build
