CC = g++
PYTHON = python2
CFLAGS = -O3 -std=c++11
LDFLAGS = -larmadillo

.PHONY: clean plot errors test

build:
	mkdir -p build

tests.cc: tridiagonal.hh

build/tridiagonal: tridiagonal.cc | build
	$(CC) $(CFLAGS) $(LDFLAGS) tridiagonal.cc -o build/tridiagonal

build/tests: tests.cc | build
	$(CC) $(CFLAGS) $(LDFLAGS) tests.cc -o build/tests

build/data.json: build/tridiagonal runner.py
	$(PYTHON) runner.py > build/data.json

plot: build/data.json
	$(PYTHON) plot.py

errors: build/data.json
	$(PYTHON) errors.py

test: build/tests
	./build/tests

clean:
	rm -f build
