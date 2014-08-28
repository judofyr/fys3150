CC = g++
PYTHON = python2
CFLAGS = -O3 -std=c++11
LDFLAGS = -larmadillo

.PHONY: clean plot

tests.cc: tridiagonal.hh

tridiagonal: tridiagonal.cc
	$(CC) $(CFLAGS) $(LDFLAGS) tridiagonal.cc -o tridiagonal

tests: tests.cc
	$(CC) $(CFLAGS) $(LDFLAGS) tests.cc -o tests

data.json: tridiagonal runner.py
	$(PYTHON) runner.py > data.json

plot: data.json
	$(PYTHON) plot.py

errors: data.json
	$(PYTHON) errors.py

clean:
	rm -f tests tridiagonal data.json
