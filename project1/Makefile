CC = g++
PYTHON = python2
CFLAGS = -O3 -std=c++11 -g
LDFLAGS = -larmadillo

.PHONY: run

tridiagonal: tridiagonal.cc
	$(CC) $(CFLAGS) $(LDFLAGS) tridiagonal.cc -o tridiagonal

run: tridiagonal
	$(PYTHON) runner.py

