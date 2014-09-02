CC = g++
PYTHON = python2
CFLAGS = -O3 -std=c++11 -Wall -pedantic
LDFLAGS = -larmadillo
DEPS = tridiagonal.hh

.PHONY: clean plot errors test

build:
	mkdir -p build

build/tridiagonal: tridiagonal.cc $(DEPS) | build
	$(CC) $(CFLAGS) $(LDFLAGS) tridiagonal.cc -o build/tridiagonal

build/tests: tests.cc $(DEPS) | build
	$(CC) $(CFLAGS) $(LDFLAGS) tests.cc -o build/tests

build/data.json: build/tridiagonal runner.py
	$(PYTHON) runner.py > build/data.json

build/perf.tex: build/data.json perf.py
	$(PYTHON) perf.py > build/perf.tex

plot: build/data.json
	$(PYTHON) plot.py

errors: build/data.json
	$(PYTHON) errors.py

test: build/tests
	./build/tests

build/article.pdf: article.tex build/perf.tex | build
	cd build; \
		pdflatex ../article.tex

article: build/article.pdf

clean:
	rm -rf build
