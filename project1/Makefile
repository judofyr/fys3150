CC = g++
PYTHON = python2
CFLAGS = -O3 -std=c++11 -Wall -pedantic
LDFLAGS = -larmadillo
DEPS = src/tridiagonal.hh

.PHONY: clean plot errors test

build:
	mkdir -p build

build/tridiagonal: src/tridiagonal.cc $(DEPS) | build
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

build/tests: test/tests.cc $(DEPS) | build
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

build/data.json: tools/runner.py build/tridiagonal
	$(PYTHON) $< > $@

build/perf.tex: tools/perf.py build/data.json
	$(PYTHON) $< > $@

plot: build/data.json
	$(PYTHON) tools/plot.py

errors: build/data.json
	$(PYTHON) tools/errors.py

test: build/tests
	./build/tests

build/article.pdf: article.tex build/perf.tex | build
	cd build; \
		pdflatex ../article.tex

article: build/article.pdf

clean:
	rm -rf build
