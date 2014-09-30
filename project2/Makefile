CC = g++
PYTHON = python2
CFLAGS = -O3 -std=c++11 -Wall -pedantic
LDFLAGS = -larmadillo
DEPS = src/jacobi.hh

.PHONY: clean test all article

all: article

build:
	mkdir -p build

build/schrodinger: src/schrodinger.cc $(DEPS) | build
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

build/tests: test/tests.cc $(DEPS) | build
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

build/iter.tex: tools/run.py build/schrodinger
	$(PYTHON) $< num_iter > $@

build/prec.tex: tools/run.py build/schrodinger
	$(PYTHON) $< prec > $@

build/perf.tex: tools/run.py build/schrodinger
	$(PYTHON) $< perf > $@

build/l1.png: tools/run.py build/schrodinger
	$(PYTHON) $< plot_w_r

plot: tools/plot.py build/wra.txt
	$(PYTHON) $<

build/article.pdf: article.tex build/iter.tex build/prec.tex build/perf.tex build/l1.png | build
	sh build.sh

build/params.tex: tools/fit.py build/schrodinger | build
	$(PYTHON) $<

article: build/article.pdf

test: build/tests
	./build/tests

clean:
	rm -rf build
