CC = g++
PYTHON = python2
CFLAGS = -O3 -std=c++11 -Wall -pedantic
LDFLAGS = -larmadillo
DEPS = src/tridiagonal.hh

.PHONY: clean test

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

build/values.png: tools/plot.py build/data.json
	$(PYTHON) $<

build/errors.tex: tools/errors.py build/data.json
	$(PYTHON) $< > $@

test: build/tests
	./build/tests

TEX_DEPS = build/perf.tex build/values.png build/errors.tex

build/article.pdf: article.tex build.sh $(TEX_DEPS)
	cd build; \
		pdflatex ../article.tex

article: build/article.pdf

clean:
	rm -rf build
