# Project 1

Source code for the second project in FYS3150.

## Getting started

Dependencies:

* C++11 compiler w/ Armadillo
* Python 2 w/ matplotlib

```
# Run tests:
make test

# Build LaTeX report:
make article
```

## Project structure

### Main algorithm files

* **src/jacobi.hh**: Class for solving tridiagonal problems.
* **src/schrodinger.cc**: Main function which solves our specific problem.
* **tools/run.py**: Runs `./schrodinger` and interprets data.
* **test/tests.cc**: Test case for JacobiEigenvalue.

