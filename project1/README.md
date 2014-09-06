# Project 1

Source code for the first project in FYS3150.

## Getting started

Dependencies:

* C++11 compiler
* Python 2 w/ matplotlib

```
# Run tests:
make test

# Plot accuracy of numerical solution:
make plot

# Build LaTeX report:
make article
```

## Project structure

### Main algorithm files

* **src/tridiagonal.hh**: Class for solving tridiagonal problems.
* **src/tridiagonal.cc**: Main function which solves a specific problem.
* **tools/runner.py**: Runs `./tridiagonal` for some specific Ns.
* **test/tests.cc**: Test case for TridiagonalProblem.

### Data analysis programs

All of these programs uses the output from `runner.py`.

* **tools/plot.py**: Plots the accuracy of the numerical solution.
* **tools/perf.py**: Outputs a LaTeX table of the performance of the solvers.
* **tools/errors.py**: Computes the errors of the numerical solution.

