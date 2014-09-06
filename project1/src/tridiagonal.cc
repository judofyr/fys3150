#include "tridiagonal.hh"
#include <iostream>

using namespace arma;

std::chrono::nanoseconds measure(std::function<void ()> f) {
  auto t1 = std::chrono::high_resolution_clock::now();
  f();
  auto t2 = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);;
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "usage: " << argv[0] << " [LU|Tri] N" << std::endl;
    return 1;
  }

  int n = atoi(argv[2]);
  double h = 1.0 / n;
  vec a(n); a.fill(-1);
  vec b(n); b.fill(2);
  vec c(n); c.fill(-1);

  auto f = [](double x) { return 100 * exp(-10*x); };
  auto problem = TridiagonalProblem(n, h, f, a, b, c);
  auto solver = &TridiagonalProblem::solve;

  if (argv[1][0] == 'L') {
    solver = &TridiagonalProblem::solve_lu;
  }

  vec v;
  auto time = measure([&v, &solver, &problem]() {
    v = (problem.*solver)();
  }).count();

  for (auto y : v) {
    std::cout << y << std::endl;
  }

  std::cerr << "elapsed: " << time << " ns" << std::endl;
}

