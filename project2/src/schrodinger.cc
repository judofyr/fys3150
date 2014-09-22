#include "jacobi.hh"
#include <armadillo>
#include <algorithm>
#include <cstdlib>

using namespace arma;

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "usage: " << argv[0]
      << " n_step rho_max (val|vec) [w_r]" << std::endl;
    return 1;
  }

  double rho_min = 0;
  int n_step = atoi(argv[1]);
  double rho_max = atof(argv[2]);
  bool print_values = (argv[3][1] == 'a');
  float w_r = (argc < 5) ? -1 : atof(argv[4]);

  double h = (rho_max - rho_min)/n_step;

  int n = n_step - 1;

  mat A; A.zeros(n, n);

  A.diag(+1).fill(-1/(h*h));
  A.diag(-1).fill(-1/(h*h));

  auto diag = A.diag();
  for (int i = 0; i < n; i++) {
    double rho = rho_min + (i+1)*h;
    double potential;

    if (w_r < 0) {
      potential = rho*rho;
    } else {
      potential = w_r*w_r*rho*rho + 1/rho;
    }

    diag(i) = 2/(h*h) + potential;
  }

  auto solver = JacobiEigenvalue(A);
  int limit = 1e8;

  auto lambda = solver.solve(1e-8, &limit);
  vec output;
  std::cerr << "iterations=" << (1e8 - limit) << std::endl;

  if (print_values) {
    output = lambda;
    std::sort(output.begin(), output.end());
  } else {
    auto min_lambda = std::min_element(lambda.begin(), lambda.end());
    auto min_i = min_lambda - lambda.begin();
    output = solver.R.col(min_i);
  }

  for (auto num : output) {
    std::cout << num << std::endl;
  }

  return 0;
}

