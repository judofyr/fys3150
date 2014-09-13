#include "jacobi.hh"
#include <armadillo>
#include <algorithm>

using namespace arma;

int main() {
  double rho_min = 0;
  double rho_max = 5;
  int n_step = 280;
  double h = (rho_max - rho_min)/n_step;

  int n = n_step - 1;

  mat A; A.zeros(n, n);

  A.diag(+1).fill(-1/(h*h));
  A.diag(-1).fill(-1/(h*h));

  auto diag = A.diag();
  for (int i = 0; i < n; i++) {
    double rho = rho_min + (i+1)*h;
    diag(i) = 2/(h*h) + rho*rho;
  }

  auto solver = JacobiEigenvalue(A);
  int limit = 1e8;
  auto lambda = solver.solve(1e-8, &limit);
  std::sort(lambda.begin(), lambda.end());

  for (auto num : lambda) {
    std::cout << num << std::endl;
  }

  return 0;
}

