#include "jacobi.hh"
#include <armadillo>
#include <algorithm>
#include <cstdlib>
#include <map>
#include <chrono>

using namespace arma;

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "usage: " << argv[0]
      << " [arma|jacobi] n_step rho_max [w_r]" << std::endl
      << std::endl
      << "Parameters:" << std::endl
      << "  [arma|jacobi]:      what algorithm to use" << std::endl
      << "  n_step/rho_max/w_r: parameters as specified for the problem" << std::endl
      << std::endl
      << "When w_r is given, it will compute the wave function for two" << std::endl
      << "electrons interecting via a repulsive Coulomb interaction." << std::endl
      << "Otherwise it will compute the wave function for a single electron." << std::endl;
    return 1;
  }

  double rho_min = 0;
  int n_step = atoi(argv[2]);
  double rho_max = atof(argv[3]);
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

  vec eigval;
  mat eigvecs;

  auto t1 = std::chrono::high_resolution_clock::now();
  if (argv[1][0] == 'a') {
    eig_sym(eigval, eigvecs, A);
  } else {
    auto solver = JacobiEigenvalue(A);
    int limit = 1e8;
    eigval = solver.solve(1e-8, &limit);
    eigvecs = solver.R;
    std::cerr << "iterations=" << (1e8 - limit) << std::endl;
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
  std::cerr << "elapsed=" << elapsed.count() << std::endl;

  for (int i = 0; i < n; i++) {
    std::cout << eigval(i) << " ";
    vec eigvec = eigvecs.col(i);
    for (auto num : eigvec) {
      std::cout << num << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}

