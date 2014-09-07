#define CATCH_CONFIG_MAIN

#include "catch.hh"
#include "../src/tridiagonal.hh"
#include <cmath>
#include <armadillo>

using namespace arma;

TEST_CASE("Tridiagonal poisson") {
  int n = 100;
  double h = 1.0 / (n+1);
  vec a(n); a.fill(-1);
  vec b(n); b.fill(2);
  vec c(n); c.fill(-1);

  SECTION("f(x) = 100 * exp(-10*x)") {
    auto f = [h](double x) { return h * h * 100 * exp(-10*x); };
    auto problem = TridiagonalProblem(n, h, h, f, a, b, c);

    SECTION("solve() and solve_lu() should return the same value") {
      auto v1 = problem.solve();
      auto v2 = problem.solve_lu();
      for (int i = 0; i < n; i++) {
        REQUIRE(v1[i] == v2[i]);
      }
    }

    SECTION("is kinda similar to the exact solution") {
      auto exact = [](double x) { return 1 - (1 - exp(-10))*x - exp(-10*x); };

      auto v = problem.solve();
      double x = 0;
      for (int i = 0; i < n; i++) {
        REQUIRE(fabs(v[i] - exact(x)) < 0.1);
        x += h;
      }
    }
  }
}

