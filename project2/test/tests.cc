#define CATCH_CONFIG_MAIN

#include "../../support/catch.hh"
#include "../src/jacobi.hh"
#include <cmath>
#include <armadillo>

using namespace arma;

#define REQUIRE_EPS(exp, act) REQUIRE(fabs((exp) - (act)) < 0.0001)

TEST_CASE("JacobiEigenvalue") {
  SECTION("A 2x2 matrix") {
    mat A;
    A << 1 << 2 << endr
      << 2 << 4 << endr;

    auto solver = JacobiEigenvalue(A);

    SECTION("step-by-step") {
      int k, l;
      double max = solver.max_off(&k, &l);
      REQUIRE(max == 2);
      REQUIRE(k == 1);
      REQUIRE(l == 0);

      double c, s;
      solver.theta(k, l, &c, &s);

      // Values taken from: http://fptchlx02.tu-graz.ac.at/cgi-bin/access.com?c1=0000&c2=0000&c3=0000&file=0638
      REQUIRE_EPS(0.8944, c);
      REQUIRE_EPS(-0.4472, s);

      solver.rotate(k, l);

      REQUIRE_EPS(5, solver.A(1, 1));
      REQUIRE_EPS(0, solver.A(0, 0));
    }

    SECTION("solving") {
      int limit = 100;
      auto lambda = solver.solve(1e-8, &limit);
      REQUIRE(limit == 99); // one iteration
      REQUIRE_EPS(0, lambda(0));
      REQUIRE_EPS(5, lambda(1));
    }
  }

  SECTION("A 3x3 matrix") {
    mat A;
    A << 1 << 2 << 3 << endr
      << 2 << 4 << 5 << endr
      << 3 << 5 << 6 << endr;

    auto solver = JacobiEigenvalue(A);

    int limit = 100;
    auto lambda = solver.solve(1e-8, &limit);

    REQUIRE(limit > 90);
    REQUIRE_EPS(-0.5157, lambda(0));
    REQUIRE_EPS(0.1709, lambda(1));
    REQUIRE_EPS(11.3448, lambda(2));
  }
}

