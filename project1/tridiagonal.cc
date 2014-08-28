#include <armadillo>
#include <iostream>
#include <chrono>
#include <functional>

using namespace arma;

typedef std::function<double (double)> doublef;

class TridiagonalProblem {
  public:
    int n;
    double h;
    vec a, b, c;
    doublef f;

    TridiagonalProblem(int n, double h, doublef f, vec a, vec b, vec c)
      : n(n), h(h), f(f), a(a), b(b), c(c) {
    }

    vec solve() {
      vec v(n);
      double x = 0;

      vec am = a;
      vec bm = b;
      vec cm = c;
      vec btwiddle(n);

      btwiddle(0) = h*h*f(x);

      for (int i = 0; i < n-1; i++) {
        x += h;
        double ratio = am(i+1)/bm(i);
        bm(i+1) = bm(i+1) - ratio*cm(i);
        btwiddle(i+1) = h*h*f(x) - ratio*btwiddle(i);
      }

      // At this point we have:
      //   b_i * v_i + c_i * v_{i+1} = btwiddle_i
      // Rearrange this into:
      //   v_i = (btwiddle_i - c_i * v_{i+1}) / b_i
      // The last element on the other hand is:
      //   v_n = btwiddle_n / b_n

      v(n-1) = btwiddle(n-1) / bm(n-1);

      for (int i = n-2; i >= 0; i--) {
        v(i) = (btwiddle(i) - cm(i)*v(i+1))/bm(i);
      }

      return v;
    }

    vec solve_lu() {
      vec bi(n);
      double x = 0;
      mat A(n, n); A.zeros();

      for (int i = 0; i < n; i++) {
        if (i > 0)   A(i, i-1) = a(i);
                     A(i, i)   = b(i);
        if (i < n-1) A(i, i+1) = c(i);

        bi(i) = h*h*f(x);
        x += h;
      }

      return arma::solve(A, bi);
    }
};

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

  std::cout << v;
  std::cerr << "elapsed: " << time << " ns" << std::endl;
}

