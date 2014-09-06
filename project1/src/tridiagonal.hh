#include <armadillo>
#include <chrono>
#include <functional>

class TridiagonalProblem {
  typedef std::function<double (double)> doublef;

  public:
    int n;
    double h;
    arma::vec a, b, c;
    doublef f;

    TridiagonalProblem(int n, double h, doublef f, arma::vec a, arma::vec b, arma::vec c)
      : n(n), h(h), a(a), b(b), c(c), f(f) {
    }

    arma::vec solve() {
      arma::vec v(n);
      double x = 0;

      arma::vec bm = b;
      arma::vec btwiddle(n);

      btwiddle(0) = h*h*f(x);

      for (int i = 0; i < n-1; i++) {
        x += h;
        double ratio = a(i+1)/bm(i);
        bm(i+1) = bm(i+1) - ratio*c(i);
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
        v(i) = (btwiddle(i) - c(i)*v(i+1))/bm(i);
      }

      return v;
    }

    arma::vec solve_lu() {
      arma::vec bi(n);
      double x = 0;
      arma::mat A(n, n); A.zeros();

      for (int i = 0; i < n; i++) {
        if (i > 0)   A(i, i-1) = a(i);
                     A(i, i)   = b(i);
        if (i < n-1) A(i, i+1) = c(i);

        bi(i) = h*h*f(x);
        x += h;
      }

      arma::mat l, u, p;
      arma::vec y;

      arma::lu(l, u, p, A);
      arma::solve(y, l, p*bi);

      return arma::solve(u, y);
    }
};

