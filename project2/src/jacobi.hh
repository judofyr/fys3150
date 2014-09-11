#include <armadillo>

class JacobiEigenvalue {
  public:
    arma::mat A;
    arma::mat R;
    int n;

    JacobiEigenvalue(arma::mat A) : A(A) {
      n = A.n_cols;
      R.eye(n, n);
    }

    arma::vec solve(double eps, int *limit) {
      int k = 0, l = 0;
      double max = max_off(&k, &l);

      while (max > eps && *limit > 0) {
        rotate(k, l);
        max = max_off(&k, &l);
        (*limit)--;
      }

      return A.diag();
    }

    double max_off(int *k, int *l) {
      int c = 0, r = 0;
      double max = 0;

      for (int c = 1; c < n; c++) {
        for (int r = 0; r < c; r++) {
          double val = fabs(A(r, c));
          if (c != r && val > max) {
            *k = c;
            *l = r;
            max = val;
          }
        }
      }

      return max;
    }

    void step_A(int k, int l, double c, double s) {
      double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
      a_kk = A(k, k);
      a_ll = A(l, l);

      A(k, k) = c*c*a_kk - 2*c*s*A(k, l) + s*s*a_ll;
      A(l, l) = s*s*a_kk + 2*c*s*A(k, l) + c*c*a_ll;
      A(k, l) = 0;
      A(l, k) = 0;

      for (int i = 0; i < n; i++) {
        if (i != k && i != l) {
          a_ik = A(i, k);
          a_il = A(i, l);
          A(i, k) = A(k, i) = c*a_ik - s*a_il;
          A(i, l) = A(l, i) = c*a_il + s*a_ik;
        }

        /*
        TODO: These eigenvectors doesn't seem to be correct?
        r_ik = R(i, k);
        r_il = R(i, l);
        R(i, k) = c*r_ik - s*r_il;
        R(i, l) = c*r_il + s*r_ik;
        */
      }
    }

    void theta(int k, int l, double *c, double *s) {
      double tau = (A(l, l) - A(k, k))/(2*A(k, l));
      double t;

      if (tau > 0) {
        t = -tau + sqrt(1 + tau*tau);
      } else {
        t = -tau - sqrt(1 + tau*tau);
      }

      *c = 1/sqrt(1+t*t);
      *s = *c * t;
    }

    void rotate(int k, int l) {
      double c, s;
      theta(k, l, &c, &s);
      step_A(k, l, c, s);
    }
};

