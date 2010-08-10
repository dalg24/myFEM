#include "../include/Polynomial.h"

Polynomial::
Polynomial() : dim(0),
	       degree(0),
	       ncoeffs(0) {
  coeffs = NULL;
}

Polynomial::
Polynomial(const int n) : dim(n),
                          degree(0),
                          ncoeffs(0) {
  coeffs = NULL;
}

Polynomial::
Polynomial(const int n,
           const int p) : dim(n),
                          degree(p) {
  ncoeffs = factorial(dim+degree)/(factorial(dim)*factorial(degree));
  coeffs = new double[ncoeffs];
  for (int i = 0; i < ncoeffs; i++) {
    coeffs[i] = 0.;
  }
}                          

Polynomial::
Polynomial(const Polynomial& poly) : dim(poly.dim),
                                     degree(poly.degree),
				     ncoeffs(poly.ncoeffs) {
  coeffs = new double[ncoeffs];
  for (int i=0; i < ncoeffs; i++)
    coeffs[i] = poly.coeffs[i];
}

Polynomial::
~Polynomial() {
  delete [] coeffs;
  coeffs = NULL;
}

Polynomial&
Polynomial::
operator=(const Polynomial& poly) {
  if (this == &poly) return *this;

  dim = poly.dim;
  degree = poly.degree;
  ncoeffs = poly.ncoeffs;

  if (coeffs != NULL) delete [] coeffs;
  coeffs = new double[ncoeffs];

  for (int i=0; i < ncoeffs; i++)
    coeffs[i] = poly.coeffs[i];

  return *this;
}

double
Polynomial::
get_value(const double x) const {
  cout << "aaa";
  return 0.;
}

double
Polynomial::
get_dx_value(const double x) const {
  cout << "aaa";
  return 0.;
}

double
Polynomial::
get_dxx_value(const double x) const {
  cout << "aaa";
  return 0.;
}

double
Polynomial::
get_value(const double x,
          const double y) const {
  cout << "bbb";
  return 0.;
}

double
Polynomial::
get_dx_value(const double x,
             const double y) const {
  cout << "bbb";
  return 0.;
}

double
Polynomial::
get_dy_value(const double x,
             const double y) const {
  cout << "bbb";
  return 0.;
}

double
Polynomial::
get_dxx_value(const double x,
              const double y) const {
  cout << "bbb";
  return 0.;
}

double
Polynomial::
get_dxy_value(const double x,
              const double y) const {
  cout << "bbb";
  return 0.;
}

double
Polynomial::
get_dyy_value(const double x,
              const double y) const {
  cout << "bbb";
  return 0.;
}

void
Polynomial::
info() const {
  cout << "reste a faire" << endl;
}
