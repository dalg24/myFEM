#include "Poly1D.h"

Poly1D::
Poly1D() : Polynomial(1) {
}

Poly1D::
Poly1D(const int p) : Polynomial(1, p) {
}

Poly1D::
Poly1D(const Poly1D& poly) : Polynomial(poly.dim,
                                        poly.degree) {
  for (int i=0; i < poly.ncoeffs; i++)
    coeffs[i] = poly.coeffs[i];
}

Poly1D::
~Poly1D() {
  delete [] coeffs;
  coeffs = NULL;
}

Poly1D&
Poly1D::
operator=(const Poly1D& poly) {
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

void
Poly1D::
set_degree(const int p) {
  degree = p;
}

void
Poly1D::
set_coeffs(...) {
  ncoeffs = factorial(dim+degree)/(factorial(dim)*factorial(degree));

  if (coeffs != NULL) delete [] coeffs;
  coeffs = new double[ncoeffs];

  va_list ap;
  va_start(ap, ncoeffs);

  for (int i = 0; i < ncoeffs; i++)
    coeffs[i] = va_arg(ap, double);

  va_end(ap);
}

double
Poly1D::
get_value(const double x) const {
  double result = 0.0;

  for (int j = 0; j < degree+1; j++) {
    result += coeffs[j]*pow(x,j);
  }

  return result;
}

double
Poly1D::
get_dx_value(const double x) const {
  double result = 0.0;

  for (int j = 1; j < degree+1; j++) {
    result += coeffs[j]*pow(x,j-1)*j;
  }

  return result;
}

double
Poly1D::
get_dxx_value(const double x) const {
  double result = 0.0;

  for (int j = 2; j < degree+1; j++) {
    result += coeffs[j]*pow(x,j-2)*(j-1)*j;
  }

  return result;
}

void
Poly1D::
info() const {
  cout << "Poly1D:\n";
  cout << "dim = " << dim << endl;
  cout << "degree = " << degree << endl;
  cout << "ncoeffs = " << ncoeffs << endl;
  if (coeffs==NULL) {cout << "coeffs not set yet!!\n";}
  else {
  cout << coeffs[0];
  for (int p = 1; p < degree+1; p++) {
    cout << " + " << coeffs[p] << " ";
    for (int l = 0; l < p; l++) {
      cout << "* x ";     
    }
  }
  }
  cout << endl;
}
