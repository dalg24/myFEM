#include "BasisFunc.h"

BasisFunc::
BasisFunc(enum element_type t) : type(t) {
  if (t == TRIANGLE_LINEAR) {
    nfunc = 3;
    basis_functions = new Poly2D[3](1);
    basis_functions[0].set_coeffs(1.0, -1.0, -1.0);
    basis_functions[1].set_coeffs(0.0, 1.0, 0.0);
    basis_functions[2].set_coeffs(0.0, 0.0, 1.0);
  }
  if (t == TRIANGLE_QUADRATIC) {
    nfunc = 6;
    basis_functions = new Poly2D[6](2);
    basis_functions[0].set_coeffs(1.0, -3.0, -3.0, 2.0, 4.0, 2.0);
    basis_functions[1].set_coeffs(0.0, -1.0, 0.0, 2.0, 0.0, 0.0);
    basis_functions[2].set_coeffs(0.0, 0.0, -1.0, 0.0, 0.0, 2.0);
    basis_functions[3].set_coeffs(0.0, 0.0, 0.0, 0.0, 4.0, 0.0);
    basis_functions[4].set_coeffs(0.0, 0.0, 4.0, 0.0, -4.0, -4.0);
    basis_functions[5].set_coeffs(0.0, 4.0, 0.0, -4.0, -4.0, 0.0);
  }
  if (t == ONED_LINEAR) {
    nfunc = 2;
    basis_functions = new Poly1D[2](1);
    basis_functions[0].set_coeffs(1.0, -1.0);
    basis_functions[1].set_coeffs(0.0, 1.0);
  }
  if (t == ONED_QUADRATIC) {
    nfunc = 3;
    basis_functions = new Poly1D[3](2);
    basis_functions[0].set_coeffs(1.0, -3.0, 2.0);
    basis_functions[1].set_coeffs(0.0, -1.0, 2.0);
    basis_functions[2].set_coeffs(0.0, 4.0, -4.0);
  }

BasisFunc::
BasisFunc(const BasisFunc& f) : type(f.type),
                                nfunc(f.nfunc) {
  basis_functions = new Polynomial[nfunc];
  for (int i = 0; i < nfunc; i++) {
    basis_functions[i] = f.basis_functions[i];
  }
}

BasisFunc::
~BasisFunc() {
  delete [] basis_functions;
}

BasisFunc&
BasisFunc::
operator=(const BasisFunc& f) {
  if (this = &m) return *this;

  type = f.type;
  nfunc = f.nfunc;
  basis_functions = new Polynomial[nfunc];
  for (int i = 0; i < nfunc; i++) {
    basis_functions[i] = f.basis_functions[i];
  }
}

int
BasisFunc::
get_number_of_basis_functions() const {
  return nfunc;
}

double
BasisFunc::
get_value(const int i,
          const double x) const {
  return basis_function[i]->get_value(x);
}

double
BasisFunc::
get_dx_value(const int i,
             const double x) const {
  return basis_function[i]->get_dx_value(x);
}

double
BasisFunc::
get_value(const int i,
          const double x,
          const double y) const {
  return basis_function[i]->get_value(x,y);
}

double
BasisFunc::
get_dx_value(const int i,
             const double x,
             const double y) const {
  return basis_function[i]->get_dx_value(x,y);
}

double
BasisFunc::
get_dy_value(const int i,
             const double x,
             const double y) const {
  return basis_function[i]->get_dy_value(x,y);
}
