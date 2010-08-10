#include "Poly2D.h"

Poly2D::
Poly2D() : dim(2),
           degree(0),
	   ncoeffs(0) {
  coeffs = NULL;
}

Poly2D::
Poly2D(const Poly2D& poly) : dim(poly.dim),
                             degree(poly.degree),
                             ncoeffs(poly.ncoeffs) {
  coeffs = new double[ncoeffs];
  for (int i=0; i < ncoeffs; i++)
    coeffs[i] = poly.coeffs[i];
}

Poly2D::
~Poly2D() {
  delete [] coeffs;
  coeffs = NULL;
}

Poly2D&
Poly2D::
operator=(const Poly2D& poly) {
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
Poly2D::
set_degree(const int p) {
  degree = p;
}

void
Poly2D::
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
Poly2D::
get_value(const double x,
          const double y) const {
  double result = 0.0;

  int i = 0;
  while (i < degree+1) {
    for (int j = 0; j < i+1; j++) {
      result += coeffs[factorial(dim+i-1)/(factorial(dim)*factorial(i-1))+j]
                * pow(x, i-j) * pow(y, j);
    }
    i++;
  }

  return result;
}

double
Poly2D::
get_dx_value(const double x,
             const double y) const {
  double result = 0.0;

  int i = 1;
  while (i < degree+1) {
    for (int j = 0; j < i; j++) {
      result += coeffs[factorial(dim+i-1)/(factorial(dim)*factorial(i-1))+j]
                * pow(x, i-j-1) * (i-j) * pow(y, j);
    }
    i++;
  }

  return result;
}
 
double
Poly2D::
get_dy_value(const double x,
             const double y) const {
  double result = 0.0;

  int i = 1;
  while (i < degree+1) {
    for (int j = 1; j < i+1; j++) {
      result += coeffs[factorial(dim+i-1)/(factorial(dim)*factorial(i-1))+j]
                * pow(x, i-j) * pow(y, j-1) * j;
    }
    i++;
  }

  return result;
}

double
Poly2D::
get_dxx_value(const double x,
              const double y) const {
  double result = 0.0;

  int i = 2;
  while (i < degree+1) {
    for (int j = 0; j < i-1; j++) {
      result += coeffs[factorial(dim+i-1)/(factorial(dim)*factorial(i-1))+j]
                * pow(x, i-j-2) * (i-j) * (i-j-1) * pow(y, j);
    }

    i++;
  }

  return result;
}

double
Poly2D::
get_dxy_value(const double x,
              const double y) const {
  double result = 0.0;

  int i = 2;
  while (i < degree+1) {
    for (int j = 1; j < i; j++) {
      result += coeffs[factorial(dim+i-1)/(factorial(dim)*factorial(i-1))+j]
                * pow(x, i-j-1) * (i-j) * pow(y, j-1) * j;
    }
    i++;
  }

  return result;
}

double
Poly2D::
get_dyy_value(const double x,
              const double y) const {
  double result = 0.0;
  
  int i = 2;
  while (i < degree+1) {
    for (int j = 2; j < i+1; j++) {
      result += coeffs[factorial(dim+i-1)/(factorial(dim)*factorial(i-1))+j]
                * pow(x, i-j) * pow(y, j-2) * j * (j-1);
    } 
    i++;
  }   
      
  return result;
}

void
Poly2D::
info() const {
  cout << "Poly2D:\n";
  cout << "(dim = " << dim;
  cout << ", degree = " << degree;
  cout << ", ncoeffs = " << ncoeffs << ")\n";
  if (coeffs==NULL) {
    cout << "Warning - coeffecients of the polynomial have not been set\n";
  }
  else {
    cout << "P(x,y) = ";
    int i = 0;
    while (i < degree+1) {
      for (int j = 0; j < i+1; j++) {
        cout << coeffs[factorial(dim+i-1)/(factorial(dim)*factorial(i-1))+j] << " ";
        for (int k = 0; k < i-j; k++) {
          cout << "* x ";
        }
        for (int l = 0; l < j; l++) {
          cout << "* y ";
        }
        cout << " + ";
      }
      i++;
    }
  cout << endl;
  }
}
