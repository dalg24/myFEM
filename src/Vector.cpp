#include "../include/Vector.h"

Vector::
Vector(const vector<double> v) : val(v),
                                 size(v.size()) {
}

Vector::
Vector(const int i) : size(i) {
  vector<double> v(i, 0.0);
  val = v;
}

Vector::
Vector(const Vector& o) : val(o.val),
                          size(o.size) {
}

Vector::
~ Vector() {
}

Vector&
Vector::
operator=(const Vector& o) {
  if (this == &o) return *this;

  val = o.val;
  size = o.size;

  return *this;
}

void
Vector::
set_zeros() {
  for (int i = 0; i < size; i++)
    val[i] = 0.0;
}

int
Vector::
set_val(const int i,
        const double v) {
  if (i > size) {
    cout << "Warning in Vector::set_val - index is out of bound" << endl;
    return 0;
  }
  val[i] = v;
  return 1;
}

int
Vector::
add_val(const int i,
        const double v) {
  if (i > size) {
    cout << "Warning in Vector::set_val - index is out of bound" << endl;
    return 0;
  }
  val[i] += v;
  return 1;
}

int
Vector::
get_size() const {
  return size;
}

int
Vector::
get_nrow() const {
  return size;
}

int
Vector::
get_ncol() const {
  return 1;
}

double
Vector::
get_val(const int i) const {
  return val[i];
}

void
Vector::
print() const {
  cout << "Vector: " << endl;
  for (int i = 0; i < size; i++) {
    cout << val[i] << endl;
  }
}

void
Vector::
info() const {
  cout << "Vector: ";
  cout << "(nrow=" << size;
  cout << ",ncol=1)";

  cout << " size=";
  cout <<  size*sizeof(double);
  cout << "bytes\n";
}

Vector
operator*(const COOMatrix& A,
          const Vector& x) {
  if (A.get_ncol() != x.get_size()) {
     cout << "Error - Cannot multiply " << A.get_nrow() << "-by-" << A.get_ncol() << "matrix ";
     cout << "by a " << x.get_size() << "x1 column vector\n";
     Vector b(x.size);
     return b;
  }

  Vector b(A.get_nrow());
  int k = 0; // iterator through triplet of the COOMatrix.
  for (int i = 0; i < A.get_nrow(); i++) {
    b.val[i] = 0.0;
    for (; A.get_row(k) == i; k++) {
	    b.val[i] += A.get_val(k) * x.val[A.get_col(k)];
      if (k == A.get_nnz() - 1) break;
    }
  }

  return b;
}

Vector
multiply(const COOMatrix& A,
         const Vector& x) {
  if (A.get_ncol() != x.get_size()) {
     cout << "Error - Cannot multiply " << A.get_nrow() << "-by-" << A.get_ncol() << "matrix ";
     cout << "by a " << x.get_size() << "x1 column vector\n";
     Vector b(x.size);
     return b;
  }

  Vector b(A.get_nrow());
  for (int i = 0; i < A.get_nrow(); i++) {
    b.val[i] = 0.0;
    for (int j = 0; j < A.get_ncol(); j++) {
      b.val[i] += A.get_val(i, j) * x.val[j];
    }
  }

  return b;
}

bool
operator==(const Vector& v,
           const Vector& w) {
  if (v.size != w.size) return false;

  for (int i = 0; i < v.size; i++)
    if (v.val[i] != w.val[i]) return false;

  return true;
}
