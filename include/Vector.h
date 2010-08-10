#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

#include "COOMatrix.h"

/** 
* @class Vector
*
* @brief This is a class for column vector intended to interface with sparse matrices in coordinate format.
*
*/

class Vector {

public:
  /// Default constructor.
  Vector(const vector<double> = vector<double>(0));

  /// Constructor.
  Vector(const int);

  /// Copy constructor.
  Vector(const Vector&);

  /// Destructor.
  ~Vector();

  ///
  Vector& operator=(const Vector&);

  /// This is a function that sets all values of the vector to zero.
  void set_zeros();

  /// This is a function that sets the ith component of the vector to a given value.
  int set_val(const int,
              const double);

  /// This is a function that adds a given value to the ith component of the vector.
  int add_val(const int,
              const double);

  /// This is a function that return the size of the vector.
  int get_size() const;

  /// This is a function that returns the number of rows of the vector (is equal to its size).
  int get_nrow() const;

  /// This is a function that returns the number of columns of the vector (is equal to 1).
  int get_ncol() const;

  /// This is a function that returns the ith component of the vector.
  double get_val(const int) const;

  /// This is a function that prints to the screen the vector.
  void print() const;

  /// This is a function that prints to the screen some informations about the vector.
  void info() const;

  friend ostream& operator<<(ostream&,
                             const Vector&);

  ///
  friend bool operator==(const Vector&,
			 const Vector&);

  /// This is a function that multiply vector by a sparse matrix in coordinate format.
  friend Vector operator*(const COOMatrix&,
                          const Vector&);

  /// This is only for testing purpose.  You might want to use operator* which is much faster.
  friend Vector multiply(const COOMatrix&,
                         const Vector&);

private:
  /// The vector.
  vector<double> val;

  /// Size of the vector.
  int size;
};

#endif // VECTOR_H

