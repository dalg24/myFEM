#ifndef COOMATRIX_H
#define COOMATRIX_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

/** 
* @class COOMatrix
*
* @brief This is a class for sparse matrices in COOrdinate format.
*
* Only non-zero elements are stored, and the coordinate of each non-zero element
* are given explicitly.
* Format is specified by three vectors: val, row, col, the size of the matrix: nrow,
* ncol, and a parameter nnz.  All three vectors have the dimension nnz.
*
* nnz: the number of non-zero elements in the matrix.
* val: a vector of double that contains the values for the non-zero elements of
*      the matrix
* row: element i of the vector of int gives the row in the matrix that contains
*      the ith value in the vector val.
* col: element i of the vector of int gives the column in the matrix that contains
*      the ith value in the vector val.
*
* Note that the triplets (row[k],col[k],val[k]) for 0 < k <= nnz are ordered such 
* that the row indices are in ascending order (ie row[k+1]>=row[k] for 0 < k < nnz)
* and that duplicate entries (row[k],col[k]) = (row[l],col[l]) for 0 < k,l <= nnz 
* with k != l) may not exist. 
*
*/

class COOMatrix {

public:
  /// Constructor.
  COOMatrix(const int n_row = 0,
	     const int n_col = 0,
	     const int nz = 0,
	     const vector<int> Ti = vector<int>(0),
	     const vector<int> Tj = vector<int>(0),
	     const vector<double> Tx = vector<double>(0));

  /// Copy constructor.
  COOMatrix(const COOMatrix&);

  /// Destructor.
  ~COOMatrix();

  ///
  COOMatrix& operator=(const COOMatrix&);

  /// This is a function that returns the number of non-zero elements in the matrix.
  int get_nnz() const;

  /// This is a function that returns the number of rows in the matrix.
  int get_nrow() const;

  /// This is a function that returns the number of columns in the matrix.
  int get_ncol() const;

  /// This is a function that returns the ith value in the vector that contains the non-zero elements of the matrix.
  double get_val(const int) const;

  /// This is a function that returns the row in the matrix of the ith value of the vector with the non-zero elements.
  int get_row(const int) const;

  /// This is a function that returns the row in the matrix of the ith value of the vector with the non-zero elements.
  int get_col(const int) const;

  /// This is a function that returns the value in the dense matrix at row i and column j.
  double get_val(const int,
		 const int) const;

  /// This is a function that adds a value in the dense matrix at row i and column j.
  int add_val(const int,
	      const int,
	      const double);

  /// This is a function to cout COOMatrix as a dense matrix.
  void print_dense() const;

  /// This is a function to cout COOMatrix as a sparse matrix (ie cout the three vectors).
  void print_sparse() const;

  /// This is a function that prints to the screen some informations about the COOMatrix.
  void info() const;

  ///
  friend ostream& operator<<(ostream&,
                             const COOMatrix&);

  ///
  friend bool operator==(const COOMatrix&,
			 const COOMatrix&);

private:
  /// Number of non-zero elements in the matrix.
  int nnz; 

  /// Number of rows.
  int nrow;

  /// Number of columns.
  int ncol;

  /// Vector that contains row indices.
  vector<int> row;

  /// Vector that contains the column indices.
  vector<int> col;

  /// The values of the non-zero elements.
  vector<double> val;
};

#endif // COOMATRIX_H

