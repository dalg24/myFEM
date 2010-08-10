#include "../include/COOMatrix.h"

COOMatrix::
COOMatrix(const int n_row,
	   const int n_col,
	   const int nz,
	   const vector<int> Ti,
	   const vector<int> Tj,
	   const vector<double> Tx) : nnz(nz),
				       nrow(n_row),
				       ncol(n_col),
				       row(Ti),
				       col(Tj),
				       val(Tx){
}

COOMatrix::
COOMatrix(const COOMatrix& m) : nnz(m.nnz),
				  nrow(m.nrow),
				  ncol(m.ncol) {
  row = m.row;
  col = m.col;
  val = m.val; 
}

COOMatrix::
~COOMatrix() {
}

COOMatrix&
COOMatrix::
operator=(const COOMatrix& m) {
  if(this==&m) return *this;

  nnz = m.nnz;
  nrow = m.nrow;
  ncol = m.ncol;
  row = m.row;
  col = m.col;
  val = m.val;

  return *this;
}

int
COOMatrix::
get_nnz() const {
  return nnz;
}

int
COOMatrix::
get_nrow() const {
  return nrow;
}

int
COOMatrix::
get_ncol() const {
  return ncol;
}

double
COOMatrix::
get_val(const int i) const {
  return val.at(i);
}

int
COOMatrix::
get_row(const int i) const {
  return row.at(i);
}

int
COOMatrix::
get_col(const int i) const {
  return col.at(i);
}

int
COOMatrix::
add_val(const int irow,
	const int jcol,
	const double v) {
  vector<int>::iterator row_it = row.begin();
  vector<int>::iterator col_it = col.begin();
  vector<double>::iterator val_it = val.begin();
  double zero_numeric = 1.0e-12;
  if (abs(v) < zero_numeric) {
    cout << "Warning in COOMatrix::add_val() - added value is a zero numeric" << endl;
    return 0;
  }

  if (irow>=nrow || jcol>=ncol) {
    cout << "Error in COOMatrix::add_val() - row or column index is out of bound" << endl;
    return 0;
  }
  else {
    while (row_it < row.end()) {
      if (*row_it > irow) {
	row.insert(row_it, irow);
	col.insert(col_it, jcol);
	val.insert(val_it, v);
	nnz++;
	return 1;
      }
      else if (*row_it == irow) {
	while (*row_it == irow) {
	  if (*col_it > jcol) {
	    row.insert(row_it, irow);
	    col.insert(col_it, jcol);
	    val.insert(val_it, v);
	    nnz++;
	    return 1;
	  }
	  else if (*col_it == jcol) {
	    if(abs(*val_it + v) < zero_numeric) {
	      row.erase(row_it);
	      col.erase(col_it);
	      val.erase(val_it);
	      nnz--;
	    }
	    else {   
	      *val_it += v;
	    }
	    return 1;
	  }
	  else {
	    row_it++;
	    col_it++;
	    val_it++; 
	  }
	} 
      }
      else {
	row_it++;
	col_it++;
	val_it++;
      }                           
    }                                 
  }
    row.push_back(irow);
    col.push_back(jcol);
    val.push_back(v);
    nnz++;
    return 1;
}

double
COOMatrix::
get_val(const int irow,
	const int jcol) const {
  vector<int>::const_iterator row_it = row.begin();
  vector<int>::const_iterator col_it = col.begin();
  vector<double>::const_iterator val_it = val.begin();
  if (irow >= nrow || jcol >= ncol) {
    cout << "Error in COOMatrix::add_val() - row or column index is out of bound" << endl;
    return 0.0;
  }   
  else {
    while (row_it < row.end()) {
      if (*row_it == irow) {
	while (*row_it == irow && col_it < col.end()) {
	  if (*col_it == jcol) {
	    return *val_it;
	  }
	  else {
	    row_it++; col_it++; val_it++;
	  }
	}
      }
      else {
	row_it++; col_it++; val_it++;
       }
    }
  }
    return 0.0;
}

void
COOMatrix::
print_dense() const {
  cout<<"print dense matrix:\n";
  for(int i = 0; i < nrow; i++) {
    for(int j = 0; j < ncol; j++) {
      cout<<get_val(i,j)<<"\t";
    }
    cout <<endl;
  }
}
		
void
COOMatrix::
print_sparse() const {
  cout<<"print sparse matrix:\n";
  cout<<"values:\n";
  for (int k =0; k < nnz; k++)
    cout<<val[k]<<"\t";
  cout<<endl;
  cout<<"rows:\n";
  for (int k = 0; k < nnz; k++)
    cout<<row[k]<<"\t";
  cout<<endl;
  cout<<"colulmns:\n";
  for(int k = 0; k< nnz; k++)
    cout<<col[k]<<"\t";
  cout<<endl;
}

void
COOMatrix::
info() const {
  cout << "COOMatrix: ";

  cout << "(nnz=" << nnz;
  cout << ",nrow=" << nrow;
  cout << ",ncol=" << ncol<<")";

  cout << " size=";
  cout << nnz*(sizeof(double)+2*sizeof(int));
  cout << "bytes\n";
}


ostream&
operator<<(ostream& os,
           const COOMatrix& m) {
  os << "COOMatrix: ";

  os << "(nnz=" << m.nnz;
  os << ",nrow=" << m.nrow;
  os << ",ncol=" << m.ncol << ")";

  os << " size=";
  os << m.nnz*(sizeof(double)+2*sizeof(int));
  os << "bytes\n";

  return os;
}

bool
operator==(const COOMatrix& a,
	   const COOMatrix& b) {
  if (a.nnz != b.nnz ||
      a.nrow != b.nrow ||
      a.ncol != b.ncol)
    //      a.row != b.row ||
      //a.col != b.col ||
    //      a.val != b.val)
    return false;
  for (int i = 0; i < a.nnz; i++)
    if (a.row[i] != b.row[i] ||
	a.col[i] != b.col[i] ||
	a.val[i] != b.val[i]) {
      return false;
    }

  return true;
}

