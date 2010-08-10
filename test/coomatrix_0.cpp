#include <iostream>

using namespace std;

#include "../include/COOMatrix.h"
#include "../include/Vector.h"

int main(int argc, char* argv[]) {

  COOMatrix a_0;
  cout << "initial values after default construction COOMatrix() /" << endl;
  cout << a_0;

  int r[] = {0,0,1,1,1};
  int c[] = {1,2,0,1,2};
  double v[] = {2,7,-0.4,1,-3};
  vector<int> row(r,r+sizeof(r)/sizeof(int));
  vector<int> col(c,c+sizeof(c)/sizeof(int));
  vector<double> val(v,v+sizeof(v)/sizeof(double));
  COOMatrix a_1(2,3,5,row,col,val);
  cout << "initial values after construction COOMatrix(2,3,5,{0,0,1,1,1},{1,2,0,1,2},{2,7,-0.4,1,-3}) /" << endl;
  cout << a_1;
  a_1.print_dense();
  a_1.print_sparse();  

  COOMatrix a_2(4,5);
  a_2.info();
  a_2.add_val(1,0,-3.1);
  a_2.add_val(1,3,17);
  a_2.add_val(2,5,2);
  a_2.add_val(0,0,0);
  a_2.add_val(2,1,0.7);
  a_2.add_val(3,2,-9);
  a_2.print_sparse();
  a_2.print_dense();
  a_2.info();

  a_0 = a_2;
  a_0.info();

  cout << "hello world!" << endl;

}
