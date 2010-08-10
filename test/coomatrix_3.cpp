#include <iostream>
#include <cmath>

using namespace std;

#include "../include/COOMatrix.h"
#include "../include/Vector.h"

int main(int argc, char* argv[]) {

  COOMatrix a_0(1000, 500);
  for (int i = 0; i < 500; i++) {
    a_0.add_val(i, i, 2);
    a_0.add_val(i+500, 499-i, -4.0);
  }
  cout << "a_1 /" << endl;
  a_0.info();
  //a_0.print_dense();

  
  vector<double> val;
  for (int i = 0; i < 500; i++) {
    val.push_back( i * pow(-1, i) ); 
  }

  Vector x_1(val);
  cout << "x_1 after construction Vector({2,7,-0.4,1,-3}) /" << endl;
  x_1.info();
  //x_1.print();

  Vector x_2;
  cout << "x_2 = a_1 * x_1 /" << endl;
  x_2 = a_0 * x_1;
  x_2.info();
  //x_2.print();

  Vector x_3;
  cout << "x_3 = multiply(a_1, x_1) /" << endl;
  x_3 = multiply(a_0, x_1);
  x_3.info();
  //x_3.print();

  cout << "x_2 == x_3 /" << endl;
  cout << (x_3 == x_2) << endl;

  cout << "hello world!" << endl;

}
