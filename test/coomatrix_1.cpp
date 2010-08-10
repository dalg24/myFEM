#include <iostream>

using namespace std;

#include "../include/COOMatrix.h"
#include "../include/Vector.h"

int main(int argc, char* argv[]) {

  Vector x_0;
  cout << "initial values after default construction Vector() /" << endl;
  x_0.print();
  x_0.info();

  double v[] = {2.0, 7.0, -0.4, 1.0, -3.0};
  vector<double> val(v,v+sizeof(v)/sizeof(double));
  Vector x_1(val);
  cout << "initial values after construction Vector({2,7,-0.4,1,-3}) /" << endl;
  x_1.info();
  x_1.print();

  Vector x_2(9);
  cout << "initial values after construction Vector(9) /" << endl;
  x_2.print();
  x_2.info();


  Vector x_3(x_1);
  x_3.info();
  x_3.print();
  x_1.set_zeros();

  cout << "x_3.print() then x_0.print() after x_0.set_zeros() /" << endl;
  x_3.print();
  x_1.print();

  cout << "hello world!" << endl;

}
