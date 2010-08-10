#include <iostream>

using namespace std;

#include "../include/Poly2D.h"

int main(int argc, char* argv[]) {

  Polynomial *p0 = new Poly2D;
  cout << "initial values after default construction polynomial() /" << endl;
  p0->info();

  p0->set_degree(1);
  p0->set_coeffs(3.0, -7.5, 0.0);
  cout << "info after set_degree(1) et set_coeffs(3.0, -7.5, 0.0) /" << endl;
  p0->info();

  cout << "polynomial value for x = 4.0 / " << p0->get_value(4.0, 1.0) << endl;
  cout << "polynomial value for x = 0.0 / " << p0->get_value(0.0, 2.0) << endl;
  cout << "polynomial value for x = -1.0 / " << p0->get_value(-1.0, 0.0) << endl;
  cout << "polynomial value for x = -10.0 / " << p0->get_value(-10.0, 32.0) << endl;

  cout << "dx value for x = 0.0, y = -0.7 / " << p0->get_dx_value(0.0, -0.7) << endl;
  cout << "dx value for x = -1.0, y = 0.0 / " << p0->get_dx_value(-1.0, 0.0) << endl;
  cout << "dx value for x = -8.5, y = 32.9 / " << p0->get_dx_value(-8.5, 32.9) << endl;
  cout << "dx value for x = 2.4, y = 18.2 / " << p0->get_dx_value(2.4, 18.2) << endl;

  cout << "dy value for x = 0.0, y = -0.7 / " << p0->get_dy_value(0.0, -0.7) << endl;
  cout << "dy value for x = -1.0, y = 0.0 / " << p0->get_dy_value(-1.0, 0.0) << endl;
  cout << "dy value for x = -8.5, y = 32.9 / " << p0->get_dy_value(-8.5, 32.9) << endl;
  cout << "dy value for x = 2.4, y = 18.2 / " << p0->get_dy_value(2.4, 18.2) << endl;

  Polynomial *p1 = new Polynomial(*p0);
  p0->set_degree(2);
  p0->set_coeffs(0.0, -0.5, 1.0, 2.0, 3.0, 4.0);

  p1->info();
  p0->info();

  cout << "polynomial value for x = 4.0, y = 0.0 / " << p0->get_value(4.0, 0.0) << endl;
  cout << "polynomial value for x = 0.0, y = 0.0 / " << p0->get_value(0.0, 0.0) << endl;
  cout << "polynomial value for x = 0.0, y = -2.0 / " << p0->get_value(0.0, -2.0) << endl;
  cout << "polynomial value for x = -10.0, y = 1.0 / " << p0->get_value(-10.0, 1.0) << endl;

  cout << "dx value for x = 0.0, y = 3.0 / " << p0->get_dx_value(0.0, 3.0) << endl;
  cout << "dx value for x = 0.0, y = 0.0 / " << p0->get_dx_value(0.0, 0.0) << endl;
  cout << "dx value for x = -8.0, y = 0.0 / " << p0->get_dx_value(-8.0, 0.0) << endl;
  cout << "dx value for x = 2.1, y = 1.0 / " << p0->get_dx_value(2.1, 1.0) << endl;

  cout << "dy value for x = 0.0, y = 3.0 / " << p0->get_dy_value(0.0, 3.0) << endl;
  cout << "dy value for x = 0.0, y = 0.0 / " << p0->get_dy_value(0.0, 0.0) << endl;
  cout << "dy value for x = -8.0, y = 0.0 / " << p0->get_dy_value(-8.0, 0.0) << endl;
  cout << "dy value for x = 2.1, y = 1.0 / " << p0->get_dy_value(2.1, 1.0) << endl;

  cout << "dxx value for x = 0.0, y = 3.0 / " << p0->get_dxx_value(0.0, 3.0) << endl;
  cout << "dxx value for x = 0.0, y = 0.0 / " << p0->get_dxx_value(0.0, 0.0) << endl;
  cout << "dxx value for x = -8.0, y = 0.0 / " << p0->get_dxx_value(-8.0, 0.0) << endl;
  cout << "dxx value for x = 2.1, y = 1.0 / " << p0->get_dxx_value(2.1, 1.0) << endl;

  cout << "dxy value for x = 0.0, y = 3.0 / " << p0->get_dxy_value(0.0, 3.0) << endl;
  cout << "dxy value for x = 0.0, y = 0.0 / " << p0->get_dxy_value(0.0, 0.0) << endl;
  cout << "dxy value for x = -8.0, y = 0.0 / " << p0->get_dxy_value(-8.0, 0.0) << endl;
  cout << "dxy value for x = 2.1, y = 1.0 / " << p0->get_dxy_value(2.1, 1.0) << endl;

  cout << "dyy value for x = 0.0, y = 3.0 / " << p0->get_dyy_value(0.0, 3.0) << endl;
  cout << "dyy value for x = 0.0, y = 0.0 / " << p0->get_dyy_value(0.0, 0.0) << endl;
  cout << "dyy value for x = -8.0, y = 0.0 / " << p0->get_dyy_value(-8.0, 0.0) << endl;
  cout << "dyy value for x = 2.1, y = 1.0 / " << p0->get_dyy_value(2.1, 1.0) << endl;

/*
  cout << "dxx value for x = -1.0 / " << p0.get_dxx_value(-1.0) << endl;
  cout << "dxx value for x = 32.8 / " << p0.get_dxx_value(32.8) << endl;
  cout << "dxx value for x = 0.0 / " << p0.get_dxx_value(0.0) << endl;


  Poly1D p1(p0);
  cout << "p1.info() after copy construction p1=polynomial(p0) /" << endl;
  p1.info();

  p1.set_coeffs(-2.3, 36.15);
  p1.info();

  p0.info();

  Poly1D p2;
  p2 = p0;
  p2.set_degree(3);
  p2.set_coeffs(1.6, 7.0, 0.0, -2.0);
  cout << "p2.info() after p2.set_degree(3) and p2.set_coeffs(1.6, 0.0, -2.0) /" << endl;
  p2.info();

  cout << "polynomial value for x = 1.0 / " << p2.get_value(1.0) << endl;
  cout << "polynomial value for x = 0.0 / " << p2.get_value(0.0) << endl;
  cout << "polynomial value for x = -1.0 / " << p2.get_value(-1.0) << endl;
  cout << "polynomial value for x = -10.0 / " << p2.get_value(-10.0) << endl;

  cout << "dx value for x = 0.0 / " << p2.get_dx_value(0.0) << endl;
  cout << "dx value for x = -1.0 / " << p2.get_dx_value(-1.0) << endl;
  cout << "dx value for x = -8.5 / " << p2.get_dx_value(-8.5) << endl;
  cout << "dx value for x = 2.4 / " << p2.get_dx_value(2.4) << endl;

  cout << "dxx value for x = -1.0 / " << p2.get_dxx_value(-1.0) << endl;
  cout << "dxx value for x = 32.8 / " << p2.get_dxx_value(32.8) << endl;
  cout << "dxx value for x = 0.0 / " << p2.get_dxx_value(0.0) << endl;
*/
}
