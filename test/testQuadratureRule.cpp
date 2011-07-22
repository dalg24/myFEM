#include <iostream>

#include "Point.hpp"
#include "QuadratureRule.hpp"

void print_quadrature_rule(QuadratureRule *quadrature_rule) {
  const unsigned short n_points = quadrature_rule->getNumberOfQuadraturePoints();
  std::cout<<quadrature_rule->getType()<<"\n";
  for (unsigned short i = 0; i < n_points; ++i) {
    std::cout<<"Point "<<quadrature_rule->getQuadraturePoint(i)<<"\t";
    std::cout<<"Weight "<<quadrature_rule->getWeight(i)<<"\n";
  }
  std::cout<<"\n";
}

int main() {

  QuadratureRule *quadrature_rule;

  for (unsigned short i = 2; i < 6; ++i) {
    quadrature_rule = makeNewPointerToQuadratureRule(i); 
    print_quadrature_rule(quadrature_rule);
    delete quadrature_rule;
  }

  std::cout<<2*sizeof(double)+sizeof(unsigned short)+sizeof(std::vector<Point>)<<"\n";

  return 0;
}
