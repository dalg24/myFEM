#include "BasisFunctions.hpp"
//#include "Exceptions.hpp"

#include <iostream>
#include <cassert>

void evaluate_basis_functions(BasisFunctions *basis_functions, unsigned int n_points) {
  const unsigned short n_dof = basis_functions->getNumberOfNodes();
  std::cout<<"==== basis functions "<<basis_functions->getType()<<" of order "<<basis_functions->getOrder()<<" ====\n";
  std::cout<<"associated nodes are\n";
  std::vector<Point> nodes = basis_functions->getNodes();
  for (unsigned int dof = 0; dof < n_dof; ++dof) {
    std::cout<<"Node "<<dof<<": "<<nodes[dof]<<"\n"; 
  }
  std::cout<<"evaluate functions at regularly spaced points\n";
  const std::vector<Point> support_points = basis_functions->getSupportPoints();
  assert(support_points.size() == 2);
  for (unsigned int p = 0; p < n_points; ++p) {
    const Point point = support_points[0]+p/double(n_points-1)*(support_points[1]-support_points[0]);
    std::cout<<"Point "<<p<<": "<<point<<"  ";
    for (unsigned int dof = 0; dof < n_dof; ++dof) {
      std::cout<<"b_"<<dof<<" = "<<basis_functions->getVal(dof, point)<<"  ";
      std::cout<<"DbDx_"<<dof<<" = "<<basis_functions->getDx(dof, point)<<"  ";
    }                  
    std::cout<<"\n";
  }
}

int main() {
  BasisFunctions *basis_functions;
  const unsigned int n_points = 11;

  for (unsigned short i = 1; i < 5; ++i) {
    std::vector<Point> support_points;
    support_points.push_back(Point(-1.0));
    support_points.push_back(Point(+1.0));
    basis_functions = new PiecewisePolynomial(i, support_points);
    evaluate_basis_functions(basis_functions, n_points);
    delete basis_functions;
  }






  return 0;
}
