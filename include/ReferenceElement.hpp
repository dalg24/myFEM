#ifndef myFEM_REFERENCE_ELEMENT_HPP
#define myFEM_REFERENCE_ELEMENT_HPP

#include "Point.hpp"
#include "BasisFunctions.hpp"

class ReferenceElement {
public:
  ReferenceElement(unsigned short u) {
    std::vector<Point> sp; 
    sp.push_back(Point(-1.0));
    sp.push_back(Point(1.0));
    bf = new PiecewisePolynomial(u, sp);
  }
  ReferenceElement(const std::string& s = "PiecewiseLinear", unsigned short u = 0) { 
    std::vector<Point> sp; 
    sp.push_back(Point(-1.0));
    sp.push_back(Point(1.0));
    if (s == "PiecewiseLinear") {
      bf = new PiecewiseLinear(sp);
    } else if (s == "PiecewiseQuadratic") {
      bf = new PiecewiseQuadratic(sp);
    } else if (s == "PiecewiseCubic") {
      bf = new PiecewiseCubic(sp);
    } else if (s == "PiecewisePolynomial") {
      assert(u != 0);
      bf = new PiecewisePolynomial(u, sp);
    } else {
      throw myFEM_INVALID_INPUT_STRING_EXCEPTION;
    }
  }
  ~ReferenceElement() { delete bf; }

  unsigned short getNumberOfNodes() const { return bf->getNumberOfNodes(); }
  Point getNode(unsigned short idof) const { return bf->getNode(idof); }
  std::vector<Point> getNodes() const { return bf->getNodes(); }
  std::vector<Point> getSupportPoints() const { return bf->getSupportPoints(); }
  double getVal(unsigned short idof, Point p) const { return bf->getVal(idof, p); }
  double getDx(unsigned short idof, Point p) const { return bf->getDx(idof, p); }
  std::string getTypeOfBasisFunctions() const { return bf->getType(); }
  unsigned short getOrder() const { return bf->getOrder(); }

protected:
  BasisFunctions *bf;
}; // end class ReferenceElement

#endif // myFEM_REFERENCE_ELEMENT_HPP
