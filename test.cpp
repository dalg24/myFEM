#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cassert>
#include <cmath>

////////////////////////// POINT //////////////////////////////////////
class Point { 
public:
  Point(double X = 0.0) : x(X) { }
  ~Point() { }
  Point(const Point& p) : x(p.x) { }
  Point& operator=(const Point& p) { if (this==&p) return *this; x = p.x; return *this; }
  friend Point operator+(const Point& l, const Point& r) { Point p; p.x = l.x + r.x; return p;}
  friend Point operator-(const Point& l, const Point& r) { Point p; p.x = l.x - r.x; return p;}
  friend Point operator/(const Point& l, const double r) { Point p; p.x = l.x / r; return p;}
  friend Point operator*(const double l, const Point& r) { Point p; p.x = l * r.x; return p;}
  friend std::ostream& operator<<(std::ostream& os, const Point& p) { os<<"( "<<p.x<<" )"; return os; }

  double x; 
}; // end class Point

////////////////////////// QUADRATURE RULE //////////////////////////////////////
class QuadratureRule {
public:
  QuadratureRule() { }
  ~QuadratureRule() { }
  QuadratureRule(const QuadratureRule& qr) : points(qr.points), weights(qr.weights) { }
  QuadratureRule& operator=(const QuadratureRule& qr) { if (this==&qr) return *this; points = qr.points; weights = qr.weights; return *this; }

  unsigned int getNumberQuadraturePoints() const { assert(points.size()==weights.size()); return points.size(); }
  Point getQuadraturePoint(unsigned int i) const { return points[i]; }
  std::vector<Point> getQuadraturePoints() const { return points; }
  double getWeight(unsigned int i) const { return weights[i]; }
  std::vector<double> getWeights() const { return weights; }

protected:
  std::vector<Point> points;
  std::vector<double> weights;
}; // end class QuadratureRule

class GaussianTwoPoints : public QuadratureRule {
public:
  GaussianTwoPoints() {
    points.push_back(Point(-1.0 / sqrt(3.0))); weights.push_back(1.0);
    points.push_back(Point(1.0 / sqrt(3.0))); weights.push_back(1.0);
  }
}; // end class GaussianTwoPoints

class GaussianThreePoints : public QuadratureRule {
public:
  GaussianThreePoints() {
    points.push_back(Point(-sqrt(3.0 / 5.0))); weights.push_back(5.0 / 9.0);
    points.push_back(Point(0.0)); weights.push_back(8.0 / 9.0);
    points.push_back(Point(sqrt(3.0 / 5.0))); weights.push_back(5.0 / 9.0);
  }
}; // end class GaussianThreePoints

////////////////////////// SHAPE FUNCTIONS //////////////////////////////////////
class BasisFunctions {
public:
  BasisFunctions() { }
  ~BasisFunctions() { }
  BasisFunctions(const BasisFunctions& bf) : nodes(bf.nodes) { }
  BasisFunctions& operator=(const BasisFunctions& bf) { if (this==&bf) return *this; nodes = bf.nodes; return *this; }


  unsigned int getNumberNodes() const { return nodes.size(); }
  Point getNode(unsigned int i) const { return nodes[i]; }
  std::vector<Point> getNodes() const { return nodes; }

  virtual double val(unsigned int, Point) const = 0;
  virtual double dx(unsigned int, Point) const = 0;

protected:
  std::vector<Point> nodes;
}; // end class BasisFunctions

class PiecewiseLinear : public BasisFunctions {
public:             
  PiecewiseLinear(std::vector<Point> sp) { nodes = sp; }

  double val(unsigned int i, Point p) const {
    assert(nodes[0].x <= p.x); assert(p.x <= nodes[1].x);
    double values;
    if (i == 0) {
      value = (nodes[1].x - p.x) / (nodes[1].x - nodes[0].x); 
    } else if (i == 1) {
      value = (p.x - nodes[0].x) / (nodes[1].x - nodes[0].x); 
    } else {
      std::cerr<<"aie"<<std::endl;
      abort();
    }
    return value;
  }

  double> dx(unsigned int i, Point p) const {
    assert(nodes[0].x <= p.x); assert(p.x <= nodes[1].x);
    double value;
    if (i == 0) {
      value = -1.0 / (nodes[1].x - nodes[0].x); 
    } else if (i == 1) {
      value = 1.0 / (nodes[1].x - nodes[0].x); 
    } else {
      std::cerr<<"aie"<<std::endl;
      abort();
    }
    return value;
  }
}; // end class PiecewiseLinear

class PiecewiseQuadratic : public BasisFunctions {
public:             
  PiecewiseQuadratic(std::vector<Point> sp) { nodes = sp; nodes.push_back((sp[0] + sp[1]) / 2.0); }

  double val(unsigned int i, Point p) const {
    assert(nodes[0].x <= p.x); assert(p.x <= nodes[1].x);
    double values(3);
    if (i == 0) {
      value = (nodes[1].x - p.x) * (nodes[2].x - p.x) / ((nodes[1].x - nodes[0].x) * (nodes[2].x - nodes[0].x)); 
    } else if (i == 1) {
      value = (nodes[2].x - p.x) * (nodes[0].x - p.x) / ((nodes[2].x - nodes[1].x) * (nodes[0].x - nodes[1].x)); 
    } else if (i == 2) {
      value = (nodes[0].x - p.x) * (nodes[1].x - p.x) / ((nodes[0].x - nodes[2].x) * (nodes[1].x - nodes[2].x)); 
    } else {
      std::cerr<<"aie"<<std::endl;
      abort();
    }
    return value;
  }

  double dx(unsigned int i, Point p) const {
    assert(nodes[0].x <= p.x); assert(p.x <= nodes[1].x);
    double values(3);
    if (i == 0) {
      value = (2.0 * p.x - nodes[1].x - nodes[2].x) / ((nodes[1].x - nodes[0].x) * (nodes[2].x - nodes[0].x)); 
    } else if (i == 1) {
      value = (2.0 * p.x - nodes[2].x - nodes[0].x) / ((nodes[2].x - nodes[1].x) * (nodes[0].x - nodes[1].x)); 
    } else if (i == 2) {
      value = (2.0 * p.x - nodes[0].x - nodes[1].x) / ((nodes[0].x - nodes[2].x) * (nodes[1].x - nodes[2].x)); 
    } else {
      std::cerr<<"aie"<<std::endl;
      abort();
    }
    return value;
  }

}; // end class PiecewiseQuadratic

////////////////////////// REFERENCE ELEMENT //////////////////////////////////////
class ReferenceElement {
public:
  ReferenceElement(std::string s = "PiecewiseLinear") { 
    sp.push_back(Point(-1.0));
    sp.push_back(Point(1.0));
    if (s == "PiecewiseLinear") {
      bf = new PiecewiseLinear(sp);
    } else if (s == "PiecewiseQuadratic") {
      bf = new PiecewiseQuadratic(sp);
    } else {
      std::cerr<<"Error: type of reference element unrecognised"<<std::endl;
      abort();
    }
  }
  ~ReferenceElement() { delete bf; }

protected:
  std::vector<Point> sp; 
  BasisFunctions *bf;
};

////////////////////////// FINITE ELEMENT //////////////////////////////////////
class FiniteElement {
public:

  Point mapGlobalToLocal(Point g) {
    Point l;
    l.x = bf->getNodes()[0].x + (bf->getNodes()[1].x - bf->getNodes()[0].x) * (g.x - supportPoints[0].x) / (supportPoints[1].x - supportPoints[0].x);
    return l;
  }

  Point mapLocalToGlobal(Point l) {

  }

protected:
  std::vector<Point> sp;
};


int main(int argc, char *argv[]) {

  Point *pp;
  pp = new Point(-3.0);
  std::cout<<*pp<<pp->x<<2.0*(*pp)<<*pp/3.0<<*pp+(*pp)<<*pp-(*pp)<<std::endl;


  QuadratureRule *qr;
  BasisFunctions *bf;

  {
  std::vector<Point> dummySupportPoints;
  dummySupportPoints.push_back(Point(-1.0));
  dummySupportPoints.push_back(Point(1.0));
  //qr = new GaussianTwoPoints;
  //bf = new PiecewiseLinear(dummySupportPoints);
  qr = new GaussianThreePoints;
  bf = new PiecewiseQuadratic(dummySupportPoints);
  unsigned int n_qp = qr->getNumberQuadraturePoints();
  std::vector<Point> qp = qr->getQuadraturePoints();

  std::cout<<"points regulary spaced\n";

  qp.clear();
  n_qp = 51;
  for (unsigned int i = 0; i < n_qp; ++i)
    qp.push_back(Point(-1.0+i*2.0/double(n_qp-1)));
  for (unsigned int i = 0; i < n_qp; ++i) {
    std::vector<double> v = bf->val(qp[i]);
    std::cout<<"x="<<qp[i].x<<"  ";
    for (unsigned int j = 0; j < v.size(); ++j) {
      std::cout<<"phi_"<<j<<"="<<v[j]<<"  ";
    }
    std::cout<<"\n";
  }

  delete bf;
  delete qr;
  }
  

  return 0;
}
