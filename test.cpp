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

  unsigned int getNumberOfQuadraturePoints() const { assert(points.size()==weights.size()); return points.size(); }
  Point getQuadraturePoint(unsigned int iqp) const { return points[iqp]; }
  std::vector<Point> getQuadraturePoints() const { return points; }
  double getWeight(unsigned int iqp) const { return weights[iqp]; }
  std::vector<double> getWeights() const { return weights; }
  std::vector<Point> getSupportPoints() const { return sp; }

protected:
  std::vector<Point> points;
  std::vector<double> weights;
  std::vector<Point> sp;
}; // end class QuadratureRule

class GaussianTwoPoints : public QuadratureRule {
public:
  GaussianTwoPoints() {
    points.push_back(Point(-1.0 / sqrt(3.0))); weights.push_back(1.0);
    points.push_back(Point(1.0 / sqrt(3.0))); weights.push_back(1.0);
    sp.push_back(Point(-1.0)); sp.push_back(Point(1.0));
  }
}; // end class GaussianTwoPoints

class GaussianThreePoints : public QuadratureRule {
public:
  GaussianThreePoints() {
    points.push_back(Point(-sqrt(3.0 / 5.0))); weights.push_back(5.0 / 9.0);
    points.push_back(Point(0.0)); weights.push_back(8.0 / 9.0);
    points.push_back(Point(sqrt(3.0 / 5.0))); weights.push_back(5.0 / 9.0);
    sp.push_back(Point(-1.0)); sp.push_back(Point(1.0));
  }
}; // end class GaussianThreePoints

////////////////////////// SHAPE FUNCTIONS //////////////////////////////////////
class BasisFunctions {
public:
  BasisFunctions() { }
  ~BasisFunctions() { }
  BasisFunctions(const BasisFunctions& bf) : nodes(bf.nodes) { }
  BasisFunctions& operator=(const BasisFunctions& bf) { if (this==&bf) return *this; nodes = bf.nodes; return *this; }


  unsigned int getNumberOfNodes() const { return nodes.size(); }
  Point getNode(unsigned int idof) const { return nodes[idof]; }
  std::vector<Point> getNodes() const { return nodes; }

  std::vector<Point> getSupportPoints() const { std::vector<Point> sp; sp.push_back(nodes[0]); sp.push_back(nodes[1]); return sp;}

  virtual double getVal(unsigned int, Point) const = 0;
  virtual double getDx(unsigned int, Point) const = 0;

protected:
  std::vector<Point> nodes;
}; // end class BasisFunctions

class PiecewiseLinear : public BasisFunctions {
public:             
  PiecewiseLinear(std::vector<Point> sp) { nodes = sp; }

  double getVal(unsigned int idof, Point p) const {
    assert(nodes[0].x <= p.x); assert(p.x <= nodes[1].x);
    double value;
    if (idof == 0) {
      value = (nodes[1].x - p.x) / (nodes[1].x - nodes[0].x); 
    } else if (idof == 1) {
      value = (p.x - nodes[0].x) / (nodes[1].x - nodes[0].x); 
    } else {
      std::cerr<<"aie"<<std::endl;
      abort();
    }
    return value;
  }

  double getDx(unsigned int idof, Point p) const {
    assert(nodes[0].x <= p.x); assert(p.x <= nodes[1].x);
    double value;
    if (idof == 0) {
      value = -1.0 / (nodes[1].x - nodes[0].x); 
    } else if (idof == 1) {
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

  double getVal(unsigned int idof, Point p) const {
    assert(nodes[0].x <= p.x); assert(p.x <= nodes[1].x);
    double value;
    if (idof == 0) {
      value = (nodes[1].x - p.x) * (nodes[2].x - p.x) / ((nodes[1].x - nodes[0].x) * (nodes[2].x - nodes[0].x)); 
    } else if (idof == 1) {
      value = (nodes[2].x - p.x) * (nodes[0].x - p.x) / ((nodes[2].x - nodes[1].x) * (nodes[0].x - nodes[1].x)); 
    } else if (idof == 2) {
      value = (nodes[0].x - p.x) * (nodes[1].x - p.x) / ((nodes[0].x - nodes[2].x) * (nodes[1].x - nodes[2].x)); 
    } else {
      std::cerr<<"aie"<<std::endl;
      abort();
    }
    return value;
  }

  double getDx(unsigned int idof, Point p) const {
    assert(nodes[0].x <= p.x); assert(p.x <= nodes[1].x);
    double value;
    if (idof == 0) {
      value = (2.0 * p.x - nodes[1].x - nodes[2].x) / ((nodes[1].x - nodes[0].x) * (nodes[2].x - nodes[0].x)); 
    } else if (idof == 1) {
      value = (2.0 * p.x - nodes[2].x - nodes[0].x) / ((nodes[2].x - nodes[1].x) * (nodes[0].x - nodes[1].x)); 
    } else if (idof == 2) {
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

  unsigned int getNumberOfNodes() const { return bf->getNumberOfNodes(); }
  Point getNode(unsigned int idof) const { return bf->getNode(idof); }
  std::vector<Point> getNodes() const { return bf->getNodes(); }
  std::vector<Point> getSupportPoints() const { return sp; }
  double getVal(unsigned int idof, Point p) const { return bf->getVal(idof, p); }
  double getDx(unsigned int idof, Point p) const { return bf->getDx(idof, p); }

protected:
  std::vector<Point> sp; 
  BasisFunctions *bf;
}; // end class ReferenceElement

////////////////////////// FINITE ELEMENT //////////////////////////////////////
class FiniteElement {
public:
  FiniteElement(const std::vector<Point>& t, ReferenceElement *b) : sp(t), re(b) { 
    if (sp.size() != b->getSupportPoints().size()) {
      std::cerr<<"pas bien"<<std::endl;
      abort();
    }
  }

  Point mapGlobalToLocal(Point g) const {
    std::vector<Point> ref_sp = re->getSupportPoints();
    assert(sp[0].x <= g.x); assert(g.x <= sp[1].x);
    Point l;
    l.x = ref_sp[0].x + (ref_sp[1].x - ref_sp[0].x) * (g.x - sp[0].x) / (sp[1].x - sp[0].x);
    return l;
  }
  Point mapLocalToGlobal(Point l) const {
    std::vector<Point> ref_sp = re->getSupportPoints();
    assert(ref_sp[0].x <= l.x); assert(l.x <= ref_sp[1].x);
    Point g;
    g.x = sp[0].x + (sp[1].x - sp[0].x) * (l.x - ref_sp[0].x) / (ref_sp[1].x - ref_sp[0].x);
    return g;
  }
  double getDeterminantOfJacobian() const { std::vector<Point> ref_sp = re->getSupportPoints(); return (sp[1].x - sp[0].x) / (ref_sp[1].x - ref_sp[0].x); }
  double getValByLocal(unsigned int idof, Point l) const { return re->getVal(idof, l); }
  double getValByGlobal(unsigned int idof, Point g) const { Point l = mapGlobalToLocal(g); return re->getVal(idof, l); }
  double getDxByLocal(unsigned int idof, Point l) const { return re->getDx(idof, l); }
  double getDxByGlobal(unsigned int idof, Point g) const { Point l = mapGlobalToLocal(g); return re->getDx(idof, l); }
  unsigned int getNumberOfNodes() const { return re->getNumberOfNodes(); }

protected:
  std::vector<Point> sp;
  ReferenceElement *re;
}; // end class FiniteElement

////////////////////////// FE Values //////////////////////////////////////
class UpdateFlags { };
class FEValues {
public:
  FEValues(FiniteElement *e, QuadratureRule *r, UpdateFlags f) : fe(e), qr(r), uf(f) { } 

  void reinit(FiniteElement *e) {
    fe = e;
    DeterminantOfJacobian = fe->getDeterminantOfJacobian();
    Weights = qr->getWeights();
    const unsigned int nqp = qr->getNumberOfQuadraturePoints();
    const unsigned int ndof = fe->getNumberOfNodes();
    ShapeValues.resize(ndof);
    ShapeDx.resize(ndof);
    for (unsigned int idof = 0; idof < ndof; ++idof) {
      ShapeValues[idof].resize(nqp);
      ShapeDx[idof].resize(nqp);
      for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
        ShapeValues[idof][iqp] = fe->getValByLocal(idof, qr->getQuadraturePoint(iqp));
        ShapeDx[idof][iqp] = fe->getDxByLocal(idof, qr->getQuadraturePoint(iqp));
      } 
    }

  }

  double getDeterminantOfJacobianTimesWeight(unsigned int iqp) const { return Weights[iqp]*DeterminantOfJacobian; }
  double getShapeValue(unsigned int idof, unsigned int iqp) const { return ShapeValues[idof][iqp]; }
  double getShapeDx(unsigned int idof, unsigned int iqp) const { return ShapeDx[idof][iqp]; }
  std::vector<Point> getQuadraturePoints() const {
    std::vector<Point> qps;
    const unsigned int nqp = qr->getNumberOfQuadraturePoints();
    for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
      qps.push_back(fe->mapLocalToGlobal(qr->getQuadraturePoint(iqp)));
    }
    return qps;
  }

protected:
  std::vector<double> Weights;
  double DeterminantOfJacobian;
  std::vector<std::vector<double> > ShapeValues;
  std::vector<std::vector<double> > ShapeDx;

  FiniteElement *fe;
  QuadratureRule *qr;
  UpdateFlags uf;
}; // end class FEValues

int main(int argc, char *argv[]) {
  ReferenceElement *referenceElement = new ReferenceElement;
  std::vector<Point> supportPoints;  supportPoints.push_back(Point(0.0));  supportPoints.push_back(Point(4.0)); 
  FiniteElement *finiteElement = new FiniteElement(supportPoints, referenceElement);
  std::cout<<finiteElement->mapGlobalToLocal(Point(1.0))<<finiteElement->mapLocalToGlobal(Point(0.0))<<std::endl;
  QuadratureRule *quadratureRule = new GaussianTwoPoints;
  UpdateFlags dummyUpdateFlags;
  FEValues *feValues = new FEValues(finiteElement, quadratureRule, dummyUpdateFlags);
  feValues->reinit(finiteElement);
  const unsigned int ndof = finiteElement->getNumberOfNodes();
  const unsigned int nqp = quadratureRule->getNumberOfQuadraturePoints();
  for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
    std::cout<<"qp="<<feValues->getQuadraturePoints()[iqp]<<"  "
             <<"JxW="<<feValues->getDeterminantOfJacobianTimesWeight(iqp)<<"  ";
    for (unsigned int idof = 0; idof < ndof; ++idof) {
      std::cout<<"phi_"<<idof<<"="<<feValues->getShapeValue(idof, iqp)<<"  "
               <<"DphiDx_"<<idof<<"="<<feValues->getShapeDx(idof, iqp)<<"  ";
    }
    std::cout<<"\n";
  }

  delete feValues;        
  delete quadratureRule;
  delete finiteElement;
  delete referenceElement;

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
  unsigned int n_qp = qr->getNumberOfQuadraturePoints();
  std::vector<Point> qp = qr->getQuadraturePoints();

  std::cout<<"points regulary spaced\n";

  qp.clear();
  n_qp = 51;
  for (unsigned int i = 0; i < n_qp; ++i)
    qp.push_back(Point(-1.0+i*2.0/double(n_qp-1)));
  for (unsigned int i = 0; i < n_qp; ++i) {
    std::cout<<"x="<<qp[i].x<<"  ";
    for (unsigned int j = 0; j < bf->getNumberOfNodes(); ++j) {
      std::cout<<"phi_"<<j<<"="<<bf->getVal(j, qp[j])<<"  ";
    }
    std::cout<<"\n";
  }

  delete bf;
  delete qr;
  }
  

  return 0;
}
