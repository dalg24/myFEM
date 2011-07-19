#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include "umfpack.h"

// TODO: templatize dimension
//       DOFHandler class
//       boundary conditions
enum Exception_t { myFEM_DOF_OUT_OF_RANGE_EXCEPTION, myFEM_INVALID_INPUT_STRING_EXCEPTION, myFEM_DIVIDE_BY_ZERO_EXCEPTION };
enum DebugLevel_t { myFEM_NO_DEBUG, myFEM_MIN_DEBUG, myFEM_MED_DEBUG, myFEM_MAX_DEBUG };
enum Object_t { myFEM_BOTH_MATRIX_AND_VECTOR, myFEM_MATRIX_ONLY, myFEM_VECTOR_ONLY };
enum Norm_t { myFEM_L1_NORM, myFEM_L2_NORM, myFEM_H1_NORM };
////////////////////////// POINT //////////////////////////////////////
class Point { 
public:
  Point(double X = 0.0) : x(X) { }
  ~Point() { }
  Point(const Point& p) : x(p.x) { }
  Point& operator=(const Point& p) { if (this==&p) return *this; x = p.x; return *this; }
  friend Point operator+(const Point& l, const Point& r) { Point p; p.x = l.x + r.x; return p;}
  friend Point operator-(const Point& l, const Point& r) { Point p; p.x = l.x - r.x; return p;}
  friend Point operator/(const Point& l, const double r) { if (r == 0) throw myFEM_DIVIDE_BY_ZERO_EXCEPTION; Point p; p.x = l.x / r; return p;}
  friend Point operator*(const double l, const Point& r) { Point p; p.x = l * r.x; return p;}
  friend Point operator*(const Point& l, const double r) { Point p; p.x = l.x * r; return p;}
  friend std::ostream& operator<<(std::ostream& os, const Point& p) { os<<"( "<<p.x<<" )"; return os; }

  double x; 
}; // end class Point
                                            
////////////////////////// CHANTIER //////////////////////////////////////
class SupportPoint { 
public:
  SupportPoint(Point p, unsigned int i, std::vector<unsigned int> l, std::vector<unsigned int> j = std::vector<unsigned int>()) : pp(p), id(i), al(l), m(j) { }
  ~SupportPoint() { }
  SupportPoint(const SupportPoint& sp) : pp(sp.pp), id(sp.id), al(sp.al), m(sp.m) { }
  SupportPoint& operator=(const SupportPoint& sp) { if (this==&sp) return *this; pp = sp.pp; id = sp.id; al = sp.al; m = sp.m; return *this; }

  void addMarker(unsigned int j) { m.push_back(j); }
  void setMarkers(std::vector<unsigned int> j) { m = j; }

  Point getPhysicalPoint() const { return pp; }
  unsigned int getID() const { return id; }
  std::vector<unsigned int> getAdjacencyList() const { return al; }
  std::vector<unsigned int> getMarkers() const { return m; }

protected:
  Point pp;
  unsigned int id;
  std::vector<unsigned int> al;
  std::vector<unsigned int> m;
}; // end class SupportPoint

class Cell { 
public:
  Cell(std::vector<unsigned int> l, unsigned int i) : id(i), t("notype"), lp(l) { }
  ~Cell() { }
  Cell(const Cell& c) : id(c.id), t(c.t), lp(c.lp) { }
  Cell& operator=(const Cell& c) { if (this==&c) return *this; id = c.id; t = c.t; lp = c.lp; return *this; }

  unsigned int getID() const { return id; }
  std::string getType() { return t; }
  std::vector<unsigned int> getSupportPointsID() const { return lp; }

protected:
  unsigned int id;
  std::string t;
  std::vector<unsigned int> lp;
}; // end class Cell

class Triangulation {
public:
  Triangulation(std::vector<SupportPoint> p = std::vector<SupportPoint>(), std::vector<Cell> c = std::vector<Cell>()) : lp(p), lc(c) { }
  ~Triangulation() { }
  Triangulation(const Triangulation& t) { }
  Triangulation& operator=(const Triangulation& t) { }

  void addSupportPoint(const SupportPoint& p) { lp.push_back(p); }
  void addCell(const Cell& c) { lc.push_back(c); }

protected:
  std::vector<SupportPoint> lp;
  std::vector<Cell> lc;
}; // end class Triangulation

////////////////////////// QUADRATURE RULE //////////////////////////////////////
class QuadratureRule {
public:
  QuadratureRule() : t("notype"), m(0) { }
  ~QuadratureRule() { }
  QuadratureRule(const QuadratureRule& qr) : qp(qr.qp), w(qr.w) { }
  QuadratureRule& operator=(const QuadratureRule& qr) { if (this==&qr) return *this; qp = qr.qp; w = qr.w; return *this; }

  unsigned int getNumberOfQuadraturePoints() const { assert(qp.size()==w.size()); return qp.size(); }
  Point getQuadraturePoint(unsigned int iqp) const { return qp[iqp]; }
  std::vector<Point> getQuadraturePoints() const { return qp; }
  double getWeight(unsigned int iqp) const { return w[iqp]; }
  std::vector<double> getWeights() const { return w; }
  std::vector<Point> getSupportPoints() const { return sp; }
  std::string getType() const { return t; }
  unsigned int getMaxOrder() const { return m; }

protected:
  std::vector<Point> qp;
  std::vector<double> w;
  std::vector<Point> sp;
  std::string t;
  unsigned int m;
}; // end class QuadratureRule

class GaussianTwoPoints : public QuadratureRule {
public:
  GaussianTwoPoints() {
    qp.push_back(Point(-1.0/sqrt(3.0))); w.push_back(1.0);
    qp.push_back(Point(+1.0/sqrt(3.0))); w.push_back(1.0);
    sp.push_back(Point(-1.0)); sp.push_back(Point(1.0));
    t = "GaussianTwoPoints";
    m = 3;
  }
}; // end class GaussianTwoPoints

class GaussianThreePoints : public QuadratureRule {
public:
  GaussianThreePoints() {
    qp.push_back(Point(0.0)); w.push_back(8.0 / 9.0);
    qp.push_back(Point(-sqrt(3.0/5.0))); w.push_back(5.0/9.0);
    qp.push_back(Point(+sqrt(3.0/5.0))); w.push_back(5.0/9.0);
    sp.push_back(Point(-1.0)); sp.push_back(Point(1.0));
    t = "GaussianThreePoints";
    m = 5;
  }
}; // end class GaussianThreePoints

class GaussianFourPoints : public QuadratureRule {
public:
  GaussianFourPoints() {
    qp.push_back(Point(-sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0))); w.push_back((18.0+sqrt(30.0))/36.0);
    qp.push_back(Point(+sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0))); w.push_back((18.0+sqrt(30.0))/36.0);
    qp.push_back(Point(-sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0))); w.push_back((18.0-sqrt(30.0))/36.0);
    qp.push_back(Point(+sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0))); w.push_back((18.0-sqrt(30.0))/36.0);
    sp.push_back(Point(-1.0)); sp.push_back(Point(1.0));
    t = "GaussianFourPoints";
    m = 7;
  }
}; // end class GaussianFourPoints

class GaussianFivePoints : public QuadratureRule {
public:
  GaussianFivePoints() {
    qp.push_back(Point(0.0)); w.push_back(128.0 / 225.0);
    qp.push_back(Point(-1.0/3.0*sqrt((5.0-2.0*sqrt(10.0/7.0))))); w.push_back((322.0+13.0*sqrt(70.))/900.0);
    qp.push_back(Point(+1.0/3.0*sqrt((5.0-2.0*sqrt(10.0/7.0))))); w.push_back((322.0+13.0*sqrt(70.))/900.0);
    qp.push_back(Point(-1.0/3.0*sqrt((5.0+2.0*sqrt(10.0/7.0))))); w.push_back((322.0-13.0*sqrt(70.))/900.0);
    qp.push_back(Point(+1.0/3.0*sqrt((5.0+2.0*sqrt(10.0/7.0))))); w.push_back((322.0-13.0*sqrt(70.))/900.0);
    sp.push_back(Point(-1.0)); sp.push_back(Point(1.0));
    t = "GaussianFivePoints";
    m = 9;
  }
}; // end class GaussianFivePoints

QuadratureRule* makeNewPointerToQuadratureRule(const std::string& t) {
  QuadratureRule* qr;
  if (t == "GaussianTwoPoints") {
    qr = new GaussianTwoPoints();
  } else if (t == "GaussianThreePoints") {
    qr = new GaussianThreePoints();
  } else if (t == "GaussianFourPoints") {
    qr = new GaussianFourPoints();
  } else if (t == "GaussianFivePoints") {
    qr = new GaussianFivePoints();
  } else {
    throw myFEM_INVALID_INPUT_STRING_EXCEPTION;
  }
  return qr;
}

QuadratureRule* makeNewPointerToQuadratureRule(unsigned int nqp) {
  QuadratureRule* qr;
  if (nqp == 2) {
    qr = new GaussianTwoPoints();
  } else if (nqp == 3) {
    qr = new GaussianThreePoints();
  } else if (nqp == 4) {
    qr = new GaussianFourPoints();
  } else if (nqp == 5) {
    qr = new GaussianFivePoints();
  } else {
    throw myFEM_INVALID_INPUT_STRING_EXCEPTION;
  }
  return qr;
}
////////////////////////// SHAPE FUNCTIONS //////////////////////////////////////
class BasisFunctions {
public:
  BasisFunctions() : t("notype"), o(0) { }
  ~BasisFunctions() { }
  BasisFunctions(const BasisFunctions& bf) : n(bf.n), t(bf.t) { }
  BasisFunctions& operator=(const BasisFunctions& bf) { if (this==&bf) return *this; n = bf.n; t = bf.t; return *this; }

  unsigned int getNumberOfNodes() const { return n.size(); }
  Point getNode(unsigned int idof) const { return n[idof]; }
  std::vector<Point> getNodes() const { return n; }
  std::string getType() const { return t; }
  unsigned int getOrder() const { return o; }

  std::vector<Point> getSupportPoints() const { std::vector<Point> sp; sp.push_back(n[0]); sp.push_back(n[1]); return sp;}

  virtual double getVal(unsigned int, const Point&) const = 0;
  virtual double getDx(unsigned int, const Point&) const = 0;

protected:
  std::vector<Point> n;
  std::string t;
  unsigned int o;
}; // end class BasisFunctions

class PiecewisePolynomial : public BasisFunctions {
public:             
  PiecewisePolynomial(unsigned int u, std::vector<Point> sp) { 
    assert(sp.size() == 2); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(0) / double(u)); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(u) / double(u)); 
    for (unsigned int a = 1; a < u; ++a) {
      n.push_back(sp[0] + (sp[1] - sp[0]) * double(a) / double(u)); 
    } // end for
    t = "PiecewisePolynomial"; 
    o = u;
    assert(o == n.size() - 1);
  }

  double getVal(unsigned int idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x); assert(idof < n.size());
    double value = 1.0;
    for (unsigned int a = 0; a < n.size(); ++a) {
      if (a != idof) {
        value *= (n[a].x - p.x) / (n[a].x - n[idof].x);
      } // end if a
    } // end for a
    return value;
  }

  double getDx(unsigned int idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x); assert(idof < n.size());
    double value = 0.0;
    static unsigned int k = 0;
//    std::cout<<">>>>"<<++k<<"<<<<  n.size()="<<n.size()<<"  o="<<o<<"  idof="<<idof<<"\n";
//    std::cout<<"parenthesis open\n";
    for (unsigned int a = 0; a < n.size(); ++a) {
      if (a != idof) {
        double tmp = 1.0;
//        std::cout<<"minus one times\n";
        for (unsigned int b = 0; b < n.size(); ++b) {
          if ((b != idof) 
              && (b != a)) {
//            std::cout<<"a="<<a<<"  b="<<b<<" times\n";
            tmp *= (n[b].x - p.x);
          } // end if b
        } // end for b
        value -= tmp;
      } //end if a
    } //end if a
//    std::cout<<"parenthesis close\n";
    for (unsigned int a = 0; a < n.size(); ++a) {
      if (a != idof) {
//         std::cout<<"divided by a="<<a<<"\n";
        value /= (n[a].x - n[idof].x);
      } //end if a
    } //end if a
//    std::cout<<std::endl;
//    if (k>=10) abort();
    return value;
  }
}; // end PiecewisePolynomial

class PiecewiseLinear : public BasisFunctions {
public:             
  PiecewiseLinear(std::vector<Point> sp) { 
    assert(sp.size() == 2); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(0) / double(1)); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(1) / double(1)); 
    t = "PiecewiseLinear"; 
    o = 1;
  }

  double getVal(unsigned int idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x);
    double value;
    if (idof == 0) {
      value = (n[1].x - p.x) / (n[1].x - n[0].x); 
    } else if (idof == 1) {
      value = (p.x - n[0].x) / (n[1].x - n[0].x); 
    } else {
      throw myFEM_DOF_OUT_OF_RANGE_EXCEPTION;
    }
    return value;
  }

  double getDx(unsigned int idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x);
    double value;
    if (idof == 0) {
      value = -1.0 / (n[1].x - n[0].x); 
    } else if (idof == 1) {
      value = 1.0 / (n[1].x - n[0].x); 
    } else {
      throw myFEM_DOF_OUT_OF_RANGE_EXCEPTION;
    }
    return value;
  }
}; // end class PiecewiseLinear

class PiecewiseQuadratic : public BasisFunctions {
public:             
  PiecewiseQuadratic(std::vector<Point> sp) { 
    assert(sp.size() == 2); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(0) / double(2)); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(2) / double(2)); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(1) / double(2)); 
    t = "PiecewiseQuadratic"; 
    o = 2;
  }

  double getVal(unsigned int idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x);
    double value;
    if (idof == 0) {
      value = (n[1].x - p.x) * (n[2].x - p.x) / ((n[1].x - n[0].x) * (n[2].x - n[0].x)); 
    } else if (idof == 1) {
      value = (n[2].x - p.x) * (n[0].x - p.x) / ((n[2].x - n[1].x) * (n[0].x - n[1].x)); 
    } else if (idof == 2) {
      value = (n[0].x - p.x) * (n[1].x - p.x) / ((n[0].x - n[2].x) * (n[1].x - n[2].x)); 
    } else {
      throw myFEM_DOF_OUT_OF_RANGE_EXCEPTION;
    }
    return value;
  }

  double getDx(unsigned int idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x);
    double value;
    if (idof == 0) {
      value = (2.0 * p.x - n[1].x - n[2].x) / ((n[1].x - n[0].x) * (n[2].x - n[0].x)); 
    } else if (idof == 1) {
      value = (2.0 * p.x - n[2].x - n[0].x) / ((n[2].x - n[1].x) * (n[0].x - n[1].x)); 
    } else if (idof == 2) {
      value = (2.0 * p.x - n[0].x - n[1].x) / ((n[0].x - n[2].x) * (n[1].x - n[2].x)); 
    } else {
      throw myFEM_DOF_OUT_OF_RANGE_EXCEPTION;
    }
    return value;
  }
}; // end class PiecewiseQuadratic

class PiecewiseCubic : public BasisFunctions {
public:             
  PiecewiseCubic(std::vector<Point> sp) { 
    assert(sp.size() == 2); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(0) / double(3)); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(3) / double(3)); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(1) / double(3)); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(2) / double(3)); 
    t = "PiecewiseCubic"; 
    o = 3;
  }

  double getVal(unsigned int idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x);
    double value;
    if (idof == 0) {
      value = (n[1].x - p.x) * (n[2].x - p.x) * (n[3].x - p.x) / ((n[1].x - n[0].x) * (n[2].x - n[0].x) * (n[3].x - n[0].x)); 
    } else if (idof == 1) {
      value = (n[2].x - p.x) * (n[3].x - p.x) * (n[0].x - p.x) / ((n[2].x - n[1].x) * (n[3].x - n[1].x) * (n[0].x - n[1].x)); 
    } else if (idof == 2) {
      value = (n[3].x - p.x) * (n[0].x - p.x) * (n[1].x - p.x) / ((n[3].x - n[2].x) * (n[0].x - n[2].x) * (n[1].x - n[2].x)); 
    } else if (idof == 3) {
      value = (n[0].x - p.x) * (n[1].x - p.x) * (n[2].x - p.x) / ((n[0].x - n[3].x) * (n[1].x - n[3].x) * (n[2].x - n[3].x)); 
    } else {
      throw myFEM_DOF_OUT_OF_RANGE_EXCEPTION;
    }
    return value;
  }

  double getDx(unsigned int idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x);
    double value;
    if (idof == 0) {
      value = ( - (n[2].x - p.x) * (n[3].x - p.x) - (n[1].x - p.x) * (n[3].x - p.x) - (n[1].x - p.x) * (n[2].x - p.x)) / ((n[1].x - n[0].x) * (n[2].x - n[0].x) * (n[3].x - n[0].x)); 
    } else if (idof == 1) {
      value = ( - (n[3].x - p.x) * (n[0].x - p.x) - (n[2].x - p.x) * (n[0].x - p.x) - (n[2].x - p.x) * (n[3].x - p.x)) / ((n[2].x - n[1].x) * (n[3].x - n[1].x) * (n[0].x - n[1].x)); 
    } else if (idof == 2) {
      value = ( - (n[0].x - p.x) * (n[1].x - p.x) - (n[3].x - p.x) * (n[1].x - p.x) - (n[3].x - p.x) * (n[0].x - p.x)) / ((n[3].x - n[2].x) * (n[0].x - n[2].x) * (n[1].x - n[2].x)); 
    } else if (idof == 3) {
      value = ( - (n[1].x - p.x) * (n[2].x - p.x) - (n[0].x - p.x) * (n[2].x - p.x) - (n[0].x - p.x) * (n[1].x - p.x)) / ((n[0].x - n[3].x) * (n[1].x - n[3].x) * (n[2].x - n[3].x)); 
    } else {
      throw myFEM_DOF_OUT_OF_RANGE_EXCEPTION;
    }
    return value;
  }
}; // end class PiecewiseCubic

////////////////////////// REFERENCE ELEMENT //////////////////////////////////////
class ReferenceElement {
public:
  ReferenceElement(unsigned int u) {
    std::vector<Point> sp; 
    sp.push_back(Point(-1.0));
    sp.push_back(Point(1.0));
    bf = new PiecewisePolynomial(u, sp);
  }
  ReferenceElement(const std::string& s = "PiecewiseLinear", unsigned int u = 0) { 
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

  unsigned int getNumberOfNodes() const { return bf->getNumberOfNodes(); }
  Point getNode(unsigned int idof) const { return bf->getNode(idof); }
  std::vector<Point> getNodes() const { return bf->getNodes(); }
  std::vector<Point> getSupportPoints() const { return bf->getSupportPoints(); }
  double getVal(unsigned int idof, Point p) const { return bf->getVal(idof, p); }
  double getDx(unsigned int idof, Point p) const { return bf->getDx(idof, p); }
  std::string getTypeOfBasisFunctions() const { return bf->getType(); }
  unsigned int getOrder() const { return bf->getOrder(); }

protected:
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

  void setDOF(const std::vector<unsigned int>& gdof) { 
    if ((dof.size() != 0) 
        &&(m.size() != 0)) {
      std::cerr<<"DOF have already been set before..."<<std::endl;
      abort();
    }
    // Get the number of DOF
    const unsigned int ndof = re->getNumberOfNodes();
    assert(gdof.size() == ndof);
    dof.resize(ndof);
    for (unsigned int idof = 0; idof < ndof; ++idof) {
      dof[idof] = gdof[idof];
      m[idof] = gdof[idof];
    } // end for idof
  }

  Point mapGlobal2Local(Point g) const {
    std::vector<Point> ref_sp = re->getSupportPoints();
    assert(sp[0].x <= g.x); assert(g.x <= sp[1].x);
    Point l;
    l.x = ref_sp[0].x + (ref_sp[1].x - ref_sp[0].x) * (g.x - sp[0].x) / (sp[1].x - sp[0].x);
    return l;
  }
  Point mapLocal2Global(Point l) const {
    std::vector<Point> ref_sp = re->getSupportPoints();
    assert(ref_sp[0].x <= l.x); assert(l.x <= ref_sp[1].x);
    Point g;
    g.x = sp[0].x + (sp[1].x - sp[0].x) * (l.x - ref_sp[0].x) / (ref_sp[1].x - ref_sp[0].x);
    return g;
  }
  double getJacobian() const { std::vector<Point> ref_sp = re->getSupportPoints(); return (sp[1].x - sp[0].x) / (ref_sp[1].x - ref_sp[0].x); }
  double getValByLocal(unsigned int idof, Point l) const { return re->getVal(idof, l); }
  double getValByGlobal(unsigned int idof, Point g) const { Point l = mapGlobal2Local(g); return re->getVal(idof, l); }
  double getDxByLocal(unsigned int idof, Point l) const { return re->getDx(idof, l); }
  double getDxByGlobal(unsigned int idof, Point g) const { Point l = mapGlobal2Local(g); return re->getDx(idof, l); }
  unsigned int getNumberOfNodes() const { return re->getNumberOfNodes(); }
  unsigned int getNumberOfDOF() const { assert(dof.size() != 0); return dof.size(); }
  std::map<unsigned int, unsigned int> getDOFMap() const { assert(m.size() != 0); return m; }
  std::string getTypeOfBasisFunctions() const { return re->getTypeOfBasisFunctions(); }
  unsigned int getOrder() const { return re->getOrder(); }

protected:
  std::vector<Point> sp;
  std::vector<unsigned int> dof;
  std::map<unsigned int, unsigned int> m;
  ReferenceElement *re;
}; // end class FiniteElement

////////////////////////// DOF HANDLER //////////////////////////////////////
class DOFHandler { 
public:
  DOFHandler(std::vector<FiniteElement*>* m) : ndof(0), fe(m) { }
  ~DOFHandler() { }
  
  void distributeDOF() {
    assert(ndof == 0);
    const unsigned int nel = fe->size();
    unsigned int prec = 0;
    for (unsigned int iel = 0; iel < nel; ++iel) {
      unsigned int iord = (*fe)[iel]->getOrder();
      unsigned int prev = (iel > 0 ? ndof - (*fe)[iel-1]->getOrder() : 0);
      std::vector<unsigned int> dof;
      if (iel > 0) { 
        dof.push_back(prev);
      } else {
        assert(ndof == 0);
        dof.push_back(ndof++);
      } // end if first element
      for (unsigned a = 0; a < iord; ++a) {
        dof.push_back(ndof++);
      } // end add nodes
      (*fe)[iel]->setDOF(dof);
    } // end for iel
    assert(ndof != 0);
  }

  unsigned int getNumberOfDOF() const { assert(ndof != 0); return ndof; }

protected:
  unsigned int ndof;
  std::vector<FiniteElement*>* fe;
}; // end class DOFHandler

////////////////////////// FE VALUES //////////////////////////////////////
class UpdateFlags { };
class FEValues {
public:
  FEValues(FiniteElement *e, QuadratureRule *r, UpdateFlags f) : fe(e), qr(r), uf(f) { } 

  void reinit(FiniteElement *e) {
    fe = e;
    DeterminantOfJacobian = fabs(fe->getJacobian());
    InverseJacobian = 1.0 / (fe->getJacobian());
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
      } // end for iqp
    } // end for idof
  }

  double getDeterminantOfJacobianTimesWeight(unsigned int iqp) const { return Weights[iqp]*DeterminantOfJacobian; }
  double getShapeValue(unsigned int idof, unsigned int iqp) const { return ShapeValues[idof][iqp]; }
  double getShapeDx(unsigned int idof, unsigned int iqp) const { return InverseJacobian * ShapeDx[idof][iqp]; }
  std::vector<Point> getQuadraturePoints() const {
    std::vector<Point> qps;
    const unsigned int nqp = qr->getNumberOfQuadraturePoints();
    for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
      qps.push_back(fe->mapLocal2Global(qr->getQuadraturePoint(iqp)));
    }
    return qps;
  }
  FiniteElement* getFiniteElement() const { return fe; }

protected:
  std::vector<double> Weights;
  double DeterminantOfJacobian;
  double InverseJacobian;
  std::vector<std::vector<double> > ShapeValues;
  std::vector<std::vector<double> > ShapeDx;

  FiniteElement *fe;
  QuadratureRule *qr;
  UpdateFlags uf;
}; // end class FEValues

// solve -d/dx(a(x) d/dx(u)) + q(x) u(x) = f(x) for u(x) and x in (0,1)
double a(double x) { return 1 + x; }
double q(double x) { return x; }
double f(double x) { return x*cos(M_PI*x) + M_PI*sin(M_PI*x) + M_PI*M_PI*(1+x)*cos(M_PI*x); }
double u(double x) { return cos(M_PI*x); }
double a(Point p) { return a(p.x); }
double q(Point p) { return q(p.x); }
double f(Point p) { return f(p.x); }
double u(Point p) { return u(p.x); }
double aMass(Point p) { return 0.0; }
double qMass(Point p) { return 1.0; }
double aStiffness(Point p) { return 1.0; }
double qStiffness(Point p) { return 0.0; }

void assembleLocal(FEValues *feValues, 
    std::vector<std::vector<double> >& localMatrix,
    std::vector<double>& localVector,
    Object_t handleObject = myFEM_BOTH_MATRIX_AND_VECTOR,
    double (*pa)(Point) = &a,
    double (*pq)(Point) = &q,
    double (*pf)(Point) = &f,
    DebugLevel_t debugLevel = myFEM_NO_DEBUG,
    std::ostream& os = std::cout) {
  // Get the quadrature points
  std::vector<Point> quadraturePoints = feValues->getQuadraturePoints();
  const unsigned int nqp = quadraturePoints.size();
  if (debugLevel > myFEM_NO_DEBUG) {
    // TODO: better debug here
    assert(localVector.size() == 0);
    assert(localMatrix.size() == 0);
    os<<"quadrature points\n";
    for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
      os<<iqp<<" = "<<quadraturePoints[iqp]<<"\n";
    } // end for iqp
  } // end if debugLevel not zero
  // Get the number of DOF
  const unsigned int ndof = feValues->getFiniteElement()->getNumberOfDOF();
  // Initialize local matrix and vector
  if (handleObject != myFEM_VECTOR_ONLY) {
    localMatrix = std::vector<std::vector<double> >(ndof, std::vector<double>(ndof, 0.0));
  } // end if handleObject not vector only
  if (handleObject != myFEM_MATRIX_ONLY) {
    localVector = std::vector<double>(ndof, 0.0);
  } // end if handleObject not matrix only
  // Compute local matrix and vector
  for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
    for (unsigned int idof = 0; idof < ndof; ++idof) {
      if (handleObject != myFEM_VECTOR_ONLY) {
        for (unsigned int jdof = 0; jdof < ndof; ++jdof) {
          localMatrix[idof][jdof] += ( pa(quadraturePoints[iqp]) * feValues->getShapeDx(idof, iqp) * feValues->getShapeDx(jdof, iqp)
                                     + pq(quadraturePoints[iqp]) * feValues->getShapeValue(idof, iqp) * feValues->getShapeValue(jdof, iqp)
                                     ) * feValues->getDeterminantOfJacobianTimesWeight(iqp);
          if (debugLevel == myFEM_MAX_DEBUG) {
            os<<"iqp="<<iqp<<quadraturePoints[iqp]<<"idof="<<idof<<"  jdof="<<jdof<<"  "
              <<"a="<<a(quadraturePoints[iqp])<<"  "
              <<"q="<<q(quadraturePoints[iqp])<<"  "
              <<"phi_"<<idof<<"="<<feValues->getShapeValue(idof, iqp)<<"  "
              <<"phi_"<<jdof<<"="<<feValues->getShapeValue(jdof, iqp)<<"  "
              <<"DphiDx_"<<idof<<"="<<feValues->getShapeDx(idof, iqp)<<"  "
              <<"DphiDx_"<<jdof<<"="<<feValues->getShapeDx(jdof, iqp)<<"  "
              <<"JxW="<<feValues->getDeterminantOfJacobianTimesWeight(iqp)<<"\n";
          } // end if debugLevel hardcore
        } // end for jdof
      } // end if handleObject not vector only
      if (handleObject != myFEM_MATRIX_ONLY) {
        localVector[idof] += pf(quadraturePoints[iqp]) * feValues->getShapeValue(idof, iqp) * feValues->getDeterminantOfJacobianTimesWeight(iqp);
        if (debugLevel == myFEM_MAX_DEBUG) {
          os<<"iqp="<<iqp<<quadraturePoints[iqp]<<"idof="<<idof<<"  "
            <<"f="<<f(quadraturePoints[iqp])<<"  "
            <<"phi_"<<idof<<"="<<feValues->getShapeValue(idof, iqp)<<"  "
            <<"JxW="<<feValues->getDeterminantOfJacobianTimesWeight(iqp)<<"\n";
        } //end if debugLevel hardcore
      } // end if handleObject not matrix only
    } // end for idof
  } //end for iqp           
}

void distributeLocal2Global(FEValues *feValues,
    const std::vector<std::vector<double> >& localMatrix,
    const std::vector<double>& localVector,
    std::vector<std::vector<double> >& globalMatrix,
    std::vector<double>& globalVector,
    Object_t handleObject = myFEM_BOTH_MATRIX_AND_VECTOR,
    DebugLevel_t debugLevel = myFEM_NO_DEBUG) {
  // Get the number of DOF
  const unsigned int nldof = feValues->getFiniteElement()->getNumberOfDOF();
  if (debugLevel > myFEM_NO_DEBUG) {
    unsigned int nldofCHECK = 0;
    unsigned int ngdof = 0;
    if (handleObject != myFEM_MATRIX_ONLY) {
      nldofCHECK = localVector.size();
      ngdof = globalVector.size();
    } else {
      nldofCHECK = localMatrix.size();
      ngdof = globalMatrix.size();
    }
    assert(nldofCHECK == nldof);
    assert(nldof != 0);
    assert(ngdof != 0);
    if (handleObject == myFEM_BOTH_MATRIX_AND_VECTOR) {
      assert(localVector.size() == localMatrix.size());
      assert(globalVector.size() == globalMatrix.size());
    } // end if handleObject not matrix only
    if (debugLevel == myFEM_MAX_DEBUG) {
      if (handleObject != myFEM_VECTOR_ONLY) {
        for (unsigned int ildof = 0; ildof < nldof; ++ildof) {
          assert(localMatrix[ildof].size() == nldof);
        } // end for ildof
        for (unsigned int igdof = 0; igdof < ngdof; ++igdof) {
          assert(globalMatrix[igdof].size() == ngdof);
        } // end for igdof
      } // end if handleObject not vector only
    } // end if debugLevel hardcore
  } // end if debugLevel not zero
  // Get the DOF map
  std::map<unsigned int, unsigned int> dofMap = feValues->getFiniteElement()->getDOFMap();
  if (debugLevel > myFEM_NO_DEBUG) {
    std::cout<<"DOF map\n";
    for (unsigned int ildof = 0; ildof < nldof; ++ildof) {
      unsigned int igdof = dofMap[ildof];
      std::cout<<ildof<<" -> "<<igdof<<"\n";
    } // end for ildof
  } // end if debugLevel not zero
  // Distribute local to global
  for (unsigned int ildof = 0; ildof < nldof; ++ildof) {
    unsigned int igdof = dofMap[ildof];
    if (handleObject != myFEM_VECTOR_ONLY) {
      for (unsigned int jldof = 0; jldof < nldof; ++jldof) {
        unsigned int jgdof = dofMap[jldof];
        globalMatrix[igdof][jgdof] += localMatrix[ildof][jldof];
      } // end for jldof
    } // end if handleObject not vector only
    if (handleObject != myFEM_MATRIX_ONLY) {
      globalVector[igdof] += localVector[ildof];
    } // end for ildof
  } // end if handleObject not matrix only
}

void printMatrixAndVector(const std::vector<std::vector<double> >& Matrix,
    const std::vector<double>& Vector,
    Object_t handleObject = myFEM_BOTH_MATRIX_AND_VECTOR,
    std::ostream& os = std::cout,
    int setwVal = 7,
    int setprecisionVal = 3,
    DebugLevel_t debugLevel = myFEM_NO_DEBUG) {
  // Get the number of DOF
  unsigned int ndof = 0;
  if (handleObject != myFEM_MATRIX_ONLY) {
    ndof = Vector.size();
  } else {
    ndof = Matrix.size();
  } // end if handleObject not matrix only
  if (debugLevel > myFEM_NO_DEBUG) {
    assert(ndof != 0);
    if (handleObject == myFEM_BOTH_MATRIX_AND_VECTOR) {
      assert(Vector.size() == Matrix.size());
    } // end if handleObject both matrix and vector
    if (debugLevel == myFEM_MAX_DEBUG) {
      if (handleObject != myFEM_VECTOR_ONLY) {
        for (unsigned int idof = 0; idof < ndof; ++idof) {
          assert(Matrix[idof].size() == ndof);
        } // end for idof
      } // end if handleObject not vector only
    } // end if debugLevel hardcore
  } // end if debugLevel not zero
  // Print the stuff out...
  for (unsigned int idof = 0; idof < ndof; ++idof) {
    if (handleObject != myFEM_VECTOR_ONLY) {
      for (unsigned int jdof = 0; jdof < ndof; ++jdof) {
        os<<std::setw(setwVal)
          <<std::setprecision(setprecisionVal)
          <<Matrix[idof][jdof]<<"  ";
      } // end for jdof
    } //end if handleObject not vector only
    if (handleObject == myFEM_BOTH_MATRIX_AND_VECTOR) {
      os<<"||  ";
    } // end if handleObject both matrix and vector
    if (handleObject != myFEM_MATRIX_ONLY) {
      os<<std::setw(setwVal)
        <<std::setprecision(setprecisionVal)
        <<Vector[idof]<<"\n";
    } else {
      os<<"\n";
    } // end if handleObject not matrix only
  } // end for idof
}
                    
// TODO: take COOMatrix as input instead of vector of vector 
const double myFEM_EPSILON = 1.0e-10;
std::vector<double> solveMatrixTimesXEqualsRHS(const std::vector<std::vector<double> >& Matrix,
    const std::vector<double>& RHS,
    bool verbose = false,
    std::ostream& os = std::cout,
    DebugLevel_t debugLevel = myFEM_NO_DEBUG) {
  // Get nRow and nCol
  const unsigned int ndof = RHS.size();
  int nRow = ndof;
  int nCol = ndof;
  if (debugLevel != myFEM_NO_DEBUG) {
    assert(Matrix.size() == ndof);
    if (debugLevel == myFEM_MAX_DEBUG) {
      for (unsigned int idof = 0; idof < ndof; ++idof) {
        assert(Matrix[idof].size() == ndof);
      } // end for idof
    } // enf if debugLevel hardcore
  } // end if debugLevel not zero
                      
  // TODO: get iRow, jCol, and ijVal directly from COOMatrix
  std::vector<double> iRowTMP, jColTMP, ijValTMP;
  // Fill iRow, jCol, and ijVal
  for (unsigned int i = 0; i < ndof; ++i) {
    for (unsigned int j = 0; j < ndof; ++j) {
      if (fabs(Matrix[i][j]) > myFEM_EPSILON) { 
        //static int k = 0;
        //std::cout<<++k<<" -> i="<<i<<" j="<<j<<" val="<<globalMatrix[i][j]<<std::endl;
        iRowTMP.push_back(i);
        jColTMP.push_back(j);
        ijValTMP.push_back(Matrix[i][j]);
      } // end if nonzero val in global matrix
    } // end for j
  } // end for i
  // Compute nNZ
  int nNZ = ijValTMP.size();
  std::cout<<"number of nonzero values = "<<nNZ<<"\n";
  // Turn vectors into arrays...
  int iRowCOO[nNZ];
  int jColCOO[nNZ];
  double ijValCOO[nNZ];
  std::copy(iRowTMP.begin(), iRowTMP.end(), iRowCOO);
  std::copy(jColTMP.begin(), jColTMP.end(), jColCOO);
  std::copy(ijValTMP.begin(), ijValTMP.end(), ijValCOO);

  double XArray[ndof];
  double RHSArray[ndof];
  std::copy(RHS.begin(), RHS.end(), RHSArray);

  int pColCC[nCol+1];
  int iRowCC[nNZ];
  double iValCC[nNZ];
  int Map[nNZ];
  int sys = UMFPACK_A;

  // TODO: print some info if verbose is true
  enum UMFPACK_STATUS_DUMMY_ENUM { TRIPLET_TO_COL, SYMBOLIC, NUMERIC, SOLVE, UMFPACK_STATUS_DUMMY_ENUM_SIZE };
  int status[UMFPACK_STATUS_DUMMY_ENUM_SIZE];

  double Info[UMFPACK_INFO];
  double Control[UMFPACK_CONTROL];
  void *Symbolic;
  void *Numeric;

  status[TRIPLET_TO_COL] = umfpack_di_triplet_to_col(nRow, nCol, nNZ, iRowCOO, jColCOO, ijValCOO, pColCC, iRowCC, iValCC, Map);
  status[SYMBOLIC] = umfpack_di_symbolic(nRow, nCol, pColCC, iRowCC, iValCC, &Symbolic, Control, Info);
  status[NUMERIC] = umfpack_di_numeric(pColCC, iRowCC, iValCC, Symbolic, &Numeric, Control, Info);
  status[SOLVE] = umfpack_di_solve(sys, pColCC, iRowCC, iValCC, XArray, RHSArray, Numeric, Control, Info);

  std::vector<double> XVector(ndof);
  std::copy(XArray, XArray + ndof, XVector.begin());

  return XVector;
}

// TODO: replace prefix myFEM_ by namespace
double computeNorm(const std::vector<std::vector<double> >& MassMatrix,
    const std::vector<double>& Vector,
    Norm_t normType = myFEM_L2_NORM,
    const std::vector<std::vector<double> >& StiffnessMatrix = std::vector<std::vector<double> >()) {
  double Norm = 0.0;
  unsigned int ndof = Vector.size();
  for (unsigned int idof = 0; idof < ndof; ++idof) {
    for (unsigned int jdof = 0; jdof < ndof; ++jdof) {
      if (normType == myFEM_L1_NORM) {
        Norm += MassMatrix[idof][jdof] * fabs(Vector[jdof]);
      } else {
        Norm += MassMatrix[idof][jdof] * Vector[jdof] * Vector[jdof];
        if (normType == myFEM_H1_NORM) {
          Norm += StiffnessMatrix[idof][jdof] * Vector[jdof] * Vector[jdof];
        } // end if H1
      } // end if not L1
    } // end for jdof
  } // end for idof
  return Norm;
}

int main(int argc, char *argv[]) {

  { /** nouveau test */
  // Default values
  unsigned int nel = 10;
  unsigned int polynomialOrder = 1;
  unsigned int numberOfGaussPoints = 2;
  //TODO: make command line argument passing smarter than this
  std::string errorMessage = std::string("Usage: argv[0] number_of_elements type_of_basis_functions type_of_quadrature_rule\n")
                           + std::string("number_of_elements -> a positive integer [default value is 10]\n")
                           + std::string("type_of_basis_functions -> [1] for PiecewiseLinear\n")
                           + std::string("                        -> 2 for PiecewiseQuadratic\n")
                           + std::string("                        -> 3 for PiecewiseCubic\n")
                           + std::string("                        -> 4 for PiecewiseQuartic\n")
                           + std::string("type_of_quadrature_rule -> [2] for GaussianTwoPoints\n")
                           + std::string("                        -> 3 for GaussianThreePoints\n")
                           + std::string("                        -> 4 for GaussianFourPoints\n")
                           + std::string("                        -> 5 for GaussianFivePoints\n");
  //for (int i = 0; i < argc; ++i) std::cout<<argv[i]<<" "; std::cout<<"\n";
  if (argc > 1) {
    nel = atoi(argv[1]);
  }
  if (argc > 2) {
    if (atoi(argv[2]) > 4) { 
      std::cerr<<"doucement garcon as tu vraiment besoin d'un ordre aussi eleve?"<<std::endl;
      abort();
    } // if order too big
      polynomialOrder = atoi(argv[2]);
  }
  if (argc > 3) {
    if (atoi(argv[3]) >  5) { 
      std::cerr<<"max quadrature rule is 5 points"<<std::endl;
      abort();
    }
    numberOfGaussPoints = atoi(argv[3]);
  }

  std::cout<<"#### BEGIN ######\n";
  // Construct nel elements
  std::cout<<"#### BUILD FINITE ELEMENTS ######\n";
  std::vector<FiniteElement*> finiteElements;
  // Create a pointer to a reference element
  ReferenceElement *referenceElement = new ReferenceElement(polynomialOrder);

  // Create triangualtion embryo
  // TODO: move towards Triangulation class
  Point startPoint(0.0); Point endPoint(1.0);
  for (unsigned int iel = 0; iel < nel; ++iel) {
    std::vector<Point> supportPoints;
    supportPoints.push_back(startPoint+double(iel+0)/double(nel)*(endPoint-startPoint));
    supportPoints.push_back(startPoint+double(iel+1)/double(nel)*(endPoint-startPoint));
    //std::cout<<iel<<"  "<<supportPoints[0]<<"  "<<supportPoints[1]<<"\n";
    finiteElements.push_back(new FiniteElement(supportPoints, referenceElement));
  } // end for iel

  // Distribute the Degrees Of Freedoms (DOF)
  std::cout<<"#### DISTRIBUTE DOF ######\n";
  std::cout<<std::endl;
  DOFHandler* dofHandler = new DOFHandler(&finiteElements);
  dofHandler->distributeDOF();
  std::cout<<"OK"<<std::endl;
  unsigned int ndof = dofHandler->getNumberOfDOF();

  // Create global matrix and global rhs
  std::vector<std::vector<double> > globalMatrix(ndof, std::vector<double>(ndof, 0.0));
  std::vector<double> globalRHS(ndof, 0.0);
  // Create mass and stiffness matrix
  std::vector<std::vector<double> > StiffnessMatrix(ndof, std::vector<double>(ndof, 0.0));
  std::vector<std::vector<double> > MassMatrix(ndof, std::vector<double>(ndof, 0.0));
  // Create null vector and null matrix
  std::vector<std::vector<double> > nullMatrix;
  std::vector<double> nullVector;

  // Create pointer to a quadrature rule
  QuadratureRule *quadratureRule = makeNewPointerToQuadratureRule(numberOfGaussPoints);

  std::cout<<"#### ASSEMBLE MATRIX AND RHS ######\n";
  UpdateFlags updateFlags;
  FEValues *feValues = new FEValues(NULL, quadratureRule, updateFlags); 
  for (unsigned int iel = 0; iel < nel; ++iel) {
    std::cout<<"element #"<<iel<<"\n";
    // compute FE values for element iel
    feValues->reinit(finiteElements[iel]);

    // compute local matrix and RHS
    std::vector<std::vector<double> > localMatrix;
    std::vector<double> localRHS; 
    assembleLocal(feValues, localMatrix, localRHS, myFEM_BOTH_MATRIX_AND_VECTOR, &a, &q, &f, myFEM_NO_DEBUG, std::cout);
    //assembleLocal(feValues, localMatrix, nullVector, myFEM_MATRIX_ONLY, &a, &q, NULL, myFEM_MAX_DEBUG, std::cout);
    //assembleLocal(feValues, nullMatrix, localRHS, myFEM_VECTOR_ONLY, NULL, NULL, &f, myFEM_MAX_DEBUG, std::cout);
    // compute local mass and stiffness matrices
    std::vector<std::vector<double> > localStiffnessMatrix;
    std::vector<std::vector<double> > localMassMatrix;
    assembleLocal(feValues, localStiffnessMatrix, nullVector, myFEM_MATRIX_ONLY, &aStiffness, &qStiffness, NULL, myFEM_NO_DEBUG, std::cout);
    assembleLocal(feValues, localMassMatrix, nullVector, myFEM_MATRIX_ONLY, &aMass, &qMass, NULL, myFEM_NO_DEBUG, std::cout);

    // distribute local to global
    distributeLocal2Global(feValues, localMatrix, localRHS, globalMatrix, globalRHS, myFEM_BOTH_MATRIX_AND_VECTOR, myFEM_NO_DEBUG);
    std::cout<<"local matrix and local rhs\n";
    printMatrixAndVector(localMatrix, localRHS, myFEM_BOTH_MATRIX_AND_VECTOR, std::cout);
    // same thing for mass and stiffness
    distributeLocal2Global(feValues, localStiffnessMatrix, nullVector, StiffnessMatrix, nullVector, myFEM_MATRIX_ONLY);
    distributeLocal2Global(feValues, localMassMatrix, nullVector, MassMatrix, nullVector, myFEM_MATRIX_ONLY);
    std::cout<<"local mass matrix\n";
    printMatrixAndVector(localMassMatrix, nullVector, myFEM_MATRIX_ONLY, std::cout);
    std::cout<<"local stiffness matrix\n";
    printMatrixAndVector(localStiffnessMatrix, nullVector, myFEM_MATRIX_ONLY, std::cout);

  } // end for iel
  std::cout<<"#### END ASSEMBLY ROUTINE ######\n";

  if (nel <= 2) {
    // Print out global matrix and global rhs
    std::cout<<"global matrix and global rhs\n";
    printMatrixAndVector(globalMatrix, globalRHS, myFEM_BOTH_MATRIX_AND_VECTOR, std::cout);
    // ... same thing with mass and stiffness
    std::cout<<"stiffness matrix\n";
    printMatrixAndVector(StiffnessMatrix, nullVector, myFEM_MATRIX_ONLY, std::cout);
    std::cout<<"mass matrix\n";
    printMatrixAndVector(MassMatrix, nullVector, myFEM_MATRIX_ONLY, std::cout);
  }
  //std::cout<<std::endl;
  //abort();

  // Print matrix and RHS to files to check with matlab
  std::fstream foutMatrix; foutMatrix.open("matrix.dat", std::fstream::out); printMatrixAndVector(globalMatrix, nullVector, myFEM_MATRIX_ONLY, foutMatrix, 9, 5, myFEM_MAX_DEBUG); foutMatrix.close();
  std::fstream foutRHS; foutRHS.open("rhs.dat", std::fstream::out); printMatrixAndVector(nullMatrix, globalRHS, myFEM_VECTOR_ONLY, foutRHS, 9, 5, myFEM_MAX_DEBUG); foutRHS.close();

#ifdef EBILE
  std::cout<<"tu te crois malin hein? gros ebile :D\n";
#endif

  // Solve the FE problem
  std::cout<<"#### SOLVE WITH UMFPACK ######\n";
  std::vector<double> solutionVector = solveMatrixTimesXEqualsRHS(globalMatrix, globalRHS);

  std::cout<<"#### POSTPROCESSING ######\n";
  // Compute L2 norm of the numerical solution
  std::cout<<"solution L2 Norm = "<<computeNorm(MassMatrix, solutionVector, myFEM_L2_NORM)<<"\n"; 
  
  // Compute exact solution
  // TODO: come up with something better than this
  //       need iterator over nodes
  std::vector<double> exactSolutionVector(ndof);
  for (unsigned int iel = 0; iel < nel; ++iel) {
    std::vector<Point> supportPoints;
    supportPoints.push_back(startPoint+double(iel)/double(nel)*(endPoint-startPoint));
    supportPoints.push_back(startPoint+double(iel+1)/double(nel)*(endPoint-startPoint));
    //std::cout<<iel<<"  "<<supportPoints[0]<<"  "<<supportPoints[1]<<"\n";
//    if (referenceElement->getTypeOfBasisFunctions() == "PiecewiseLinear") {
    if (referenceElement->getOrder() == 1) {
      if (iel > 0) {
        exactSolutionVector[1*iel-0] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(0) / double(1));
      } else {
        exactSolutionVector[1*iel+0] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(0) / double(1));
      } // end if iel
      exactSolutionVector[1*iel+1] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(1) / double(1));
//    } else if (referenceElement->getTypeOfBasisFunctions() == "PiecewiseQuadratic") {
    } else if (referenceElement->getOrder() == 2) {
      if (iel > 0) {
        exactSolutionVector[2*iel-1] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(0) / double(2));
      } else {
        exactSolutionVector[2*iel+0] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(0) / double(2));
      } // end if iel
      exactSolutionVector[2*iel+1] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(2) / double(2));
      exactSolutionVector[2*iel+2] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(1) / double(2));
//    } else if (referenceElement->getTypeOfBasisFunctions() == "PiecewiseCubic") {
    } else if (referenceElement->getOrder() == 3) {
      if (iel > 0) {
        exactSolutionVector[3*iel-2] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(0) / double(3));
      } else {
        exactSolutionVector[3*iel+0] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(0) / double(3));
      } // end if iel
      exactSolutionVector[3*iel+1] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(3) / double(3));
      exactSolutionVector[3*iel+2] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(1) / double(3));
      exactSolutionVector[3*iel+3] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(2) / double(3));
//    } else if (referenceElement->getTypeOfBasisFunctions() == "PiecewiseQuartic") {
    } else if (referenceElement->getOrder() == 4) {
      if (iel > 0) {
        exactSolutionVector[4*iel-3] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(0) / double(4));
      } else {
        exactSolutionVector[4*iel+0] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(0) / double(4));
      } // end if iel
      exactSolutionVector[4*iel+1] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(4) / double(4));
      exactSolutionVector[4*iel+2] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(1) / double(4));
      exactSolutionVector[4*iel+3] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(2) / double(4));
      exactSolutionVector[4*iel+4] = u(supportPoints[0] + (supportPoints[1] - supportPoints[0]) * double(3) / double(4));
    } else {
      throw myFEM_INVALID_INPUT_STRING_EXCEPTION;
    }
  } // end for iel
  // Get its L2 norm
  std::cout<<"exact solution L2 Norm = "<<computeNorm(MassMatrix, exactSolutionVector, myFEM_L2_NORM)<<"\n"; 
  if (nel <= 2) {
    std::cout<<"solution\n";
    printMatrixAndVector(nullMatrix, solutionVector, myFEM_VECTOR_ONLY);
    std::cout<<"exact solution\n";
    printMatrixAndVector(nullMatrix, exactSolutionVector, myFEM_VECTOR_ONLY);
  }
    
  // Compute exact error
  std::vector<double> exactErrorVector;
  for (unsigned int idof = 0; idof < ndof; ++idof) {
    exactErrorVector.push_back(solutionVector[idof] - exactSolutionVector[idof]);
  } // end for idof
  //std::cout<<"exact error\n";
  //printMatrixAndVector(std::vector<std::vector<double> >(), exactErrorVector, myFEM_VECTOR_ONLY, std::cout);

  // ... and get its L1, L2, and H1 norms
  std::cout<<"error L1 Norm = "<<computeNorm(MassMatrix, exactErrorVector, myFEM_L1_NORM)<<"\n"; 
  std::cout<<"error L2 Norm = "<<computeNorm(MassMatrix, exactErrorVector, myFEM_L2_NORM)<<"\n"; 
  std::cout<<"error H1 Norm = "<<computeNorm(MassMatrix, exactErrorVector, myFEM_H1_NORM, StiffnessMatrix)<<"\n"; 

  std::cout<<"number of elements = "<<nel<<"\n";
  std::cout<<"number of DOF = "<<ndof<<"\n";
  std::cout<<"basis functions = "<<referenceElement->getTypeOfBasisFunctions()<<"\n";
  std::cout<<"quadrature rule = "<<quadratureRule->getType()<<"\n";

  // Delete allocated memory
  delete feValues;
  for (unsigned int iel = 0; iel < nel; ++iel) { delete finiteElements[iel]; }
  delete referenceElement;
  delete quadratureRule;
  std::cout<<"#### END ######\n";
  }


  if (false)
  { /** test finite element */
  ReferenceElement *referenceElement = new ReferenceElement;
  std::vector<Point> supportPoints;  supportPoints.push_back(Point(0.0));  supportPoints.push_back(Point(4.0)); 
  FiniteElement *finiteElement = new FiniteElement(supportPoints, referenceElement);
  std::cout<<finiteElement->mapGlobal2Local(Point(1.0))<<finiteElement->mapLocal2Global(Point(0.0))<<std::endl;
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
  }


  if (false)
  { /** test basis functions */
  // create a basis of shape functions
  std::vector<Point> dummySupportPoints;
  dummySupportPoints.push_back(Point(0.0));
  dummySupportPoints.push_back(Point(0.5));
  //BasisFunctions *bf = new PiecewiseLinear(dummySupportPoints);
  //BasisFunctions *bf = new PiecewiseQuadratic(dummySupportPoints);
  BasisFunctions *bf = new PiecewiseCubic(dummySupportPoints);

  // fill vector of points
  std::vector<Point> p;
  const unsigned int np = 151;
  const unsigned int ndof = bf->getNumberOfNodes();
  for (unsigned int ip = 0; ip < np; ++ip)
    p.push_back(dummySupportPoints[0]+ip/double(np-1)*(dummySupportPoints[1]-dummySupportPoints[0]));

  // evaluate shape functions at points
  for (unsigned int ip = 0; ip < np; ++ip) {
    std::cout<<"x="<<p[ip].x<<"  ";
    for (unsigned int idof = 0; idof < ndof; ++idof) {
      std::cout<<"phi_"<<idof<<"="<<bf->getVal(idof, p[ip])<<"  ";
    }
    for (unsigned int idof = 0; idof < ndof; ++idof) {
      std::cout<<"DphiDx_"<<idof<<"="<<bf->getDx(idof, p[ip])<<"  ";
    }
    std::cout<<"\n";
  }

  std::vector<Point> nodes = bf->getNodes();
  unsigned int nn = bf->getNumberOfNodes();
  for (unsigned int i = 0; i < nn; ++i) std::cout<<i<<nodes[i]<<"\n"; 
  delete bf;
  }
  
  return 0;
}
