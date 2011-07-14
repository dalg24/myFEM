#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
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
  QuadratureRule() : type("notype") { }
  ~QuadratureRule() { }
  QuadratureRule(const QuadratureRule& qr) : points(qr.points), weights(qr.weights) { }
  QuadratureRule& operator=(const QuadratureRule& qr) { if (this==&qr) return *this; points = qr.points; weights = qr.weights; return *this; }

  unsigned int getNumberOfQuadraturePoints() const { assert(points.size()==weights.size()); return points.size(); }
  Point getQuadraturePoint(unsigned int iqp) const { return points[iqp]; }
  std::vector<Point> getQuadraturePoints() const { return points; }
  double getWeight(unsigned int iqp) const { return weights[iqp]; }
  std::vector<double> getWeights() const { return weights; }
  std::vector<Point> getSupportPoints() const { return sp; }
  std::string getType() const { return type; }

protected:
  std::vector<Point> points;
  std::vector<double> weights;
  std::vector<Point> sp;
  std::string type;
}; // end class QuadratureRule

class GaussianTwoPoints : public QuadratureRule {
public:
  GaussianTwoPoints() {
    points.push_back(Point(-1.0 / sqrt(3.0))); weights.push_back(1.0);
    points.push_back(Point(1.0 / sqrt(3.0))); weights.push_back(1.0);
    sp.push_back(Point(-1.0)); sp.push_back(Point(1.0));
    type = "GaussianTwoPoints";
  }
}; // end class GaussianTwoPoints

class GaussianThreePoints : public QuadratureRule {
public:
  GaussianThreePoints() {
    points.push_back(Point(-sqrt(3.0 / 5.0))); weights.push_back(5.0 / 9.0);
    points.push_back(Point(0.0)); weights.push_back(8.0 / 9.0);
    points.push_back(Point(sqrt(3.0 / 5.0))); weights.push_back(5.0 / 9.0);
    sp.push_back(Point(-1.0)); sp.push_back(Point(1.0));
    type = "GaussianThreePoints";
  }
}; // end class GaussianThreePoints

////////////////////////// SHAPE FUNCTIONS //////////////////////////////////////
class BasisFunctions {
public:
  BasisFunctions() : type("notype") { }
  ~BasisFunctions() { }
  BasisFunctions(const BasisFunctions& bf) : nodes(bf.nodes) { }
  BasisFunctions& operator=(const BasisFunctions& bf) { if (this==&bf) return *this; nodes = bf.nodes; return *this; }

  unsigned int getNumberOfNodes() const { return nodes.size(); }
  Point getNode(unsigned int idof) const { return nodes[idof]; }
  std::vector<Point> getNodes() const { return nodes; }
  std::string getType() const { return type; }

  std::vector<Point> getSupportPoints() const { std::vector<Point> sp; sp.push_back(nodes[0]); sp.push_back(nodes[1]); return sp;}

  virtual double getVal(unsigned int, Point) const = 0;
  virtual double getDx(unsigned int, Point) const = 0;

protected:
  std::vector<Point> nodes;
  std::string type;
}; // end class BasisFunctions

class PiecewiseLinear : public BasisFunctions {
public:             
  PiecewiseLinear(std::vector<Point> sp) { nodes = sp; type = "PiecewiseLinear"; }

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
  PiecewiseQuadratic(std::vector<Point> sp) { nodes = sp; nodes.push_back((sp[0] + sp[1]) / 2.0); type = "PiecewiseQuadratic"; }

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
  std::string getTypeOfBasisFunction() const { return bf->getType(); }

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

  void setDOF(unsigned int& dc) { 
    if ((dof.size() != 0) 
        &&(m.size() != 0)) {
      std::cerr<<"DOF have already been set before..."<<std::endl;
      abort();
    }
    unsigned int ndof = re->getNumberOfNodes();
    for (unsigned int idof = 0; idof < ndof; ++idof, ++dc) {
      dof.push_back(dc);
      m[idof] = dc;
    }
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
  double getDeterminantOfJacobian() const { std::vector<Point> ref_sp = re->getSupportPoints(); return (sp[1].x - sp[0].x) / (ref_sp[1].x - ref_sp[0].x); }
  double getValByLocal(unsigned int idof, Point l) const { return re->getVal(idof, l); }
  double getValByGlobal(unsigned int idof, Point g) const { Point l = mapGlobal2Local(g); return re->getVal(idof, l); }
  double getDxByLocal(unsigned int idof, Point l) const { return re->getDx(idof, l); }
  double getDxByGlobal(unsigned int idof, Point g) const { Point l = mapGlobal2Local(g); return re->getDx(idof, l); }
  unsigned int getNumberOfNodes() const { return re->getNumberOfNodes(); }
  unsigned int getNumberOfDOF() const { assert(dof.size() != 0); return dof.size(); }
  std::map<unsigned int, unsigned int> getDOFMap() const { assert(m.size() != 0); return m; }

protected:
  std::vector<Point> sp;
  std::vector<unsigned int> dof;
  std::map<unsigned int, unsigned int> m;
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
      } // end for iqp
    } // end for idof
  }

  double getDeterminantOfJacobianTimesWeight(unsigned int iqp) const { return Weights[iqp]*DeterminantOfJacobian; }
  double getShapeValue(unsigned int idof, unsigned int iqp) const { return ShapeValues[idof][iqp]; }
  double getShapeDx(unsigned int idof, unsigned int iqp) const { return ShapeDx[idof][iqp]; }
  std::vector<Point> getQuadraturePoints() const {
    std::vector<Point> qps;
    const unsigned int nqp = qr->getNumberOfQuadraturePoints();
    for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
      qps.push_back(fe->mapLocal2Global(qr->getQuadraturePoint(iqp)));
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

// solve -d/dx(a(x) d/dx(u)) + q(x) u(x) = f(x) for u(x) and x in (0,1)
double a(double x) { return 1 + x; };
double q(double x) { return x; };
double f(double x) { return x*cos(M_PI*x) + M_PI*sin(M_PI*x) + M_PI*M_PI*(1+x)*cos(M_PI*x); };
double a(Point p) { return a(p.x); };
double q(Point p) { return q(p.x); };
double f(Point p) { return f(p.x); };

int main(int argc, char *argv[]) {

  { /** nouveau test */
  unsigned int nel = 10;
  // construct nel elements
  ReferenceElement referenceElement("PiecewiseLinear");
  //ReferenceElement referenceElement("PiecewiseQuadratic");
  std::vector<FiniteElement> finiteElements;
  Point startPoint(0.0); Point endPoint(1.0);
  for (unsigned int iel; iel < nel; ++iel) {
    std::vector<Point> supportPoints;
    supportPoints.push_back(startPoint+double(iel)/double(nel)*(endPoint-startPoint));
    supportPoints.push_back(startPoint+double(iel+1)/double(nel)*(endPoint-startPoint));
    //std::cout<<iel<<"  "<<supportPoints[0]<<"  "<<supportPoints[1]<<"\n";
    finiteElements.push_back(FiniteElement(supportPoints, &referenceElement));
  } // end for iel

  // distribute the Degrees Of Freedoms
  // TODO: need to come up with something better than this
  unsigned int dofCounter = 0;
  for (unsigned int iel; iel < nel; ++iel) {
    finiteElements[iel].setDOF(dofCounter);
    dofCounter--;
  } // end for iel
  dofCounter++;

  std::vector<std::vector<double> > globalMatrix(dofCounter, std::vector<double>(dofCounter, 0.0));
  std::vector<double> globalRHS(dofCounter, 0.0);

  GaussianTwoPoints quadratureRule;
  UpdateFlags updateFlags;

  FEValues feValues(NULL, &quadratureRule, updateFlags); 
  for (unsigned int iel = 0; iel < nel; ++iel) {
    feValues.reinit(&finiteElements[iel]);

    std::vector<Point> quadraturePoints = feValues.getQuadraturePoints();

    const unsigned int nqp = quadraturePoints.size();
    std::cout<<iel<<"  ";
    for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
      std::cout<<"  "<<quadraturePoints[iqp]<<"  ";
    } // end for iqp
    std::cout<<"\n";

    std::map<unsigned int, unsigned int> dofMap = finiteElements[iel].getDOFMap();
    const unsigned int ndof = finiteElements[iel].getNumberOfDOF();
    /*for (unsigned int idof = 0; idof < ndof; ++idof) {
      std::cout<<idof<<" -> "<<dofMap[idof]<<"\n";
    } // end for idof*/

    std::vector<std::vector<double> > localMatrix(ndof, std::vector<double>(ndof, 0.0));
    for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
      for (unsigned int idof = 0; idof < ndof; ++idof) {
        for (unsigned int jdof = 0; jdof < ndof; ++jdof) {
          localMatrix[idof][jdof] += ( a(quadraturePoints[iqp]) * feValues.getShapeDx(idof, iqp) * feValues.getShapeDx(jdof, iqp)
                                       + q(quadraturePoints[iqp]) * feValues.getShapeValue(idof, iqp) * feValues.getShapeValue(jdof, iqp)
                                     ) * feValues.getDeterminantOfJacobianTimesWeight(iqp);
        } // end for jdof
      } // end for idof
    } //end for iqp

    std::vector<double> localRHS(ndof, 0.0);
    for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
      for (unsigned int idof = 0; idof < ndof; ++idof) {
        localRHS[idof] += f(quadraturePoints[iqp]) * feValues.getShapeValue(idof, iqp) * feValues.getDeterminantOfJacobianTimesWeight(iqp);
        /*
        std::cout<<iqp<<quadraturePoints[iqp]<<idof<<"  "
          <<"f="<<f(quadraturePoints[iqp])<<"  "
          <<"phi_"<<idof<<"="<<feValues.getShapeValue(idof, iqp)<<"  "
          <<"JxW="<<feValues.getDeterminantOfJacobianTimesWeight(iqp)<<"\n";
        */
      } // end for idof
    } //end for iqp

    for (unsigned int idof = 0; idof < ndof; ++idof) {
      for (unsigned int jdof = 0; jdof < ndof; ++jdof) {
        globalMatrix[dofMap[idof]][dofMap[jdof]] += localMatrix[idof][jdof];
        std::cout<<localMatrix[idof][jdof]<<"  ";
      } // end for jdof
      globalRHS[dofMap[idof]] += localRHS[idof];
      std::cout<<"||  "<<localRHS[idof]<<"\n";
    } // end for idof

  } // end for iel

  std::cout<<dofCounter<<"\n";

  std::fstream foutMatrix;
  std::fstream foutRHS;
  foutMatrix.open("matrix.dat", std::fstream::out);
  foutRHS.open("rhs.dat", std::fstream::out);
  for (unsigned int idof = 0; idof < dofCounter; ++idof) {
    for (unsigned int jdof = 0; jdof < dofCounter; ++jdof) {
      std::cout<<std::setw(7)<<std::setprecision(3)<<globalMatrix[idof][jdof]<<"  ";
      foutMatrix<<std::setw(9)<<std::setprecision(5)<<globalMatrix[idof][jdof]<<"  ";
    } // end for jdof
    std::cout<<"||  "<<std::setprecision(3)<<globalRHS[idof]<<"\n";
    foutMatrix<<"\n";
    foutRHS<<std::setprecision(5)<<globalRHS[idof]<<"\n";
  } // end for idof
    foutMatrix.close();
    foutRHS.close();

  std::cout<<"just checking\n";
  }


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

  { /** test point operations */
  Point *pp;
  pp = new Point(-3.0);
  std::cout<<*pp<<pp->x<<2.0*(*pp)<<*pp/3.0<<*pp+(*pp)<<*pp-(*pp)<<std::endl;
  delete pp;
  }

  { /** test basis functions */
  // create a basis of shape functions
  std::vector<Point> dummySupportPoints;
  dummySupportPoints.push_back(Point(-1.0));
  dummySupportPoints.push_back(Point(1.0));
  BasisFunctions *bf = new PiecewiseLinear(dummySupportPoints);
  //BasisFunctions *bf = new PiecewiseQuadratic(dummySupportPoints);

  // fill vector of points
  std::vector<Point> p;
  const unsigned int np = 51;
  const unsigned int ndof = bf->getNumberOfNodes();
  for (unsigned int ip = 0; ip < np; ++ip)
    p.push_back(Point(-1.0+ip*2.0/double(np-1)));

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

  delete bf;
  }
  
  return 0;
}
