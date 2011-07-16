#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include "umfpack.h"

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
  std::string getTypeOfBasisFunctions() const { return bf->getType(); }

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
  double getJacobian() const { std::vector<Point> ref_sp = re->getSupportPoints(); return (sp[1].x - sp[0].x) / (ref_sp[1].x - ref_sp[0].x); }
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

std::vector<std::vector<double> > computeLocalMatrix(FEValues *feValues, 
    bool verbose = false, 
    std::ostream& os = std::cout,
    unsigned int debugLevel = 1,
    double (*pa)(Point) = &a,
    double (*pq)(Point) = &q) {
  /** get quadrature points */
  std::vector<Point> quadraturePoints = feValues->getQuadraturePoints();
  const unsigned int nqp = quadraturePoints.size();
  if ((verbose)
      && (debugLevel > 3)) {
    os<<"quadrature points\n";
    for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
      os<<iqp<<" = "<<quadraturePoints[iqp]<<"\n";
    } // end for iqp
  } // end if verbose

  /** compute local matrix */
  const unsigned int ndof = feValues->getFiniteElement()->getNumberOfDOF();
  std::vector<std::vector<double> > localMatrix(ndof, std::vector<double>(ndof, 0.0));
  for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
    for (unsigned int idof = 0; idof < ndof; ++idof) {
      for (unsigned int jdof = 0; jdof < ndof; ++jdof) {
        localMatrix[idof][jdof] += ( pa(quadraturePoints[iqp]) * feValues->getShapeDx(idof, iqp) * feValues->getShapeDx(jdof, iqp)
                                     + pq(quadraturePoints[iqp]) * feValues->getShapeValue(idof, iqp) * feValues->getShapeValue(jdof, iqp)
                                   ) * feValues->getDeterminantOfJacobianTimesWeight(iqp);
        if ((verbose)
            && (debugLevel > 6)) {
          os<<"iqp="<<iqp<<quadraturePoints[iqp]<<"idof="<<idof<<"  jdof="<<jdof<<"  "
            <<"a="<<a(quadraturePoints[iqp])<<"  "
            <<"q="<<q(quadraturePoints[iqp])<<"  "
            <<"phi_"<<idof<<"="<<feValues->getShapeValue(idof, iqp)<<"  "
            <<"phi_"<<jdof<<"="<<feValues->getShapeValue(jdof, iqp)<<"  "
            <<"DphiDx_"<<idof<<"="<<feValues->getShapeDx(idof, iqp)<<"  "
            <<"DphiDx_"<<jdof<<"="<<feValues->getShapeDx(jdof, iqp)<<"  "
            <<"JxW="<<feValues->getDeterminantOfJacobianTimesWeight(iqp)<<"\n";
        } // end if verbose
      } // end for jdof
    } // end for idof
  } //end for iqp           

  if (verbose) {
    for (unsigned int idof = 0; idof < ndof; ++idof) {
      for (unsigned int jdof = 0; jdof < ndof; ++jdof) {
        os<<std::setw(7)<<std::setprecision(3)<<localMatrix[idof][jdof]<<"  ";
      } // end for jdof
      std::cout<<"\n";
    } // end for idof
  } // end if verbose

  return localMatrix;
}

std::vector<double> computeLocalRHS(FEValues *feValues, 
    bool verbose = false, 
    std::ostream& os = std::cout,
    unsigned int debugLevel = 1,
    double (*pf)(Point) = &f) {
  /** get quadrature points */
  std::vector<Point> quadraturePoints = feValues->getQuadraturePoints();
  const unsigned int nqp = quadraturePoints.size();
  if ((verbose)
      && (debugLevel > 3)) {
    os<<"quadrature points\n";
    for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
      os<<iqp<<" = "<<quadraturePoints[iqp]<<"\n";
    } // end for iqp
  } // end if verbose

  /** compute local RHS */
  const unsigned int ndof = feValues->getFiniteElement()->getNumberOfDOF();
  std::vector<double> localRHS(ndof, 0.0);
  for (unsigned int iqp = 0; iqp < nqp; ++iqp) {
    for (unsigned int idof = 0; idof < ndof; ++idof) {
      localRHS[idof] += pf(quadraturePoints[iqp]) * feValues->getShapeValue(idof, iqp) * feValues->getDeterminantOfJacobianTimesWeight(iqp);
      if ((verbose)
          && (debugLevel > 6)) {
        os<<"iqp="<<iqp<<quadraturePoints[iqp]<<"idof="<<idof<<"  "
          <<"f="<<f(quadraturePoints[iqp])<<"  "
          <<"phi_"<<idof<<"="<<feValues->getShapeValue(idof, iqp)<<"  "
          <<"JxW="<<feValues->getDeterminantOfJacobianTimesWeight(iqp)<<"\n";
      } //end if verbose
    } // end for idof
  } //end for iqp
  if (verbose) {
    for (unsigned int idof = 0; idof < ndof; ++idof) {
      os<<std::setprecision(3)<<localRHS[idof]<<"\n";
    } // end for idof
  } //end if verbose

  return localRHS;
}

enum DebugLevel_t { myFEM_NO_DEBUG, myFEM_MIN_DEBUG, myFEM_MED_DEBUG, myFEM_MAX_DEBUG };
enum Object_t { myFEM_BOTH_MATRIX_AND_VECTOR, myFEM_MATRIX_ONLY, myFEM_VECTOR_ONLY };
void distributeLocal2Global(FEValues *feValues,
    const std::vector<std::vector<double> >& localMatrix,
    const std::vector<double>& localVector,
    std::vector<std::vector<double> >& globalMatrix,
    std::vector<double>& globalVector,
    Object_t handleObject = myFEM_BOTH_MATRIX_AND_VECTOR,
    DebugLevel_t debugLevel = myFEM_NO_DEBUG) {
  // get number of DOF
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
    } // end if not matrix only
    if (debugLevel == myFEM_MAX_DEBUG) {
      if (handleObject != myFEM_VECTOR_ONLY) {
        for (unsigned int ildof = 0; ildof < nldof; ++ildof) {
          assert(localMatrix[ildof].size() == nldof);
        } // end for ildof
        for (unsigned int igdof = 0; igdof < ngdof; ++igdof) {
          assert(globalMatrix[igdof].size() == ngdof);
        } // end for igdof
      } // end if not vector only
    } // end if debugLevel hardcore
  } // end if debugLevel not zero
  // get DOF map
  std::map<unsigned int, unsigned int> dofMap = feValues->getFiniteElement()->getDOFMap();
  if (debugLevel > myFEM_NO_DEBUG) {
    std::cout<<"DOF map\n";
    for (unsigned int ildof = 0; ildof < nldof; ++ildof) {
      unsigned int igdof = dofMap[ildof];
      std::cout<<ildof<<" -> "<<igdof<<"\n";
    } // end for ildof
  } // end if debugLevel not zero
  // distribute local to global
  for (unsigned int ildof = 0; ildof < nldof; ++ildof) {
    unsigned int igdof = dofMap[ildof];
    if (handleObject != myFEM_VECTOR_ONLY) {
      for (unsigned int jldof = 0; jldof < nldof; ++jldof) {
        unsigned int jgdof = dofMap[jldof];
        globalMatrix[igdof][jgdof] += localMatrix[ildof][jldof];
      } // end for jldof
    } // end if not vector only
    if (handleObject != myFEM_MATRIX_ONLY) {
      globalVector[igdof] += localVector[ildof];
    } // end for ildof
  } // end if not matrix only
}

void printMatrixAndVector(const std::vector<std::vector<double> >& Matrix,
    const std::vector<double>& Vector,
    Object_t handleObject = myFEM_BOTH_MATRIX_AND_VECTOR,
    std::ostream& os = std::cout,
    int setwVal = 7,
    int setprecisionVal = 3,
    DebugLevel_t debugLevel = myFEM_NO_DEBUG) {
  // get number of DOF
  unsigned int ndof = 0;
  if (handleObject != myFEM_MATRIX_ONLY) {
    ndof = Vector.size();
  } else {
    ndof = Matrix.size();
  } // end if not matrix only
  if (debugLevel > myFEM_NO_DEBUG) {
    assert(ndof != 0);
    if (handleObject == myFEM_BOTH_MATRIX_AND_VECTOR) {
      assert(Vector.size() == Matrix.size());
    } // end if both matrix and vector
    if (debugLevel == myFEM_MAX_DEBUG) {
      if (handleObject != myFEM_VECTOR_ONLY) {
        for (unsigned int idof = 0; idof < ndof; ++idof) {
          assert(Matrix[idof].size() == ndof);
        } // end for idof
      } // end if not vector only
    } // end if debugLevel hardcore
  } // end if debugLevel not zero
  // print the stuff
  for (unsigned int idof = 0; idof < ndof; ++idof) {
    if (handleObject != myFEM_VECTOR_ONLY) {
      for (unsigned int jdof = 0; jdof < ndof; ++jdof) {
        os<<std::setw(setwVal)
          <<std::setprecision(setprecisionVal)
          <<Matrix[idof][jdof]<<"  ";
      } // end for jdof
    } //end if not vector only
    if (handleObject == myFEM_BOTH_MATRIX_AND_VECTOR) {
      os<<"||  ";
    } // end if both matrix and vector
    if (handleObject != myFEM_MATRIX_ONLY) {
      os<<std::setw(setwVal)
        <<std::setprecision(setprecisionVal)
        <<Vector[idof]<<"\n";
    } else {
      os<<"\n";
    } // end if not matrix only
  } // end for idof
}
                    
// TODO: replace prefix myFEM_ by namespace
enum Norm_t { myFEM_L1, myFEM_L2, myFEM_H1 };
double computeNorm(const std::vector<std::vector<double> >& MassMatrix,
    const std::vector<double>& Vector,
    Norm_t normType = myFEM_L2,
    const std::vector<std::vector<double> >& StiffnessMatrix = std::vector<std::vector<double> >()) {
  double Norm = 0.0;
  unsigned int ndof = Vector.size();
  for (unsigned int idof = 0; idof < ndof; ++idof) {
    for (unsigned int jdof = 0; jdof < ndof; ++jdof) {
      if (normType == myFEM_L1) {
        Norm += MassMatrix[idof][jdof] * fabs(Vector[jdof]);
      } else {
        Norm += MassMatrix[idof][jdof] * Vector[jdof] * Vector[jdof];
        if (normType == myFEM_H1) {
          Norm += StiffnessMatrix[idof][jdof] * Vector[jdof] * Vector[jdof];
        } // end if H1
      } // end if not L1
    } // end for jdof
  } // end for idof
  return Norm;
}

int main(int argc, char *argv[]) {

  { /** nouveau test */
  unsigned int nel = 10;
  std::string basisType = "PiecewiseLinear";
  std::string quadType = "GaussianThreePoints";
  //TODO: make command line argument passing smarter than this
  std::string errorMessage = std::string("Error: Bad command line arguments\n")
                           + std::string("Usage: argv[0] number_of_elements type_of_basis_functions type_of_quadrature_rule\n")
                           + std::string("Where: number_of_elements must be a positive integer\n")
                           + std::string("       type_of_basis_functions must be 1 for PiecewiseLinear or 2 for PiecewiseQuadratici\n")
                           + std::string("       type_of_quadrature_rule must be 2 for GaussianTwoPoints or 3 for GaussianThreePoints\n");
  for (int i = 0; i < argc; ++i) {
    std::cout<<argv[i]<<" ";
  }
  std::cout<<"\n";

  if (argc > 1) {
    nel = atoi(argv[1]);
  }
  if (argc > 2) {
    if (atoi(argv[2]) == 1) { 
      basisType = "PiecewiseLinear"; 
    } else if (atoi(argv[2]) == 2) {
      basisType = "PiecewiseQuadratic"; 
    } else {
      std::cerr<<errorMessage<<std::endl;
      abort();
    }
  }
  if (argc > 3) {
    if (atoi(argv[3]) == 2) { 
      quadType = "GaussianTwoPoints"; 
    } else if (atoi(argv[3]) == 3) {
      quadType = "GaussianThreePoints"; 
    } else {
      std::cerr<<errorMessage<<std::endl;
      abort();
    }
  }

  //std::cout<<basisType<<std::endl;
  //std::cout<<quadType<<std::endl;
  //abort();

  std::cout<<"#### BEGIN ######\n";
  std::cout<<"number of elements = "<<nel<<"\n";
  // construct nel elements
  ReferenceElement *referenceElement = new ReferenceElement(basisType);
  std::cout<<"basis functions = "<<referenceElement->getTypeOfBasisFunctions()<<"\n";

  std::vector<FiniteElement> finiteElements;
  Point startPoint(0.0); Point endPoint(1.0);
  for (unsigned int iel; iel < nel; ++iel) {
    std::vector<Point> supportPoints;
    supportPoints.push_back(startPoint+double(iel)/double(nel)*(endPoint-startPoint));
    supportPoints.push_back(startPoint+double(iel+1)/double(nel)*(endPoint-startPoint));
    //std::cout<<iel<<"  "<<supportPoints[0]<<"  "<<supportPoints[1]<<"\n";
    finiteElements.push_back(FiniteElement(supportPoints, referenceElement));
  } // end for iel

  // distribute the Degrees Of Freedoms (DOF)
  // TODO: need to come up with something better than this
  unsigned int dofCounter = 0;
  for (unsigned int iel; iel < nel; ++iel) {
    finiteElements[iel].setDOF(dofCounter);
    dofCounter--;
  } // end for iel
  dofCounter++;

  // create global matrix and global rhs
  std::vector<std::vector<double> > globalMatrix(dofCounter, std::vector<double>(dofCounter, 0.0));
  std::vector<std::vector<double> > gStiffnessMatrix(dofCounter, std::vector<double>(dofCounter, 0.0));
  std::vector<std::vector<double> > gMassMatrix(dofCounter, std::vector<double>(dofCounter, 0.0));
  std::vector<double> globalRHS(dofCounter, 0.0);
  std::vector<std::vector<double> > nullMatrix;
  std::vector<double> nullVector;

  QuadratureRule *quadratureRule;
  if (quadType == "GaussianTwoPoints") {
    //std::cout<<"Iwasthere"<<std::endl; 
    //abort();
    quadratureRule = new GaussianTwoPoints;
  } else if (quadType == "GaussianThreePoints") {
    quadratureRule = new GaussianThreePoints;
  } else {
    std::cerr<<"pb with qr"<<std::endl;
    abort();
  }
  std::cout<<"quadrature rule = "<<quadratureRule->getType()<<"\n";
  UpdateFlags updateFlags;

  FEValues feValues(NULL, quadratureRule, updateFlags); 
  for (unsigned int iel = 0; iel < nel; ++iel) {
  std::cout<<"##########\n";
    std::cout<<"element #"<<iel<<"\n";
    // compute FE values for element iel
    feValues.reinit(&finiteElements[iel]);

    // compute local matrix
    std::vector<std::vector<double> > localMatrix = computeLocalMatrix(&feValues, false, std::cout, 5);
    std::vector<std::vector<double> > lStiffnessMatrix = computeLocalMatrix(&feValues, false, std::cout, 5, &aStiffness, &qStiffness);
    std::vector<std::vector<double> > lMassMatrix = computeLocalMatrix(&feValues, false, std::cout, 5, &aMass, &aStiffness);

    // compute local rhs
    std::vector<double> localRHS = computeLocalRHS(&feValues, false, std::cout, 5);

    /** get DOF map */
    std::map<unsigned int, unsigned int> dofMap = feValues.getFiniteElement()->getDOFMap();
    const unsigned int ndof = feValues.getFiniteElement()->getNumberOfDOF();
    std::cout<<"DOF map\n";
    for (unsigned int idof = 0; idof < ndof; ++idof) {
      std::cout<<idof<<" -> "<<dofMap[idof]<<"\n";
    } // end for idof

    // distribute local to global
    distributeLocal2Global(&feValues, localMatrix, localRHS, globalMatrix, globalRHS, myFEM_BOTH_MATRIX_AND_VECTOR, myFEM_MAX_DEBUG);
    distributeLocal2Global(&feValues, lStiffnessMatrix, nullVector, gStiffnessMatrix, nullVector, myFEM_MATRIX_ONLY);
    distributeLocal2Global(&feValues, lMassMatrix, nullVector, gMassMatrix, nullVector, myFEM_MATRIX_ONLY);
    /*
    for (unsigned int idof = 0; idof < ndof; ++idof) {
      for (unsigned int jdof = 0; jdof < ndof; ++jdof) {
        globalMatrix[dofMap[idof]][dofMap[jdof]] += localMatrix[idof][jdof];
        gStiffnessMatrix[dofMap[idof]][dofMap[jdof]] += lStiffnessMatrix[idof][jdof];
        gMassMatrix[dofMap[idof]][dofMap[jdof]] += lMassMatrix[idof][jdof];
      } // end for jdof
      globalRHS[dofMap[idof]] += localRHS[idof];
    } // end for idof
    */
    std::cout<<"local matrix and local rhs\n";
    printMatrixAndVector(localMatrix, localRHS, myFEM_BOTH_MATRIX_AND_VECTOR, std::cout);

  } // end for iel

  std::cout<<"##########\n";
  std::cout<<"global number of DOF = "<<dofCounter<<"\n";

  // print out global matrix and global rhs
  std::cout<<"global matrix and global rhs\n";
  printMatrixAndVector(globalMatrix, globalRHS, myFEM_BOTH_MATRIX_AND_VECTOR, std::cout);

  // print matrix and RHS to files to check with matlab
  std::fstream foutMatrix;
  std::fstream foutRHS;
  foutMatrix.open("matrix.dat", std::fstream::out);
  foutRHS.open("rhs.dat", std::fstream::out);
  printMatrixAndVector(globalMatrix, nullVector, myFEM_MATRIX_ONLY, foutMatrix, 9, 5, myFEM_MAX_DEBUG);
  printMatrixAndVector(nullMatrix, globalRHS, myFEM_VECTOR_ONLY, foutRHS, 9, 5, myFEM_MAX_DEBUG);
  foutMatrix.close();
  foutRHS.close();

  //std::cout<<"stiffness matrix\n";
  //printMatrixAndVector(gStiffnessMatrix, nullVector, myFEM_MATRIX_ONLY, std::cout);
  //std::cout<<"mass matrix\n";
  //printMatrixAndVector(gMassMatrix, nullVector, myFEM_MATRIX_ONLY, std::cout);

#ifdef EBILE
  std::cout<<"tu te crois malin hein? gros ebile :D\n";
#endif

  std::cout<<"##########\n";
  std::cout<<"solve using umfpack\n";
  enum UMFPACK_STATUS_DUMMY_ENUM { TRIPLET_TO_COL, SYMBOLIC, NUMERIC, SOLVE, UMFPACK_STATUS_DUMMY_ENUM_SIZE };
  int status[UMFPACK_STATUS_DUMMY_ENUM_SIZE];

  double Info[UMFPACK_INFO];
  double Control[UMFPACK_CONTROL];
  void *Symbolic;
  void *Numeric;

  std::vector<double> iRowTMP, jColTMP, ijValTMP;
  int nRow = dofCounter;
  int nCol = dofCounter;
  double EPSILON = 1.0e-10;
  for (unsigned int i = 0; i < nRow; ++i) {
    for (unsigned int j = 0; j < nCol; ++j) {
      if (fabs(globalMatrix[i][j]) > EPSILON) { 
        //static int k = 0;
        //std::cout<<++k<<" -> i="<<i<<" j="<<j<<" val="<<globalMatrix[i][j]<<std::endl;
        iRowTMP.push_back(i);
        jColTMP.push_back(j);
        ijValTMP.push_back(globalMatrix[i][j]);
      } // end if nonzero val in global matrix
    } // end for j
  } // end for i
  int nNZ = ijValTMP.size();
  std::cout<<"number of nonzero values = "<<nNZ<<"\n";
  std::cout<<std::endl;
  int iRowCOO[nNZ];
  int jColCOO[nNZ];
  double ijValCOO[nNZ];
  std::copy(iRowTMP.begin(), iRowTMP.end(), iRowCOO);
  std::copy(jColTMP.begin(), jColTMP.end(), jColCOO);
  std::copy(ijValTMP.begin(), ijValTMP.end(), ijValCOO);

  double solVec[dofCounter];
  double rhsVec[dofCounter];
  std::copy(globalRHS.begin(), globalRHS.end(), rhsVec);

  int pColCC[nCol+1];
  int iRowCC[nNZ];
  double iValCC[nNZ];
  int Map[nNZ];
  int sys = UMFPACK_A;

  status[TRIPLET_TO_COL] = umfpack_di_triplet_to_col(nRow, nCol, nNZ, iRowCOO, jColCOO, ijValCOO, pColCC, iRowCC, iValCC, Map);
  status[SYMBOLIC] = umfpack_di_symbolic(nRow, nCol, pColCC, iRowCC, iValCC, &Symbolic, Control, Info);
  status[NUMERIC] = umfpack_di_numeric(pColCC, iRowCC, iValCC, Symbolic, &Numeric, Control, Info);
  status[SOLVE] = umfpack_di_solve(sys, pColCC, iRowCC, iValCC, solVec, rhsVec, Numeric, Control, Info);

  std::vector<double> solVecTMP(dofCounter);
  std::copy(solVec, solVec + dofCounter, solVecTMP.begin());
  //std::cout<<"solution\n";
  //printMatrixAndVector(std::vector<std::vector<double> >(), solVecTMP, VECTOR_ONLY, std::cout);
  std::cout<<"solution L2 Norm = "<<computeNorm(gMassMatrix, solVecTMP, myFEM_L2)<<"\n"; 
  
  /** compute exact error */
  // TODO: come up with something better than this
  std::vector<double> exactSolVec(dofCounter);
  for (unsigned int iel = 0; iel < nel; ++iel) {
    std::vector<Point> supportPoints;
    supportPoints.push_back(startPoint+double(iel)/double(nel)*(endPoint-startPoint));
    supportPoints.push_back(startPoint+double(iel+1)/double(nel)*(endPoint-startPoint));
    //std::cout<<iel<<"  "<<supportPoints[0]<<"  "<<supportPoints[1]<<"\n";
    if (referenceElement->getTypeOfBasisFunctions() == "PiecewiseLinear") {
      exactSolVec[iel] = u(supportPoints[0]);
      exactSolVec[iel+1] = u(supportPoints[1]);
    } else if (referenceElement->getTypeOfBasisFunctions() == "PiecewiseQuadratic") {
      exactSolVec[iel] = u(supportPoints[0]);
      exactSolVec[iel+1] = u(supportPoints[1]);
      exactSolVec[iel+2] = u((supportPoints[0] + supportPoints[1]) / 2.0);
    } else {
      std::cerr<<"mais qu'est-ce que tu fais la?\n";
      abort();
    }
  } // end for iel
  std::cout<<"exact solution L2 Norm = "<<computeNorm(gMassMatrix, exactSolVec, myFEM_L2)<<"\n"; 

  std::vector<double> errVec;
  for (unsigned int idof = 0; idof < dofCounter; ++idof) {
    errVec.push_back(solVec[idof] - exactSolVec[idof]);
  } // end for idof
  //std::cout<<"exact error\n";
  //printMatrixAndVector(std::vector<std::vector<double> >(), errVec, VECTOR_ONLY, std::cout);

  double errL2Norm = computeNorm(gMassMatrix, errVec, myFEM_L2);
  std::cout<<"error L1 Norm = "<<computeNorm(gMassMatrix, errVec, myFEM_L1)<<"\n"; 
  std::cout<<"error L2 Norm = "<<computeNorm(gMassMatrix, errVec, myFEM_L2)<<"\n"; 
  std::cout<<"error H1 Norm = "<<computeNorm(gMassMatrix, errVec, myFEM_H1, gStiffnessMatrix)<<"\n"; 

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
  { /** test point operations */
  Point *pp;
  pp = new Point(-3.0);
  std::cout<<*pp<<pp->x<<2.0*(*pp)<<*pp/3.0<<*pp+(*pp)<<*pp-(*pp)<<std::endl;
  delete pp;
  }

  if (false)
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
