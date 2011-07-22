#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <cstring>
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
#include "Point.hpp"
////////////////////////// QUADRATURE RULE //////////////////////////////////////
#include "QuadratureRule.hpp"
////////////////////////// SHAPE FUNCTIONS //////////////////////////////////////
#include "BasisFunctions.hpp"
////////////////////////// REFERENCE ELEMENT //////////////////////////////////////
#include "ReferenceElement.hpp"
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

  std::vector<Point> getNodes() const { return std::vector<Point>(); }

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
  if (normType == myFEM_H1_NORM) { assert(StiffnessMatrix != std::vector<std::vector<double> >()); }
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
  if (normType != myFEM_L1_NORM) { Norm = sqrt(Norm); }
  return Norm;
}

void print_version() {
  std::cout<<"myFEM - Version 0.0\nCopyright 2011 (c) Damien L-G\n\n";
}

void print_help() {
  std::cout<<"myFEM - Version 0.0\nCopyright 2011 (c) Damien L-G\n\n";

  std::cout<<"Usage: myFEM [--help] [--version] <options> [<args>] \n\n";
  std::cout<<"The following options are available:\n";
  std::cout<<"-h = this help\n";
  std::cout<<"-v = verbose output\n";
  std::cout<<"-n = number of elements  default value is 10.\n";
  std::cout<<"-q = quadrature rule     default value is 2 for Gaussian two points.\n";
  std::cout<<"                         for now 5 is the maximum value allowed.\n";
  std::cout<<"                         yields exact integration for polynomials of\n";
  std::cout<<"                         degree 2n - 1 or less.\n";
  std::cout<<"-p = polynomial order    default value is 1 for piecewise linear\n";
  std::cout<<"                         basis shape functions.\n";
  std::cout<<"                         set to 2 for quadratics, 3 for cubics, etc.\n";
  std::cout<<"                         quadrature rule must be set accordingly so that\n";
  std::cout<<"                         numerical integration is exact.\n";
  std::cout<<"\n";
}

bool check_input_arguments(int n_elements, int n_quadrature_points, int polynomial_order) {
  if (n_elements < 1) {
    return false;
  }
  if ((n_quadrature_points < 2)
      || (n_quadrature_points > 5)) {
    return false;
  }
  if ((polynomial_order < 1)
      || (polynomial_order > 6)) {
    return false;
  }
  if (2 * n_quadrature_points - 1 < 2 * polynomial_order) {
    std::cerr<<"myFEM: the order of the quadrature rule is too low!\n";
    abort();
  }
  if (n_elements > 100000) {
    std::cerr<<"myFEM: "<<n_elements<<" is a bit ambitious.\n";
    abort();
  }
  return true;
}

int main(int argc, char *argv[]) {
  std::cout<<"*********************************************\n";
  std::cout<<"* Running: "; for (int i = 0; i < argc; ++i) std::cout<<argv[i]<<" "; std::cout<<"\n";
  std::cout<<"*********************************************\n";
  { /** nouveau test */
  // Default values
  int n_elements = 10;
  int n_quadrature_points = 2;
  int polynomial_order = 1;
  bool verbose = false;

  for (int i = 1; i < argc; ++i) {
    if ((strcmp(argv[i], "-h") == 0)
        || (strcmp(argv[i], "--help") == 0)) {
      print_help();
      return 0;
    } else if (strcmp(argv[i], "--version") == 0) {
      print_version();
      return 0;
    } else if (strcmp(argv[i], "-n") == 0) {
      if (++i >= argc) { 
        print_help(); 
        return -1;
      } else {
        n_elements = atoi(argv[i]);
      }
    } else if (strcmp(argv[i], "-q") == 0) {
      if (++i >= argc) { 
        print_help(); 
        return -1;
      } else {
        n_quadrature_points = atoi(argv[i]);
      }
    } else if (strcmp(argv[i], "-p") == 0) {
      if (++i >= argc) { 
        print_help(); 
        return -1;
      } else {
        polynomial_order = atoi(argv[i]);
      }
    } else if (strcmp(argv[i], "-v") == 0) {
      verbose = true;
    } else {
      std::cerr<<"myFEM: '"<<argv[i]<<"' is not a valid option. See './myFEM --help'."<<std::endl;
      return -1;
    }
  } // end for i

  if (!check_input_arguments(n_elements, 
                             n_quadrature_points, 
                             polynomial_order)) {
    print_help();
    return -1;
  }

  try {
  // run problem with input arguments
  }
  catch(Exception_t &e) {
    // handle exception
  }
  unsigned int nel = n_elements;
  unsigned int numberOfGaussPoints = n_quadrature_points;
  unsigned int polynomialOrder = polynomial_order;


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
  std::vector<double> solutionVector = solveMatrixTimesXEqualsRHS(globalMatrix, globalRHS);

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
  std::cout<<"error (absolute) L1 Norm = "<<computeNorm(MassMatrix, exactErrorVector, myFEM_L1_NORM)<<"  relative percent = "<<100.0*computeNorm(MassMatrix, exactErrorVector, myFEM_L1_NORM)/computeNorm(MassMatrix, exactSolutionVector, myFEM_L1_NORM)<<"\n"; 
  std::cout<<"error L2 Norm = "<<computeNorm(MassMatrix, exactErrorVector, myFEM_L2_NORM)<<"  relative percent = "<<100.0*computeNorm(MassMatrix, exactErrorVector, myFEM_L2_NORM)/computeNorm(MassMatrix, exactSolutionVector, myFEM_L2_NORM)<<"\n";  
  std::cout<<"error H1 Norm = "<<computeNorm(MassMatrix, exactErrorVector, myFEM_H1_NORM, StiffnessMatrix)<<"  relative percent = "<<100.0*computeNorm(MassMatrix, exactErrorVector, myFEM_H1_NORM, StiffnessMatrix)/computeNorm(MassMatrix, exactSolutionVector, myFEM_H1_NORM, StiffnessMatrix)<<"\n";

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
  
  std::cout<<"*********************************************\n";
  std::cout<<"* Done running: "; for (int i = 0; i < argc; ++i) std::cout<<argv[i]<<" "; std::cout<<"\n";
  std::cout<<"*********************************************\n";
  return 0;
}
