#ifndef myFEM_BASIS_FUNCTIONS_HPP
#define myFEM_BASIS_FUNCTIONS_HPP

#include <string>
#include <vector>

#include "Point.hpp"

class BasisFunctions {
public:
  BasisFunctions() : o(0) { }
  ~BasisFunctions() { }
  BasisFunctions(const BasisFunctions& bf) : n(bf.n), o(bf.o) { }
  BasisFunctions& operator=(const BasisFunctions& bf) { if (this==&bf) return *this; n = bf.n; o = bf.o; return *this; }

  unsigned short getNumberOfNodes() const { return n.size(); }
  const Point& getNode(unsigned short idof) const { return n[idof]; }
  const std::vector<Point>& getNodes() const { return n; }
  const unsigned short& getOrder() const { return o; }

  std::vector<Point> getSupportPoints() const { std::vector<Point> sp; sp.push_back(n[0]); sp.push_back(n[1]); return sp;}

  virtual double getVal(unsigned short, const Point&) const = 0;
  virtual double getDx(unsigned short, const Point&) const = 0;
  virtual std::string getType() const = 0;

protected:
  std::vector<Point> n;
  unsigned short o;
}; // end class BasisFunctions

class PiecewisePolynomial : public BasisFunctions {
public:             
  PiecewisePolynomial(unsigned short u, std::vector<Point> sp) { 
    assert(sp.size() == 2); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(0) / double(u)); 
    n.push_back(sp[0] + (sp[1] - sp[0]) * double(u) / double(u)); 
    for (unsigned int a = 1; a < u; ++a) {
      n.push_back(sp[0] + (sp[1] - sp[0]) * double(a) / double(u)); 
    } // end for
    o = u;
    assert(o == n.size() - 1);
  }
  std::string getType() const { return "PiecewisePolynomial"; }

  double getVal(unsigned short idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x); assert(idof < n.size());
    double value = 1.0;
    for (unsigned int a = 0; a < n.size(); ++a) {
      if (a != idof) {
        value *= (n[a].x - p.x) / (n[a].x - n[idof].x);
      } // end if a
    } // end for a
    return value;
  }

  double getDx(unsigned short idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x); assert(idof < n.size());
    double value = 0.0;
//    static unsigned int k = 0;
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
    o = 1;
  }
  std::string getType() const { return "PiecewiseLinear"; }

  double getVal(unsigned short idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x);
    double value;
    if (idof == 0) {
      value = (n[1].x - p.x) / (n[1].x - n[0].x); 
    } else if (idof == 1) {
      value = (p.x - n[0].x) / (n[1].x - n[0].x); 
    } else {
      assert(false);
    }
    return value;
  }

  double getDx(unsigned short idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x);
    double value;
    if (idof == 0) {
      value = -1.0 / (n[1].x - n[0].x); 
    } else if (idof == 1) {
      value = 1.0 / (n[1].x - n[0].x); 
    } else {
      assert(false);
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
    o = 2;
  }
  std::string getType() const { return "PiecewiseQuadratic"; }

  double getVal(unsigned short idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x);
    double value;
    if (idof == 0) {
      value = (n[1].x - p.x) * (n[2].x - p.x) / ((n[1].x - n[0].x) * (n[2].x - n[0].x)); 
    } else if (idof == 1) {
      value = (n[2].x - p.x) * (n[0].x - p.x) / ((n[2].x - n[1].x) * (n[0].x - n[1].x)); 
    } else if (idof == 2) {
      value = (n[0].x - p.x) * (n[1].x - p.x) / ((n[0].x - n[2].x) * (n[1].x - n[2].x)); 
    } else {
      assert(false);
    }
    return value;
  }

  double getDx(unsigned short idof, const Point& p) const {
    assert(n[0].x <= p.x); assert(p.x <= n[1].x);
    double value;
    if (idof == 0) {
      value = (2.0 * p.x - n[1].x - n[2].x) / ((n[1].x - n[0].x) * (n[2].x - n[0].x)); 
    } else if (idof == 1) {
      value = (2.0 * p.x - n[2].x - n[0].x) / ((n[2].x - n[1].x) * (n[0].x - n[1].x)); 
    } else if (idof == 2) {
      value = (2.0 * p.x - n[0].x - n[1].x) / ((n[0].x - n[2].x) * (n[1].x - n[2].x)); 
    } else {
      assert(false);
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
    o = 3;
  }
  std::string getType() const { return "PiecewiseCubic"; }

  double getVal(unsigned short idof, const Point& p) const {
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
      assert(false);
    }
    return value;
  }

  double getDx(unsigned short idof, const Point& p) const {
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
      assert(false);
    }
    return value;
  }
}; // end class PiecewiseCubic
#endif // myFEM_BASIS_FUNCTIONS_HPP
