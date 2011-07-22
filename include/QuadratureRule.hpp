#ifndef myFEM_QUADRATURE_RULE_HPP
#define myFEM_QUADRATURE_RULE_HPP

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#include <cmath>

#include "Point.hpp"

class QuadratureRule {
public:
  QuadratureRule() : points(NULL), weights(NULL), n_points(0) { }
  ~QuadratureRule() { 
    if (points != NULL) { delete [] points; points = NULL; } 
    if (weights != NULL) { delete [] weights; weights = NULL; } 
  }
  QuadratureRule(const QuadratureRule& quadrature_rule) {
    n_points = quadrature_rule.n_points;
    points = new double[n_points];
    weights = new double[n_points];
    std::copy(quadrature_rule.points, quadrature_rule.points + n_points, points);
    std::copy(quadrature_rule.weights, quadrature_rule.weights + n_points, weights);
sp = quadrature_rule.sp;
  }
  QuadratureRule& operator=(const QuadratureRule& quadrature_rule) { 
    if (this == &quadrature_rule) { return *this; }
    if (points != NULL) { delete [] points; }; 
    if (weights != NULL) { delete [] weights; }; 
    n_points = quadrature_rule.n_points;
    points = new double[n_points];
    weights = new double[n_points];
    std::copy(quadrature_rule.points, quadrature_rule.points + n_points, points);
    std::copy(quadrature_rule.weights, quadrature_rule.weights + n_points, weights);
sp = quadrature_rule.sp;
    return *this; 
  }

  unsigned short getNumberOfQuadraturePoints() const { return n_points; }
  Point getQuadraturePoint(unsigned short i) const { assert(i < n_points); return Point(points[i]); }
  std::vector<Point> getQuadraturePoints() const { 
    std::vector<Point> p(n_points);
for (unsigned int i = 0; i < n_points; ++i) {
  p[i] = Point(points[i]);
}
    return p;
  }
  double getWeight(unsigned short i) const { assert(i < n_points); return weights[i]; }
  std::vector<double> getWeights() const {
    std::vector<double> w(weights, weights + n_points);
    return w;
  }

std::vector<Point> getSupportPoints() const { return sp; }

  virtual std::string getType() const = 0;
  virtual unsigned short getMaxOrder() const = 0;

protected:
  double *points;
  double *weights;
  unsigned short n_points;
std::vector<Point> sp;
}; // end class QuadratureRule

class GaussianTwoPoints : public QuadratureRule {
public:
  GaussianTwoPoints() {
    n_points = 2;
    points = new double[2];
    weights = new double[2];
    points[0] = - 1.0 / sqrt(3.0);
    weights[0] = 1.0;
    points[1] = + 1.0 / sqrt(3.0); 
    weights[1] = 1.0;
sp.push_back(Point(-1.0)); sp.push_back(Point(1.0));
  }
  std::string getType() const { return "GaussianTwoPoints"; }
  unsigned short getMaxOrder() const { return 3; };
}; // end class GaussianTwoPoints

class GaussianThreePoints : public QuadratureRule {
public:
  GaussianThreePoints() {
    n_points = 3;
    points = new double[3];
    weights = new double[3];
    points[0] = 0.0; 
    weights[0] = 8.0 / 9.0;
    points[1] = - sqrt(3.0 / 5.0); 
    weights[1] = 5.0 / 9.0;
    points[2] = + sqrt(3.0 / 5.0); 
    weights[2] = 5.0 / 9.0;
sp.push_back(Point(-1.0)); sp.push_back(Point(1.0));
  }
  std::string getType() const { return "GaussianThreePoints"; }
  unsigned short getMaxOrder() const { return 5; };
}; // end class GaussianThreePoints

class GaussianFourPoints : public QuadratureRule {
public:
  GaussianFourPoints() {
    n_points = 4;
    points = new double[4];
    weights = new double[4];
    points[0] = - sqrt((3.0 - 2.0 * sqrt(6.0/5.0)) / 7.0); 
    weights[0] = (18.0 + sqrt(30.0)) / 36.0;
    points[1] = + sqrt((3.0 - 2.0 * sqrt(6.0/5.0)) / 7.0); 
    weights[1] = (18.0 + sqrt(30.0)) / 36.0;
    points[2] = - sqrt((3.0 + 2.0 * sqrt(6.0/5.0)) / 7.0); 
    weights[2] = (18.0 - sqrt(30.0)) / 36.0;
    points[3] = + sqrt((3.0 + 2.0 * sqrt(6.0/5.0)) / 7.0); 
    weights[3] = (18.0 - sqrt(30.0)) / 36.0;
sp.push_back(Point(-1.0)); sp.push_back(Point(1.0));
  }
  std::string getType() const { return "GaussianFourPoints"; }
  unsigned short getMaxOrder() const { return 7; };
}; // end class GaussianFourPoints

class GaussianFivePoints : public QuadratureRule {
public:
  GaussianFivePoints() {
    n_points = 5;
    points = new double[5];
    weights = new double[5];
    points[0] = 0.0; 
    weights[0] = 128.0 / 225.0;
    points[1] = - 1.0 / 3.0 * sqrt((5.0 - 2.0 * sqrt(10.0/7.0))); 
    weights[1] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
    points[2] = + 1.0 / 3.0 * sqrt((5.0 - 2.0 * sqrt(10.0/7.0))); 
    weights[2] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
    points[3] = - 1.0 / 3.0 * sqrt((5.0 + 2.0 * sqrt(10.0/7.0))); 
    weights[3] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
    points[4] = + 1.0 / 3.0 * sqrt((5.0 + 2.0 * sqrt(10.0/7.0))); 
    weights[4] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
sp.push_back(Point(-1.0)); sp.push_back(Point(1.0));
  }
  std::string getType() const { return "GaussianFivePoints"; }
  unsigned short getMaxOrder() const { return 9; };
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
    assert(false);
  }
  return qr;
}

QuadratureRule* makeNewPointerToQuadratureRule(unsigned short nqp) {
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
    assert(false);
  }
  return qr;
}

#endif // myFEM_QUADRATURE_RULE_HPP
