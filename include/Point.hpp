#ifndef myFEM_POINT_HPP
#define myFEM_POINT_HPP

#include <iostream>
#include <cassert>

class Point { 
public:
  Point(double X = 0.0) : x(X) { }
  ~Point() { }
  Point(const Point& p) : x(p.x) { }
  Point& operator=(const Point& p) { if (this==&p) return *this; x = p.x; return *this; }
  friend Point operator+(const Point& l, const Point& r) { Point p; p.x = l.x + r.x; return p;}
  friend Point operator-(const Point& l, const Point& r) { Point p; p.x = l.x - r.x; return p;}
  friend Point operator/(const Point& l, const double r) { assert(r != 0); Point p; p.x = l.x / r; return p;}
  friend Point operator*(const double l, const Point& r) { Point p; p.x = l * r.x; return p;}
  friend Point operator*(const Point& l, const double r) { Point p; p.x = l.x * r; return p;}
  friend std::ostream& operator<<(std::ostream& os, const Point& p) { os<<"( "<<p.x<<" )"; return os; }
  double x; 
}; // end class Point
                                            
#endif // myFEM_POINT_HPP
