#ifndef REFELEMENT_H
#define REFELEMENT_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

/** 
* @class RefElement
*
* @brief This is a class for reference finite elements.
*
*
* @author DaLG <lebrungd@neo.tamu.edu>
*/

class RefElement {

public:
  /// Constructor.
  RefElement(const element_type);

  /// Copy constructor.
  RefElement(const RefElement&);

  /// Destructor.
  ~RefElement();

  RefElement& operator=(const RefElement&);

  int get_nv() const;

  elem_type get_type() const;

  point get_lv(const int) const;

  shapefunction get_psi(const int) const;

  /// This is a function to
  void info() const;

  friend ostream& operator<<(ostream&,
                             const RefElement&);

private:
  /// Type of element.
  enum element_type {
    /// Triangles with piecewise linear basis functions.
    TRIANGLE_LINEAR,
    /// Triangles with piecewise quadratic basis functions.
    TRIANGLE_QUADRATIC,
    /// One-dimensional element with piecewise piecewise linear basis functions.
    ONED_LINEAR,
    /// One-dimensional element with piecewise piecewise quadratic basis functions.
    ONED_QUADRATIC
  } type;
  /// Number of vertices per reference element.
  int nv; 
  /// Vertices (local coordinate).
  point *lv;
  /// Basis functions.
  shapefunction *psi;

};

#endif // REFELEMENT_H

