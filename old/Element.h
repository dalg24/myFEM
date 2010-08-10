#ifndef ELEMENT_H
#define ELEMENT_H

#include <iostream>
#include "Node.h"
#include "Quad.h"
#include "RefElement.h"

using namespace std;

/** 
* @class Element
*
* @brief This is a class for finite elements.
*
*
* @author DaLG <lebrungd@neo.tamu.edu>
*/

class Element {

public:
  /// Constructor.
  Element(enum element_type,
          const int,
          ...);

  /// Copy constructor.
  Element(const Element&);

  /// Destructor.
  ~Element();

  ///
  Element& operator=(const Element&);

  ///
  int get_element_number() const;

  ///
  int get_node_per_element() const;

  ///
  Node get_node(const int) const;

  ///
  void info() const;

  ///
  friend ostream& operator<<(ostream&,
                             const Element&);

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

  /// Element #
  int element_number;

  /// Nodes per element
  int node_per_element;

  /// Nodes
  Node* node;

};

#endif // ELEMENT_H

