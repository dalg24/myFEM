#ifndef FINITEELEMENT_H
#define FINITEELEMENT_H

#include <iostream>
#include "Node.h"
#include "Quad.h"
#include "RefElement.h"
#include "Element"

using namespace std;

/** 
* @class FiniteElement
*
* @brief This is a class for finite elements.
*
*
* @author DaLG <lebrungd@neo.tamu.edu>
*/

class FiniteElement {

public:
  /// Constructor.
  FiniteElement(const RefElement*,
                const Quad*,
                const int,
                const int,
                ...);

  /// Copy constructor.
  FiniteElement(const FiniteElement&);

  /// Destructor.
  ~FiniteElement();

  ///
  FiniteElement& operator=(const FiniteElement&);

  ///
  int get_element_number() const;

  ///
  int get_node_per_element() const;

  ///
  Node get_node(const int) const;

  /// This is a function that prints to the sreen some information about the finite elements.
  void info() const;

private:
  /// Number of elements.
  int nelem;

  /// Pointer to pointer to elements.
  Element** elements;

  /// Pointer to reference element.
  RefElement* reference_element;

  /// Pointer to quadrature for integration.
  Quad* quadrature;

  /// Pointer to basis function.
  ShapeFunc* basis_func;

  /// Number of external functions.
  int nfunc;

  /// External functions.
  double (**ext_func)(double x, double y);
};

#endif // FINITEELEMENT_H

