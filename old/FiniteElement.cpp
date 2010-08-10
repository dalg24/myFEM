#include "FiniteElement.h"

FiniteElement::
FiniteElement(const RefElement* ref_elem,
              const Quad* quad,
              const BasisFunc* basis_func) : reference_element(ref_elem),
                                             quadrature(quad),
                                             basis_function(basis_func) {
}

FiniteElement::
FiniteElement(const FiniteElement& fe) : reference_element(fe.reference_element),
                                         quadrature(fe.quadrature),
                                         basis_function(fe.basis_function) {
}

FiniteElement::
~FiniteElement() {
}

FiniteElement::
FiniteElement&
operator=(const FiniteElement& fe) {
  if (this == &fe) return *this;

  reference_element = fe.reference_element;
  quadrature = fe.quadrature_element;
  basis_function = fe.basis_function;

  return *this;
}

int
FiniteElement::
get_number_of_elements() const {
  return nelem;
}

void
FiniteElement::
info() const {
  cout << "print some info about the finite elements " << endl;
}
