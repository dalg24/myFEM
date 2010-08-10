#include "Element.h"

Element::
Element(enum element_type t,
        const int elem_numb,
        ...) : type(t),
               element_number(elem_numb) {
  if (t == TRIANGLE_LINEAR) {
    node_per_element = 3;
  }
  if (t == TRIANGLE_QUADRATIC) {
    node_per_element = 6;
  }
  if (t == ONED_LINEAR) {
    node_per_element = 2;
  }
  if (t == ONED_QUADRATIC) {
    node_per_element = 3;
  }
  node = new Node[node_per_element];
  va_list vl;
  va_start(vl, node_per_element);
  for (int i = 0; i < node_per_element, i++) {
    node[i] = va_arg(vl, Node);
  }
  va_end(vl);
}

Element::
Element(const Element& elem) : type(elem.type),
                               reference_element(elem.reference_element),
                               quadrature(elem.quadrature),
                               element_number(elem.element_number),
                               node_per_element(elem.node_per_element) {
  node =  new Node[node_per_element];
  for (int i = 0; i < node_per_element, i++) {
    node[i] = elem.node[i];
  }
}

Element::
~Element() {
delete [] node;
}

Element::
Element&
operator=(const Element& elem) {
  if (this == &elem) return *this;

  type = elem.type;
  reference_element = elem.reference_element;
  quadrature = elem.quadrature_element;
  element_number = elem.element_number;
  node_per_element = elem.node_per_element;
  node = new Node[node_per_element];
  for (int i = 0; i < node_per_element, i++) {
    node[i] = elem.node[i];
  }

  return *this;
}

Element::
int
get_element_number() const {
  return element_number;
}

Element::
int
get_node_per_element() const {
  return node_per_element;
}

Element::
Node
get_node(const int i) const {
  return node[i];
}

Element::
void
info() const {
  cout << "print some info about the element" << endl;
}
