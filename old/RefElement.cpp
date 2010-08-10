#include "RefElement.h"

RefElement::
RefElement(const element_type t) : type(t) {
}

RefElement::
RefElement(const RefElement& re) : type(re.type),
				   nv(re.nv) {
  for (int i = 0, i < nv, i++) {
    lv[i] = re.lv[i];
    psi[i] = re.psi[i];
  }
}

RefElement::
~RefElement() {
}

RefElement&
RefElement::
operator=(const RefElement& re) {
  if(this==&re) return *this;

  type = re.type;
  nv = re.nv;
  for (int i = 0, i < nv, i++) {
    lv[i] = re.lv[i];
    psi[i] = re.psi[i];
  }

  return *this;
}

int
RefElement::
get_nv() const {
  return nv;
}

point
RefElement::
get_lv(const int i) const {
  return lv[i];
}

shapefunction
RefElement::
get_psi(const int i) const {
  return psi[i];
}

void
RefElement::
info() const {
  cout << "RefElement: ";

  cout << "(nv=" << nv;
  cout << ",...=" << nv <<")";
}


ostream&
operator<<(ostream& os,
           const RefElement& re) {
  os << "RefElement: ";

  os << "(nv=" << re.nv;
  os << ",...=" << re.nv << ")";

  return os;
}
