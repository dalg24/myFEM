#include "Quad.h"

Quad::
Quad(enum quad_type t) {
  if (t == TWOPT_GAUSS_ONED) {
    type = TWOPT_GAUSS_ONED;
    np = 2;
    w = new double[2];
    w[0] = 1.0 / 2.0;
    w[1] = 1.0 / 2.0;
    p = new Vertex[2];
    p[0].x = 0.5 * (1.0 - sqrt(1.0/3.0));
    p[1].x = 0.5 * (1.0 + sqrt(1.0/3.0));
  }

  if (t == MIDPOINT_TWOD) {
    type = MIDPOINT_TWOD;
    np = 3;
    w = new double[3];
    w[0] = 1.0 / 6.0;
    w[1] = 1.0 / 6.0;
    w[2] = 1.0 / 6.0;
    p = new Vertex[3];
    p[0].x = 0.0;  p[0].y = 0.5;
    p[1].x = 0.5;  p[1].y = 0.0;
    p[2].x = 0.5;  p[2].y = 0.5;
  }
}

Quad::
Quad(const Quad& q) : type(q.type),
                      np(q.np) {
  w = new double[np];
  p = new Vertex[np];
  for (int i = 0; i < np; i++) {
    w[i] = q.w[i];
    p[i] = q.p[i];
  }
}

Quad::
~Quad() {
  delete [] p;
  delete [] w;
}

Quad&
Quad::
operator=(const Quad& q) {
  if (this == &q) return *this;

  type = q.type;
  np = q.np;
  delete [] p;
  delete [] w;
  w = new double[np];
  p = new Vertex[np];
  for (int i = 0; i < np; i++) {
    w[i] = q.w[i];
    p[i] = q.p[i];
  }

  return *this;
}

int
Quad::
get_np() const {
  return np;
}

double
Quad::
get_w(const int i) const {
  return w[i];
}

Vertex
Quad::
get_p(const int i) const {
  return p[i];
}
