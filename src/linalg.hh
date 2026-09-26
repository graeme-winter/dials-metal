// Three-vectors, three-by-three matrices, and a small dense solver.
//
// Hand-rolled rather than pulled in, for the same reason the spot finder
// hand-rolled its msgpack writer: this has to build on a laptop and on a
// beamline machine with nothing installed, and the total content here is
// about three hundred lines. Eigen would be better if there were more of it.
//
// Everything is double. The reciprocal-space work sets positions to a
// millipixel over a 4000-pixel detector, which is one part in 4e6, and float
// has 7 digits. That is not enough margin to be casual about.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <ostream>

namespace mxi {

struct Vec3 {
  double x = 0.0, y = 0.0, z = 0.0;

  constexpr Vec3() = default;
  constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  double operator[](std::size_t i) const { return (&x)[i]; }
  double &operator[](std::size_t i) { return (&x)[i]; }

  Vec3 operator+(const Vec3 &o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3 operator-(const Vec3 &o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3 operator-() const { return {-x, -y, -z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
  Vec3 &operator+=(const Vec3 &o) {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  Vec3 &operator-=(const Vec3 &o) {
    x -= o.x;
    y -= o.y;
    z -= o.z;
    return *this;
  }

  double dot(const Vec3 &o) const { return x * o.x + y * o.y + z * o.z; }
  Vec3 cross(const Vec3 &o) const {
    return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
  }
  double norm_squared() const { return dot(*this); }
  double norm() const { return std::sqrt(norm_squared()); }
  Vec3 normalized() const {
    const double n = norm();
    return n > 0.0 ? *this / n : Vec3{};
  }
};

inline Vec3 operator*(double s, const Vec3 &v) { return v * s; }

inline std::ostream &operator<<(std::ostream &os, const Vec3 &v) {
  return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

// Row-major 3x3. Row-major because that is how DIALS serialises a matrix into
// a flat array of nine, and a transpose hidden in the I/O layer is the single
// easiest way to lose a day.
struct Mat3 {
  std::array<double, 9> m{};

  constexpr Mat3() = default;
  Mat3(std::initializer_list<double> values) {
    std::size_t i = 0;
    for (double v : values) {
      if (i < 9)
        m[i++] = v;
    }
  }

  static Mat3 identity() { return {1, 0, 0, 0, 1, 0, 0, 0, 1}; }

  // Basis vectors as the *rows*, which is the real-space convention: a crystal
  // whose real_space_a/b/c are the rows has A = inverse of that matrix, with
  // a*, b*, c* as the columns.
  static Mat3 from_rows(const Vec3 &r0, const Vec3 &r1, const Vec3 &r2) {
    return {r0.x, r0.y, r0.z, r1.x, r1.y, r1.z, r2.x, r2.y, r2.z};
  }
  static Mat3 from_columns(const Vec3 &c0, const Vec3 &c1, const Vec3 &c2) {
    return {c0.x, c1.x, c2.x, c0.y, c1.y, c2.y, c0.z, c1.z, c2.z};
  }

  double operator()(std::size_t r, std::size_t c) const { return m[3 * r + c]; }
  double &operator()(std::size_t r, std::size_t c) { return m[3 * r + c]; }

  Vec3 row(std::size_t r) const {
    return {m[3 * r], m[3 * r + 1], m[3 * r + 2]};
  }
  Vec3 column(std::size_t c) const { return {m[c], m[3 + c], m[6 + c]}; }

  Vec3 operator*(const Vec3 &v) const {
    return {m[0] * v.x + m[1] * v.y + m[2] * v.z,
            m[3] * v.x + m[4] * v.y + m[5] * v.z,
            m[6] * v.x + m[7] * v.y + m[8] * v.z};
  }

  Mat3 operator*(const Mat3 &o) const {
    Mat3 r;
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j) {
        double s = 0.0;
        for (std::size_t k = 0; k < 3; ++k)
          s += (*this)(i, k) * o(k, j);
        r(i, j) = s;
      }
    return r;
  }

  Mat3 operator*(double s) const {
    Mat3 r = *this;
    for (double &v : r.m)
      v *= s;
    return r;
  }
  Mat3 operator-(const Mat3 &o) const {
    Mat3 r;
    for (std::size_t i = 0; i < 9; ++i)
      r.m[i] = m[i] - o.m[i];
    return r;
  }

  Mat3 transpose() const {
    return {m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]};
  }

  double determinant() const {
    return m[0] * (m[4] * m[8] - m[5] * m[7]) -
           m[1] * (m[3] * m[8] - m[5] * m[6]) +
           m[2] * (m[3] * m[7] - m[4] * m[6]);
  }

  // Returns the identity and sets `ok` false on a singular matrix rather than
  // returning infinities: a caller that ignores the flag then gets a wrong
  // answer that looks wrong, instead of NaNs that propagate silently.
  Mat3 inverse(bool *ok = nullptr) const {
    const double d = determinant();
    if (ok)
      *ok = std::abs(d) > 1e-30;
    if (std::abs(d) <= 1e-30)
      return identity();
    const double i = 1.0 / d;
    return {(m[4] * m[8] - m[5] * m[7]) * i, (m[2] * m[7] - m[1] * m[8]) * i,
            (m[1] * m[5] - m[2] * m[4]) * i, (m[5] * m[6] - m[3] * m[8]) * i,
            (m[0] * m[8] - m[2] * m[6]) * i, (m[2] * m[3] - m[0] * m[5]) * i,
            (m[3] * m[7] - m[4] * m[6]) * i, (m[1] * m[6] - m[0] * m[7]) * i,
            (m[0] * m[4] - m[1] * m[3]) * i};
  }
};

// Rotation by `angle` radians about `axis`, right-handed. Rodrigues' formula
// written out rather than built from a quaternion: it is used in the inner
// loop of prediction and this form avoids the round trip.
inline Mat3 rotation(const Vec3 &axis, double angle) {
  const Vec3 u = axis.normalized();
  const double c = std::cos(angle), s = std::sin(angle), t = 1.0 - c;
  return {
      t * u.x * u.x + c,       t * u.x * u.y - s * u.z, t * u.x * u.z + s * u.y,
      t * u.x * u.y + s * u.z, t * u.y * u.y + c,       t * u.y * u.z - s * u.x,
      t * u.x * u.z - s * u.y, t * u.y * u.z + s * u.x, t * u.z * u.z + c};
}

// Angle of a rotation matrix, in radians, in [0, pi]. Clamped before acos
// because accumulated round-off puts the trace a few ULP outside the valid
// range often enough to matter, and acos of 1+1e-16 is NaN.
inline double rotation_angle(const Mat3 &r) {
  const double cosine = 0.5 * (r(0, 0) + r(1, 1) + r(2, 2) - 1.0);
  return std::acos(std::fmax(-1.0, std::fmin(1.0, cosine)));
}

// Solve a symmetric positive-definite system by Cholesky, in place.
// Returns false if the matrix is not positive definite, which for a normal
// matrix means the problem is rank deficient and the caller must not proceed.
bool solve_spd(double *a, double *b, std::size_t n);

} // namespace mxi
