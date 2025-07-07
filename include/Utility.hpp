#pragma once
#include <cmath>
#include <type_traits>
#include <unordered_map>
#pragma pack(push, 1)
template <typename T> struct Vector3x {
  T x, y, z;

  // Default constructor (required for std::vector)
  Vector3x() : x(0), y(0), z(0) {}

  // Parameterized constructor
  Vector3x(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

  template <typename P> explicit Vector3x(const Vector3x<P> &l) {
    x = static_cast<T>(l.x);
    y = static_cast<T>(l.y);
    z = static_cast<T>(l.z);
  }
  template <typename P> Vector3x<T> &operator+=(const Vector3x<P> &rhs) {
    x += static_cast<T>(rhs.x);
    y += static_cast<T>(rhs.y);
    z += static_cast<T>(rhs.z);
    return *this;
  }
  template <typename P> Vector3x<T> &operator-=(const Vector3x<P> &rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
  }
  T &operator[](int idx) { return *(&x + idx); }
  T &operator*=(double rhs) {
    x *= rhs;
    y *= rhs;
    z *= rhs;
    return *this;
  }

  // Copy constructor (auto-generated is fine)
};
#pragma pack(pop) // Restore default packing
using Vector3d = Vector3x<double>;
using Vector3ld = Vector3x<long double>;
using Vector3i = Vector3x<int>;

template <typename T> struct Particle {

  // HINT: global represent type family
  //  0 - metal
  //  1- gas
  //  local for indexing in that families
  int m_type_id_global, m_type_id_local;
  Vector3x<T> pos, vel, force;
  Particle(int _type, int _id)
      : m_type_id_global(_type), m_type_id_local(_id), pos({0., 0., 0.}),
        vel({0., 0., 0.}), force({0., 0., 0.}) {}
};

template <typename T, typename P>
auto operator/=(Vector3x<P> &v,
                T t) -> std::enable_if_t<std::is_arithmetic_v<T>, void> {
  v.x /= t;
  v.y /= t;
  v.z /= t;
}
template <typename T, typename P>
auto operator*=(Vector3x<P> &v,
                T t) -> std::enable_if_t<std::is_arithmetic_v<T>, void> {
  v.x *= t;
  v.y *= t;
  v.z *= t;
}

template <typename T, typename P,
          std::enable_if_t<std::is_arithmetic<P>::value, int> = 0>
auto operator*(const Vector3x<T> &lhs,
               P val) -> Vector3x<std::common_type_t<T, P>> {
  return {lhs.x * val, lhs.y * val, lhs.z * val};
}

template <typename T, typename P,
          std::enable_if_t<std::is_arithmetic<P>::value, int> = 0>
auto operator/(const Vector3x<T> &lhs,
               P val) -> Vector3x<std::common_type_t<T, P>> {
  return {lhs.x / val, lhs.y / val, lhs.z / val};
}

template <typename T, typename P,
          std::enable_if_t<std::is_arithmetic<P>::value, double> = 0>
auto operator+(const Vector3x<T> &lhs,
               P val) -> Vector3x<std::common_type_t<T, P>> {
  return {lhs.x + val, lhs.y + val, lhs.z + val};
}

template <typename T, typename P>
auto operator-(const Vector3x<T> &lhs,
               const Vector3x<P> &rhs) -> Vector3x<std::common_type_t<T, P>> {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

template <typename T, typename P>
auto operator+(const Vector3x<T> &lhs,
               const Vector3x<P> &rhs) -> Vector3x<std::common_type_t<T, P>> {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

template <typename T>
Vector3x<T> operator*(const Vector3x<T> &lhs, const Vector3x<T> &rhs) {
  return {lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z};
}
template <typename T, typename P>
auto operator*(P multiplyier, const Vector3x<T> &lhs)
    -> std::enable_if_t<std::is_arithmetic_v<P>, Vector3x<T>> {
  return {lhs.x * multiplyier, lhs.y * multiplyier, lhs.z * multiplyier};
}

template <typename T> T length(const Vector3x<T> &t) {
  return sqrt(pow(t.x, 2) + pow(t.y, 2) + pow(t.z, 2));
}
// a[1]*b[2] - a[2]*b[1],   // x‐component: a_y*b_z − a_z*b_y
// a[2]*b[0] - a[0]*b[2],   // y‐component: a_z*b_x − a_x*b_z
// a[0]*b[1] - a[1]*b[0]

template <typename T, typename P>
auto crossProduct(const Vector3x<T> &f,
                  const Vector3x<P> &s) -> Vector3x<std::common_type_t<T, P>> {
  using common_type = std::common_type_t<T, P>;
  return Vector3x<common_type>{static_cast<common_type>(f.y * s.z - f.z * s.y),
                               static_cast<common_type>(f.z * s.x - f.x * s.z),
                               static_cast<common_type>(f.x * s.y - f.y * s.x)};
}
template <typename T, typename P>
auto dotProduct(const Vector3x<T> &f,
                const Vector3x<P> &s) -> std::common_type_t<T, P> {
  return f.x * s.x + f.y * s.y + f.z * s.z;
}
template <typename T,
          std::enable_if_t<std::is_arithmetic<T>::value, double> = 0>
struct ProcSubGrid3D {
  Vector3x<T> lo, hi;
  int rank;

  bool contains(const Vector3x<T> &pos) {
    return pos.x > lo.x && pos.x < hi.x && pos.y > lo.y && pos.y < hi.y &&
           pos.z > lo.z && pos.z < hi.z;
  }
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
struct SystemTraits {
  Vector3x<T> m_centr_mass_rad;
  Vector3x<T> m_centr_mass_vel;
  Vector3x<T> m_angular_momentum;
  Vector3x<T> diagonal_pressure_tensor;
  Vector3x<T> non_diagonal_pressure_tensor;

  /*
  due to storage potentials in vector, types, with whose coefficients and data
  they will work, will be changed to next_counter suppose for types 6 and 7
  setPotential has been called, then, type 6 gets index 0 and type 7 gets
  index 1.
  then, when system needs to be dumped in loggs, initial types's
  indexes have to be restored
  */
  std::unordered_map<int, int> m_types_mapping;
  T m_average_lattice;
  T m_cohesive_energy;
  T m_temperature;
  T m_kinetic_energy;
  T m_potential_energy;
  T m_volume;
  T m_pressure;
  T m_system_mass;
  T m_mass_of_the_free;

  T m_r_cut;
  T m_r_on;
  T m_timestep;
  Vector3i cluster_structure;
  int type_of_atomic_constaints;
  // TODO: for sake of accuracy, let it be double[3]
  double gravity[3];
};