#include <Utility.hpp>
#include <cmath>
#include <optional>
#include <type_traits>

template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
class SimulationBox {

#if 0
  // HINT: lower corner
  Vector3x<T> origin;

  // HINT: Lattice vectors spanning the cell:
  //    the three edges of your (possibly triclinic) parallelepiped.
  Vector3x<T> a, b, c;
#endif
  bool pbc[3];

public:
  SimulationBox(const Vector3d &origin_, const Vector3d &a_, const Vector3d &b_,
                const Vector3d &c_);
  //: origin(origin_), a(a_), b(b_), c(c_) {}
  /// Convert a Cartesian position `x` into fractional coords f in [0,1)³,
  /// solving   x = origin + f.x*a1 + f.y*a2 + f.z*a3
  Vector3x<T> toFrac(const Vector3ld &x) const;

  /// Convert fractional coords `f` back to Cartesian:
  ///   x = origin + f.x*a1 + f.y*a2 + f.z*a3
  Vector3x<T> toCart(const Vector3ld &x) const;

  Vector3x<T> applyPBC(const Vector3ld &x) const;

  /// Compute cell volume = |a1·(a2×a3)|
  // T volume() const { return std::abs(dotProduct(a, crossProduct(b, c))); }

  /// Given a desired cutoff+skin, compute how many cell‐list bins
  /// to split each axis into (floor only if you want integer counts).
  /// You can precompute these once per build.
  std::array<int, 3> cellListGrid(double cutoff_skin) const;

  /// Compute the minimum image displacement between two points
  /// under PBC (wraps Δ = x2−x1 into the [−0.5,0.5)³ box).
  Vector3x<T> minImage(const Vector3x<T> &x1, const Vector3x<T> &x2) const;

  std::optional<Vector3x<T>> getFractionalCoordinates();
};
