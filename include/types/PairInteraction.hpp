#pragma once
#include <cmath>
#include <type_traits>
namespace MD {
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
struct alignas(64) PairInteraction {

  int i, j;
  T r_sq, dx, dy, dz;
  PairInteraction() = default;

  PairInteraction(int particle_i, int particle_j, T distance_sq, T delta_x,
                  T delta_y, T delta_z)
      : i(particle_i), j(particle_j), r_sq(distance_sq), dx(delta_x),
        dy(delta_y), dz(delta_z) {}

  T distance() const { return std::sqrt(r_sq); }
  T distance_squared() const { return r_sq; }
};

} // namespace MD