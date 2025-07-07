#pragma once
#include <Utility.hpp>
#include <random>
#include <type_traits>
#include <vector>
namespace Position {

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
std::vector<Vector3x<T>> genRandomPositions(int count,
                                            const Vector3x<T> &box_size) {
  std::vector<Vector3x<T>> __res;
  __res.reserve(count);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::array<std::uniform_real_distribution<T>, 3> dist(
      {std::uniform_real_distribution<T>(0.0, box_size.x),
       std::uniform_real_distribution<T>(0.0, box_size.y),
       std::uniform_real_distribution<T>(0.0, box_size.z)});
  for (int t_ = 0; t_ != count; ++t_) {
    __res.emplace_back(dist[0](gen), dist[1](gen), dist[2](gen));
  }
  return __res;
}

} // namespace Position