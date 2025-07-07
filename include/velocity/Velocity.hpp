#pragma once
#include <Utility.hpp>
#include <random>
#include <type_traits>
#include <vector>
namespace Velocity {
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
T get_uniform(void) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<T> dis(1e-10, 1.0 - 1e-10); // Avoid 0
  return dis(gen);
}
template <typename T>
std::vector<Vector3x<T>> genMaxwellVelocity(T temp, T mass, int count) {
  constexpr T KB = 1.38064852;
  constexpr long double PI = 3.14159265358979323846264338327950280;
  std::vector<Vector3x<T>> __res(count);
  long double variance = sqrt(KB * temp / mass);

  long double _r1 = 0.0, _r2 = 0.0;
  int c_ = 0;
  for (; count > 0; --count) {
    _r1 = get_uniform<T>(), _r2 = get_uniform<T>();
    __res.at(c_).x = variance * (sqrt(-2 * log(_r1)) * cos(2 * PI * _r2));
    __res.at(c_).y = variance * (sqrt(-2 * log(_r2)) * cos(2 * PI * _r1));
    __res.at(c_).z = variance * (sqrt(-2 * log(get_uniform<T>())) *
                                 cos(2 * PI * get_uniform<T>()));
  }

  return __res;
}

} // namespace Velocity