#pragma once
#include <Utility.hpp>
#include <type_traits>
#include <vector>
template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class IMDManager {
public:
  T getPressure() const noexcept = 0;
  T getTemperature() const noexcept = 0;
  std::vector<Particle<T>> &getParticles() const noexcept = 0;
};