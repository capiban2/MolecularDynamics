#pragma once
#include <Utility.hpp>
#include <type_traits>
#include <vector>
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class IMDManager {
public:
  virtual T getCurrentPressure() const noexcept = 0;
  virtual T getCurrentTemperature() const noexcept = 0;
  virtual std::vector<Particle<T>> &getParticles() const noexcept = 0;

  virtual T getTimeStep() const noexcept = 0;
};