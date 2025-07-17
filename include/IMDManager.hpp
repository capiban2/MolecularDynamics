#pragma once
#include <Utility.hpp>
#include <type_traits>
#include <vector>
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class IMDManager {
public:
  virtual T getCurrentPressure() const noexcept = 0;

  // HINT: mainly for thermostats
  virtual T getThermostatGroupTemperature(int) const noexcept = 0;
  virtual const Vector3x<T> &getThermostatGroupVelCom(int) const noexcept = 0;
  virtual int getThermostatGroupSize(int) const noexcept = 0;

  virtual std::vector<Particle<T>> &getParticles() const noexcept = 0;

  virtual T getTimeStep() const noexcept = 0;
  virtual void updateForces() = 0;

  virtual int getFirstGhostIndex() const noexcept = 0;
  virtual Vector3x<T> getExternalForce() const noexcept = 0;

  // HINT: this function returns true for first iteration after particles
  // redistribution for afterwards using that information in potentials
  virtual bool isFirstAfterRebuilt() const noexcept = 0;
};