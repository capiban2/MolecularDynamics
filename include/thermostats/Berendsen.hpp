#pragma once
#include "Thermostat.hpp"

template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class BerendsenThermostat final : public Thermostat<T> {
  T m_ratio;
  T m_tar_temp;

public:
  BerendsenThermostat(IMDManager<T> *mdm, T ratio, T target_temp)
      : Thermostat<T>(mdm), m_ratio(ratio), m_tar_temp(target_temp) {}
  virtual void applyThermostat() override;
};