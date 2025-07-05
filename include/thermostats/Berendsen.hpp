#pragma once
#include "Thermostat.hpp"

template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class BerendsenThermostat final : public Thermostat<T> {
  T m_ratio;

public:
  BerendsenThermostat(IMDManager<T> *mdm, T ratio,
                      const std::vector<std::pair<int, T>> &_preserve_temp)
      : Thermostat<T>(mdm, _preserve_temp), m_ratio(ratio) {}
  virtual void applyThermostat() override;
};