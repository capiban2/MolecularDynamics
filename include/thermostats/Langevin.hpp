#pragma once
#include "Thermostat.hpp"
template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class LangevinThermostat final : public Thermostat<T> {

public:
  // TODO: add necessary parameters
  LangevinThermostat(IMDManager<T> *mdm,
                     const std::vector<std::pair<int, T>> &_preserve_temp)
      : Thermostat<T>(mdm, _preserve_temp) {}
  virtual void applyThermostat() override;
};