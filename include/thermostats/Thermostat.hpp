#pragma once
#include <IMDManager.hpp>
#include <Utility.hpp>
#include <memory>
#include <type_traits>
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class Thermostat {
protected:
  std::unique_ptr<IMDManager<T>> m_parent;

public:
  enum class Type { Unitary, Andersen, Berendsen, Langevin, NoseHoover };
  Thermostat(IMDManager<T> *mdm) : m_parent(mdm) {}
  virtual void applyThermostat() = 0;
};
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class UnitaryThermostat final : public Thermostat<T> {

public:
  UnitaryThermostat(IMDManager<T> *mdm) : Thermostat<T>(mdm) {}
  virtual void applyThermostat() override {}
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>,
          typename... Args>
std::unique_ptr<Thermostat<T>> createThermostat(typename Thermostat<T>::Type &,
                                                Args &&...args);