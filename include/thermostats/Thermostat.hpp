#pragma once
#include <IMDManager.hpp>
#include <Utility.hpp>
#include <memory>
#include <type_traits>
#include <utility>
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class Thermostat {
protected:
  std::unique_ptr<IMDManager<T>> m_parent;

  /*
HINT: this vector contains pairs (int, T),
where int is a count if iterations, and
T is the temperature that needs to be preserve.
 Every iteration this counter decreases, then, if it is zero
next pair is going to be in work.
 If no pair left, then last pair's temperature is a temperature, that will
 be preserved all the way to the end of the modelling
For doing so, if no pairs left and counter is a zero, then counter gets
increase to some big arbitary value.
*/
  std::vector<std::pair<int, T>> m_time_nd_temp_to_preserve;

  T m_time_step;

public:
  enum class Type { Unitary, Andersen, Berendsen, Langevin, NoseHoover };
  Thermostat(IMDManager<T> *mdm,
             const std::vector<std::pair<int, T>> &_preserve_temp)
      : m_parent(mdm), m_time_nd_temp_to_preserve(_preserve_temp) {
    m_time_step = m_parent->getTimeStep();
  }
  virtual void applyThermostat() = 0;

  // HINT: becauase dynamic target temperature has been supported, then, it may
  // return different target temperatures, that depend on current timestemp
  virtual T getTargetTemperature() const noexcept;
};
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class UnitaryThermostat final : public Thermostat<T> {

public:
  UnitaryThermostat(IMDManager<T> *mdm) : Thermostat<T>(mdm) {}
  virtual void applyThermostat() override {}
  virtual T getTargetTemperature() const noexcept override {
    return this->m_parent->getCurrentTemperature();
  }
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>,
          typename... Args>
std::unique_ptr<Thermostat<T>> createThermostat(typename Thermostat<T>::Type &,
                                                Args &&...args);