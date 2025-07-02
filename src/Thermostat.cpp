#include <Thermostat.hpp>

template <> void LangevinThermostat<double>::applyThermostat() {}
template <> void LangevinThermostat<long double>::applyThermostat() {}