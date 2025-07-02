#pragma once
#include <barostats/Andersen.hpp>
#include <barostats/Berendsen.hpp>
#include <barostats/Langevin.hpp>
#include <barostats/MTK.hpp>
#include <barostats/ParinelloRahman.hpp>
#include <memory>
#include <stdexcept>
#include <type_traits>

#include <thermostats/Andersen.hpp>
#include <thermostats/Berendsen.hpp>
#include <thermostats/Langevin.hpp>
#include <thermostats/NoseHoover.hpp>

#include <potentials/EAM.hpp>
#include <potentials/LJ.hpp>

template <typename T, typename... Args>
std::unique_ptr<Barostat<T>> createBarostat(typename Barostat<T>::Type &type,
                                            Args &&...args) {

  switch (type) {

  case (Barostat<T>::Type::Berendsen):
    return std::make_unique<BerendsenBarostat<T>>(std::forward<Args>(args)...);

  case (Barostat<T>::Type::Andersen):
    return std::make_unique<AndersenBarostat<T>>(std::forward<Args>(args)...);

  case (Barostat<T>::Type::Langevin):
    return std::make_unique<LangevinBarostat<T>>(std::forward<Args>(args)...);

  case (Barostat<T>::Type::MTK):
    return std::make_unique<MTKBarostat<T>>(std::forward<Args>(args)...);

  case (Barostat<T>::Type::PR):
    return std::make_unique<ParinelloRahmanBarostat<T>>(
        std::forward<Args>(args)...);
  default:
    throw std::runtime_error("There's no such barostat!");
  }
}

template <typename T, typename... Args>
std::unique_ptr<Barostat<T>>
createThermostat(typename Thermostat<T>::Type &type, Args &&...args) {

  switch (type) {

  case (Thermostat<T>::Type::Berendsen):
    return std::make_unique<BerendsenThermostat<T>>(
        std::forward<Args>(args)...);

  case (Thermostat<T>::Type::Andersen):
    return std::make_unique<AndersenThermostat<T>>(std::forward<Args>(args)...);

  case (Thermostat<T>::Type::Langevin):
    return std::make_unique<LangevinThermostat<T>>(std::forward<Args>(args)...);

  case (Thermostat<T>::Type::NoseHoover):
    return std::make_unique<NoseHooverThermostat<T>>(
        std::forward<Args>(args)...);

  default:
    throw std::runtime_error("There's no such thermostat!");
  }
}

template <typename T, typename... Args>
std::shared_ptr<Potential<T>> createPotential(typename Potential<T>::Type &type,
                                              Args &&...args) {

  switch (type) {
  case (Potential<T>::Type::EAM):
    return std::make_shared<EAMPotential<T>>(std::forward<Args>(args)...);
  case (Potential<T>::Type::LJ):
    return std::make_shared<LJPotential<T>>(std::forward<Args>(args)...);
  default:
    throw std::runtime_error("Given potential is not provided yet!");
  }
}