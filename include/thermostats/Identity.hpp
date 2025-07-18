#pragma once

#include "Thermostat.hpp"
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class IdentityThermostat final : public Thermostat<T> {

  void __updateKinetic(ParticleData<T> &_data,
                       const std::vector<double> &mass) {
    this->m_full_step_kinetic = this->__calculateKinetic(
        _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
        this->m_half_vel_com, mass);

    this->m_full_vel_com = this->__calculateVelCOM(_data.vel_x, _data.vel_y,
                                                   _data.vel_z, mass.data());
  }
  void __updateKinetic(ParticleData<T> &_data, const std::vector<int> &_idx,
                       const std::vector<double> &mass) {
    this->m_full_step_kinetic = this->__calculateKinetic(
        _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
        this->m_half_vel_com, mass, _idx);

    this->m_full_vel_com = this->__calculateVelCOM(
        _data.vel_x, _data.vel_y, _data.vel_z, mass.data(), _idx);
  }

public:
  IdentityThermostat(IMDManager<T> *mdm,
                     const std::vector<std::pair<int, T>> &_preserve_temp)
      : Thermostat<T>(mdm, _preserve_temp) {}

  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<double> &_mass,
                               bool half_step) {
    if (half_step)
      return;
    __updateKinetic(_data, _mass);
  }

  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<int> &_idx,
                               const std::vector<double> &_mass,
                               bool half_step) override {
    if (half_step)
      return;
    __updateKinetic(_data, _idx, _mass);
  }
};