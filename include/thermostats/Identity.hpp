#pragma once

#include "Thermostat.hpp"
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class IdentityThermostat final : public Thermostat<T> {
  T KB = 1.380649;

  void __updateKinetic(ParticleData<T> &_data,
                       const std::vector<double> &mass) {
    int first_ghost = this->m_mdmanager->getFirstGhostIndex();

    Vector3x<T> vcom = this->m_mdmanager->getThermostatGroupVelCom(
        this->m_thermostat_group_id);
    this->m_full_step_kinetic =
        this->calculateKinetic(_data.vel_x.data(), _data.vel_y.data(),
                               _data.vel_z.data(), vcom, mass, first_ghost);
  }
  void __updateKinetic(ParticleData<T> &_data, const std::vector<int> &_idx,
                       const std::vector<double> &mass) {
    Vector3x<T> vcom = this->m_mdmanager->getThermostatGroupVelCom(
        this->m_thermostat_group_id);
    this->m_full_step_kinetic =
        this->calculateKinetic(_data.vel_x.data(), _data.vel_y.data(),
                               _data.vel_z.data(), vcom, mass, _idx);
  }

public:
  IdentityThermostat(IMDManager<T> *mdm,
                     const std::vector<std::pair<int, T>> &_preserve_temp,
                     int group_id)
      : Thermostat<T>(mdm, _preserve_temp, group_id) {}
  virtual T getTargetTemperature() const noexcept override {}
  virtual T getHalfStepTemperature() const noexcept override {}
  virtual T getFullStepTemperature() const noexcept override {}
  // HINT: for global thermostat; for calculating current system's
  // temperature(kinetic temperature) to work with it afterwards
  virtual void initThermostat(ParticleData<T> &_data,
                              const std::vector<double> &mass) override {

    __updateKinetic(_data, mass);
  }
  // HINT: for local thermostat; for calculating current system's
  // temperature(kinetic temperature) to work with it afterwards
  virtual void initThermostat(ParticleData<T> &_data,
                              const std::vector<int> &_idx,
                              const std::vector<double> &mass) override {

    __updateKinetic(_data, _idx, mass);
  }
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