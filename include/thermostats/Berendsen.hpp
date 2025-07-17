#pragma once
#include "Thermostat.hpp"
#include <Utility.hpp>
#include <algorithm>
#include <cmath>
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class BerendsenThermostat final : public Thermostat<T> {
  T m_ratio;
  int m_current_temp_index = 0;
  T KB = 1.380649;

public:
  BerendsenThermostat(IMDManager<T> *mdm, T ratio,
                      const std::vector<std::pair<int, T>> &_preserve_temp,
                      int group_id)
      : Thermostat<T>(mdm, _preserve_temp, group_id), m_ratio(ratio) {}

  BerendsenThermostat() = delete;
  ~BerendsenThermostat() = default;

  virtual T getTargetTemperature() const noexcept override {

    return this->m_time_nd_temp_to_preserve.at(m_current_temp_index).second;
  }

  virtual T getHalfStepTemperature() const noexcept override {
    int thermostat_group_size =
        this->m_mdmanager->getThermostatGroupSize(this->m_thermostat_group_id);
    return 2. / (3 * thermostat_group_size * KB) * this->m_half_step_kinetic;
  }

  virtual T getFullStepTemperature() const noexcept override {
    int thermostat_group_size =
        this->m_mdmanager->getThermostatGroupSize(this->m_thermostat_group_id);
    return 2. / (3 * thermostat_group_size * KB) * this->m_full_step_kinetic;
  }
  // HINT: for global thermostat; for calculating current system's
  // temperature(kinetic temperature) to work with it afterwards
  virtual void initThermostat(ParticleData<T> &_data,
                              const std::vector<double> &mass) override {
    int first_ghost = this->m_mdmanager->getFirstGhostIndex();

    Vector3x<T> vcom = this->m_mdmanager->getThermostatGroupVelCom(
        this->m_thermostat_group_id);
    this->m_full_step_kinetic =
        this->calculateKinetic(_data.vel_x.data(), _data.vel_y.data(),
                               _data.vel_z.data(), vcom, mass, first_ghost);
  }
  // HINT: for local thermostat; for calculating current system's
  // temperature(kinetic temperature) to work with it afterwards
  virtual void initThermostat(ParticleData<T> &_data,
                              const std::vector<int> &_idx,
                              const std::vector<double> &mass) override {
    Vector3x<T> vcom = this->m_mdmanager->getThermostatGroupVelCom(
        this->m_thermostat_group_id);
    this->m_full_step_kinetic =
        this->calculateKinetic(_data.vel_x.data(), _data.vel_y.data(),
                               _data.vel_z.data(), vcom, mass, _idx);
  }
  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<double> &mass,
                               bool half_step) override {
    int first_ghost = this->m_mdmanager->getFirstGhostIndex();

    T current_temp = getFullStepTemperature();

    Vector3x<T> external_force = this->m_mdmanager->getExternalForce();
    T time_stamp = this->m_mdmanager->getTimeStep();
    // HINT: update anything only on full step
    if (!half_step) {
      if (this->m_time_nd_temp_to_preserve.at(m_current_temp_index).first <=
          0) {
        if (this->m_time_nd_temp_to_preserve.size() == m_current_temp_index + 1)
          // HINT: just add big value to preserve given temperature as long as
          // modelling goes
          this->m_time_nd_temp_to_preserve.at(m_current_temp_index).first +=
              1e+05;
        else
          m_current_temp_index++;
      }
    }
    T target_temp =
        this->m_time_nd_temp_to_preserve.at(m_current_temp_index).second;
    T lambda;
    if (half_step)
      lambda = sqrt(1 + m_ratio * (target_temp - current_temp) / current_temp);
    if (half_step) {
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd
#endif
      for (int t_ = 0; t_ < first_ghost; ++t_) {

        _data.vel_x[t_] =
            lambda * (_data.vel_x[t_] +
                      (_data.force_x[t_] / mass[t_] + external_force.x) / 2 *
                          this->m_time_step);
        _data.vel_y[t_] =
            lambda * (_data.vel_y[t_] +
                      (_data.force_y[t_] / mass[t_] + external_force.y) / 2 *
                          this->m_time_step);
        _data.vel_z[t_] =
            lambda * (_data.vel_z[t_] +
                      (_data.force_z[t_] / mass[t_] + external_force.z) / 2 *
                          this->m_time_step);
      }

    } else {
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd

#endif
      for (int t_ = 0; t_ < first_ghost; ++t_) {

        _data.vel_x[t_] = (_data.vel_x[t_] +
                           (_data.force_x[t_] / mass[t_] + external_force.x) /
                               2 * this->m_time_step);
        _data.vel_y[t_] = (_data.vel_y[t_] +
                           (_data.force_y[t_] / mass[t_] + external_force.y) /
                               2 * this->m_time_step);
        _data.vel_z[t_] = (_data.vel_z[t_] +
                           (_data.force_z[t_] / mass[t_] + external_force.z) /
                               2 * this->m_time_step);
      }
    }
  }
  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<int> &_idx,
                               const std::vector<double> &mass,
                               bool half_step) override {
    // HINT: think about it; may be thermostat itself responsible for providing
    // this stuff
    T current_temp = this->m_mdmanager->getThermostatGroupTemperature(
        this->m_thermostat_group_id);

    Vector3x<T> external_force = this->m_mdmanager->getExternalForce();
    T time_stamp = this->m_mdmanager->getTimeStep();
    // HINT: update anything only on full step
    if (!half_step) {
      if (this->m_time_nd_temp_to_preserve.at(m_current_temp_index).first <=
          0) {
        if (this->m_time_nd_temp_to_preserve.size() == m_current_temp_index + 1)
          // HINT: just add big value to preserve given temperature as long as
          // modelling goes
          this->m_time_nd_temp_to_preserve.at(m_current_temp_index).first +=
              1e+05;
        else
          m_current_temp_index++;
      }
    }
    T target_temp =
        this->m_time_nd_temp_to_preserve.at(m_current_temp_index).second;
    T lambda;
    if (half_step)
      lambda = sqrt(1 + m_ratio * (target_temp - current_temp) / current_temp);

    if (half_step) {
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd
#endif
      for (int t_ = 0; t_ < _idx.size(); ++t_) {

        _data.vel_x[_idx[_idx[t_]]] =
            lambda *
            (_data.vel_x[_idx[t_]] +
             (_data.force_x[_idx[t_]] / mass[_idx[t_]] + external_force.x) / 2 *
                 this->m_time_step);
        _data.vel_y[_idx[t_]] =
            lambda *
            (_data.vel_y[_idx[t_]] +
             (_data.force_y[_idx[t_]] / mass[_idx[t_]] + external_force.y) / 2 *
                 this->m_time_step);
        _data.vel_z[_idx[t_]] =
            lambda *
            (_data.vel_z[_idx[t_]] +
             (_data.force_z[_idx[t_]] / mass[_idx[t_]] + external_force.z) / 2 *
                 this->m_time_step);
      }

    } else {
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd

#endif
      for (int t_ = 0; t_ < _idx.size(); ++t_) {

        _data.vel_x[_idx[t_]] =
            (_data.vel_x[_idx[t_]] +
             (_data.force_x[_idx[t_]] / mass[_idx[t_]] + external_force.x) / 2 *
                 this->m_time_step);
        _data.vel_y[_idx[t_]] =
            (_data.vel_y[_idx[t_]] +
             (_data.force_y[_idx[t_]] / mass[_idx[t_]] + external_force.y) / 2 *
                 this->m_time_step);
        _data.vel_z[_idx[t_]] =
            (_data.vel_z[_idx[t_]] +
             (_data.force_z[_idx[t_]] / mass[_idx[t_]] + external_force.z) / 2 *
                 this->m_time_step);
      }
    }
  }
};