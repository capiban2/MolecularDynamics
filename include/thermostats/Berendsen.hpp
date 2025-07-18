#pragma once
#include "Thermostat.hpp"
#include <Utility.hpp>
#include <algorithm>
#include <cmath>
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class BerendsenThermostat final : public Thermostat<T> {
  T m_ratio;

public:
  BerendsenThermostat(IMDManager<T> *mdm, T ratio,
                      const std::vector<std::pair<int, T>> &_preserve_temp)
      : Thermostat<T>(mdm, _preserve_temp), m_ratio(ratio) {}

  BerendsenThermostat() = delete;
  ~BerendsenThermostat() = default;

  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<double> &mass,
                               bool half_step) override {

    T current_temp = this->getFullStepTemperature();

    Vector3x<T> external_force = this->m_mdmanager->getExternalForce();

    // HINT: update anything only on full step
    if (!half_step) {
      this->__updateTargetTemperature();
    }
    T target_temp =
        this->m_time_nd_temp_to_preserve.at(this->m_current_temp_index).second;
    T lambda;
    if (half_step)
      lambda = sqrt(1 + m_ratio * (target_temp - current_temp) / current_temp);
    if (half_step) {
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd
#endif
      for (int t_ = 0; t_ < this->m_first_ghost_index; ++t_) {

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
      this->m_half_step_kinetic = this->__calculateKinetic(
          _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
          this->m_full_vel_com, mass.data());

      this->m_half_vel_com = this->__calculateVelCOM(_data.vel_x, _data.vel_y,
                                                     _data.vel_z, mass.data());

    } else {
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd

#endif
      for (int t_ = 0; t_ < this->m_first_ghost_index; ++t_) {

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
      this->m_full_step_kinetic = this->__calculateKinetic(
          _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
          this->m_half_vel_com, mass.data());

      this->m_full_vel_com = this->__calculateVelCOM(_data.vel_x, _data.vel_y,
                                                     _data.vel_z, mass.data());
    }
  }
  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<int> &_idx,
                               const std::vector<double> &mass,
                               bool half_step) override {

    T current_temp = this->getFullStepTemperature();
    Vector3x<T> external_force = this->m_mdmanager->getExternalForce();

    // HINT: update anything only on full step
    if (!half_step)
      this->__updateTargetTemperature();

    T target_temp =
        this->m_time_nd_temp_to_preserve.at(this->m_current_temp_index).second;
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
      this->m_half_step_kinetic = this->__calculateKinetic(
          _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
          this->m_full_vel_com, mass.data(), _idx);

      this->m_half_vel_com = this->__calculateVelCOM(
          _data.vel_x, _data.vel_y, _data.vel_z, mass.data(), _idx);

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
      this->m_full_step_kinetic = this->__calculateKinetic(
          _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
          this->m_half_vel_com, mass.data(), _idx);

      this->m_full_vel_com = this->__calculateVelCOM(
          _data.vel_x, _data.vel_y, _data.vel_z, mass.data(), _idx);
    }
  }
};