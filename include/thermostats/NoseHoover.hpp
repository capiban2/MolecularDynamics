#pragma once
#include "Thermostat.hpp"

template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class NoseHooverThermostat final : public Thermostat<T> {

  // HINT: initially xi set to zero to reflect no tension at start
  T m_xi = 0;
  const T m_Q;

public:
  NoseHooverThermostat(IMDManager<T> *mdm,
                       const std::vector<std::pair<int, T>> &_preserve_temp,
                       T _Q, T _xi = 0)
      : Thermostat<T>(mdm, _preserve_temp), m_Q(_Q), m_xi(_xi) {}

  NoseHooverThermostat() = delete;

  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<double> &mass,
                               bool half_step) override {

    Vector3x<T> external_force = this->m_mdmanager->getExternalForce();
    T time_stamp = this->m_mdmanager->getTimeStep();
    // HINT: update anything only on full step
    if (!half_step)
      this->__updateTargetTemperature();

    T target_temp =
        this->m_time_nd_temp_to_preserve.at(this->m_current_temp_index).second;
    if (half_step) {
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd
#endif
      for (int t_ = 0; t_ < this->m_first_ghost_index; ++t_) {
        _data.vel_x[t_] +=

            (this->m_time_step / 2.) *
            (_data.force_x[t_] / mass[t_] + external_force.x -
             m_xi * _data.vel_x[t_]);
        _data.vel_y[t_] +=

            (this->m_time_step / 2.) *
            (_data.force_y[t_] / mass[t_] + external_force.y -
             m_xi * _data.vel_y[t_]);
        _data.vel_z[t_] +=

            (this->m_time_step / 2.) *
            (_data.force_z[t_] / mass[t_] + external_force.z -
             m_xi * _data.vel_z[t_]);
      }

      this->m_half_step_kinetic = this->__calculateKinetic(
          _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
          this->m_full_vel_com, mass.data());

      this->m_half_vel_com = this->__calculateVelCOM(_data.vel_x, _data.vel_y,
                                                     _data.vel_z, mass.data());

      m_xi += this->m_time_step / (2 * m_Q) *
              (this->m_half_step_kinetic -
               3 * this->m_group_size * this->KB * target_temp / 2);
    } else {
      T denom = (1 + (this->m_time_step / 2.) * m_xi);
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd
#endif
      for (int t_ = 0; t_ < this->m_first_ghost_index; ++t_) {
        _data.vel_x[t_] = (_data.vel_x[t_] + (this->m_time_step / 2.) *
                                                 (_data.force_x[t_] / mass[t_] +
                                                  external_force.x)) /
                          denom;
        _data.vel_y[t_] = (_data.vel_y[t_] + (this->m_time_step / 2.) *
                                                 (_data.force_y[t_] / mass[t_] +
                                                  external_force.y)) /
                          denom;
        _data.vel_z[t_] = (_data.vel_z[t_] + (this->m_time_step / 2.) *
                                                 (_data.force_z[t_] / mass[t_] +
                                                  external_force.z)) /
                          denom;
      }
      this->m_full_step_kinetic = this->__calculateKinetic(
          _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
          this->m_half_vel_com, mass.data());

      this->m_full_vel_com = this->__calculateVelCOM(_data.vel_x, _data.vel_y,
                                                     _data.vel_z, mass.data());

      m_xi += this->m_time_step / (2 * m_Q) *
              (this->m_full_step_kinetic -
               3 * this->m_group_size * this->KB * target_temp / 2);
    }
  }
  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<int> &_idx,
                               const std::vector<double> &mass,
                               bool half_step) override {

    Vector3x<T> external_force = this->m_mdmanager->getExternalForce();

    if (!half_step)
      this->__updateTargetTemperature();

    T target_temp =
        this->m_time_nd_temp_to_preserve.at(this->m_current_temp_index).second;
    if (half_step) {
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd
#endif
      for (int t_ = 0; t_ < _idx.size(); ++t_) {
        _data.vel_x[_idx[t_]] +=

            (this->m_time_step / 2.) *
            (_data.force_x[_idx[t_]] / mass[_idx[t_]] + external_force.x -
             m_xi * _data.vel_x[_idx[t_]]);
        _data.vel_y[_idx[t_]] +=

            (this->m_time_step / 2.) *
            (_data.force_y[_idx[t_]] / mass[_idx[t_]] + external_force.y -
             m_xi * _data.vel_y[_idx[t_]]);
        _data.vel_z[_idx[t_]] +=

            (this->m_time_step / 2.) *
            (_data.force_z[_idx[t_]] / mass[_idx[t_]] + external_force.z -
             m_xi * _data.vel_z[_idx[t_]]);
      }

      this->m_half_step_kinetic = this->__calculateKinetic(
          _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
          this->m_full_vel_com, mass.data(), _idx);

      this->m_half_vel_com = this->__calculateVelCOM(
          _data.vel_x, _data.vel_y, _data.vel_z, mass.data(), _idx);

      m_xi += this->m_time_step / (2 * m_Q) *
              (this->m_half_step_kinetic -
               3 * this->m_group_size * this->KB * target_temp / 2);
    } else {
      T denom = (1 + (this->m_time_step / 2.) * m_xi);
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd
#endif
      for (int t_ = 0; t_ < _idx.size(); ++t_) {
        _data.vel_x[_idx[t_]] = (_data.vel_x[_idx[t_]] +
                                 (this->m_time_step / 2.) *
                                     (_data.force_x[_idx[t_]] / mass[_idx[t_]] +
                                      external_force.x)) /
                                denom;
        _data.vel_y[_idx[t_]] = (_data.vel_y[_idx[t_]] +
                                 (this->m_time_step / 2.) *
                                     (_data.force_y[_idx[t_]] / mass[_idx[t_]] +
                                      external_force.y)) /
                                denom;
        _data.vel_z[_idx[t_]] = (_data.vel_z[_idx[t_]] +
                                 (this->m_time_step / 2.) *
                                     (_data.force_z[_idx[t_]] / mass[_idx[t_]] +
                                      external_force.z)) /
                                denom;
      }
      this->m_full_step_kinetic = this->__calculateKinetic(
          _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
          this->m_half_vel_com, mass.data(), _idx);

      this->m_full_vel_com = this->__calculateVelCOM(
          _data.vel_x, _data.vel_y, _data.vel_z, mass.data(), _idx);

      m_xi += this->m_time_step / (2 * m_Q) *
              (this->m_full_step_kinetic -
               3 * this->m_group_size * this->KB * target_temp / 2);
    }
  }
};