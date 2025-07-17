#pragma once
#include "Thermostat.hpp"

template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class NoseHooverThermostat final : public Thermostat<T> {
  int m_current_temp_index = 0;

  T KB = 1.380649;
  // HINT: initially xi set  to zero to reflect no tension at start
  T m_xi = 0;
  const T m_Q;

public:
  // TODO: add necessary parameters
  NoseHooverThermostat(IMDManager<T> *mdm, T _Q,
                       const std::vector<std::pair<int, T>> &_preserve_temp,
                       int group_id, T _xi = 0)
      : Thermostat<T>(mdm, _preserve_temp, group_id), m_Q(_Q), m_xi(_xi) {}

  NoseHooverThermostat() = delete;
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
#if 0
    int first_ghost = this->m_mdmanager->getFirstGhostIndex();

    Vector3x<T> vcom = this->m_mdmanager->getThermostatGroupVelCom(
        this->m_thermostat_group_id);
    this->m_full_step_kinetic =
        this->calculateKinetic(_data.vel_x.data(), _data.vel_y.data(),
                               _data.vel_z.data(), vcom, mass, first_ghost);
#endif
  }
  // HINT: for local thermostat; for calculating current system's
  // temperature(kinetic temperature) to work with it afterwards
  virtual void initThermostat(ParticleData<T> &_data,
                              const std::vector<int> &_idx,
                              const std::vector<double> &mass) override {
#if 0
    Vector3x<T> vcom = this->m_mdmanager->getThermostatGroupVelCom(
        this->m_thermostat_group_id);
    this->m_full_step_kinetic =
        this->calculateKinetic(_data.vel_x.data(), _data.vel_y.data(),
                               _data.vel_z.data(), vcom, mass, _idx);
#endif
  }
  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<double> &mass,
                               bool half_step) override {
    int first_ghost = this->m_mdmanager->getFirstGhostIndex();

    Vector3x<T> external_force = this->m_mdmanager->getExternalForce();
    T time_stamp = this->m_mdmanager->getTimeStep();
    int thermostat_group_size =
        this->m_mdmanager->getThermostatGroupSize(this->m_thermostat_group_id);

    Vector3x<T> vcom = this->m_mdmanager->getThermostatGroupVelCom(
        this->m_thermostat_group_id);

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
    if (half_step) {
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd
#endif
      for (int t_ = 0; t_ < first_ghost; ++t_) {
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

      this->m_half_step_kinetic =
          this->calculateKinetic(_data.vel_x.data(), _data.vel_y.data(),
                                 _data.vel_z.data(), vcom, mass, first_ghost);

      m_xi += this->m_time_step / (2 * m_Q) *
              (this->m_half_step_kinetic -
               3 * thermostat_group_size * KB * target_temp / 2);
    } else {
      T denom = (1 + (this->m_time_step / 2.) * m_xi);
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd
#endif
      for (int t_ = 0; t_ < first_ghost; ++t_) {
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
      this->m_full_step_kinetic =
          this->calculateKinetic(_data.vel_x.data(), _data.vel_y.data(),
                                 _data.vel_z.data(), vcom, mass, first_ghost);

      m_xi += this->m_time_step / (2 * m_Q) *
              (this->m_full_step_kinetic -
               3 * thermostat_group_size * KB * target_temp / 2);
    }
  }
  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<int> &_idx,
                               const std::vector<double> &mass,
                               bool half_step) override {
    Vector3x<T> external_force = this->m_mdmanager->getExternalForce();
    T time_stamp = this->m_mdmanager->getTimeStep();
    int thermostat_group_size =
        this->m_mdmanager->getThermostatGroupSize(this->m_thermostat_group_id);

    Vector3x<T> vcom = this->m_mdmanager->getThermostatGroupVelCom(
        this->m_thermostat_group_id);

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

      this->m_half_step_kinetic =
          this->calculateKinetic(_data.vel_x.data(), _data.vel_y.data(),
                                 _data.vel_z.data(), vcom, mass, _idx);

      m_xi += this->m_time_step / (2 * m_Q) *
              (this->m_half_step_kinetic -
               3 * thermostat_group_size * KB * target_temp / 2);
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
      this->m_full_step_kinetic =
          this->calculateKinetic(_data.vel_x.data(), _data.vel_y.data(),
                                 _data.vel_z.data(), vcom, mass, _idx);

      m_xi += this->m_time_step / (2 * m_Q) *
              (this->m_full_step_kinetic -
               3 * thermostat_group_size * KB * target_temp / 2);
    }
  }
};