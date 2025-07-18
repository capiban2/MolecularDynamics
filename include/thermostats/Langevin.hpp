#pragma once
#include "Thermostat.hpp"
#include <random>

template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class LangevinThermostat final : public Thermostat<T> {
  T m_gamma; // Friction/damping coefficient
  std::mt19937 m_rng;
  std::normal_distribution<T> m_norm;

public:
  LangevinThermostat(IMDManager<T> *mdm,
                     const std::vector<std::pair<int, T>> &_preserve_temp,
                     T _gamma, unsigned int seed = 42)
      : Thermostat<T>(mdm, _preserve_temp), m_gamma(_gamma), m_rng(seed),
        m_norm(0., 1.) {}

  LangevinThermostat() = delete;
  // HINT: for global thermostat
  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<double> &mass,
                               bool half_step) override {

    Vector3x<T> external_force = this->m_mdmanager->getExternalForce();
    T dt = this->m_mdmanager->getTimeStep();
    T kB = this->KB;
    // HINT: update anything only on full step
    if (!half_step)
      this->__updateTargetTemperature();

    T target_temp =
        this->m_time_nd_temp_to_preserve.at(this->m_current_temp_index).second;
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd
#endif
    for (int t_ = 0; t_ < this->m_first_ghost_index; ++t_) {
      T sqrt_coeff = std::sqrt(2. * m_gamma * kB * target_temp / dt / mass[t_]);

      // Deterministic force
      T fx = _data.force_x[t_] / mass[t_] + external_force.x;
      T fy = _data.force_y[t_] / mass[t_] + external_force.y;
      T fz = _data.force_z[t_] / mass[t_] + external_force.z;

      // For half-step, apply (dt/2); for full-step, apply dt
      T scale = dt / 2.;

      _data.vel_x[t_] += scale * (fx - m_gamma * _data.vel_x[t_]) +
                         sqrt_coeff * m_norm(m_rng) * std::sqrt(scale);
      _data.vel_y[t_] += scale * (fy - m_gamma * _data.vel_y[t_]) +
                         sqrt_coeff * m_norm(m_rng) * std::sqrt(scale);
      _data.vel_z[t_] += scale * (fz - m_gamma * _data.vel_z[t_]) +
                         sqrt_coeff * m_norm(m_rng) * std::sqrt(scale);
    }

    // Update kinetic energy and velocity center of mass
    if (half_step) {
      this->m_half_step_kinetic = this->__calculateKinetic(
          _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
          this->m_full_vel_com, mass.data());
      this->m_half_vel_com = this->__calculateVelCOM(_data.vel_x, _data.vel_y,
                                                     _data.vel_z, mass.data());
    } else {
      this->m_full_step_kinetic = this->__calculateKinetic(
          _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
          this->m_half_vel_com, mass.data());
      this->m_full_vel_com = this->__calculateVelCOM(_data.vel_x, _data.vel_y,
                                                     _data.vel_z, mass.data());
    }
  }
  // HINT: for local thermostat
  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<int> &_idx,
                               const std::vector<double> &mass,
                               bool half_step) override {
    Vector3x<T> external_force = this->m_mdmanager->getExternalForce();
    T dt = this->m_mdmanager->getTimeStep();
    T kB = this->KB;
    // HINT: update anything only on full step
    if (!half_step)
      this->__updateTargetTemperature();

    T target_temp =
        this->m_time_nd_temp_to_preserve.at(this->m_current_temp_index).second;
#ifdef OPENMP_ENABLED
#pragma omp parallel for simd
#endif
    for (int t_ = 0; t_ < _idx.size(); ++t_) {
      T sqrt_coeff = std::sqrt(2. * m_gamma * kB * target_temp / dt / mass[t_]);

      // Deterministic force
      T fx = _data.force_x[_idx[t_]] / mass[_idx[t_]] + external_force.x;
      T fy = _data.force_y[_idx[t_]] / mass[_idx[t_]] + external_force.y;
      T fz = _data.force_z[_idx[t_]] / mass[_idx[t_]] + external_force.z;

      // For half-step, apply (dt/2); for full-step, apply dt
      T scale = dt / 2.;

      _data.vel_x[_idx[t_]] += scale * (fx - m_gamma * _data.vel_x[_idx[t_]]) +
                               sqrt_coeff * m_norm(m_rng) * std::sqrt(scale);
      _data.vel_y[_idx[t_]] += scale * (fy - m_gamma * _data.vel_y[_idx[t_]]) +
                               sqrt_coeff * m_norm(m_rng) * std::sqrt(scale);
      _data.vel_z[_idx[t_]] += scale * (fz - m_gamma * _data.vel_z[_idx[t_]]) +
                               sqrt_coeff * m_norm(m_rng) * std::sqrt(scale);
    }

    // Update kinetic energy and velocity center of mass
    if (half_step) {
      this->m_half_step_kinetic = this->__calculateKinetic(
          _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
          this->m_full_vel_com, mass.data(), _idx);
      this->m_half_vel_com = this->__calculateVelCOM(
          _data.vel_x, _data.vel_y, _data.vel_z, mass.data(), _idx);
    } else {
      this->m_full_step_kinetic = this->__calculateKinetic(
          _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
          this->m_half_vel_com, mass.data());
      this->m_full_vel_com = this->__calculateVelCOM(
          _data.vel_x, _data.vel_y, _data.vel_z, mass.data(), _idx);
    }
  }
};