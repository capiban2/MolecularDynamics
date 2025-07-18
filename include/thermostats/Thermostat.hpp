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
  int m_current_temp_index = 0;
  T KB = 1.380649;
  T m_half_step_kinetic = 0, m_full_step_kinetic = 0;
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

  Vector3x<T> m_half_vel_com;
  Vector3x<T> m_full_vel_com;
  T m_time_step;

  double m_group_mass = 0;
  int m_group_size = 0;
  int m_first_ghost_index = -1;

  T __calculateKinetic(T *velx, T *vely, T *velz, const Vector3x<T> &vcom,
                       double *mass) const noexcept {
    T kinetic = 0.;
    T dx, dy, dz;
#ifdef OPENMP_ENABLED
#pragma omp parallel for reduction(+ : kinetic)
#endif
    for (int t_ = 0; t_ != m_first_ghost_index; ++t_) {
      dx = velx[t_] - vcom.x;
      dy = vely[t_] - vcom.y;
      dz = velz[t_] - vcom.z;
      kinetic += mass[t_] * 0.5 * (dx * dx + dy * dy + dz * dz);
    }
    return kinetic;
  }
  T __calculateKinetic(T *velx, T *vely, T *velz, const Vector3x<T> &vcom,
                       double *mass,
                       const std::vector<int> &idx) const noexcept {
    T kinetic = 0.;
    T dx, dy, dz;
#ifdef OPENMP_ENABLED
#pragma omp parallel for reduction(+ : kinetic)
#endif
    for (int t_ = 0; t_ != idx.size(); ++t_) {
      dx = velx[idx[t_]] - vcom.x;
      dy = vely[idx[t_]] - vcom.y;
      dz = velz[idx[t_]] - vcom.z;
      kinetic += mass[idx[t_]] * 0.5 * (dx * dx + dy * dy + dz * dz);
    }
    return kinetic;
  }

  Vector3x<T> __calculateVelCOM(T *velx, T *vely, T *velz,
                                double *mass) const noexcept {

    Vector3x<T> com{0., 0., 0.};
#ifdef OPENMP_ENABLED
#pragma omp simd
#endif
    for (int t_ = 0; t_ != this->m_first_ghost_index; ++t_) {
      com.x += velx[t_] * mass[t_];
      com.y += vely[t_] * mass[t_];
      com.z += velz[t_] * mass[t_];
    }
    com /= m_group_mass;
    return com;
  }
  T __calculateVelCOM(T *velx, T *vely, T *velz, double *mass,
                      const std::vector<int> &idx) const noexcept {
    Vector3x<T> com{0., 0., 0.};
#ifdef OPENMP_ENABLED
#pragma omp simd
#endif
    for (int t_ = 0; t_ != idx.size(); ++t_) {
      com.x += velx[idx[t_]] * mass[idx[t_]];
      com.y += vely[idx[t_]] * mass[idx[t_]];
      com.z += velz[idx[t_]] * mass[idx[t_]];
    }
    com /= m_group_mass;
    return com;
  }

  void __updateTargetTemperature() {
    if (m_time_nd_temp_to_preserve.at(m_current_temp_index).first <= 0) {
      if (m_time_nd_temp_to_preserve.size() == m_current_temp_index + 1)
        // HINT: just add big value to preserve given temperature as long as
        // modelling goes
        m_time_nd_temp_to_preserve.at(m_current_temp_index).first += 1e+05;
      else
        m_current_temp_index++;
    }
  }

public:
  enum class Type { Identity, Andersen, Berendsen, Langevin, NoseHoover };
  Thermostat() = delete;
  virtual ~Thermostat() = default;
  Thermostat(IMDManager<T> *mdm,
             const std::vector<std::pair<int, T>> &_preserve_temp)
      : m_parent(mdm), m_time_nd_temp_to_preserve(_preserve_temp) {
    m_time_step = m_parent->getTimeStep();
  }

  // HINT: call it whenever particles got redistributed and new group have been
  // created
  void initThermostatGroup(int size, double mass, int first_ghost) {
    m_group_mass = mass;
    m_group_size = size;
    m_first_ghost_index = first_ghost;
  }

  // HINT: becauase dynamic target temperature has been supported, then, it may
  // return different target temperatures, that depend on current timestemp
  T getTargetTemperature() const noexcept {

    return m_time_nd_temp_to_preserve.at(m_current_temp_index).second;
  }

  T getHalfStepTemperature() const noexcept {

    return 2. / (3 * m_group_size * KB) * m_half_step_kinetic;
  }

  T getFullStepTemperature() const noexcept override {

    return 2. / (3 * m_group_size * KB) * m_full_step_kinetic;
  }

  const Vector3x<T> &getHalfVelocityCOM() const noexcept {
    return m_half_vel_com;
  }
  const Vector3x<T> &getFullVelocityCOM() const noexcept {
    return m_full_vel_com;
  }
  // HINT: for global thermostat; for calculating current system's
  // temperature(kinetic temperature) to work with it afterwards
  void initThermostat(ParticleData<T> &_data, const std::vector<double> &mass) {

    m_full_vel_com =
        __calculateVelCOM(_data.vel_x, _data.vel_y, _data.vel_z, mass.data());

    m_full_step_kinetic =
        __calculateKinetic(_data.vel_x.data(), _data.vel_y.data(),
                           _data.vel_z.data(), m_full_vel_com, mass);
  }
  // HINT: for local thermostat; for calculating current system's
  // temperature(kinetic temperature) to work with it afterwards
  void initThermostat(ParticleData<T> &_data, const std::vector<int> &_idx,
                      const std::vector<double> &mass) {

    m_full_vel_com = this->__calculateVelCOM(_data.vel_x, _data.vel_y,
                                             _data.vel_z, mass.data(), _idx);

    m_full_step_kinetic = this->__calculateKinetic(
        _data.vel_x.data(), _data.vel_y.data(), _data.vel_z.data(),
        m_full_vel_com, mass, _idx);
  }

  // HINT: for global thermostat
  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<double> &mass,
                               bool half_step) = 0;
  // HINT: for local thermostat
  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<int> &_idx,
                               const std::vector<double> &mass,
                               bool half_step) = 0;
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>,
          typename... Args>
std::unique_ptr<Thermostat<T>> createThermostat(typename Thermostat<T>::Type &,
                                                Args &&...args);
