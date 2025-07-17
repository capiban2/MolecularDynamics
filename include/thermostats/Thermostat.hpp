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

  T m_time_step;

  // HINT: necessary for systems with hybrid thermostating, to get temperature
  // of proper group
  // this id should be installed by MDManager
  int m_thermostat_group_id;

  T calculateKinetic(T *velx, T *vely, T *velz, const Vector3x<T> &vcom,
                     double *mass, int first_ghost) {
    T kinetic = 0.;
    T dx, dy, dz;
#ifdef OPENMP_ENABLED
#pragma omp parallel for reduction(+ : kinetic)
#endif
    for (int t_ = 0; t_ != first_ghost; ++t_) {
      dx = velx[t_] - vcom.x;
      dy = vely[t_] - vcom.y;
      dz = velz[t_] - vcom.z;
      kinetic += mass[t_] * 0.5 * (dx * dx + dy * dy + dz * dz);
    }
    return kinetic;
  }
  T calculateKinetic(T *velx, T *vely, T *velz, const Vector3x<T> &vcom,
                     double *mass, const std::vector<int> &idx) {
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

public:
  enum class Type { Identity, Andersen, Berendsen, Langevin, NoseHoover };
  Thermostat() = delete;
  virtual ~Thermostat() = default;
  Thermostat(IMDManager<T> *mdm,
             const std::vector<std::pair<int, T>> &_preserve_temp, int group_id)
      : m_parent(mdm), m_time_nd_temp_to_preserve(_preserve_temp),
        m_thermostat_group_id(group_id) {
    m_time_step = m_parent->getTimeStep();
  }

  void setData(IMDManager<T> *mdm,
               const std::vector<std::pair<int, T>> &_preserve_temp) noexcept {
    m_parent = mdm;
    m_time_nd_temp_to_preserve = _preserve_temp;
    m_time_step = m_parent->getTimeStep();
  }
  // HINT: becauase dynamic target temperature has been supported, then, it may
  // return different target temperatures, that depend on current timestemp
  virtual T getTargetTemperature() const noexcept;
  virtual T getHalfStepTemperature() const noexcept;
  virtual T getFullStepTemperature() const noexcept;

  // HINT: for global thermostat
  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<double> &mass,
                               bool half_step) = 0;
  // HINT: for local thermostat
  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<int> &_idx,
                               const std::vector<double> &mass,
                               bool half_step) = 0;

  // HINT: for global thermostat; for calculating current system's
  // temperature(kinetic temperature) to work with it afterwards
  virtual void initThermostat(ParticleData<T> &_data,
                              const std::vector<double> &mass) = 0;
  // HINT: for local thermostat; for calculating current system's
  // temperature(kinetic temperature) to work with it afterwards
  virtual void initThermostat(ParticleData<T> &_data,
                              const std::vector<int> &_idx,
                              const std::vector<double> &mass) = 0;
};
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class IdentityThermostat final : public Thermostat<T> {

public:
  IdentityThermostat(IMDManager<T> *mdm,
                     const std::vector<std::pair<int, T>> &_preserve_temp,
                     int group_id)
      : Thermostat<T>(mdm, _preserve_temp, group_id) {}

  virtual void applyThermostat(ParticleData<T> &_data, bool to_half_velocity) {}

  virtual void applyThermostat(ParticleData<T> &_data,
                               const std::vector<int> &_idx,
                               bool to_half_velocity) override {}

  virtual T getTargetTemperature() const noexcept override {
    return this->m_parent->getCurrentTemperature();
  }
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>,
          typename... Args>
std::unique_ptr<Thermostat<T>> createThermostat(typename Thermostat<T>::Type &,
                                                Args &&...args);
