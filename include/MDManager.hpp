#include <CellList.hpp>
#include <DomainDecomposition.hpp>
#include <Eigen/Dense>
#include <IMDManager.hpp>
#include <Utility.hpp>
#include <VerletList.hpp>
#include <barostats/Barostat.hpp>
#include <memory>
#include <nlohmann/json.hpp>
#include <omp.h>
#include <potentials/EAM.hpp>
#include <potentials/LJ.hpp>
#include <thermostats/Thermostat.hpp>
#include <type_traits>
#include <vector>
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class MDManager : public IMDManager<T> {

  CellList<T> cell_list;
  std::unique_ptr<DomainDecomposition<T>> m_domain_decomposition;
  std::unique_ptr<Barostat<T>> m_barostat;
  std::unique_ptr<Thermostat<T>> m_thermostat;
  std::unique_ptr<VerletList<T>> m_verlet_list;
  std::unique_ptr<CellList<T>> m_cell_list;
  std::vector<std::vector<std::shared_ptr<Potential<T>>>> m_potentials;

  std::vector<std::vector<T>> m_cutoff_radiuses;
  SystemTraits<T> m_system_traits;
  Vector3x<T> m_global_lo, m_global_hi;
  int m_ghost_index_first;
  std::vector<Particle<T>> m_particles;
  // HINT: for metals only
  std::vector<T> m_electron_density;

  Vector3d m_cell_size;
  T m_max_cutoff, m_skin;
#if 0
  Eigen::Matrix3<T> m_distances;
#endif
  std::string m_path_to_config;

  void calculateInternalForces();

  void __calculateTotalElectronDensity() {
    std::fill(m_electron_density.begin(), m_electron_density.end(), 0.);
#ifdef OPENMP_ENABLED
#pragma omp parallel
#endif
    for (int t_ = 0; t_ != m_particles.size(); ++t_) {
      const Particle<T> &particle = m_particles.at(t_);
      if (particle.m_global_id != 0)
        continue;
      int f_type = particle.m_type_id_local;
      const std::vector<int> &nghbrs = m_verlet_list->getNeighbors(t_);

      for (int ngh_idx : nghbrs) {
        const Particle<T> &nghb = m_particles.at(ngh_idx);
        if (nghb.m_type_id_global != 0)
          continue;
        int ngh_type = nghb.m_type_id_local;
        T dist = length(particle.pos - nghb.pos);
        if (dist > m_cutoff_radiuses[f_type][ngh_type])
          continue;

        int min_idx = std::min(f_type, ngh_type);
        T contrib_in_me =
            m_potentials[f_type][ngh_type]->calculateElectronDensity(
                dist, ngh_type % min_idx);

#ifdef OPENMP_ENABLED
#pragma omp atomic
#endif
        m_electron_density.at(t_) += contrib_in_me;
        if (f_type != ngh_type)
          m_electron_density.at(ngh_idx) +=
              m_potentials.at(f_type).at(ngh_type)->calculateElectronDensity(
                  dist, f_type % min_idx);
        else
#ifdef OPENMP_ENABLED
#pragma omp atomic
#endif
          m_electron_density[ngh_idx] += contrib_in_me;
      }
    }
  }
  void redistributeParticles(std::vector<Particle<T>> &all_particles) {

    m_particles.clear();
    m_particles = std::vector<Particle<T>>();
    const ProcSubGrid3D<T> &local_simbox =
        m_domain_decomposition->m_local_simbox;
    for (auto &part : all_particles) {
      if (local_simbox.contains(part)) {
        m_particles.emplace(std::move(part));
      }
    }
  }

public:
  MDManager(int types_count, const std::string &_config_path)
      : m_potentials(std::vector<std::vector<std::unique_ptr<Potential<T>>>>(
            types_count, std::vector<std::unique_ptr<Potential<T>>>(
                             types_count, std::unique_ptr<Potential<T>>()))),
        m_path_to_config(_config_path),
        m_max_cutoff(std::numeric_limits<T>::min()) {}

  void initDomainDecomposition(Vector3x<T> global_lo, Vector3x<T> global_hi,
                               MPI_Comm comm) {
    m_domain_decomposition = std::make_unique<DomainDecomposition<T>>();
    m_domain_decomposition->initCart(comm, global_lo, global_hi);
  }

  void setPotential(int f_type, typename Potential<T>::Type &type,
                    int s_type = -1) {
    static int potential_set_count = 0;
    Potential<T> new_potential;
    switch (type) {
    case (Potential<T>::Type::EAM_ALLOY): {
      new_potential = std::make_unique<EAMPotentialAlloy<T>>();
      break;
    }
    case (Potential<T>::Type::EAM_PURE): {
      new_potential = std::make_unique<EAMPotentialPure<T>>();
      break;
    }

    case (Potential<T>::Type::LJ): {
      new_potential = std::make_unique<LJPotential<T>>();
      break;
    }
    }

    if (s_type != -1) {
      f_type = std::min(f_type, s_type);
      s_type = std::max(f_type, s_type);
    }
    new_potential->loadParameters(f_type, m_path_to_config, s_type);
    T current_cutoff = new_potential.getCutoffRadius();
    m_max_cutoff = std::max(m_max_cutoff, current_cutoff);
    m_skin = 0.1 * m_max_cutoff;

    /*
    //HINT this thing seems redundant
    in place where potential is obtaining sorting can be done for types
    and then it can be called
    m_potentials[std::min(f_type,s_type)][std::max(f_type,s_type)]
    */

    if (s_type != -1) {
      m_potentials[potential_set_count][potential_set_count + 1] =
          new_potential;
      m_potentials[potential_set_count + 1][potential_set_count] =
          new_potential;
      m_cutoff_radiuses[potential_set_count][potential_set_count + 1] =
          current_cutoff;
      m_cutoff_radiuses[potential_set_count + 1][potential_set_count] =
          current_cutoff;

      m_system_traits.m_types_mapping[f_type] = potential_set_count;
      m_system_traits.m_types_mapping[s_type] = potential_set_count + 1;

      potential_set_count += 2;
    } else {
      if (auto el = m_system_traits.m_types_mapping.find(f_type);
          el != m_system_traits.m_types_mapping.end()) {
        int potential_count = el->second;
        m_potentials[potential_count][potential_count] = new_potential;
        m_cutoff_radiuses[potential_count][potential_count] = current_cutoff;
        return;
      }

      m_system_traits.m_types_mapping[f_type] = potential_set_count;

      m_potentials[potential_set_count][potential_set_count] = new_potential;
      m_cutoff_radiuses[potential_set_count][potential_set_count] =
          current_cutoff;
      potential_set_count++;
    }
  }

  void createCellLists(std::vector<Particle<T>> &all_particles) {
    redistributeParticles(all_particles);

    remapParticleTypes();

    if (m_electron_density.size() < m_particles.size())
      m_electron_density.resize(m_particles.size());

    m_ghost_index_first = m_particles.size();
    m_domain_decomposition->m_ghost_index_first = m_particles.size();

    m_domain_decomposition->exchangeGhosts(m_particles, m_max_cutoff + m_skin,
                                           &m_ghost_index_first);
    m_cell_list = std::make_unique<CellList<T>>();
    m_cell_list->build(m_particles, m_domain_decomposition->m_local_simbox,
                       m_max_cutoff + m_skin);
    m_verlet_list = std::make_unique<VerletList<T>>(m_max_cutoff, m_skin,
                                                    m_cutoff_radiuses);
    m_verlet_list->build(m_particles, m_cell_list, m_ghost_index_first);
  }

  void evolve();

  template <template <class> class CustomBarostat,
            typename = std::enable_if_t<
                std::is_base_of_v<Barostat<T>, CustomBarostat<T>>>,
            typename... Args>
  void setBarostat(Args &&...args) {
    m_barostat =
        std::make_unique<CustomBarostat<T>>(std::forward<Args>(args)...);
  }

  template <typename... Args>
  void setBarostat(typename Barostat<T>::Type &type, Args &&...args) {
    m_barostat = createBarostat<T>(type, std::forward<Args>(args)...);
  }

  // HINT: this function have to be called after all setPotential's calls have
  // been made to conform actual indexes of materials with their representation,
  // that will be used in computational process

  void remapParticleTypes() noexcept(false) {
    for (Particle<T> &particle : m_particles) {
      particle.m_type_id_local =
          m_system_traits.m_types_mapping.at(particle.m_type_id_local);
    }
  }

  // HINT: Implementation of IMDManager interface
  virtual T getCurrentPressure() const noexcept override {
    return m_system_traits.m_pressure;
  }
  virtual T getCurrentTemperature() const noexcept override {
    return m_system_traits.m_temperature;
  }
  virtual std::vector<Particle<T>> &getParticles() const noexcept override;
  virtual T getTimeStep() const noexcept override {
    return m_system_traits.m_timestep;
  }

  Barostat<T> &getBarostat();

  template <template <class> class CustomThermostat,
            typename = std::enable_if_t<
                std::is_base_of_v<Thermostat<T>, CustomThermostat<T>>>,
            typename... Args>

  void setThermostat(Args &&...args) {
    m_thermostat =
        std::make_unique<CustomThermostat<T>>(std::forward<Args>(args)...);
  }
  template <typename... Args>
  void setThermostat(typename Thermostat<T>::Type &type, Args &&...args) {
    m_thermostat = createThermostat<T>(type, std::forward<Args>(args)...);
  }

  // template <typename Integrator> void setIntegrator();

  // void setPBC(const Vector3x<bool> &pbc, const Vector3d &lo,
  //             const Vector3d &hi);
};
#if 0
template <typename T> void MDManager<T>::evolve() { calculateInternalForces(); }
template <typename T> void MDManager<T>::calculateInternalForces() {

  // HINT: for metal atoms only
  calculateElectronDensity();

  for (int i : m_cell_list) {
    for (int j : m_cell_list[i]) {
      for (int k : m_verlet_list[j]) {
        if (radvec[flat_idx] > r_cut + 1e-08) {
          flat_idx++;
          continue;
        }
        double partial_force = calculatePartialForce(j, k);
        m_pInternalAccForce.at(j) -= partial_force;
        m_pInternalAccForce.at(k) += partial_force;
      }
    }
  }
}
#endif