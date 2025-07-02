#include <CellList.hpp>
#include <DomainDecomposition.hpp>
#include <IMDManager.hpp>
#include <Utility.hpp>
#include <VerletList.hpp>
#include <barostats/Barostat.hpp>
#include <memory>
#include <nlohmann/json.hpp>
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
  Vector3x<T> m_global_lo, m_global_hi;
  int m_ghost_index_first;
  std::vector<Particle<T>> m_particles;
  std::vector<T> m_electron_density;
  Vector3d m_cell_size;
  T m_cutoff, m_skin;

  std::string m_path_to_config;
  void calculateInternalForces();

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
  MDManager(Vector3x<T> global_lo, Vector3x<T> global_hi, Vector3d cell_size,
            T cutoff, T skin, int types_count, const std::string &_config_path)
      : m_potentials(std::vector<std::vector<std::unique_ptr<Potential<T>>>>(
            types_count, std::vector<std::unique_ptr<Potential<T>>>(
                             types_count, std::unique_ptr<Potential<T>>()))),
        m_path_to_config(_config_path) {}

  void initDomainDecomposition(Vector3x<T> global_lo, Vector3x<T> global_hi,
                               MPI_Comm comm) {
    m_domain_decomposition = std::make_unique<DomainDecomposition<T>>();
    m_domain_decomposition->initCart(comm, global_lo, global_hi);
  }

  void setPotential(int f_type, int s_type, typename Potential<T>::Type &type) {
    static int potential_set_count = 0;
    Potential<T> new_potential;
    switch (type) {
    case (Potential<T>::Type::EAM): {
      new_potential = std::make_unique<EAMPotential<T>>();
      break;
    }

    case (Potential<T>::Type::LJ): {
      new_potential = std::make_unique<LJPotential<T>>();
      break;
    }
    }

    new_potential->loadParameters(m_path_to_config);

    /*
    //HINT this thing seems redundant
    in place where potential is obtaining sorting can be done for types
    and then it can be called
    m_potentials[std::min(f_type,s_type)][std::max(f_type,s_type)]
    */

    if (f_type != s_type) {
      m_potentials[potential_set_count][potential_set_count + 1] =
          new_potential;
      m_potentials[potential_set_count + 1][potential_set_count] =
          new_potential;
      potential_set_count += 2;
    } else {
      m_potentials[potential_set_count][potential_set_count] = new_potential;
      potential_set_count++;
    }
  }
  void createCellLists(std::vector<Particle<T>> &all_particles,
                       Vector3d cell_size) {
    redistributeParticles(all_particles);
    m_ghost_index_first = m_particles.size();
    m_domain_decomposition->m_ghost_index_first = m_particles.size();

    m_domain_decomposition->exchangeGhosts(m_particles, m_cutoff + m_skin,
                                           &m_ghost_index_first);
    m_cell_list = std::make_unique<CellList<T>>();
    m_cell_list->build(m_particles, m_domain_decomposition->m_local_simbox,
                       m_cutoff + m_skin);
    m_verlet_list = std::make_unique<VerletList<T>>(m_cutoff, m_skin);
    m_verlet_list->build(m_particles, m_cell_list, m_ghost_index_first);
    m_electron_density.clear();
    m_electron_density.resize(m_ghost_index_first);
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

  T getPressure() const noexcept override;
  T getTemperature() const noexcept override;
  std::vector<Particle<T>> &getParticles() const noexcept override;
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

  template <typename Integrator> void setIntergrator();

  // void setPBC(const Vector3x<bool> &pbc, const Vector3d &lo,
  //             const Vector3d &hi);
};

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
