#pragma once
#include "Potential.hpp"
#include <cmath>
#include <fstream>
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class LJPotential : public Potential<T> {

  T m_sigma, m_epsilon;

#ifdef WITH_SMOOTHING
  T m_pot_at_rcut;
  T m_derpot_at_rcut;
#endif
private:
  inline T __calculateLJ(const T &dist) {

    return 4 * m_epsilon * (pow(m_sigma / dist, 12) - pow(m_sigma / dist, 6));
  }
  inline T __calculateLJDer(const T &dist) {
    return 24 * m_epsilon / dist *
           (2 * pow(m_sigma / dist, 12) - pow(m_sigma / dist, 6));
  }

public:
  ~LJPotential() = default;
  virtual void loadParameters(int f_type, const std::string &path,
                              int s_type) noexcept(false) override {

    using json = nlohmann::json;
    if (!std::filesystem::exists(path))
      throw std::runtime_error(std::string("Provided file ") + path +
                               std::string(" does not exist!"));
    std::ifstream _file(path);
    json data = json::parse(_file)["LJ"][std::to_string(f_type)];
    m_sigma = data["sigma"].get<T>();
    m_epsilon = data["epsilon"].get<T>();

    // HINT: here not actual rcut being passed, but multiplyier for sigma
    this->m_rcut = data["r_cut"].get<T>() * m_epsilon;
#ifdef WITH_SMOOTHING
    m_pot_at_rcut = __calculateLJ(r_cut);
    m_derpot_at_rcut = __calculateLJDer(r_cut);
#endif
  }

  virtual void initStep(const std::vector<Particle<T>> &particles,
                        std::unique_ptr<VerletList<T>> verlet_list) override {}

  virtual T computeEnergy(const Particle<T> &p1,
                          const Particle<T> &p2) const noexcept override {

    T dist = length(p1.pos - p2.pos);
    if (dist > this->m_rcut)
      return 0.;
#ifdef WITH_SMOOTHING

    return __calculateLJ(dist) - m_pot_at_rcut - m_derpot_at_rcut;
#else

    return __calculateLJ(dist);
#endif
  }
  virtual Vector3x<T>
  computeForce(const Particle<T> &p1, const Particle<T> &p2,
               const std::pair<int, int> &indx = {}) const noexcept override {
    T dist = length(p1.pos - p2.pos);
    if (dist > this->m_rcut)
      return Vector3x<T>{-0., 0., 0.};
#ifdef WITH_SMOOTHING
    return (p1.pos - p2.pos) * (__calculateLJDer(dist) - m_derpot_at_rcut);
#else
    return (p1.pos - p2.pos) * -__calculateLJDer(dist);
#endif
  }
};