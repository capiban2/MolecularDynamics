#pragma once
#include "Potential.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <utility>
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class LJPotential final : public Potential<T> {

  std::vector<std::vector<T>> m_sigma, m_epsilon;

  std::vector<std::pair<int, int>> m_types;
  std::string m_config_path;
  int m_n_types = 0;
  int m_max_type = 0;

  std::vector<std::vector<int>> m_types2map;
  T m_max_r = 0, t_min_r = 0;

#ifdef WITH_SMOOTHING
  T m_pot_at_rcut;
  T m_derpot_at_rcut;
#endif
private:
  inline T __calculateLJ(const T &dist, const T &sigm, const T &eps) {

    return 4 * eps * (pow(sigm / dist, 12) - pow(sigm / dist, 6));
  }
  inline T __calculateLJDer(const T &dist, const T &sigm, const T &eps) {
    return 24 * eps / dist * (2 * pow(sigm / dist, 12) - pow(sigm / dist, 6));
  }

  void allocate() {
    auto max_pair = *std::max_element(m_types.begin(), m_types.end());
    m_max_type = std::max(max_pair.first, max_pair.second);

    m_sigma = std::vector<std::vector<T>>(m_types.size(),
                                          std::vector<T>(m_types.size(), T{}));
    m_epsilon = std::vector<std::vector<T>>(
        m_types.size(), std::vector<T>(m_types.size(), T{}));
    m_types2map = std::vector<std::vector<int>>(
        m_max_type, std::vector<int>(m_max_type, -1));
  }

  void file2array() {
    using json = nlohmann::json;
    if (!std::filesystem::exists(m_config_path))
      throw std::runtime_error(std::string("Provided file ") + m_config_path +
                               std::string(" does not exist!"));
    std::ifstream _file(m_config_path);
    json lj_data = json::parse(_file)["LJ"];
    json type2name = json::parse(_file)["Type2Name"];
    for (auto const &type_pair : m_types) {

      int f_type = type_pair.first, s_type = type_pair.second;
      std::string f_type_name = type2name[std::to_string(f_type)],
                  s_type_name = type2name[std::to_string(s_type)];
      json pair_data = lj_data[f_type_name][s_type_name];
      T sigma = pair_data["sigma"].get<T>(),
        epsilon = pair_data["eps"].get<T>();

      // HINT: here not actual rcut being passed, but multiplyier for sigma
      m_max_r = std::max(m_max_r, pair_data["r_cut"].get<T>() * epsilon);

      m_epsilon[m_types2map[f_type][s_type]] = epsilon;
      m_epsilon[m_types2map[s_type][f_type]] = epsilon;
      m_sigma[m_types2map[f_type][s_type]] = sigma;
      m_sigma[m_types2map[s_type][f_type]] = sigma;
    }
  }

  void array2spline() {}

public:
  virtual T getCutoffRadius() const noexcept override { return m_max_r; }
  // HINT: these types should be unique
  LJPotential(const std::vector<std::pair<int, int>> &types,
              const std::string &path)
      : m_types(types), m_config_path(path) {}
  void configure() override {

    allocate();

    for (int t_ = 0; t_ != m_types.size(); ++t_) {
      m_types2map[m_types[t_].first][m_types[t_].second] = t_;
      m_types2map[m_types[t_].second][m_types[t_].first] = t_;
    }

    file2array();
  }
  ~LJPotential() = default;
#if 0
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
#endif
#if 0
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
#endif
  virtual enum PotentialType rtti() override { return PotentialType::LJ; }
  virtual void
  computeBulkForces(const std::vector<MD::PairInteraction<T>> &_interactions,
                    std::vector<T> &_f_x, std::vector<T> &_f_y,
                    std::vector<T> &_f_z, const std::vector<int> types,
                    int first_ghost) override {

    T x, y, z, delx, dely, delz, dist, fpair;
    int i_, j_;
#ifdef OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int t_ = 0; t_ < _interactions.size(); ++t_) {
      auto &interact = _interactions[t_];

      if (interact.r >= m_max_r + 1e-08)
        continue;
      i_ = interact.i, j_ = interact.j;
      fpair =
          -__calculateLJDer(dist, m_sigma[m_types2map[types[i_]][types[j_]]],
                            m_epsilon[m_types2map[types[i_]][types[j_]]]);
      if (i_ < first_ghost) {
#ifdef OPENMP_ENABLED
#pragma openmp critical {
#endif
        _f_x[i_] += delx * fpair;
        _f_y[i_] += dely * fpair;
        _f_z[i_] += delz * fpair;

#ifdef OPENMP_ENABLED
      }
#endif
    }
    if (j_ < first_ghost) {
#ifdef OPENMP_ENABLED
#pragma openmp critical {
#endif

      _f_x[j_] -= delx * fpair;
      _f_y[j_] -= dely * fpair;
      _f_z[j_] -= delz * fpair;
#ifdef OPENMP_ENABLED
    }
#endif
  }
}
}

virtual void computeBulkForces(const VerletList<T> &neigh_list,
                               ParticleData<T> &_data,
                               int first_ghost) override {

  T *fx = _data.force_x, *fy = _data.force_y, *fz = _data.force_z;
  const std::vector<int> &types = _data.m_type_id;
  T x, y, z, delx, dely, delz, dist, fpair;
#ifdef OPENMP_ENABLED
#pragma omp parallel for
#endif
  for (int t_ = 0; t_ < first_ghost; ++t_) {
    int f_type = types[t_];
    x = _data.pos_x[t_];
    y = _data.pos_y[t_];
    z = _data.pos_z[t_];

    for (int j_ : neigh_list.neighbours(t_)) {

      // TODO: if j_ > t_ here necessary
      if (j_ >= first_ghost && j_ > t_)
        continue;
      int s_type = types[j_];
      delx = x - _data.pos_x[j_];
      dely = y - _data.pos_y[j_];
      delz = z - _data.pos_z[j_];
      dist = sqrt(delx * delx + dely * dely + delz * delz);
      // TODO: think about 1e-08
      if (dist >= m_max_r + 1e-08)
        continue;

      fpair = -__calculateLJDer(dist, m_sigma[m_types2map[f_type][s_type]],
                                m_epsilon[m_types2map[f_type][s_type]]);
#ifdef OPENMP_ENABLED
#pragma openmp critical {
#endif
      fx[t_] += delx * fpair;
      fy[t_] += dely * fpair;
      fz[t_] += delz * fpair;
      fx[j_] -= delx * fpair;
      fy[j_] -= dely * fpair;
      fz[j_] -= delz * fpair;
#ifdef OPENMP_ENABLED
    }
#endif
  }
}
}
}
;