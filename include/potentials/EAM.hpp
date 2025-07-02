#pragma once
#include "Potential.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stdexcept>
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class EAMPotential : public Potential<T> {

protected:
  struct Constants {
    T r_e;
    T f_e;
    T rho_e;
    T rho_n;
    T rho_zero;
    T rho_s;
    T alpha;
    T betta;
    T kappa;
    T lambda;
    T A;
    T B;
    T omega_n[4];
    T omega_e[4];
    T omega_s;
    T etta;
    int idx;
  };
  T r_cut, r_on;
  // HINT: C2 smoothing
  inline void taperQuinticDer(const T &dist, T &S, T &dSdr) {
    const auto dr = r_cut - r_on;
    if (dist < r_on) {
      dSdr = 0.0;
      S = 1.0;
      return;
    }
#if 0
    if (r >= r_c) {
        dSdr = 0.0;
        return 0.0;
    }
#endif
    auto xi = (dist - r_on) / dr;
    auto xi2 = xi * xi;
    auto xi3 = xi2 * xi;
    auto xi4 = xi3 * xi;
    auto xi5 = xi4 * xi;
    // Switch
    S = 1.0 - 10.0 * xi3 + 15.0 * xi4 - 6.0 * xi5;
    // Derivative: dS/dξ = -30ξ² + 60ξ³ - 30ξ⁴
    // dS/dr = (dS/dξ) * (dξ/dr) = (dS/dξ) / dr
    auto dSdxi = -30.0 * xi2 + 60.0 * xi3 - 30.0 * xi4;
    dSdr = dSdxi / dr;
  }
  inline void taperQuintic(const T &dist, T &S) {
    const auto dr = r_cut - r_on;
    if (dist < r_on) {
      S = 1.0;
      return;
    }
    auto xi = (dist - r_on) / dr;
    auto xi2 = xi * xi;
    auto xi3 = xi2 * xi;
    auto xi4 = xi3 * xi;
    auto xi5 = xi4 * xi;
    // Switch
    S = 1.0 - 10.0 * xi3 + 15.0 * xi4 - 6.0 * xi5;
  }
  // Eigen::Matrix3<T> m_bfunc;
  std::vector<Constants> m_constants;
  std::vector<T> m_total_electron_density;
  // std::vector<T> m_bfunc_precomputed;
  virtual void __loadParameters(int type_idx, int storage_idx,
                                const std::string &path) noexcept(false) {
    using json = nlohmann::json;
    if (!std::filesystem::exists(path))
      throw std::runtime_error(std::string("Provided file ") + path +
                               std::string(" is not exist!"));
    std::ifstream _file(path);
    json data = json::parse(_file)[std::to_string(type_idx)];
    m_constants.at(storage_idx).r_e = data["r_e"].get<T>();
    m_constants.at(storage_idx).f_e = data["f_e"].get<T>();
    m_constants.at(storage_idx).rho_e = data["rho_e"].get<T>();
    m_constants.at(storage_idx).rho_zero = data["rho_zero"].get<T>();
    m_constants.at(storage_idx).rho_s = data["rho_s"].get<T>();
    m_constants.at(storage_idx).alpha = data["alpha"].get<T>();
    m_constants.at(storage_idx).betta = data["betta"].get<T>();
    m_constants.at(storage_idx).kappa = data["kappa"].get<T>();
    m_constants.at(storage_idx).lambda = data["lambda"].get<T>();
    m_constants.at(storage_idx).A = data["A"].get<T>();
    m_constants.at(storage_idx).B = data["B"].get<T>();
    for (const auto &om : data["omega_n"]) {
      m_constants.at(storage_idx).omega_n++ = om.get<T>();
    }
    m_constants.at(storage_idx).omega_n -= 4;
    for (const auto &om : data["omega_e"]) {
      m_constants.at(storage_idx).omega_e++ = om.get<T>();
    }
    m_constants.at(storage_idx).omega_e -= 4;
    m_constants.at(storage_idx).omega_s = data["omega_s"].get<T>();
    m_constants.at(storage_idx).etta = data["etta"].get<T>();
    m_constants.at(storage_idx).idx = type_idx;
  }

  // std::vector<T> m_bfunc_der_precomputed;
  T calculateElectronDensity(const T &dist, const T &betta, T &r_e,
                             const T &lambda) const noexcept {
    return exp(-betta * (dist / r_e)) / (1 + pow(dist / r_e), 20);
  }

  T calculateElectronDensityDer(const T &precomp_exp, const T &dist,
                                const T &betta, T &r_e,
                                const T &lambda) const noexcept {
    double ratio = dist / r_e;
    double ratio_kappa = ratio - lambda;
    return (-1. / r_e) * precomp_exp *
           (betta + (20 * pow(ratio_kappa, 19)) / (1 + pow(ratio_kappa, 20)));
  }

  T calculateEmbedded(int atom_idx, int type_local_idx) const noexcept {
    const struct Constants &eam_const = m_constants.at(type_local_idx);
    T rho_zero = eam_const.rho_zero, rho_n = eam_const.rho_n,
      rho_e = eam_const.rho_e, rho_s = eam_const.rho_s;
    T total_atom_density = m_total_electron_density.at(atom_idx);
    T result = 0;
    if (total_atom_density < rho_n) {

      for (size_t k_ = 0; k_ != 4; ++k_) {

        result += ((eam_const.omega_n[k_] * k_) / rho_n) *
                  pow(total_atom_density / rho_n - 1, k_);
      }
      return result;
    }
    if (total_atom_density < rho_zero && total_atom_density >= rho_n) {

      for (size_t k_ = 0; k_ != 4; ++k_) {
        result += ((eam_const.omega_e[k_] * k_) / rho_e) *
                  pow(total_atom_density / rho_e - 1, k_);
      }
      return result;
    }
    return -eam_const.omega_s *
           (1 - eam_const.etta * log(total_atom_density / rho_s)) *
           (pow(total_atom_density / rho_s, eam_const.etta));
  }

  T calculateEmbeddedDer(int atom_idx, int type_local_idx) const noexcept {
    const struct Constants &eam_const = m_constants.at(type_local_idx);
    T rho_zero = eam_const.rho_zero, rho_n = eam_const.rho_n,
      rho_e = eam_const.rho_e, rho_s = eam_const.rho_s;
    T total_atom_density = m_total_electron_density.at(atom_idx);
    T result = 0;
    if (total_atom_density < rho_n) {
      // printf("[%d]At case total_atom_density < rho_n, and rho_n =
      // %e,total_atom_density = %e\n", m_rank,
      //        rho_n, total_atom_density);
      for (size_t k_ = 1; k_ != 4; ++k_) {

        result += ((eam_const.omega_n[k_] * k_) / rho_n) *
                  pow(total_atom_density / rho_n - 1, k_ - 1);
      }
      return result;
    }
    if (total_atom_density < rho_zero && total_atom_density >= rho_n) {

      for (size_t k_ = 1; k_ != 4; ++k_) {
        result += ((eam_const.omega_e[k_] * k_) / rho_e) *
                  pow(total_atom_density / rho_e - 1, k_ - 1);
      }
      return result;
    }

    return -eam_const.omega_s * eam_const.etta / rho_s *
           pow(total_atom_density / rho_s, eam_const.etta - 1) *
           log(total_atom_density / rho_s) * eam_const.etta;
  }

  T calculatePairHomogen(const T &dist, int type_local_idx) {
    const struct Costants &eam_const = m_constants.at(type_local_idx);

    return eam_const.A * calculateElectronDensity(dist, eam_const.r_e,
                                                  eam_const.alpha,
                                                  eam_const.kappa) -
           // HINT: here can be used b_precomputed
           eam_const.B * calculateElectronDensity(dist, eam_const.r_e,
                                                  eam_const.betta,
                                                  eam_const.lambda);
  }
  T calculatePairHomogenDer(const T &dist, int type_local_idx) {
    const struct Costants &eam_const = m_constants.at(type_local_idx);
    // HINT: here can be used b_precomputed
    T b_func_der = calculateElectronDensityDer(
        calculateElectronDensity(dist, eam_const.r_e, eam_const.betta,
                                 eam_const.lambda),
        dist, eam_const.r_e, eam_const.betta, eam_const.lambda);
    T a_func_der = calculateElectronDensityDer(
        calculateElectronDensity(dist, eam_const.r_e, eam_const.alpha,
                                 eam_const.kappa),
        dist, eam_const.r_e, eam_const.betta, eam_const.lambda);
    return eam_const.A * a_func_der - eam_const.B * b_func_der;
  }

  virtual T calculatePairPotential(const T &dist) const noexcept = 0;
  virtual T calculatePairPotentialDer(const T &dist) const noexcept = 0;

public:
  virtual ~EAMPotential() = default;
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class EAMPotentialPure : public EAMPotential<T> {
private:
  virtual T calculatePairPotential(const T &dist) const noexcept override {

    return this->calculatePairHomogen(dist, 0);
  }
  virtual T calculatePairPotentialDer(const T &dist) const noexcept override {
    return this->calculatePairHomogenDer(dist, 0);
  }

public:
  EAMPotentialPure() { this->m_constants.resize(1); }
  virtual void loadParameters(int f_type, const std::string &path,
                              int s_type) override {
    this->__loadParameters(f_type, 0, path);
  }
  virtual T computeEnergy(const Particle<T> &p1,
                          const Particle<T> &p2) const noexcept override;
  virtual Vector3x<T>
  computeForce(const Particle<T> &p1,
               const Particle<T> &p2) const noexcept override;
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class EAMPotentialAlloy : public EAMPotential<T> {
private:
  // const T &dist, const T &betta, T &r_e,
  //                              const T &lambda
  virtual T calculatePairPotential(const T &dist) const noexcept override {

    T f_e_first =
          this->m_constants.at(0).f_e *
          this->calculateElectronDensity(dist, this->m_constants.at(0).betta,
                                         this->m_constants.at(0).r_e,
                                         this->m_constants.at(0).lambda),
      f_e_second =
          this->m_constants.at(1).f_e *
          this->calculateElectronDensity(dist, this->m_constants.at(1).betta,
                                         this->m_constants.at(1).r_e,
                                         this->m_constants.at(1).lambda);

    return 0.5 * ((f_e_second / f_e_first) * this->calculatePairHomogen(0) +
                  (f_e_first / f_e_second) * this->calculatePairHomogen(1));
  }
  virtual T calculatePairPotentialDer(const T &dist,
                                      int atom_idx) const noexcept override {

    T f_e_first = this->m_constants.at(0).f_e,
      f_e_second = this->m_constants.at(0).f_e,
      dens_first = this->calculateElectronDensity(
          dist, this->m_constants.at(0).betta, this->m_constants.at(0).r_e,
          this->m_constants.at(0).lambda),
      dens_second = this->calculateElectronDensity(
          dist, this->m_constants.at(1).betta, this->m_constants.at(1).r_e,
          this->m_constants.at(1).lambda),
      dens_der_first =

          this->calculateElectronDensityDer(
              dens_first, dist, this->m_constants.at(0).betta,
              this->m_constants.at(0).r_e, this->m_constants.at(0).lambda),
      dens_der_second =

          this->calculateElectronDensityDer(
              dens_second, dist, this->m_constants.at(1).betta,
              this->m_constants.at(1).r_e, this->m_constants.at(1).lambda),
      f_pair_pot_homogen = this->calculatePairHomogen(0),
      s_pair_pot_homogen = this->calculatePairHomogen(1),
      f_pair_pot_der_homogen = this->calculatePairHomogenDer(0),
      s_pair_pot_der_homogen = this->calculatePairHomogenDer(1);

    return (f_e_second / f_e_first) *
               (f_pair_pot_homogen / (dens_first) *
                    (dens_der_second -
                     (dens_second * dens_der_first) / dens_first) +
                dens_second / dens_first * f_pair_pot_der_homogen) +
           (f_e_first / f_e_second) *
               (s_pair_pot_homogen / (dens_second) *
                    (dens_der_first -
                     dens_first * dens_der_second / dens_second) +
                dens_der_first / dens_second * s_pair_pot_der_homogen);
  }

public:
  EAMPotentialAlloy() { this->m_constants.resize(2); }
  virtual void loadParameters(int f_type, int s_type,
                              const std::string &path) override {
    this->__loadParameters(f_type, 0, path);
    this->__loadParameters(s_type, 1, path);
  }
  virtual T computeEnergy(const Particle<T> &p1,
                          const Particle<T> &p2) const noexcept override;
  virtual Vector3x<T>
  computeForce(const Particle<T> &p1,
               const Particle<T> &p2) const noexcept override;
};