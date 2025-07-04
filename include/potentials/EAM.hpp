#pragma once
#include "Potential.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <numbers>
#include <ranges>
#include <stdexcept>
#ifndef TABULATED

constexpr int chebyshev_nodes_count = 10000;
#endif
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class EAMPotential : public Potential<T> {

protected:
#ifndef TABULATED

  std::vector<T> m_chebyshev_nodes;
  std::vector<T> m_precomputed_pair;
  std::vector<T> m_precomputed_pair_der;

#endif
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
                               std::string(" does not exist!"));
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

#ifndef TABULATED
  virtual void precomputePairPotential(T r_cut) {

    m_chebyshev_nodes = std::views::iota(chebyshev_nodes_count) |
                        std::views::transform([&](int _v) {
                          T x_k = cos((2 * _v + 1) * std::numbers::pi /
                                      (2 * chebyshev_nodes_count));
                          return 0.5 * r_cut - 0.5 * r_cut * x_k;
                        }) |
                        to<std::vector<T>>();
    m_precomputed_pair =
        m_chebyshev_nodes | std::views::transform([this](const T &_node) {
          return this->calculatePairPotential(_node);
        });
    m_precomputed_pair_der =
        m_chebyshev_nodes | std::views::transform([this](const T &_node) {
          return this->calculatePairPotentialDer(_node);
        });
  }

  // TODO: verify all these cubd sqrd
  T __calculateHermiteInterpolantPairPot(const T &h, const T &t,
                                         int idx) const noexcept {
    T t_sqrd = pow(t, 2), t_cubd = pow(t, 3);
    T h00 = 2 * t_cubd - 3 * t_sqrd + 1, h10 = t_cubd - 2 * t_sqrd + t,
      h01 = -2 * t_cubd + 3 * t_sqrd, h11 = t_cubd - t_sqrd;
    T result = h00 * m_precomputed_pair[idx] +
               h10 * h * m_precomputed_pair_der[idx] +
               h01 * m_precomputed_pair[idx + 1] +
               h11 * h * m_precomputed_pair_der[idx + 1];
    return result;
  }
  T __calculateHermiteInterpolantPairPotDer(const T &h, const T &t,
                                            int idx) const noexcept {
    T t_sqrd = pow(t, 2);
    T h00 = 6 * t_sqrd - 6 * t;
    T h10 = 3 * t_sqrd - 4 * t + 1;
    T h01 = -6 * t_sqrd + 6 * t;
    T h11 = 3 * t_sqrd - 2 * t;
    T result = h00 * m_precomputed_pair[idx] +
               h10 * h * m_precomputed_pair_der[idx] +
               h01 * m_precomputed_pair[idx + 1] +
               h11 * h * m_precomputed_pair_der[idx + 1];
    return result;
  }
  // TODO: make r_cut attribute of the class or smth like that for not asking it
  // everywhere without a reason

  T interpolatePairPotential(const T &dist, T r_cut) const noexcept {
    // HINT: map arbitary value from (0, r_cut) to [-1,1] that chebyshev's node
    // can work with
    T mapped = (2 * dist - r_cut) / r_cut;

    // HINT: dont do any assertion because as it seems any value will be inside
    // given range, because they all in (0, r_cut)
    int first_greater_idx =
        std::distance(m_chebyshev_nodes.begin(),
                      std::upper_bound(m_chebyshev_nodes.begin(),
                                       m_chebyshev_nodes.end(), mapped));
    T step = m_chebyshev_nodes[first_greater_idx] -
             m_chebyshev_nodes[first_greater_idx - 1],
      step2 = mapped - m_chebyshev_nodes[first_greater_idx - 1];
    return __calculateHermiteInterpolantPairPot(step, step2, first_greater_idx);
  }
  T interpolatePairPotentialDer(const T &dist, T r_cut) const noexcept {
    // HINT: map arbitary value from (0, r_cut) to [-1,1] that chebyshev's node
    // can work with
    T mapped = (2 * dist - r_cut) / r_cut;

    // HINT: dont do any assertion because as it seems any value will be inside
    // given range, because they all in (0, r_cut)
    int first_greater_idx =
        std::distance(m_chebyshev_nodes.begin(),
                      std::upper_bound(m_chebyshev_nodes.begin(),
                                       m_chebyshev_nodes.end(), mapped));
    T step = m_chebyshev_nodes[first_greater_idx] -
             m_chebyshev_nodes[first_greater_idx - 1],
      step2 = mapped - m_chebyshev_nodes[first_greater_idx - 1];
    return __calculateHermiteInterpolantPairPotDer(step, step2,
                                                   first_greater_idx);
  }
#endif

public:
  virtual ~EAMPotential() = default;
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class EAMPotentialPure : public EAMPotential<T> {
private:
#ifndef TABULATED
  virtual void precomputePairPotential(T r_cut) override {

    EAMPotential<T>::precomputePairPotential(r_cut);
    this->m_precomputed_pair =
        this->m_chebyshev_nodes | std::views::transform([this](const T &_node) {
          return this->calculatePairHomogen(_node, 0);
        });
    this->m_precomputed_pair_der =
        this->m_chebyshev_nodes | std::views::transform([this](const T &_node) {
          return this->calculatePairHomogenDer(_node, 0);
        });
  }
#endif
  virtual T calculatePairPotential(const T &dist) const noexcept override {

#ifdef TABULATED
    return this->interpolatePairPotential(dist, this->m_constants.at(0).r_cut);
#else

    return this->calculatePairHomogen(dist, 0);
#endif
  }
  virtual T calculatePairPotentialDer(const T &dist) const noexcept override {
#ifdef TABULATED
    return this->interpolatePairPotentialDer(dist,
                                             this->m_constants.at(0).r_cut);
#else

    return this->calculatePairHomogenDer(dist, 0);
#endif
  }

public:
  EAMPotentialPure() { this->m_constants.resize(1); }
  virtual void loadParameters(int f_type, const std::string &path,
                              int s_type) noexcept(false) override {
    this->__loadParameters(f_type, 0, path);
  }
  virtual T computeEnergy(const Particle<T> &p1,
                          const Particle<T> &p2) const noexcept override {
    return T{};
  }
  virtual Vector3x<T>
  computeForce(const Particle<T> &p1, const Particle<T> &p2,
               const std::pair<int, int> &indx) const noexcept override {
    T dist = length(p1.pos - p2.pos);
    const struct Constants &eam_const = this->m_constants.at(0);

#ifdef WITH_SMOOTHING
    T smooth, der_smooth;
    this->taperQuinticDer(dist, smooth, der_smooth);

    T pair_part = this->calculatePairPotential(dist);
    T elect_density = this->calculateElectronDensity(
        dist, eam_const.betta, eam_const.r_e, eam_const.lambda);
    T embedded_part_der =
        (this->calculateEmbeddedDer(indx.first, 0) +
         this->calculateEmbeddedDer(indx.second, 1)) *
        (this->calculateElectronDensityDer(elect_density, dist, eam_const.betta,
                                           eam_const.r_e, eam_const.lambda) *
             smooth +
         elect_density * der_smooth

        );
    T pair_part_der =
        this->calculatePairPotentialDer(dist) * smooth + pair_part * der_smooth;
#else

    T embedded_part_der =
        (this->calculateEmbeddedDer(indx.first, 0) +
         this->calculateEmbeddedDer(indx.second, 1)) *
        this->calculateElectronDensityDer(
            this->calculateElectronDensity(dist, eam_const.betta, eam_const.r_e,
                                           eam_const.lambda),
            dist, eam_const.betta, eam_const.r_e, eam_const.lambda);
    T pair_part_der = this->calculatePairPotentialDer(dist);

    return (p1.pos - p2.pos) * -(embedded_part_der + pair_part_der);
#endif
  }
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class EAMPotentialAlloy : public EAMPotential<T> {
private:
  T __calculatePairHeterogen(const T &dist) const noexcept {
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

    return 0.5 *
           ((f_e_second / f_e_first) * this->calculatePairHomogen(dist, 0) +
            (f_e_first / f_e_second) * this->calculatePairHomogen(dist, 1));
  }
  T __calculatePairHeterogenDer(const T &dist) const noexcept {
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
      f_pair_pot_homogen = this->calculatePairHomogen(dist, 0),
      s_pair_pot_homogen = this->calculatePairHomogen(dist, 1),
      f_pair_pot_der_homogen = this->calculatePairHomogenDer(dist, 0),
      s_pair_pot_der_homogen = this->calculatePairHomogenDer(dist, 1);

    return 0.5 * (f_e_second / f_e_first) *
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
#ifndef TABULATED
  virtual void precomputePairPotential(T r_cut) override {

    EAMPotential<T>::precomputePairPotential(r_cut);
#ifdef WITH_SMOOTHING
    this->m_precomputed_pair =
        this->m_chebyshev_nodes | std::views::transform([this](const T &_node) {
          return this->calculatePairHeterogen(_node);
        });
#endif
    this->m_precomputed_pair_der =
        this->m_chebyshev_nodes | std::views::transform([this](const T &_node) {
          return this->calculatePairHeterogenDer(_node);
        });
  }
#endif

  // TODO: check this r_cut here, or utilize it, because cannot be 2 different
  // r_cuts
  virtual T calculatePairPotential(const T &dist) const noexcept override {
#ifdef TABULATED
    return this->interpolatePairPotential(dist, this->m_constants.at(0).r_cut);
#else

    return __calculatePairHeterogen(dist);
#endif
  }
  virtual T calculatePairPotentialDer(const T &dist) const noexcept override {
#ifdef TABULATED
    return this->interpolatePairPotentialDer(dist,
                                             this->m_constants.at(0).r_cut);
#else
    return __calculatePairHeterogenDer(dist);
#endif
  }

public:
  EAMPotentialAlloy() { this->m_constants.resize(2); }
  virtual void
  loadParameters(int f_type, int s_type,
                 const std::string &path) noexcept(false) override {
    this->__loadParameters(f_type, 0, path);
    this->__loadParameters(s_type, 1, path);
  }
  virtual T computeEnergy(const Particle<T> &p1,
                          const Particle<T> &p2) const noexcept override {
    return T{};
  }
  virtual Vector3x<T>
  computeForce(const Particle<T> &p1, const Particle<T> &p2,
               const std::pair<int, int> &indx) const noexcept override {
    T dist = length(p1.pos - p2.pos);
    const struct Constants &eam_const_ftype = this->m_constants.at(0);
    const struct Constants &eam_const_stype = this->m_constants.at(1);
#ifdef WITH_SMOOTHING
    T smooth, der_smooth;
    this->taperQuinticDer(dist, smooth, der_smooth);

    T pair_part = this->calculatePairPotential(dist);
    T elect_density_f = this->calculateElectronDensity(
          dist, eam_const_ftype.betta, eam_const_ftype.r_e,
          eam_const_ftype.lambda),
      elect_density_s = this->calculateElectronDensity(
          dist, eam_const_stype.betta, eam_const_stype.r_e,
          eam_const_stype.lambda);

    T embedded_part_der = (this->calculateEmbeddedDer(indx.first, 0) *
                               (this->calculateElectronDensityDer(
                                    elect_density_s, dist, eam_const_f.betta,
                                    eam_const_f.r_e, eam_const_f.lambda) *
                                    smooth +
                                elect_density_f * der_smooth) +
                           this->calculateEmbeddedDer(indx.second, 1)) *
                          (this->calculateElectronDensityDer(
                               elect_density_f, dist, eam_const_f.betta,
                               eam_const_f.r_e, eam_const_f.lambda) *
                               smooth +
                           elect_density_f * der_smooth

                          );
    T pair_part_der =
        this->calculatePairPotentialDer(dist) * smooth + pair_part * der_smooth;
#else
    T embedded_part_der =
        (this->calculateEmbeddedDer(indx.first, 0) *
             this->calculateElectronDensityDer(
                 this->calculateElectronDensity(dist, eam_const_stype.betta,
                                                eam_const_stype.r_e,
                                                eam_const_stype.lambda),
                 dist, eam_const_stype.betta, eam_const_stype.r_e,
                 eam_const_stype.lambda)

         + this->calculateEmbeddedDer(indx.second, 1) *
               this->calculateElectronDensityDer(
                   this->calculateElectronDensity(dist, eam_const_ftype.betta,
                                                  eam_const_ftype.r_e,
                                                  eam_const_ftype.lambda),
                   dist, eam_const_ftype.betta, eam_const_ftype.r_e,
                   eam_const_ftype.lambda)

        );
    T pair_part_der = this->calculatePairPotentialDer(dist);

    return (p1.pos - p2.pos) * -(embedded_part_der + pair_part_der);
#endif
  }
};