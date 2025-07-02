#pragma once
#include "Potential.hpp"
#include <cmath>
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
  std::vector<Constants> m_constants;

  std::vector<T> m_total_electron_density;
  std::vector<T> m_bfunc_precomputed;
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

  T calculatePairHomogen(const T &dist, int atom_idx, int type_local_idx) {
    const struct Costants &eam_const = m_constants.at(type_local_idx);

    return eam_const.A * calculateElectronDensity(dist, eam_const.r_e,
                                                  eam_const.alpha,
                                                  eam_const.betta) -
           eam_const.B * m_bfunc_precomputed.at(atom_idx);
  }
  T calculatePairHomogenDer(const T &dist, int atom_idx, int type_local_idx) {
    const struct Costants &eam_const = m_constants.at(type_local_idx);

    T b_func_der = calculateElectronDensityDer(
        m_bfunc_precomputed.at(atom_idx), dist, eam_const.r_e, eam_const.betta,
        eam_const.lambda);
    T a_func_der = calculateElectronDensityDer(
        calculateElectronDensity(dist, eam_const.r_e, eam_const.alpha,
                                 eam_const.kappa),
        dist, eam_const.r_e, eam_const.betta, eam_const.lambda);
    return eam_const.A * a_func_der - eam_const.B * b_func_der;
  }

  virtual T calculatePairPotential(const T &dist,
                                   int atom_idx) const noexcept = 0;
  virtual T calculatePairPotentialDer(const T &dist,
                                      int atom_idx) const noexcept = 0;

public:
  virtual ~EAMPotential() = default;
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class EAMPotentialPure : public EAMPotential<T> {
private:
  virtual T calculatePairPotential(const T &dist,
                                   int atom_idx) const noexcept override {
    return this->calculatePairHomogen(dist, atom_idx, 0);
  }
  virtual T calculatePairPotentialDer(const T &dist,
                                      int atom_idx) const noexcept override {
    return this->calculatePairHomogenDer(dist, atom_idx, 0);
  }

public:
  EAMPotentialPure() = default;
  virtual void loadParameters(int f_type, int s_type,
                              const std::string &path) override;
  virtual T computeEnergy(const Particle<T> &p1,
                          const Particle<T> &p2) const noexcept override;
  virtual Vector3x<T>
  computeForce(const Particle<T> &p1,
               const Particle<T> &p2) const noexcept override;
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class EAMPotentialAlloy : public EAMPotential<T> {
private:
  virtual T calculatePairPotential(const T &dist,
                                   int atom_idx) const noexcept override {
                                    
                                   }
  virtual T calculatePairPotentialDer(const T &dist,
                                      int atom_idx) const noexcept override {}

public:
  EAMPotentialAlloy() = default;
  virtual void loadParameters(int f_type, int s_type,
                              const std::string &path) override;
  virtual T computeEnergy(const Particle<T> &p1,
                          const Particle<T> &p2) const noexcept override;
  virtual Vector3x<T>
  computeForce(const Particle<T> &p1,
               const Particle<T> &p2) const noexcept override;
};