#pragma once

#include "Utility.hpp"
#include <VerletList.hpp>
#include <array>
#include <cmath>
#include <memory>
#include <nlohmann/json.hpp>
#include <types/PairInteraction.hpp>

enum class PotentialType { LJ, EAM };

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class Potential {

public:
  Potential() = default;
  virtual ~Potential() = default;

  virtual T getCutoffRadius() const noexcept = 0;

  virtual std::string getDescription() const noexcept = 0;

  virtual void configure() = 0;

  // HINT: for multipotential modelling
  virtual void
  computeBulkForces(const std::vector<MD::PairInteraction<T>> &_interactions,
                    std::vector<T> &_f_x, std::vector<T> &_f_y,
                    std::vector<T> &_f_z, const std::vector<int> types,
                    int _first_ghost_idx) = 0;

  // HINT: for monopotential modelling
  virtual void computeBulkForces(const VerletList<T> &t, ParticleData<T> &_data,
                                 int _first_ghost_idx) = 0;

  virtual enum PotentialType rtti() = 0;
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>,
          typename... Args>
std::shared_ptr<Potential<T>>
createPotential(const typename Potential<T>::Type &type, Args &&...args);

#if 0
template <typename T, int types_count> class PairPotential {
  using PoFunc = T (*)(double);
  static T computeLennardJones(double);
  template <int f_type, int s_type> T static computeEAM(double);

  std::array<std::array<PoFunc, types_count>, types_count> potentials;

public:
  enum PotentialName { EAM = 0, LJ };
  virtual ~PairPotential(){};
  template <PotentialName potential_name, int f_type, int s_type>
  void set_potential();

  virtual T compute(std::vector<Particle<T>> &particles,
                    const VerletList<T> &vlist) = 0;
};

template <typename T, int types_count>
template <typename PairPotential<T, types_count>::PotentialName potential_name,
          int f_type, int s_type>
void PairPotential<T, types_count>::set_potential() {

  if constexpr (potential_name == PairPotential::EAM) {
    potentials[f_type][s_type] = &PairPotential::computeEAM<f_type, s_type>;
  }

  if constexpr (potential_name == PairPotential::LJ) {
    potentials[f_type][s_type] = &PairPotential::computeLennardJones;
  }
}
#endif
#if 0
template <typename T> class LennardJones : public PairPotential<T> {
  double m_eps, m_sigma6;

public:
  LennardJones(double eps, double sigma)
      : m_eps(eps), m_sigma6(std::pow(sigma, 6.0)) {}

  T compute(std::vector<Particle<T>> &particles,
            const VerletList<T> &vlist) override;
};

template <typename T> class EAM : public PairPotential<T> {
public:
  EAM() {}

  T compute(std::vector<Particle<T>> &particles,
            const VerletList<T> &vlist) override;
};
#endif
