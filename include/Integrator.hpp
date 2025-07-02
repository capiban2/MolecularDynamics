#pragma once
#include "PairPotential.hpp"
#include "Utility.hpp"
#include "VerletList.hpp"

template <typename T> class Integrator {
public:
  virtual ~Integrator() {}

  virtual void step(std::vector<Particle<T>> &particles,
                    PairPotential<T, 0> &potential, VerletList<T> &vlist,
                    double dt) = 0;
};

template <typename T> class VelocityVerlet : Integrator<T> {
public:
  virtual void step(std::vector<Particle<T>> &particles,
                    PairPotential<T, 0> &potential, VerletList<T> &vlist,
                    double dt) override;
};
