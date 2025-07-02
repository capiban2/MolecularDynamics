#pragma once
#include "Potential.hpp"
template <typename T> class LJPotential : public Potential<T> {

  T sigma, epsilon;

public:
  ~LJPotential() = default;
  virtual void loadParameters(int f_type, const std::string &path,
                              int s_type) override;
};