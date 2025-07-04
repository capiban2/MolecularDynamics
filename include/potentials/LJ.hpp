#pragma once
#include "Potential.hpp"
#include <cmath>
#include <fstream>
template <typename T> class LJPotential : public Potential<T> {

  T m_sigma, m_epsilon;

public:
  ~LJPotential() = default;
  virtual void loadParameters(int f_type, const std::string &path,
                              int s_type) noexcept(false) override {

    using json = nlohmann::json;
    if (!std::filesystem::exists(path))
      throw std::runtime_error(std::string("Provided file ") + path +
                               std::string(" does not exist!"));
    std::ifstream _file(path);
    json data = json::parse(_file)[std::to_string(f_type)];
    m_sigma = data["sigma"].get<T>();
    m_epsilon = data["epsilon"].get<T>();
  }
  virtual T computeEnergy(const Particle<T> &p1,
                          const Particle<T> &p2) const noexcept override {

    T dist = length(p1.pos - p2.pos);
    return 4 * m_epsilon * (pow(m_sigma / dist, 12) - pow(m_sigma / dist, 6));
  }
  virtual Vector3x<T>
  computeForce(const Particle<T> &p1, const Particle<T> &p2,
               const std::pair<int, int> &indx = {}) const noexcept override {
    T dist = length(p1.pos - p2.pos);

    return (p1.pos - p2.pos) *
           (-24 * m_epsilon / dist *
            (2 * pow(m_sigma / dist, 12) - pow(m_sigma / dist, 6)));
  }
};