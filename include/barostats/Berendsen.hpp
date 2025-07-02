#pragma once
#include "Barostat.hpp"

template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class BerendsenBarostat final : public Barostat<T> {
  T m_ratio;
  T m_tar_pressure;

public:
  BerendsenBarostat(IMDManager<T> *mdm, T ratio, T target_pressure)
      : Barostat<T>(mdm), m_ratio(ratio), m_tar_pressure(target_pressure) {}
  void applyBarostat() override;
};