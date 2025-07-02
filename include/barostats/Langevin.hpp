#pragma once
#include "Barostat.hpp"

template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class LangevinBarostat final : public Barostat<T> {

public:
  // TODO: add necessary parameters
  LangevinBarostat(IMDManager<T> *mdm) : Barostat<T>(mdm) {}
  void applyBarostat() override;
};