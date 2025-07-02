#pragma once
#include "Barostat.hpp"

template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class MTKBarostat final : public Barostat<T> {

public:
  // TODO: add necessary parameters
  MTKBarostat(IMDManager<T> *mdm) : Barostat<T>(mdm) {}
  void applyBarostat() override;
};