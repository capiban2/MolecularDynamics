#pragma once
#include "Barostat.hpp"
template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
class ParinelloRahmanBarostat final : public Barostat<T> {

  T m_p;
  T m_bulk_modulus;

public:
  ParinelloRahmanBarostat(IMDManager<T> *mdm, T p, T bulk_mod)
      : Barostat<T>(mdm), m_p(p), m_bulk_modulus(bulk_mod) {}

  void applyBarostat() override;
};