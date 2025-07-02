#pragma once
#include <IMDManager.hpp>
#include <Utility.hpp>
#include <memory>
#include <type_traits>
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class Barostat {
protected:
  std::unique_ptr<IMDManager<T>> m_parent;

public:
  enum class Type { Unitary, Andersen, Berendsen, Langevin, MTK, PR };
  Barostat(IMDManager<T> *mdm) : m_parent(mdm) {}
  virtual void applyBarostat() = 0;
};
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class UnitaryBarostat final : public Barostat<T> {

  UnitaryBarostat(IMDManager<T> *mdm) : Barostat<T>(mdm) {}
  virtual void applyBarostat() override {}
};

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>,
          typename... Args>
std::unique_ptr<Barostat<T>> createBarostat(typename Barostat<T>::Type &,
                                            Args &&...args);
