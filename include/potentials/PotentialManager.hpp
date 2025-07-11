#pragma once
#include <IMDManager.hpp>
#include <VerletList.hpp>
#include <memory>
#include <potentials/Potential.hpp>
#include <type_traits>
#include <types/PairInteraction.hpp>
#include <utility>
#include <vector>
namespace MD {
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>

class PotentialManager {

  struct PotentialGroup {
    std::unique_ptr<Potential<T>> potential;
    std::vector<MD::PairInteraction<T>> interactions;
    std::vector<std::pair<int, int>> type_pairs;
  };

  std::vector<PotentialGroup> m_potential_groups;
  std::vector<std::vector<int>> m_t2t_map;

  std::unique_ptr<IMDManager<T>> m_mdmanager;
  void buildTypeMapping() {

    int max_types = 0;
    for (const auto &group : m_potential_groups) {
      for (const auto &type_pair : group.type_pairs) {
        max_types = std::max(max_types,
                             std::max(type_pair.first, type_pair.second) + 1);
      }
    }

    // Initialize the 2D mapping array
    m_t2t_map.resize(max_types);
    for (int i = 0; i < max_types; ++i) {
      m_t2t_map[i].resize(max_types, -1); // -1 means no potential assigned
    }
    // Fill the mapping
    for (size_t group_idx = 0; group_idx < m_potential_groups.size();
         ++group_idx) {
      for (const auto &type_pair : m_potential_groups[group_idx].type_pairs) {
        m_t2t_map[type_pair.first][type_pair.second] = group_idx;
        m_t2t_map[type_pair.second][type_pair.first] = group_idx; // Symmetric
      }
    }
  }

public:
  PotentialManager(std::unique_ptr<IMDManager<T>> parent)
      : m_mdmanager(parent) {}
  PotentialManager() = delete;
  void addPotential(std::unique_ptr<Potential<T>> potential,
                    const std::vector<std::pair<int, int>> &type_pairs) {

    PotentialGroup group;
    group.potential = std::move(potential);
    group.type_pairs = type_pairs;
    m_potential_groups.emplace_back(std::move(group));
  }
  void computeForces(const VerletList<T> &neigh_list, ParticleData<T> &_data) {
#ifdef OPENM_ENABLED
#pragma omp simd
#endif
    for (int t_ = 0; t_ != _data.force_x.size(); ++t_) {
      _data.force_x[t_] = 0.0;
      _data.force_y[t_] = 0.0;
      _data.force_z[t_] = 0.0;
    }

    // HINT: mono-potential case
    if (m_potential_groups.size() == 1) {
      m_potential_groups[0].potential->computeBulkForces(
          _data.force_x, _data.force_y, _data.force_z, _data.m_type_id);
      return;
    }

    // HINT: multi-potential case

    // HINT: if first iteration after particles redistribution
    if (m_mdmanager->isFirstAfterRebuild) {
    }
  }
  void setup() { buildTypeMapping(); }
};

} // namespace MD
