#pragma once
#include "CellList.hpp"
#include "Utility.hpp"
#include <cmath>
#include <vector>

template <typename T> class VerletList {
  double m_act_cutoff, m_cutoff2;

  std::vector<std::vector<int>> m_nbrs;

public:
  VerletList(double cutoff, double skin)
      : m_act_cutoff(cutoff + skin),
        m_cutoff2((cutoff + skin) * (cutoff + skin)) {}
  void build(const std::vector<Particle<T>> &particles,
             const CellList<T> &clist, int ghost_index_first);

  const std::vector<int> &neighbours(int idx) const { return m_nbrs.at(idx); }
};

template <typename T>
void VerletList<T>::build(const std::vector<Particle<T>> &particles,
                          const CellList<T> &clist, int ghost_index_first) {

  m_nbrs.clear();
  m_nbrs.resize(ghost_index_first);

  for (int i = 0; i < ghost_index_first; ++i) {
    const auto &pi = particles[i];
    int cell_id = clist.getCellIndex(pi.pos);
    auto neighbors = clist.getNeighborCells(cell_id);

    for (int cid : neighbors) {
      for (int j : clist.getCell(cid)) {
        if (j <= i)
          continue;

        const auto &pj = particles[j];

        T dx = pi.pos[0] - pj.pos[0];
        T dy = pi.pos[1] - pj.pos[1];
        T dz = pi.pos[2] - pj.pos[2];
        T r2 = dx * dx + dy * dy + dz * dz;

        if (r2 < m_cutoff2) {
          m_nbrs[i].push_back(j);
        }
      }
    }
  }
}