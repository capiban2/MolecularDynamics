#pragma once
#include "CellList.hpp"
#include "Utility.hpp"
#include <cmath>
#include <iterator>
#include <vector>

template <typename T> class VerletList {
  T m_act_cutoff, m_cutoff2;

  T m_skin_half_sq;
  T m_skin;
  std::vector<std::vector<int>> m_nbrs;

  std::vector<std::vector<T>> m_all_cutoffs;

  std::vector<T> m_last_x, m_last_y, m_last_z;

  int m_calls_since_build = 0;

  T dt;

  // HINT: this thing got called every 5-10 calls
  bool __needsUpdateVelocityBased(const std::vector<T> &x,
                                  const std::vector<T> &y,
                                  const std::vector<T> &z,
                                  const std::vector<T> &vx,
                                  const std::vector<T> &vy,
                                  const std::vector<T> &vz) const;

  bool __needsUpdatePositionBased(const std::vector<T> &x,
                                  const std::vector<T> &y,
                                  const std::vector<T> &z) const;

public:
  VerletList(double cutoff, double skin,
             const std::vector<std::vector<T>> &all_cutoffs)
      : m_act_cutoff(cutoff + skin),
        m_cutoff2((cutoff + skin) * (cutoff + skin)),
        m_all_cutoffs(all_cutoffs), m_skin_half_sq((skin * 0.5) * (skin * 0.5)),
        m_skin(skin) {}
  void build(const ParticleData<T> &particles, const CellList<T> &clist,
             int ghost_index_first);

  const std::vector<int> &neighbours(int idx) const { return m_nbrs.at(idx); }

  const std::vector<std::vector<T>> &getCutoffs() const noexcept {
    return m_all_cutoffs;
  }

  bool needsUpdate(const std::vector<T> &x, const std::vector<T> &y,
                   const std::vector<T> &z, const std::vector<T> &vx,
                   const std::vector<T> &vy, const std::vector<T> &vz);
};

template <typename T>
void VerletList<T>::build(const ParticleData<T> &particles,
                          const CellList<T> &clist, int ghost_index_first) {

  m_nbrs.clear();
  m_nbrs.resize(ghost_index_first);
  m_last_x.resize(ghost_index_first);
  m_last_y.resize(ghost_index_first);
  m_last_z.resize(ghost_index_first);
  std::copy(particles.pos_x.begin(),
            particles.pos_x.begin() + ghost_index_first, m_last_x);
  std::copy(particles.pos_y.begin(),
            particles.pos_y.begin() + ghost_index_first, m_last_y);
  std::copy(particles.pos_z.begin(),
            particles.pos_z.begin() + ghost_index_first, m_last_z);

  for (int i = 0; i < ghost_index_first; ++i) {
    const auto &pi = particles[i];
    int cell_id = clist.getCellIndex(
        {particles.pos_x[i], particles.pos_y[i], particles.pos_z[i]});
    auto neighbors = clist.getNeighborCells(cell_id);

    for (int cid : neighbors) {
      for (int j : clist.getCell(cid)) {
        if (j <= i)
          continue;

        T dx = particles.pos_x[i] - particles.pos_x[j];
        T dy = particles.pos_y[i] - particles.pos_y[j];
        T dz = particles.pos_z[i] - particles.pos_z[j];
        T r2 = dx * dx + dy * dy + dz * dz;

        if (r2 < m_cutoff2) {
          m_nbrs[i].push_back(j);

          // HINT: think about it a great deal
          m_nbrs[j].push_back(i);
        }
      }
    }
  }
  m_calls_since_build = 0;
}
template <typename T>
bool VerletList<T>::needsUpdate(const std::vector<T> &x,
                                const std::vector<T> &y,
                                const std::vector<T> &z,
                                const std::vector<T> &vx,
                                const std::vector<T> &vy,
                                const std::vector<T> &vz) {
  if ((m_calls_since_build + 1) % 5 == 0) {
    if (__needsUpdateVelocityBased(x, y, z, vx, vy, vz))
      return true;
    m_calls_since_build++;
    return false;
  }
  if (__needsUpdatePositionBased(x, y, z))
    return true;
  m_calls_since_build++;
  return false;
}

template <typename T>
bool VerletList<T>::__needsUpdateVelocityBased(const std::vector<T> &x,
                                               const std::vector<T> &y,
                                               const std::vector<T> &z,
                                               const std::vector<T> &vx,
                                               const std::vector<T> &vy,
                                               const std::vector<T> &vz) const {
  auto len = x.size();
  for (size_t i = 0; i < len; ++i) {
    // Current displacement vector
    T dx = x[i] - m_last_x[i];
    T dy = y[i] - m_last_y[i];
    T dz = z[i] - m_last_z[i];

    T current_disp = std::sqrt(dx * dx + dy * dy + dz * dz);

    // Velocity vector
    T par_vx = vx[i];
    T par_vy = vy[i];
    T par_vz = vz[i];

    // Project velocity onto displacement direction
    if (current_disp > 1e-10) { // Avoid division by zero
      T dot_product = (dx * vx + dy * vy + dz * vz) / current_disp;

      // If moving in same direction as displacement, predict faster growth
      if (dot_product > 0) {
        T predicted_disp = current_disp + dot_product * dt;
        if (predicted_disp > m_skin * 0.5) {
          return true;
        }
      }
    }

    // Fallback to magnitude-based prediction
    T v_mag = std::sqrt(vx * vx + vy * vy + vz * vz);
    T total_predicted = current_disp + v_mag * dt;
    if (total_predicted > m_skin * 0.5) {
      return true;
    }
  }

  return false;
}
template <typename T>
bool VerletList<T>::__needsUpdatePositionBased(const std::vector<T> &x,
                                               const std::vector<T> &y,
                                               const std::vector<T> &z) const {

  // HINT: here iteration only to firstGhostIndex, were not intrested in ghosts
  // displacements
#ifdef OPENMP_ENABLED
#pragma omp simd
#endif
  for (int t_ = 0; t_ != m_last_x.size(); ++t_) {
    T dx = x[t_] - m_last_x[t_];
    T dy = y[t_] - m_last_y[t_];
    T dz = z[t_] - m_last_z[t_];
    T disp_sq = dx * dx + dy * dy + dz * dz;
    if (disp_sq > m_skin_half_sq)
      return true;
  }
  return false;
}