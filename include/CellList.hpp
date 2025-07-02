#pragma once
#include "Utility.hpp"
#include <algorithm>
#include <vector>

template <typename T> class CellList {

  Vector3i m_ncell;
  Vector3d m_cell_size;

  // HINT: if other tasks will arise in future where it is necessary to have
  // anisotropic hence triclinic space organization, nowadays it would just
  // unnecessary increase complexity
#if 0
  // HINT: lower corner
  Vector3x<T> origin;

  // HINT: Lattice vectors spanning the cell:
  //    the three edges of your (possibly triclinic) parallelepiped.
  Vector3x<T> a, b, c;
#endif

  std::vector<std::vector<int>> m_cell_content;
  int flatten(int i, int j, int k);

public:
  CellList(const Vector3x<T> &lo, const Vector3x<T> &hi,

           Vector3d cell_size, T ghost_width);

  void build(const std::vector<Particle<T>> &particles,
             const ProcSubGrid3D<T> &simbox, T ghost_width);

  const std::vector<int> &cellContents(int cell_index) const;

  std::vector<int> neighboringCells(int cell_index) const;

  int cellIndexOf(const Vector3x<T> &pos, const ProcSubGrid3D<T> &simbox) const;
};

template <typename T>
CellList<T>::CellList(const Vector3x<T> &lo, const Vector3x<T> &hi,

                      Vector3d cell_size, T ghost_width)
    : m_cell_size(cell_size) {

  auto boxsize = hi - lo;
  int n_ghost_x = std::ceil(ghost_width / m_cell_size.x);
  int n_ghost_y = std::ceil(ghost_width / m_cell_size.y);
  int n_ghost_z = std::ceil(ghost_width / m_cell_size.z);

  // TODO: check that
  m_ncell[0] = boxsize.x / cell_size.x + 2 * n_ghost_x;
  m_ncell[1] = boxsize.y / cell_size.y + 2 * n_ghost_y;
  m_ncell[2] = boxsize.z / cell_size.z + +2 * n_ghost_z;

  int total_cells = m_ncell[0] * m_ncell[1] * m_ncell[2];
  m_cell_content.resize(total_cells);
}

// TODO: fix this, need to multiply simbox.lo.whatever to
// cell_size.whatever
template <typename T>
void CellList<T>::build(const std::vector<Particle<T>> &particles,
                        const ProcSubGrid3D<T> &simbox, T ghost_width) {
  for (auto &cell : m_cell_content)
    cell.clear();

  int n_ghost_x = std::ceil(ghost_width / m_cell_size.x);
  int n_ghost_y = std::ceil(ghost_width / m_cell_size.y);
  int n_ghost_z = std::ceil(ghost_width / m_cell_size.z);

  Vector3x<T> extended_lo = simbox.lo - Vector3x<T>(n_ghost_x * m_cell_size.x,
                                                    n_ghost_y * m_cell_size.y,
                                                    n_ghost_z * m_cell_size.z);

  for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
    const auto &p = particles[i];
    int cx = (p.pos[0] - extended_lo.x) / m_cell_size[0];
    int cy = (p.pos[1] - extended_lo.y) / m_cell_size[1];
    int cz = (p.pos[2] - extended_lo.z) / m_cell_size[2];

    cx = std::clamp(cx, 0, m_ncell[0] - 1);
    cy = std::clamp(cy, 0, m_ncell[1] - 1);
    cz = std::clamp(cz, 0, m_ncell[2] - 1);
    int idx = flatten(cx, cy, cz);
    m_cell_content[idx].push_back(i);
  }
}

template <typename T> int CellList<T>::flatten(int i, int j, int k) {
  return (k * m_ncell[1] + j) * m_ncell[0] + i;
}

template <typename T>
const std::vector<int> &CellList<T>::cellContents(int index) const {
  return m_cell_content[index];
}
// TODO: fix this, need to multiply simbox.lo.whatever to
// cell_size.whatever

template <typename T>
int CellList<T>::cellIndexOf(const Vector3x<T> &pos,
                             const ProcSubGrid3D<T> &simbox) const {
  int cx = (pos.x - simbox.lo.x) / m_cell_size.x;
  int cy = (pos.y - simbox.lo.y) / m_cell_size.y;
  int cz = (pos.z - simbox.lo.z) / m_cell_size.z;

  cx = std::clamp(cx, 0, m_ncell.x - 1);
  cy = std::clamp(cy, 0, m_ncell.y - 1);
  cz = std::clamp(cz, 0, m_ncell.z - 1);

  return flatten(cx, cy, cz);
}

template <typename T>
std::vector<int> CellList<T>::neighboringCells(int cell_index) const {
  std::vector<int> neighbors;

  int k = cell_index / (m_ncell.x * m_ncell.y);
  int j = (cell_index / m_ncell.x) % m_ncell.y;
  int i = cell_index % m_ncell.x;

  for (int dk = -1; dk <= 1; ++dk)
    for (int dj = -1; dj <= 1; ++dj)
      for (int di = -1; di <= 1; ++di) {
        int ni = i + di;
        int nj = j + dj;
        int nk = k + dk;

        if (ni < 0 || ni >= m_ncell.x)
          continue;
        if (nj < 0 || nj >= m_ncell.y)
          continue;
        if (nk < 0 || nk >= m_ncell.z)
          continue;

        neighbors.push_back(flatten(ni, nj, nk));
      }

  return neighbors;
}
