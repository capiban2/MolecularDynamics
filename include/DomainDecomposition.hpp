#pragma once
#include "Utility.hpp"
#include <array>
#include <memory>
#include <mpi.h>
#include <numeric>
#include <vector>

template <typename T> class DomainDecomposition {

  MPI_Comm m_comm;
  int m_mpi_size{}, m_mpi_rank{};
  MPI_Comm m_cart_comm{};
  Vector3x<T> m_global_lo, m_global_hi;
  std::array<int, 3> mpi_procs_dims{}, m_periods_{}, m_mpi_cart_coords{};
  std::array<int, 6> m_nbrs_{}; // -X,+X, -Y,+Y, -Z,+Z
  std::vector<Particle<T>> sendBuffer(const std::vector<Particle<T>> &,
                                      const std::array<double, 3> &,
                                      const std::array<double, 3> &) const;
  ProcSubGrid3D<T> m_local_simbox;
  void receiveBuffer(std::vector<Particle<T>> &recvbuf, int neighbor) const;

  int findSplitIndex(const std::vector<int> &weights,
                     int total_weight) const noexcept;
  void rcb3D(const std::vector<std::vector<std::vector<int>>> &grid, int i0,
             int i1, int j0, int j1, int k0, int k1, int first_rank,
             int n_ranks, std::vector<ProcSubGrid3D<T>> &output, int axis = 0);

public:
  DomainDecomposition() = default;

#if 0
  // Use external SimulationBox (e.g., after restart)
  template <typename P>
  void setFromExistingBox(const SimulationBox<P> &externalBox) {
    m_SimBox = externalBox;
  }
#endif

  void initCart(MPI_Comm communicator, const Vector3x<T> &global_lo,
                const Vector3x<T> &global_hi);

  // TODO: i think ghost_wide should be obtained from somewhere else but for now
  // let it be like that
  void exchangeGhosts(std::vector<Particle<T>> &particles, T ghost_wide,
                      int *ghost_index_first) const noexcept;

  const Vector3x<T> &localLo() const { return m_local_simbox.lo; };
  const Vector3x<T> &localHi() const { return m_local_simbox.hi; };
  void rearrangeWithRCB(const std::vector<Particle<T>> &_particles);
  std::array<int, 6>
  findRCBNeighbors(const ProcSubGrid3D<T> &my_procinfo,
                   const std::vector<ProcSubGrid3D<T>> &all_procs) const;
  std::pair<std::array<int, 6>, bool>
  getNewNeighbors(const std::array<int, 6> &prev_nghbrs,
                  const ProcSubGrid3D<T> &proc_info, MPI_Comm comm, int m_size);

  const int *getNeighbors() const { return m_nbrs_.data(); }
  const int *getCoords() const { return m_mpi_cart_coords.data(); }
  MPI_Comm getCartComm() const { return m_cart_comm; }
  int rank() const { return m_mpi_rank; }
  int size() const { return m_mpi_size; }
};

template <typename T>
void DomainDecomposition<T>::initCart(MPI_Comm communicator,
                                      const Vector3x<T> &global_lo,
                                      const Vector3x<T> &global_hi) {
  m_comm = communicator;

  m_periods_[0] = m_periods_[1] = m_periods_[2] = 1;

  MPI_Comm_size(m_comm, &m_mpi_size);
  MPI_Comm_rank(m_comm, &m_mpi_rank);
  MPI_Dims_create(m_mpi_size, 3, mpi_procs_dims.data());
  MPI_Cart_create(m_comm, 3, mpi_procs_dims.data(), m_periods_.data(), 1,
                  &m_cart_comm);
  MPI_Cart_coords(m_cart_comm, m_mpi_rank, 3, m_mpi_cart_coords.data());

  // Determine neighbors
  for (int i = 0; i < 3; ++i) {
    MPI_Cart_shift(m_cart_comm, i, 1, &m_nbrs_[2 * i], &m_nbrs_[2 * i + 1]);
  }
  m_local_simbox = ProcSubGrid3D<T>{
      .lo = {-1, -1, -1}, .hi = {-1, -1, -1}, .rank = m_mpi_rank};
  // Векторы локального поддомена
  Vector3x<T> box_size = global_hi - global_lo;
  Vector3x<T> a_local = box_size.x / mpi_procs_dims[0];
  Vector3x<T> b_local = box_size.y / mpi_procs_dims[1];
  Vector3x<T> c_local = box_size.z / mpi_procs_dims[2];

  Vector3x<T> local_lo = global_lo + m_mpi_cart_coords[0] * a_local +
                         m_mpi_cart_coords[1] * b_local +
                         m_mpi_cart_coords[2] * c_local;
  Vector3x<T> local_hi = local_lo + Vector3x<T>(a_local, b_local, c_local);

  m_local_simbox.lo = local_lo;
  m_local_simbox.hi = local_hi;
}

template <typename T>
int DomainDecomposition<T>::findSplitIndex(const std::vector<int> &weights,
                                           int total_weight) const noexcept {
  int cumulative = 0;
  for (int i = 0; i < weights.size(); ++i) {
    cumulative += weights[i];
    if (cumulative >= total_weight / 2)
      return i;
  }
  return static_cast<int>(weights.size()) - 1;
}
template <typename T>
void DomainDecomposition<T>::rcb3D(
    const std::vector<std::vector<std::vector<int>>> &grid, int i0, int i1,
    int j0, int j1, int k0, int k1, int first_rank, int n_ranks,
    std::vector<ProcSubGrid3D<T>> &output, int axis) {
  if (n_ranks == 1) {
    output[first_rank] = ProcSubGrid3D<T>{
        .lo = {i0, j0, k0}, .hi = {i1, j1, k1}, .rank = first_rank};
    return;
  }

  int next_axis = (axis + 1) % 3;
  std::vector<int> weights;

  if (axis == 0) {
    weights.resize(i1 - i0 + 1, 0);
    for (int i = i0; i <= i1; ++i)
      for (int j = j0; j <= j1; ++j)
        for (int k = k0; k <= k1; ++k)
          weights[i - i0] += grid[i][j][k];
  } else if (axis == 1) {
    weights.resize(j1 - j0 + 1, 0);
    for (int j = j0; j <= j1; ++j)
      for (int i = i0; i <= i1; ++i)
        for (int k = k0; k <= k1; ++k)
          weights[j - j0] += grid[i][j][k];
  } else {
    weights.resize(k1 - k0 + 1, 0);
    for (int k = k0; k <= k1; ++k)
      for (int i = i0; i <= i1; ++i)
        for (int j = j0; j <= j1; ++j)
          weights[k - k0] += grid[i][j][k];
  }

  int total = std::accumulate(weights.begin(), weights.end(), 0);

  // Handle empty region (no particles): just split spatially in the middle
  int cut;
  if (total == 0) {
    if (axis == 0)
      cut = (i0 + i1) / 2;
    else if (axis == 1)
      cut = (j0 + j1) / 2;
    else
      cut = (k0 + k1) / 2;
  } else {
    // Find split index such that both sides have nearly equal weight
    int cumulative = 0;
    int split_idx = 0;
    for (; split_idx < weights.size(); ++split_idx) {
      cumulative += weights[split_idx];
      if (cumulative >= total / 2)
        break;
    }
    // Map split_idx to grid index
    if (axis == 0)
      cut = i0 + split_idx;
    else if (axis == 1)
      cut = j0 + split_idx;
    else
      cut = k0 + split_idx;
  }

  int n1 = n_ranks / 2;
  int n2 = n_ranks - n1;

  if (axis == 0) {
    rcb3D(grid, i0, cut, j0, j1, k0, k1, first_rank, n1, output, next_axis);
    rcb3D(grid, cut + 1, i1, j0, j1, k0, k1, first_rank + n1, n2, output,
          next_axis);
  } else if (axis == 1) {
    rcb3D(grid, i0, i1, j0, cut, k0, k1, first_rank, n1, output, next_axis);
    rcb3D(grid, i0, i1, cut + 1, j1, k0, k1, first_rank + n1, n2, output,
          next_axis);
  } else {
    rcb3D(grid, i0, i1, j0, j1, k0, cut, first_rank, n1, output, next_axis);
    rcb3D(grid, i0, i1, j0, j1, cut + 1, k1, first_rank + n1, n2, output,
          next_axis);
  }
}
template <typename T>
void DomainDecomposition<T>::rearrangeWithRCB(
    const std::vector<Particle<T>> &_particles) {
  std::vector<ProcSubGrid3D<T>> partitions(m_mpi_size);
  if (m_mpi_rank == 0) {
    std::array<int, 3> even_dims;
    MPI_Dims_create(m_mpi_size, 3, even_dims.data());
    Vector3x<T> bin_size{m_global_hi.x / even_dims.at(0),
                         m_global_hi.x / even_dims.at(1),
                         m_global_hi.x / even_dims.at(2)};
    int ni = even_dims.at(0), nj = even_dims.at(1), nk = even_dims.at(2);
    int size = ni * nj * nk;

    std::vector<std::vector<std::vector<int>>> global_grid(
        ni, std::vector<std::vector<int>>(nj, std::vector<int>(nk, 0)));
    for (auto const &part : _particles) {
      const auto &pos = part.pos;

      int i = pos.x / bin_size.x;
      int j = pos.y / bin_size.y;
      int k = pos.z / bin_size.z;
      global_grid[i][j][k]++;
    }
    rcb3D(global_grid, 0, ni - 1, 0, nj - 1, 0, nk - 1, size, partitions);
  }
  MPI_Allgather(MPI_IN_PLACE, m_mpi_size * sizeof(ProcSubGrid3D<T>), MPI_BYTE,
                partitions.data(), m_mpi_size * sizeof(ProcSubGrid3D<T>),
                MPI_BYTE, m_cart_comm);
  m_local_simbox = partitions.at(m_mpi_rank);
  m_nbrs_ = findRCBNeighbors(m_local_simbox, partitions);
}

template <typename T>
std::array<int, 6> DomainDecomposition<T>::findRCBNeighbors(
    const ProcSubGrid3D<T> &my_procinfo,
    const std::vector<ProcSubGrid3D<T>> &all_procs) const {
  auto const &my_lo = my_procinfo.lo;
  auto const &my_hi = my_procinfo.hi;
  int m_rank = my_procinfo.rank;
  std::array<int, 6> nbrs = {-1, -1, -1, -1, -1, -1}; // -x,+x,-y,+y,-z,+z

  for (const auto &other : all_procs) {
    if (other.rank == m_rank)
      continue;
    const auto &lo = other.lo;
    const auto &hi = other.hi;
    // shared region in YZ plane, touching in -X
    if (other.hi.x + 1 == my_lo.x &&
        !(other.hi.y < my_lo.y || other.lo.y > my_hi.y) &&
        !(other.hi.z < my_lo.z || lo.z > my_hi.z)) {
      nbrs[0] = other.rank;
    }

    if (lo.x - 1 == my_hi.x && !(hi.y < my_lo.y || lo.y > my_hi.y) &&
        !(hi.z < my_lo.z || lo.z > my_hi.z)) {
      nbrs[1] = other.rank;
    }

    // аналогично ±Y:
    if (hi.y + 1 == my_lo.y && !(hi.x < my_lo.x || lo.x > my_hi.x) &&
        !(hi.z < my_lo.z || lo.z > my_hi.z)) {
      nbrs[2] = other.rank;
    }

    if (lo.y - 1 == my_hi.y && !(hi.x < my_lo.x || lo.x > my_hi.x) &&
        !(hi.z < my_lo.z || lo.z > my_hi.z)) {
      nbrs[3] = other.rank;
    }

    // ±Z:
    if (hi.z + 1 == my_lo.z && !(hi.x < my_lo.x || lo.x > my_hi.x) &&
        !(hi.y < my_lo.y || lo.y > my_hi.y)) {
      nbrs[4] = other.rank;
    }

    if (lo.z - 1 == my_hi.z && !(hi.x < my_lo.x || lo.x > my_hi.x) &&
        !(hi.y < my_lo.y || lo.y > my_hi.y)) {
      nbrs[5] = other.rank;
    }
  }

  return nbrs;
}

// HINT: call it just after making new decomposition from rcb
// to find out is it necessary to exchange all data to from and all over again
// or we can just send/receive data with previous neighbors, like shift borders
// a little
template <typename T>
std::pair<std::array<int, 6>, bool>
DomainDecomposition<T>::getNewNeighbors(const std::array<int, 6> &prev_nghbrs,
                                        const ProcSubGrid3D<T> &proc_info,
                                        MPI_Comm comm, int m_size) {

  int m_rank = proc_info.rank;
  std::vector<ProcSubGrid3D<T>> whole_info(m_size);
  MPI_Allgather(&proc_info, sizeof(ProcSubGrid3D<T>), MPI_BYTE,
                whole_info.data(), sizeof(ProcSubGrid3D<T>) * m_size, MPI_BYTE,
                comm);
  auto new_neigh = findRCBNeighbors(proc_info, whole_info);
  bool is_neighbors_preserved = true;
  for (int t_ = 0; t_ < 6; ++t_)
    is_neighbors_preserved &= prev_nghbrs[t_] == new_neigh[t_];
  bool going_hard_way;
  int local_int = is_neighbors_preserved ? 1 : 0;
  int global_int = 0;
  MPI_Allreduce(&local_int, &global_int, 1, MPI_INT, MPI_LAND, comm);

  going_hard_way = global_int != 1;

  return std::make_pair(new_neigh, going_hard_way);
}
template <typename T>
void DomainDecomposition<T>::exchangeGhosts(
    std::vector<Particle<T>> &particles, T ghost_width,
    int *ghost_index_first) const noexcept {
  // Step 1: Remove existing ghost particles
  particles.erase(particles.begin() + *ghost_index_first, particles.end());

  // Step 2: Prepare communication buffers
  std::array<std::vector<Particle<T>>, 6> send_buffers;
  std::array<std::vector<Particle<T>>, 6> recv_buffers;
  std::array<MPI_Request, 12> requests; // 6 sends + 6 receives
  int req_count = 0;

  // Step 3: Identify particles in ghost regions
  const auto &lo = m_local_simbox.lo;
  const auto &hi = m_local_simbox.hi;

  for (const auto &p : particles) {
    // -X direction (left)
    if (p.pos.x <= lo.x + ghost_width)
      send_buffers[0].push_back(p);
    // +X direction (right)
    if (p.pos.x >= hi.x - ghost_width)
      send_buffers[1].push_back(p);
    // -Y direction (bottom)
    if (p.pos.y <= lo.y + ghost_width)
      send_buffers[2].push_back(p);
    // +Y direction (top)
    if (p.pos.y >= hi.y - ghost_width)
      send_buffers[3].push_back(p);
    // -Z direction (back)
    if (p.pos.z <= lo.z + ghost_width)
      send_buffers[4].push_back(p);
    // +Z direction (front)
    if (p.pos.z >= hi.z - ghost_width)
      send_buffers[5].push_back(p);
  }

  // Step 4: Exchange particle counts
  std::array<int, 6> send_counts, recv_counts;
  for (int dir = 0; dir < 6; dir++) {
    send_counts[dir] = send_buffers[dir].size();
    MPI_Isend(&send_counts[dir], 1, MPI_INT, m_nbrs_[dir], 0, m_cart_comm,
              &requests[req_count++]);
    MPI_Irecv(&recv_counts[dir], 1, MPI_INT, m_nbrs_[dir], 0, m_cart_comm,
              &requests[req_count++]);
  }
  MPI_Waitall(req_count, requests.data(), MPI_STATUSES_IGNORE);
  req_count = 0;

  // Step 5: Exchange particle data
  for (int dir = 0; dir < 6; dir++) {
    if (recv_counts[dir] > 0) {
      recv_buffers[dir].resize(recv_counts[dir]);
      MPI_Irecv(recv_buffers[dir].data(),
                recv_counts[dir] * sizeof(Particle<T>), MPI_BYTE, m_nbrs_[dir],
                1, m_cart_comm, &requests[req_count++]);
    }
    if (send_counts[dir] > 0) {
      MPI_Isend(send_buffers[dir].data(),
                send_counts[dir] * sizeof(Particle<T>), MPI_BYTE, m_nbrs_[dir],
                1, m_cart_comm, &requests[req_count++]);
    }
  }
  MPI_Waitall(req_count, requests.data(), MPI_STATUSES_IGNORE);

  // Step 6: Append received ghost particles
  for (int dir = 0; dir < 6; dir++) {
    particles.insert(particles.end(), recv_buffers[dir].begin(),
                     recv_buffers[dir].end());
  }
  *ghost_index_first = particles.size() - std::accumulate(recv_counts.begin(),
                                                          recv_counts.end(), 0);
}