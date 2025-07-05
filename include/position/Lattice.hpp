#pragma once
#include <Utility.hpp>
#include <array>
#include <limits>
#include <mpi.h>
#include <type_traits>
#include <vector>
namespace Position {
int countInFCC(const Vector3i &unit_sizes) {
  int x = unit_sizes.x, y = unit_sizes.y, z = unit_sizes.z;
  return

      (x + 1) * (y + 1) * (z + 1) +

      (x + 1) * (y) * (z) +

      (x) * (y + 1) * (z) +

      (x) * (y) * (z + 1);
}
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
std::vector<Vector3x<T>> genFCCLattice(const Vector3i &size,
                                       const Vector3x<T> &side_lengths) {

  int count_in_lattice = countInFCC(size);
  std::vector<Vector3x<T>> __result(count_in_lattice);
  int t = 0;
  // Generate cubic lattice points
  for (unsigned int _z = 0; _z <= size.z; _z++) {
    for (unsigned int _y = 0; _y <= size.y; _y++) {
      for (unsigned int _x = 0; _x <= size.x; _x++) {
        __result.at(t) = _x * side_lengths.x;
        __result.at(t) = _y * side_lengths.y;
        __result.at(t++) = _z * side_lengths.y;
      }
    }
  }

  // Add face-centered points on yz planes (front and back)
  for (unsigned int _z = 0; _z <= size.z; _z++) {
    for (unsigned int _y = 0; _y < size.y; _y++) {
      for (unsigned int _x = 0; _x < size.x; _x++) {
        /*printf("Loop[1], t = %d\n", t);*/
        __result.at(t) = (_x + 0.5) * side_lengths.x;
        __result.at(t) = _y * side_lengths.y + side_lengths.y / 2;
        __result.at(t++) = _z * side_lengths.z;
      }
    }
  }

  for (double _y = 0; _y <= size.y; _y++) {
    for (double _z = 0; _z < size.z; _z++) {
      for (double _x = 0; _x < size.x; ++_x) {
        /*printf("Loop[2], t = %d\n", t);*/
        __result.at(t) = _x * side_lengths.x + side_lengths.x / 2;
        __result.at(t) = _y * side_lengths.y;
        __result.at(t++) = _z * side_lengths.z + side_lengths.z / 2;
        // Add face-centered points on xz planes (left and right)
      }
    }
  }

  // Add face-centered points on xy planes (bottom and top)
  //

  for (double _x = 0; _x <= size.x; _x++) {
    for (double _z = 0; _z < size.z; _z++) {
      for (double _y = 0; _y < size.y; ++_y) {
        /*printf("Loop[3], t = %d\n", t);*/
        __result.at(t) = _x * side_lengths.x;
        __result.at(t) = _y * side_lengths.y + side_lengths.y / 2;
        __result.at(t++) = _z * side_lengths.z + side_lengths.z / 2;
      }
    }
  }
}
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
std::pair<std::vector<Vector3x<T>>, Vector3x<T>>
generateFCCCompleted(const Vector3x<T> &leng, const Vector3i &sizes) {

  std::array<Vector3x<T>, 4> basis{
      Vector3x<T>{0., 0., 0.},
      Vector3x<T>{0.5, 0.5, 0.},
      Vector3x<T>{0.5, 0., 0.5},
      Vector3x<T>{0., 0.5, 0.5},
  };
  std::vector<Vector3x<T>> result;
  for (unsigned int i = 0; i != sizes.x; ++i) {
    for (unsigned int j = 0; j != sizes.y; ++j) {
      for (unsigned int k = 0; k != sizes.z; ++k) {
        for (auto &atom : basis) {
          result.emplace_back(Vector3x<T>{
              (atom.x + i) * leng.x,
              (atom.y + j) * leng.y,
              (atom.z + k) * leng.z,
          });
        }
      }
    }
  }
  return std::make_pair(result, Vector3x<T>{leng.x * sizes.x + 1e-08,
                                            leng.y * sizes.y + 1e-08,
                                            leng.z * sizes.z + 1e-08});
}
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
std::vector<Vector3x<T>>
generateComplexFCCStructure(const std::vector<Vector3i> &sizes,
                            const std::vector<Vector3x<T>> &side_length,
                            const std::vector<double> &distances) {
  int count = sizes.size();

  if (sizes.size() != side_length.size()) {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  std::vector<std::pair<double, double>> xz_lattice_shift(count);

  for (int t_ = 0; t_ != count; ++t_) {
    xz_lattice_shift[t_].first = sizes[t_].x * side_length[t_].x / 2;
    xz_lattice_shift[t_].second = sizes[t_].z * side_length[t_].y / 2;
  }
#if 0
  std::pair<double, double> min_lat_size;
  min_lat_size.first =
      std::min_element(xz_lattice_shift.begin(), xz_lattice_shift.end(),
                       [](auto t, auto v) { return t.first < t.first; })
          ->first;
  min_lat_size.second =
      std::min_element(xz_lattice_shift.begin(), xz_lattice_shift.end(),
                       [](auto t, auto v) { return t.second < t.second; })
          ->second;
  for (int t_ = 0; t_ != count; ++t_) {
    xz_lattice_shift[t_].first =
        (xz_lattice_shift[t_].first - min_lat_size.first) / 2;
    xz_lattice_shift[t_].second =
        (xz_lattice_shift[t_].second - min_lat_size.second) / 2;
  }
#endif
  int total_atoms_count = 0;
  for (int t_ = 0; t_ != count; ++t_) {
    total_atoms_count += countInFCC(sizes[t_]);
  }

  std::vector<Vector3x<T>> res(total_atoms_count);
  // centralize acros y-axis

  double current_start = 0;
  double new_start, new_end;
  double total_lattice_ylen_half = sizes[0].y * side_length[0].y;
  for (int t_ = 0; t_ != sizes.size(); ++t_) {
    new_start = current_start;
    new_end = new_start + sizes[t_].y * side_length[t_].y;
    current_start = new_end + (t_ < sizes.size() - 2 ? distances[t_] : 0);
  }
  total_lattice_ylen_half = new_end / 2;

  current_start = -total_lattice_ylen_half;
  int counter = 0;
  for (int t_ = 0; t_ != count; ++t_) {
    std::vector<Vector3x<T>> lattice =
        genFCCLattice(sizes[t_], side_length[t_]);

    double new_start = current_start;
    double new_end = new_start + sizes[t_].y * side_length[t_].y;
    move_cubic(lattice, -xz_lattice_shift[t_].first, new_start,
               -xz_lattice_shift[t_].second);
    // for (int t_ = 0; t_ < 10; ++t_) {
    //   std::cout << "Lattice[" << t_ << "] = (" << lattice[t_].x << ","
    //             << lattice[t_].y << "," << lattice[t_].z << ")" << std::endl;
    // }
    current_start = new_end + (t_ < sizes.size() - 1 ? distances[t_] : 0);
    for (int t_ = 0; t_ != lattice.size(); ++t_) {
      res[counter++] = lattice[t_];
    }
  }

  return res;
}

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
std::vector<int> getIndexesOfBorder(bool positive, int idx, T distance,
                                    const std::vector<Vector3x<T>> &positions) {
  std::vector<int> result;

  // Validate input index
  if (idx < 0 || idx > 2) {
    return result; // Return empty vector for invalid index
  }

  // Early return for empty input
  if (positions.empty()) {
    return result;
  }

  // Find the border value (max for positive, min for negative)
  const auto compare = [positive](const T &a, const T &b) {
    return positive ? (a < b) : (a > b);
  };

  T border_value = positions[0][idx]; // Initialize with first element
  for (const auto &atom : positions) {
    if (compare(border_value, atom[idx])) {
      border_value = atom[idx];
    }
  }

  distance += distance / 10;
  const T eps =
      std::max(std::numeric_limits<T>::epsilon(),
               std::abs(border_value) * std::numeric_limits<T>::epsilon());

  for (int i = 0; i < positions.size(); ++i) {
    if (std::abs(border_value - positions[i][idx]) <= eps + distance) {
      result.push_back(i);
    }
  }

  return result;
}
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
std::pair<Vector3x<T>, Vector3x<T>>
getBoundingRect(const std::vector<Vector3x<T>> &positions, T margin = 0.35196) {
  Vector3x<T> min_val = {std::numeric_limits<T>::max(),
                         std::numeric_limits<T>::max(),
                         std::numeric_limits<T>::max()},
              max_val = {std::numeric_limits<T>::min(),
                         std::numeric_limits<T>::min(),
                         std::numeric_limits<T>::min()};
  for (const auto &particle : positions) {
    min_val.x = std::min(min_val.x, particle.x);
    min_val.y = std::min(min_val.y, particle.y);
    min_val.z = std::min(min_val.z, particle.z);
    max_val.x = std::max(max_val.x, particle.x);
    max_val.y = std::max(max_val.y, particle.y);
    max_val.z = std::max(max_val.z, particle.z);
  }
  min_val -= Vector3x<T>{margin / 2., margin / 2., margin / 2.};
  max_val += Vector3x<T>{margin / 2., margin / 2., margin / 2.};
  return std::make_pair(min_val, max_val);
}

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
void movePositions(std::vector<Vector3x<T>> &holder, double x_padding,
                   double y_padding, double z_padding = 0.) {
  for (size_t t_ = 0; t_ != holder.size(); ++t_) {
    holder[t_].x += x_padding;
    holder[t_].y += y_padding;
    holder[t_].z += z_padding;
  }
}

}; // namespace Position