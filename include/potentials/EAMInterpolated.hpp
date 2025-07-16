#pragma once
#include "VerletList.hpp"
#include <IMDManager.hpp>
#include <algorithm>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cmath>
#include <fstream>
#include <numeric>
#include <potentials/Potential.hpp>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <types/PairInteraction.hpp>
#include <unordered_map>
#include <vector>
template <typename T> T string_to_float(const std::string &str) {
  if constexpr (std::is_same_v<T, float>) {
    return std::stof(str);
  } else if constexpr (std::is_same_v<T, double>) {
    return std::stod(str);
  } else if constexpr (std::is_same_v<T, long double>) {
    return std::stold(str);
  } else {
    static_assert(sizeof(T) == 0, "Unsupported floating-point type");
  }
}

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
class EAMInterpolated final : public Potential<T> {

  void delete2d(T **__p, size_t t) {
    if (__p == nullptr)
      return;

    for (int t_ = 0; t_ != t; ++t_) {
      delete[] __p[t_];
    }
    delete[] __p;
  }

  void delete3d(T ***__p, size_t t, size_t t2) {
    if (__p == nullptr)
      return;

    for (size_t t_ = 0; t_ != t; ++t_)
      delete2d(__p[t_], t2);
    delete[] __p;
  }
  std::vector<std::string> m_interp_files;
  // HINT: to provide neighboring algoritm, max across all types
  T m_max_cutoff = 0, m_min_r = 0;
  T m_max_rho = 0, m_min_rho = 0;
  T m_dr = 0, m_drho = 0;
  int n_type_max = 0, n_types = 0;
  int m_nr = 0, m_nrho = 0;
  std::vector<T> m_current_rho, m_dembed_energy;
  T m_virial = 0;
  int *m_type2map = nullptr;
  int *type2frho = nullptr, **type2rho = nullptr, **type2z2r = nullptr;
  T **frho = nullptr, **rho = nullptr, **z2r = nullptr;

  T ***rho_spline = nullptr, ***frho_spline = nullptr, ***z2r_spline = nullptr;

  void file2array() {

    T r_min, r_max, rho_max, rho_min, dr_rho, frho_max, frho_min, drho, z2r_min,
        z2r_max, dr_z2r, dr;
    std::unordered_map<int, bool> seen;

    // TODO: implement unification of all grids; now expect only grids of equal
    // range; otherwise will work incorrectly
    T r, p, cof1, cof2, cof3, cof4;

    for (const auto &filename : m_interp_files) {
      std::ifstream f;
      f.open(filename, std::ios::in);
      if (!f.is_open()) {
        throw std::runtime_error("Couldnt open " + filename + " file");
      }
      // HINT: skip header and comments

      for (int i = 0; i < 3; ++i)
        f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      std::string row;
      std::vector<std::string> data;
      boost::split(std::getline(f, row), data, boost::is_any_of(" \t"));
      int count = std::stoi(data.at(0));
      std::array<int, 2> types;
      for (int t_ = 0; t_ < count; ++t_) {
        boost::split(std::getline(f, row), data, boost::is_any_of(" \t"));
        types[t_] = std::stoi(data.at(0));
        if (seen[types[t_]])
          continue;
        seen[types[t_]] = true;
        r_min = string_to_float<T>(data.at(3));

        // HINT: this usually do nothing for to this point all grids should be
        // equal
        m_max_cutoff = std::max(string_to_float<T>(data.at(4)), m_max_cutoff);
        dr = string_to_float<T>(data.at(5));
        rho_min = string_to_float<T>(data.at(6));
        drho = string_to_float<T>(data.at(8));

        // HINT: the same as for max_cutoff
        m_max_rho = std::max(m_max_rho, string_to_float<T>(data.at(7)));

        // HINT: there should be unifying of the grid, but now already unifyied
        // grids are expected
        if (m_min_r > 0 && fabs(m_min_r - r_min) > 1e-02)
          throw std::runtime_error("m_r_min mismatch!! Fix your grids");
        else
          m_min_r = r_min;
        if (m_dr > 0 && fabs(m_dr - dr) > 1e-02)
          throw std::runtime_error("dr mismatch!! Fix your grids");
        else
          m_dr = dr;
        if (m_min_rho > 0 && fabs(m_min_rho - rho_min) > 1e-02)
          throw std::runtime_error("min_rho mismatch!! Fix your grids");
        else
          m_min_rho = rho_min;
        if (m_drho > 0 && fabs(m_drho - drho) > 1e-02)
          throw std::runtime_error("drho mismatch!! Fix your grids");
        else
          m_drho = drho;
      }
    }
    m_nrho = std::floor((m_max_rho - m_min_rho) / m_drho + 0.5);
    m_nr = std::floor((m_max_cutoff - m_min_r) / m_dr + 0.5);

    for (int i = 1; i <= n_type_max; ++i)
      type2frho[i] = m_type2map[i];
    for (int i = 1; i <= n_type_max; ++i)
      for (int j = 1; j <= n_type_max; ++j)
        type2rho[i][j] = m_type2map[j];

    int irow, icol, n = 0;
    for (int i = 1; i <= n_type_max; i++) {
      for (int j = 1; j <= n_type_max; j++) {
        irow = m_type2map[i];
        icol = m_type2map[j];
        if (irow == -1 || icol == -1) {
          type2z2r[i][j] = 0;
          continue;
        }
        if (irow < icol) {
          irow = m_type2map[j];
          icol = m_type2map[i];
        }
        n = 0;
        for (int m = 0; m < irow; m++)
          n += m + 1;
        n += icol;
        type2z2r[i][j] = n;
      }
    }
  }
  void interpolate(int n, T delta, T *f, T **spline) {
    for (int m = 1; m <= n; m++)
      spline[m][6] = f[m];

    spline[1][5] = spline[2][6] - spline[1][6];
    spline[2][5] = 0.5 * (spline[3][6] - spline[1][6]);
    spline[n - 1][5] = 0.5 * (spline[n][6] - spline[n - 2][6]);
    spline[n][5] = spline[n][6] - spline[n - 1][6];

    for (int m = 3; m <= n - 2; m++)
      spline[m][5] = ((spline[m - 2][6] - spline[m + 2][6]) +
                      8.0 * (spline[m + 1][6] - spline[m - 1][6])) /
                     12.0;

    for (int m = 1; m <= n - 1; m++) {
      spline[m][4] = 3.0 * (spline[m + 1][6] - spline[m][6]) -
                     2.0 * spline[m][5] - spline[m + 1][5];
      spline[m][3] = spline[m][5] + spline[m + 1][5] -
                     2.0 * (spline[m + 1][6] - spline[m][6]);
    }

    spline[n][4] = 0.0;
    spline[n][3] = 0.0;

    for (int m = 1; m <= n; m++) {
      spline[m][2] = spline[m][5] / delta;
      spline[m][1] = 2.0 * spline[m][4] / delta;
      spline[m][0] = 3.0 * spline[m][3] / delta;
    }
  }
  void array2spline() {

    rho_spline = new T **[n_types];
    frho_spline = new T **[n_types];
    z2r_spline = new T **[n_types * (n_types + 1) / 2];
    for (int t_ = 0; t_ != n_types; ++t_) {
      rho_spline[t_] = new T *[n_types];
      for (int j_ = 0; j_ != n_types; ++j_)
        rho_spline[t_][j_] = new T[7];
      frho_spline[t_] = new T *[n_types];
      for (int j_ = 0; j_ != n_types; ++j_)
        frho_spline[t_][j_] = new T[7];
      z2r_spline[t_] = new T *[n_types];
      for (int j_ = 0; j_ != n_types; ++j_)
        z2r_spline[t_][j_] = new T[7];
    }
    std::unordered_map<int, bool> seen;
    for (const auto &filename : m_interp_files) {
      std::ifstream f;
      f.open(filename, std::ios::in);
      if (!f.is_open()) {
        throw std::runtime_error("Couldnt open " + filename + " file");
      }
      // HINT: skip header and comments

      for (int i = 0; i < 3; ++i)
        f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      std::string row;
      std::vector<std::string> data;
      boost::split(std::getline(f, row), data, boost::is_any_of(" \t"));
      int count = std::stoi(data.at(0));
      std::array<int, 2> types;
      for (int t_ = 0; t_ < count; ++t_) {
        boost::split(std::getline(f, row), data, boost::is_any_of(" \t"));
        int type = std::stoi(data.at(0));
        types[t_] = type;
        if (seen[type])
          continue;
        seen[type] = true;
      }
      // HINT: if it is hybrid file, then skip alloy's constants, 'case already
      // handled previously
      if (count > 1)
        f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      std::vector<T> func;
      for (int t_ = 0; t_ < count; ++t_) {
        if (seen[types[t_]])
          continue;
        int current_type_mapped = m_type2map[types[t_]];
        T delta = m_dr;
        func.resize(m_nr);
        for (int t_ = 0; t_ != m_nr; ++t_) {
          boost::split(std::getline(f, row), data, boost::is_any_of(" \t"));
          func[t_] = string_to_float<T>(data.at(0));
        }
        interpolate(m_nr, delta, func.data(),
                    rho_spline[type2rho[types[t_]][types[t_]]]);

        delta = m_drho;
        func.resize(m_nrho);
        for (int t_ = 0; t_ != m_nrho; ++t_) {
          boost::split(std::getline(f, row), data, boost::is_any_of(" \t"));
          func[t_] = string_to_float<T>(data.at(0));
        }
        interpolate(m_nrho, delta, func.data(),
                    frho_spline[type2frho[types[t_]]]);

        delta = m_dr;

        func.resize(m_nr);
        for (int t_ = 0; t_ != m_nr; ++t_) {
          boost::split(std::getline(f, row), data, boost::is_any_of(" \t"));
          func[t_] = string_to_float<T>(data.at(0));
        }
        interpolate(m_nr, delta, func.data(),
                    z2r[type2z2r[types[t_]][types[t_]]]);
      }
      // HINT: pure case, no need in handling hybrid pair
      if (count <= 1)
        continue;
      int ftype_mapped = m_type2map[types[0]],
          stype_mapped = m_type2map[types[1]];
      T delta = m_dr;
      int nz2r = m_nr;
      func.resize(nz2r);
      for (int t_ = 0; t_ != nz2r; ++t_) {
        boost::split(std::getline(f, row), data, boost::is_any_of(" \t"));
        func[t_] = string_to_float<T>(data.at(0));
      }
      interpolate(nz2r, delta, func.data(), z2r[type2z2r[types[0]][types[1]]]);
    }
  }

  void allocate() {

    type2frho = new int[n_type_max + 1];
    type2rho = new int *[n_type_max + 1];
    type2z2r = new int *[n_type_max + 1];
    for (int t_ = 1; t_ <= n_type_max; ++t_) {
      type2rho[t_] = new int[n_type_max + 1];
      type2z2r[t_] = new int[n_type_max + 1];
    }
  }

public:
  EAMInterpolated(const std::vector<std::string> &intepolation_files)
      : m_interp_files(intepolation_files) {}

  // HINT: for dumping into logs and then using it like the leverage to
  // recreate computation
  std::string getDescription() const noexcept override;

  virtual void configure() override {

    for (const auto &filename : m_interp_files) {
      std::ifstream f;
      f.open(filename, std::ios::in);
      if (!f.is_open()) {
        throw std::runtime_error("Couldnt open " + filename + " file");
      }
      // HINT: skip header and comments

      std::set<int> types;
      for (int i = 0; i < 3; ++i)
        f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      std::string row;
      std::vector<std::string> data;
      boost::split(std::getline(f, row), data, boost::is_any_of(" \t"));
      int count = std::stoi(data.at(0));

      for (int t_ = 0; t_ < count; ++t_) {
        boost::split(std::getline(f, row), data, boost::is_any_of(" \t"));
        int element_id = std::stoi(data.at(0));
        types.insert(element_id);
        n_type_max = std::max(n_type_max, element_id);
      }
      n_types = types.size();
      std::vector<int> vec_types(types.begin(), types.end());
      m_type2map = new int[n_type_max + 1];
      for (int t_ = 0; t_ != vec_types.size(); ++t_)
        m_type2map[vec_types[t_]] = t_;
    }

    allocate();
    file2array();
  }
  ~EAMInterpolated() {

    if (type2frho)
      delete[] type2frho;
    if (type2rho) {
      for (int t_ = 1; t_ <= n_type_max; ++t_) {
        if (type2rho[t_])
          delete[] type2rho[t_];
      }
      delete[] type2rho;
    }
    if (type2z2r) {
      for (int t_ = 1; t_ <= n_type_max; ++t_) {
        if (type2z2r[t_])
          delete[] type2z2r[t_];
      }
      delete[] type2z2r;
    }
  }

  virtual void
  computeBulkForces(const std::vector<MD::PairInteraction<T>> &_interactions,
                    std::vector<T> &_f_x, std::vector<T> &_f_y,
                    std::vector<T> &_f_z,
                    const std::vector<int> types) override {

    T *coeff;
    T x, y, z, delx, dely, delz, dist, p, rhotp, rhojp, z2p, z2, psip, fpair;
    // TODO: define it later; this value can be vastly greater than actual atoms
    // , that interact via EAM, but for simplicity of indexing will use global
    // first_ghost
    int first_ghost;
    m_current_rho.resize(first_ghost);
    m_dembed_energy.resize(first_ghost);
    std::fill(m_current_rho.begin(), m_current_rho.end(), 0.);
    // calculate densities
#ifdef OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int t_ = 0; t_ != _interactions.size(); ++t_) {
      auto &interact = _interactions[t_];

      if (interact.r >= m_max_cutoff + 1e-08)
        continue;

      int f_type = types[interact.i], s_type = types[interact.i],
          i_ = interact.i, j_ = interact.j;
      int m = static_cast<int>((dist - m_min_r) / m_dr);
      m = std::clamp(m, 0, static_cast<int>((m_max_cutoff - m_min_r) / 2) - 2);
      p = std::clamp(((dist - m_dr * m) / m_dr), 0., 1.);
      coeff = rho_spline[type2rho[f_type][s_type]][m];
      T res = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

#ifdef OPENMP_ENABLED
#pragma omp atomic
#endif
      m_current_rho[i_] += res;
      if (j_ >= first_ghost)
        continue;

      if (f_type != s_type) {
        coeff = rho_spline[type2rho[s_type][f_type]][m];
#ifdef OPENMP_ENABLED
#pragma omp atomic
#endif
        rho[j_] += ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
      } else
#ifdef OPENMP_ENABLED
#pragma omp atomic
#endif
        rho[j_] += res;
    }
// HINT: this thing better be enabled
#ifdef OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int t_ = 0; t_ < first_ghost; ++t_) {
      T density = m_current_rho[t_];
      if (density > m_max_rho || density < m_min_rho) {
        printf("Handle it somehow if energy is calculating");
      }
      int type = types[t_];
      // HINT: skip particles that is not interact via EAM
      if (m_type2map[type] == -1)
        continue;

      int m = static_cast<int>((density - m_min_rho) / m_drho);
      m = std::clamp(m, 0, static_cast<int>((m_max_rho - m_min_rho) / 2) - 2);
      p = std::clamp(((density - m_drho * m) / m_drho), 0., 1.);

      coeff = frho_spline[type2frho[type]][m];

      m_dembed_energy[t_] =
          ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
    }
#ifdef OPENMP_ENABLED
#pragma omp parallel for
#endif
    for (int t_ = 0; t_ < _interactions.size(); ++t_) {
      auto &interact = _interactions[t_];

      if (interact.r >= m_max_cutoff + 1e-08)
        continue;

      int f_type = types[interact.i], s_type = types[interact.i],
          i_ = interact.i, j_ = interact.j;

      if (i_ >= first_ghost || j_ >= first_ghost)
        continue;

      // TODO: may be use here typ2rho and then m_drho will be 1d
      int m = static_cast<int>((dist - m_min_r) / m_dr);
      m = std::clamp(m, 0, static_cast<int>((m_max_cutoff - m_min_r) / 2) - 2);
      p = std::clamp(((dist - m_dr * m) / m_dr), 0., 1.);

      coeff = rho_spline[type2rho[f_type][s_type]][m];
      rhotp = (coeff[0] * p + coeff[1]) * p + coeff[2];
      coeff = rho_spline[type2rho[s_type][f_type]][m];
      rhojp = (coeff[0] * p + coeff[1]) * p + coeff[2];

      coeff = z2r_spline[type2z2r[f_type][s_type]][m];
      // HINT: this is phi'
      z2p = (coeff[0] * p + coeff[1]) * p + coeff[2];
      // HINT: this is phi
      z2 = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

      psip = m_dembed_energy[t_] * rhojp + m_dembed_energy[j_] * rhotp;

      // HINT: here can be some scalling
      fpair = -psip;
#ifdef OPNEMP_ENABLED
#pragma omp criticial {
#endif
      _f_x[t_] += delx * fpair;
      _f_y[t_] += dely * fpair;
      _f_z[t_] += delz * fpair;
      _f_x[j_] -= delx * fpair;
      _f_y[j_] -= dely * fpair;
      _f_z[j_] -= delz * fpair;
#ifdef OPNEMP_ENABLED
    }
#endif
  }
}

virtual void
computeBulkForces(const VerletList<T> &neigh_list,
                  ParticleData<T> &_data) override {

  T *fx = _data.force_x, *fy = _data.force_y, *fz = _data.force_z;
  const std::vector<int> &types = _data.m_type_id;
  T *coeff;
  T x, y, z, delx, dely, delz, dist, p, rhotp, rhojp, z2p, z2, psip, fpair;
  // TODO: define it later
  int first_ghost;

  int total = _data.force_x.size();
  m_current_rho.resize(first_ghost);
  m_dembed_energy.resize(first_ghost);
  std::fill(m_current_rho.begin(), m_current_rho.end(), 0.0);
  // HINT: calculate densities

#ifdef OPENMP_ENABLED
#pragma omp parallel for
#endif
  for (int t_ = 0; t_ != total; ++t_) {
    int f_type = types[t_];
    x = _data.pos_x[t_];
    y = _data.pos_y[t_];
    z = _data.pos_z[t_];

    for (int j_ : neigh_list.neighbours(t_)) {
      int s_type = types[j_];
      delx = x - _data.pos_x[j_];
      dely = y - _data.pos_y[j_];
      delz = z - _data.pos_z[j_];

      // TODO: boys from lammps for here with squared, and r_max squared too,
      // but i dunno
      dist = sqrt(delx * delx + dely * dely + delz * delz);
      // TODO: think about 1e-08
      if (dist >= m_max_cutoff + 1e-08)
        continue;

      int m = static_cast<int>((dist - m_min_r) / m_dr);
      m = std::clamp(m, 0, static_cast<int>((m_max_cutoff - m_min_r) / 2) - 2);
      p = std::clamp(((dist - m_dr * m) / m_dr), 0., 1.);
      coeff = rho_spline[type2rho[f_type][s_type]][m];
      T res = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

      m_current_rho[t_] += res;
      if (j_ >= first_ghost)
        continue;

      if (f_type != s_type) {
        coeff = rho_spline[type2rho[s_type][f_type]][m];
        rho[j_] += ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
      } else
        rho[j_] += res;
    }
  }

// HINT: this thing better be enabled
#ifdef OPENMP_ENABLED
#pragma omp parallel for
#endif
  for (int t_ = 0; t_ < first_ghost; ++t_) {
    T density = m_current_rho[t_];
    if (density > m_max_rho || density < m_min_rho) {
      printf("Handle it somehow if energy is calculating");
    }
    int type = types[t_];

    // TODO: may be use here typ2rho and then m_drho will be 1d
    int m = static_cast<int>((density - m_min_rho) / m_drho);
    m = std::clamp(m, 0, static_cast<int>((m_max_rho - m_min_rho) / 2) - 2);
    p = std::clamp(((density - m_drho * m) / m_drho), 0., 1.);
    coeff = frho_spline[type2frho[type]][m];
    m_dembed_energy[t_] =
        ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
  }
#ifdef OPENMP_ENABLED
#pragma omp parallel for
#endif
  for (int t_ = 0; t_ < first_ghost; ++t_) {
    int f_type = types[t_];
    x = _data.pos_x[t_];
    y = _data.pos_y[t_];
    z = _data.pos_z[t_];

    for (int j_ : neigh_list.neighbours(t_)) {
      if (j_ >= first_ghost)
        continue;
      int s_type = types[j_];
      delx = x - _data.pos_x[j_];
      dely = y - _data.pos_y[j_];
      delz = z - _data.pos_z[j_];
      dist = sqrt(delx * delx + dely * dely + delz * delz);
      // TODO: think about 1e-08
      if (dist >= m_max_cutoff + 1e-08)
        continue;

      // TODO: may be use here typ2rho and then m_drho will be 1d
      int m = static_cast<int>((dist - m_min_r) / m_dr);
      m = std::clamp(m, 0, static_cast<int>((m_max_cutoff - m_min_r) / 2) - 2);
      p = std::clamp(((dist - m_dr * m) / m_dr), 0., 1.);

      coeff = rho_spline[type2rho[f_type][s_type]][m];
      rhotp = (coeff[0] * p + coeff[1]) * p + coeff[2];
      coeff = rho_spline[type2rho[s_type][f_type]][m];
      rhojp = (coeff[0] * p + coeff[1]) * p + coeff[2];

      coeff = z2r_spline[type2z2r[f_type][s_type]][m];
      // HINT: this is phi'
      z2p = (coeff[0] * p + coeff[1]) * p + coeff[2];
      // HINT: this is phi
      z2 = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

      psip = m_dembed_energy[t_] * rhojp + m_dembed_energy[j_] * rhotp;

      // HINT: here can be some scalling
      fpair = -psip;
#ifdef OPENMP_ENABLED
#pragma openmp critical {
#endif
      fx[t_] += delx * fpair;
      fy[t_] += dely * fpair;
      fz[t_] += delz * fpair;
      fx[j_] -= delx * fpair;
      fy[j_] -= dely * fpair;
      fz[j_] -= delz * fpair;
#ifdef OPENMP_ENABLED
    }
#endif
  }
}
}
}
;