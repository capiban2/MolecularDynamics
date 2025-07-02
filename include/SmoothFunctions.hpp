#include <cmath>

// HINT: c1 smooth
inline void taper_cosine(double dist, double r_on, double r_cut, double &S,
                         double &dSdr) {
  const double dr = r_cut - r_on;
  if (dist < r_on) {
    dSdr = 0.0;
    S = 1.0;
  }
// TODO: i dunno, i guess doesnt come to here if that satisfyies
#if 0
  if (r >= r_c) {
    dSdr = 0.0;
    return 0.0;
  }
#endif
  double xi = (dist - r_on) / dr;
  double theta = M_PI * xi;
  double c = std::cos(theta);
  double s = std::sin(theta);
  // Switch
  S = 0.5 * (1.0 + c);
  // Derivative: dS/dr = -½ * sin(θ) * dθ/dr
  // dθ/dr = π / dr
  dSdr = -0.5 * (M_PI / dr) * s;
}

// HINT: c2 smooth
inline void taper_quintic(double dist, double r_on, double r_cut, double &S,
                          double &dSdr) {
  const double dr = r_cut - r_on;
  if (dist < r_on) {
    dSdr = 0.0;
    S = 1.0;
  }
#if 0
    if (r >= r_c) {
        dSdr = 0.0;
        return 0.0;
    }
#endif
  double xi = (dist - r_on) / dr;
  double xi2 = xi * xi;
  double xi3 = xi2 * xi;
  double xi4 = xi3 * xi;
  double xi5 = xi4 * xi;
  // Switch
  S = 1.0 - 10.0 * xi3 + 15.0 * xi4 - 6.0 * xi5;
  // Derivative: dS/dξ = -30ξ² + 60ξ³ - 30ξ⁴
  // dS/dr = (dS/dξ) * (dξ/dr) = (dS/dξ) / dr
  double dSdxi = -30.0 * xi2 + 60.0 * xi3 - 30.0 * xi4;
  dSdr = dSdxi / dr;
}
