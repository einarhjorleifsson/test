// whack_sda.cpp — compiled core of whack_sda(): a speed-distance-angle track filter.
//
// Derived from the algorithm of argosfilter::sdafilter() (Freitas et al. 2008, Marine
// Ecology Progress Series 367:31-42). Re-implemented, not copied: the R package is a
// set of row-by-row for-loops (about 120 microseconds per ping), this is one compiled
// pass per vessel. Checked against the R function on real tracks: identical removals.
//
// Steps, per vessel (input sorted by group, then time):
//   1. vmask: the rms speed of each ping to its 2 neighbours on either side; local maxima
//      above vmax are removed; repeat on the thinned track until no survivor exceeds vmax.
//   2. a vmask removal stands only if the ping is more than vmask_min_dist from its
//      predecessor in the ORIGINAL track (0 = every removal stands).
//   3. angle/distance spikes: remove a ping whose two legs enclose an angle <= ang[k] and
//      are both longer than distlim[k] metres, or, where speedlim[k] > 0, both faster than
//      speedlim[k] m/s (the time-aware form); repeat on the thinned track until none.
//
// Differences from the R package, all deliberate: haversine rather than acos (which gives
// NaN or 0 under ~1 m); bearings from atan2; the vmask loop stops if a pass removes nothing
// (the R loop can spin on a plateau); no Argos location-class argument; a time step shorter
// than min_dt counts as min_dt (the R package adds 1 s to every step; 2026-10-08).
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

namespace {

const double R_EARTH = 6371000.0;
inline double rad(double d) { return d * M_PI / 180.0; }

// Haversine between two pings from precomputed radians and cos(lat): 2 sin + 1 asin per pair
// instead of 2 sin + 2 cos + 1 asin, and no degree-to-radian conversion in the inner loops.
struct Rad {
  std::vector<double> phi, lam, cphi;
  Rad(const double* lat, const double* lon, int a, int b) : phi(b - a), lam(b - a), cphi(b - a) {
    for (int i = a; i < b; ++i) {
      phi[i - a] = rad(lat[i]); lam[i - a] = rad(lon[i]); cphi[i - a] = std::cos(phi[i - a]);
    }
  }
  inline double dist(int i, int j) const {            // i, j are indices relative to the track start
    double sp = std::sin((phi[j] - phi[i]) / 2), sl = std::sin((lam[j] - lam[i]) / 2);
    double a = sp * sp + cphi[i] * cphi[j] * sl * sl;
    return 2 * R_EARTH * std::asin(std::min(1.0, std::sqrt(a)));
  }
};

inline double bearing_deg(double lat1, double lon1, double lat2, double lon2) {
  double p1 = rad(lat1), p2 = rad(lat2), dl = rad(lon2 - lon1);
  double y = std::sin(dl) * std::cos(p2);
  double x = std::cos(p1) * std::sin(p2) - std::sin(p1) * std::cos(p2) * std::cos(dl);
  double b = std::atan2(y, x) * 180.0 / M_PI;
  return b < 0 ? b + 360.0 : b;
}

// One track, indices [a, b) into the full vectors. Writes a reason code into `removed`:
// 1 = vmask peak, 2 = angle/distance (or speed) spike. Nothing is dropped here: the caller labels.
void sda_track(const double* lat, const double* lon, const double* t, int a, int b,
               double vmax, const std::vector<double>& ang, const std::vector<double>& distlim,
               const std::vector<double>& speedlim, double vmask_min_dist, double min_dt,
               int* removed) {
  const int N = b - a;
  if (N < 5) return;

  // 1. vmask --------------------------------------------------------------------
  std::vector<int> cur(N);
  for (int i = 0; i < N; ++i) cur[i] = a + i;
  std::vector<char> vm(N, 0);
  std::vector<double> v, s1, s2;
  const Rad rd(lat, lon, a, b);
  while (true) {
    const int n = cur.size();
    if (n < 5) break;
    // Leg speeds once per leg, not once per ping: s1[i] = i to i+1, s2[i] = i to i+2.
    // v[i] is the rms of s2[i-2], s1[i-1], s1[i], s2[i]  (the four legs the R code recomputes).
    s1.assign(n, 0.0); s2.assign(n, 0.0); v.assign(n, 0.0);
    for (int i = 0; i + 1 < n; ++i) {
      int p = cur[i], q = cur[i + 1];
      s1[i] = rd.dist(p - a, q - a) / std::max(std::fabs(t[p] - t[q]), min_dt);
      if (i + 2 < n) { int r2 = cur[i + 2]; s2[i] = rd.dist(p - a, r2 - a) / std::max(std::fabs(t[p] - t[r2]), min_dt); }
    }
    for (int i = 2; i <= n - 3; ++i)
      v[i] = std::sqrt((s2[i - 2] * s2[i - 2] + s1[i - 1] * s1[i - 1] + s1[i] * s1[i] + s2[i] * s2[i]) / 4.0);
    std::vector<int> peaks;
    bool ascending = true; double curr_peak = 0, curr_null = 0;
    for (int i = 2; i <= n - 3; ++i) {
      if (ascending) {
        if (v[i] > curr_peak) curr_peak = v[i];
        else { ascending = false; curr_null = v[i]; peaks.push_back(i - 1); }
      } else {
        if (v[i] < curr_null) curr_null = v[i];
        else { ascending = true; curr_peak = v[i]; }
      }
    }
    if (ascending) peaks.push_back(n - 3);
    std::vector<char> drop(n, 0);
    int ndrop = 0;
    for (int p : peaks) if (v[p] > vmax) { drop[p] = 1; ++ndrop; }
    if (ndrop == 0) break;
    std::vector<int> nxt; nxt.reserve(n - ndrop);
    double maxi = 0;
    for (int i = 0; i < n; ++i) {
      if (drop[i]) vm[cur[i] - a] = 1;
      else { nxt.push_back(cur[i]); maxi = std::max(maxi, v[i]); }
    }
    cur.swap(nxt);
    if (maxi <= vmax) break;
  }

  // 2. a vmask removal stands only > vmask_min_dist from the original predecessor -----
  std::vector<int> keep; keep.reserve(N);
  for (int i = 0; i < N; ++i) {
    bool rem = false;
    if (vm[i] && i > 0)
      rem = rd.dist(i - 1, i) > vmask_min_dist;
    if (rem) removed[a + i] = 1; else keep.push_back(a + i);
  }

  // 3. angle / distance (or speed) spikes, repeated -------------------------------------
  const int K = ang.size();
  while (K > 0) {
    const int n = keep.size();
    if (n < 3) break;
    std::vector<char> drop(n, 0); int ndrop = 0;
    for (int i = 1; i <= n - 2; ++i) {
      int p = keep[i - 1], c = keep[i], q = keep[i + 1];
      double dprev = rd.dist(p - a, c - a);
      double dnext = rd.dist(c - a, q - a);
      double an;
      if ((lat[c] == lat[q] && lon[c] == lon[q]) || (lat[c] == lat[p] && lon[c] == lon[p])) an = 180;
      else {
        an = std::fabs(bearing_deg(lat[c], lon[c], lat[p], lon[p]) -
                       bearing_deg(lat[c], lon[c], lat[q], lon[q]));
        if (an > 180) an = 360 - an;
      }
      double sprev = dprev / std::max(std::fabs(t[c] - t[p]), min_dt);
      double snext = dnext / std::max(std::fabs(t[q] - t[c]), min_dt);
      for (int k = 0; k < K; ++k) {
        bool far = (k < (int)speedlim.size() && speedlim[k] > 0)
                     ? (sprev > speedlim[k] && snext > speedlim[k])
                     : (dprev > distlim[k] && dnext > distlim[k]);
        if (an <= ang[k] && far) { drop[i] = 1; break; }
      }
      if (drop[i]) ++ndrop;
    }
    if (ndrop == 0) break;
    std::vector<int> nxt; nxt.reserve(n - ndrop);
    for (int i = 0; i < n; ++i) { if (drop[i]) removed[keep[i]] = 2; else nxt.push_back(keep[i]); }
    keep.swap(nxt);
  }
}

} // namespace

// All vessels in one call. `grp` is an integer code per ping, sorted so each vessel is a
// contiguous run, and `time` is sorted within it. Returns 0 (fine), 1 (vmask) or 2 (spike).
// [[Rcpp::export]]
IntegerVector whack_sda_cpp(NumericVector lat, NumericVector lon, NumericVector time,
                            IntegerVector grp, double vmax, NumericVector ang,
                            NumericVector distlim, NumericVector speedlim,
                            double vmask_min_dist, double min_dt) {
  const int n = lat.size();
  std::vector<int> removed(n, 0);
  std::vector<double> A(ang.begin(), ang.end()), D(distlim.begin(), distlim.end()),
                      S(speedlim.begin(), speedlim.end());
  int start = 0;
  for (int i = 1; i <= n; ++i) {
    if (i == n || grp[i] != grp[start]) {
      sda_track(&lat[0], &lon[0], &time[0], start, i, vmax, A, D, S, vmask_min_dist, min_dt, removed.data());
      start = i;
    }
  }
  return IntegerVector(removed.begin(), removed.end());
}
