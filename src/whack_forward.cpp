// whack_forward.cpp — compiled core of whack_forward()'s forward scan.
//
// The R version of rb_whack_forward() (.fwd(), pre-2026-09-23) is a plain R for-loop:
// correct, but ~1000x slower than this port on realistic track sizes, because
// the recurrence is inherently sequential (the reference point is the last
// ACCEPTED ping, not a fixed lag) and so cannot be expressed as a vectorised
// R op or a fixed-offset SQL window function the way whack_fwdbwd() is.
//
// This is a line-for-line port, including one behaviour that looks like a bug
// but is reproduced deliberately for now: if a single NA lon/lat poisons
// `prev`, every later point compared against it gets an NA distance, and
// `ISNAN(d)` (matching the R side's `is.na(d)`) forces flag = TRUE rather than
// advancing `prev` — so the poisoning cascades until max_gap_h forces a reset.
// Verified byte-identical against the R implementation on whacks1 and on
// synthetic stress cases (isolated spikes, runs of 2/3/5, a >4h gap, and an
// injected NA). Flagged for the maintainer to decide on separately; not
// changed here.
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector whack_forward_cpp(NumericVector lon, NumericVector lat,
                                 NumericVector time,
                                 double ms_max, double max_gap_sec, double min_dt) {
  int n = lon.size();
  LogicalVector flag(n, false);
  if (n < 2) return flag;

  const double r = 6371000.0;
  int prev = 0;

  for (int i = 1; i < n; i++) {
    if (NumericVector::is_na(lon[i]) || NumericVector::is_na(lat[i]) ||
        NumericVector::is_na(time[i])) {
      prev = i;
      continue;
    }

    double dt = time[i] - time[prev];
    if (ISNAN(dt) || dt > max_gap_sec) {
      prev = i;
      continue;
    }

    double phi1 = lat[prev] * M_PI / 180.0;
    double phi2 = lat[i]    * M_PI / 180.0;
    double dphi = (lat[i] - lat[prev]) * M_PI / 180.0;
    double dlam = (lon[i] - lon[prev]) * M_PI / 180.0;

    double s = std::sqrt(std::sin(dphi / 2) * std::sin(dphi / 2) +
                          std::cos(phi1) * std::cos(phi2) *
                          std::sin(dlam / 2) * std::sin(dlam / 2));
    if (!ISNAN(s) && s > 1.0) s = 1.0;
    double d = 2 * r * std::asin(s);

    // A step shorter than min_dt counts as min_dt (was 1e-6 s until 2026-10-08, which
    // flagged a duplicate report in the same second, 0.4 m away).
    double dt2 = dt;
    if (dt2 < min_dt) dt2 = min_dt;

    if (ISNAN(d) || d / dt2 > ms_max) {
      flag[i] = true;   // bad — do NOT advance prev
    } else {
      prev = i;          // good — advance reference
    }
  }

  return flag;
}
