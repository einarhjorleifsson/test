# rb_whack_clean.R — label whacky (implausible) vessel positions.

#' Label whacky positions: vmask, leg-speed spike, forward scan
#'
#' The recommended filter for vessel positions. It labels, it does not remove: every
#' ping is returned, so the track can be rebuilt with or without the flagged ones and the
#' reason for each flag is kept.
#'
#' Two stages, because neither is enough on its own:
#' 1. [rb_whack_sda()] - judges each ping from **both sides** (rms speed to four neighbours,
#'    then both legs at the ping), so it can flag a bad ping that the forward scan would
#'    adopt as its anchor, e.g. the first ping after a long gap.
#' 2. [rb_whack_forward()] on what stage 1 left - compares with the last **kept** ping, so it
#'    resolves runs and interleaved streams that stage 1 cannot judge from the neighbours.
#'
#' Measured on 543.6 M Icelandic AIS pings (2007-2026): 1.45 M flagged (0.27 %), against 2.87 M with the
#' former 1-microsecond time step. The difference is almost all near-duplicate reports from two
#' merged feeds (96 % between feeds, 92 % within 50 m of the ping before); the kept track is
#' 0.003 % shorter and keeps 706 legs over 1 km at over 25 kn, against 716.
#'
#' The time-step floor `min_dt_s = 10` was chosen on June 2016 and June 2023 (7.2 M pings):
#' with a floor of up to 10 s no kept leg jumps more than 200 m at over 25 kn; at 20 s, 575 do,
#' and at 60 s, 6,881. Same-second reports come mostly from two feeds merged.
#'
#' @section Data frame and DuckDB:
#' The labelling always runs in compiled code, in R; DuckDB never runs the filter itself.
#' * **data.frame** - the whole frame is sorted by `vid`, `time` and labelled in one
#'   compiled call. A table too big for memory is split by the caller into batches of whole
#'   vessels, which is exact for the reason given below.
#' * **lazy DuckDB table** (`tbl_lazy`) - the algorithm is sequential per vessel (a ping's
#'   fate depends on the pings kept before it), which SQL window functions cannot express,
#'   so the table is processed in **batches of whole vessels** (a vessel's result depends on
#'   no other vessel, so batching is exact). Each batch is collected with only
#'   `vid`, `time`, `lon`, `lat`, labelled in compiled code, and only the **flagged** rows
#'   (well under 1 percent) are written back to a temporary table on the same connection. The
#'   result is a lazy table: `x` left-joined to those labels, so the next verbs, a filter
#'   or a `write_dataset()`, run in DuckDB and the full table is never held in R. Row order
#'   of the lazy result is not guaranteed.
#'
#' @param x A data frame or a lazy DuckDB table with `vid`, `lon`, `lat` and `time`.
#' @param kn_max Speed threshold in knots (default 25).
#' @param max_gap_h Passed to [rb_whack_forward()] (default 4).
#' @param min_dt_s Shortest time step for a speed, in seconds, used by both stages
#'   (default 10); see [rb_whack_forward()].
#' @param batch_pings Lazy tables only: approximate pings per batch of whole vessels
#'   (default 5 million, which keeps R below a few GB).
#' @param ... Passed to [rb_whack_sda()].
#'
#' @return `x` with `whack` (logical) and `whack_stage` (`NA`, `"vmask"`, `"spike"` or
#'   `"forward"`). A data frame comes back sorted by `vid` and `time`; a lazy table comes back lazy.
#'
#' @examples
#' # whacky_sim: seven simulated tracks, one scenario each; `truth` marks the injected errors
#' d <- rb_whack_clean(whacky_sim)
#' table(d$scenario, d$whack_stage, useNA = "ifany")
#' track <- d[!d$whack, ]
#'
#' # DuckDB / parquet: the pings pass through R one batch of whole vessels at a time;
#' # the result is a lazy table
#' if (requireNamespace("duckdb", quietly = TRUE) && requireNamespace("dbplyr", quietly = TRUE)) {
#'   con <- DBI::dbConnect(duckdb::duckdb())
#'   print(dplyr::copy_to(con, whacky_sim, "pings") |> rb_whack_clean() |> dplyr::count(whack_stage))
#'   DBI::dbDisconnect(con, shutdown = TRUE)
#' }
#' @export
rb_whack_clean <- function(x, kn_max = 25, max_gap_h = 4, min_dt_s = 10, batch_pings = 5e6, ...) {
  if (inherits(x, "tbl_lazy")) return(.whack_clean_lazy(x, kn_max, max_gap_h, min_dt_s, batch_pings, ...))
  .whack_clean_df(x, kn_max, max_gap_h, min_dt_s, ...)
}

.whack_clean_df <- function(x, kn_max, max_gap_h, min_dt_s, ...) {
  grp_vars <- dplyr::group_vars(x)
  out <- rb_whack_sda(x, kn_max = kn_max, min_dt_s = min_dt_s, ...)
  s1 <- !is.na(out$whack_sda)
  fwd <- rb_whack_forward(out[!s1, ], kn_max = kn_max, max_gap_h = max_gap_h, min_dt_s = min_dt_s)
  # rb_whack_forward() re-sorts by vid, time; out[!s1, ] already is, so the order is the same
  stage <- out$whack_sda
  stage[!s1] <- ifelse(fwd$whack2, "forward", NA_character_)
  out$whack_sda <- NULL
  out$whack_stage <- stage
  out$whack <- !is.na(stage)
  if (length(grp_vars) > 0) out <- dplyr::group_by(out, dplyr::across(dplyr::all_of(grp_vars)))
  out
}

.whack_clean_lazy <- function(x, kn_max, max_gap_h, min_dt_s, batch_pings, ...) {
  con <- dbplyr::remote_con(x)
  # A deterministic row number per vessel, computed in DuckDB and identically in the batch
  # (the batch is sorted by time, lon, lat before it is labelled), is the join key back.
  keyed <- x |>
    dplyr::group_by(vid) |>
    dbplyr::window_order(time, lon, lat) |>
    dplyr::mutate(.whack_rn = dplyr::row_number()) |>
    dplyr::ungroup()
  vn <- x |> dplyr::count(vid) |> dplyr::collect() |> dplyr::arrange(dplyr::desc(n))
  batches <- split(vn$vid, cumsum(as.numeric(vn$n)) %/% batch_pings)
  tmp <- paste0("whack_labels_", format(as.hexmode(sample.int(.Machine$integer.max, 1)), width = 8))
  first <- TRUE
  for (g in batches) {
    # filter to the batch BEFORE numbering, so DuckDB sorts only this batch, not the table
    d <- x |>
      dplyr::filter(vid %in% g) |>
      dplyr::select(vid, time, lon, lat) |>
      dplyr::group_by(vid) |>
      dbplyr::window_order(time, lon, lat) |>
      dplyr::mutate(.whack_rn = dplyr::row_number()) |>
      dplyr::ungroup() |>
      dplyr::collect() |>
      dplyr::arrange(vid, time, lon, lat)
    lab <- .whack_clean_df(d, kn_max, max_gap_h, min_dt_s, ...)
    lab <- lab[lab$whack, c("vid", ".whack_rn", "whack_stage")]
    lab$.whack_rn <- as.numeric(lab$.whack_rn)
    if (first) {
      DBI::dbWriteTable(con, tmp, lab, temporary = TRUE, overwrite = TRUE)
      first <- FALSE
    } else if (nrow(lab) > 0) {
      DBI::dbAppendTable(con, tmp, lab)
    }
  }
  if (first) stop("rb_whack_clean(): the table has no vessels.", call. = FALSE)
  # window_order() clears the numbering's order: left on the result it is added to a later
  # arrange(), which fails after a count() or summarise() ("time" must appear in GROUP BY)
  keyed |>
    dplyr::left_join(dplyr::tbl(con, tmp), by = c("vid", ".whack_rn")) |>
    dplyr::mutate(whack = !is.na(whack_stage)) |>
    dplyr::select(-.whack_rn) |>
    dbplyr::window_order()
}


#' Speed-distance-angle filter for whacky positions (compiled)
#'
#' Flags implausible positions with the speed-distance-angle algorithm of
#' `argosfilter::sdafilter()` (Freitas et al. 2008), re-implemented as one compiled pass
#' over all vessels. Nothing is removed: every row is returned, with a label. Stage 1 of
#' [rb_whack_clean()].
#'
#' @section Algorithm:
#' Per vessel, on the track sorted by time:
#' 1. **vmask** - the root-mean-square speed from each ping to its two neighbours on either
#'    side. Local maxima of that speed above `kn_max` are flagged and taken out of the
#'    working track; this repeats until no remaining ping exceeds `kn_max`.
#' 2. A vmask flag stands only if the ping is more than `vmask_min_dist` metres from its
#'    predecessor (default 0: every flag stands; `argosfilter` uses 5000).
#' 3. **Spikes** - a ping is flagged when the two legs at it enclose an angle of at most
#'    `ang[k]` degrees and either both legs are longer than `distlim[k]` metres or, where
#'    `speedlim_kn[k] > 0`, both legs are faster than that. Repeats on the thinned track
#'    until none is left. The speed form is the time-aware one: the same leg is long or
#'    short depending on the time it took.
#'
#' Every speed uses a time step of at least `min_dt_s`; `argosfilter` adds 1 s to every
#' step instead.
#'
#' The package defaults tuned for Argos seal tracks (15/25 degrees, 2.5/5 km, 5 km) leave
#' tens of thousands of speed spikes in vessel AIS. The defaults here were tuned on two
#' years of Icelandic AIS (2012 and 2016, 32 variants): a speed test, `vmask_min_dist = 0`.
#'
#' @param x A data frame with `vid`, `lon`, `lat` and `time` (`POSIXct`).
#' @param kn_max Speed threshold in knots (default 25).
#' @param ang Angle limits in degrees (default 25).
#' @param distlim Distance limits in metres, one per `ang`; used where `speedlim_kn` is 0.
#' @param speedlim_kn Speed limits in knots, one per `ang` (default `kn_max`); 0 = use `distlim`.
#' @param vmask_min_dist See above (default 0).
#' @param min_dt_s Shortest time step, in seconds, used for a speed (default 10). A shorter
#'   step counts as `min_dt_s`: over a few seconds the error in position and the whole-second
#'   timestamps dominate the speed. Within `min_dt_s` a ping is therefore flagged only if it
#'   is more than `kn_max * min_dt_s` away (about 130 m at the defaults).
#'
#' @return `x`, sorted by `vid`, `time` (ties by `lon`, `lat`), with `whack_sda`: `NA` = not flagged,
#'   `"vmask"` or `"spike"` = the step that flagged it.
#'
#' @references Freitas, C., Lydersen, C., Fedak, M.A. and Kovacs, K.M. (2008). A simple
#'   new algorithm to filter marine mammal Argos locations. Marine Mammal Science 24:315-325.
#'
#' @seealso [rb_whack_clean()], the full recipe; [rb_whack_forward()], its second stage.
#'
#' @examples
#' d <- rb_whack_sda(whacky_sim)
#' table(d$scenario, d$whack_sda, useNA = "ifany")
#'
#' # argosfilter's own settings
#' d <- rb_whack_sda(whacky_sim, ang = c(15, 25), distlim = c(2500, 5000),
#'                   speedlim_kn = c(0, 0), vmask_min_dist = 5000)
#' @export
rb_whack_sda <- function(x, kn_max = 25, ang = 25, distlim = rep(0, length(ang)),
                      speedlim_kn = rep(kn_max, length(ang)), vmask_min_dist = 0,
                      min_dt_s = 10) {
  grp_vars <- dplyr::group_vars(x)
  # lon and lat break ties in time, so the result does not depend on the row order of the input
  # (parquet reads come back in a different order from run to run)
  x <- dplyr::arrange(dplyr::ungroup(x), vid, time, lon, lat)
  code <- match(x$vid, unique(x$vid))
  r <- whack_sda_cpp(x$lat, x$lon, as.numeric(x$time), code, kn_max * 0.514444,
                     as.numeric(ang), as.numeric(distlim), as.numeric(speedlim_kn) * 0.514444,
                     vmask_min_dist, min_dt_s)
  x$whack_sda <- c(NA_character_, "vmask", "spike")[r + 1L]
  if (length(grp_vars) > 0) x <- dplyr::group_by(x, dplyr::across(dplyr::all_of(grp_vars)))
  x
}


#' Forward-only speed filter for whacky positions (compiled)
#'
#' Stage 2 of [rb_whack_clean()], run on the pings [rb_whack_sda()] left. Where a
#' neighbour-based test misses consecutive clusters (because the bad points form an
#' internally slow but geographically wrong track), `rb_whack_forward` catches them by always
#' comparing against the **last retained** ping rather than the immediate predecessor.
#'
#' @section Algorithm:
#' Iterates forward through the track (sorted by `vid`, `time`). For each ping,
#' the implied speed from the last *retained* ping is computed. If that speed
#' exceeds `kn_max` **and** the time gap is below `max_gap_h`, the ping is
#' flagged and skipped — the reference ping does not advance. The reference
#' resets unconditionally when the time gap exceeds `max_gap_h` (legitimate
#' long gaps, e.g. vessel off-air overnight, should not trigger the filter).
#'
#' A time step shorter than `min_dt_s` counts as `min_dt_s`.
#'
#' On its own it is aggressive: a bad ping that it adopts as the reference, the first
#' ping of a track or the first after a gap, makes it flag the good pings that follow.
#' That is why [rb_whack_clean()] runs [rb_whack_sda()] first.
#'
#' @param x A data frame (or tibble). Must contain columns `vid`, `lon`, `lat`,
#'   `time` (`POSIXct`).
#' @param kn_max Speed threshold in knots (default: 25).
#' @param max_gap_h Time gaps larger than this (hours) reset the reference
#'   without flagging (default 4). Prevents legitimate port-to-port transits
#'   from being erroneously removed.
#' @param min_dt_s Shortest time step, in seconds, used for a speed (default 10). A shorter
#'   step counts as `min_dt_s`: over a few seconds the error in position and the whole-second
#'   timestamps dominate the speed. Within `min_dt_s` a ping is therefore flagged only if it
#'   is more than `kn_max * min_dt_s` away (about 130 m at the defaults).
#'
#' @return The input data frame, sorted by `vid`, `time`, with an added logical
#'   column `whack2`: `TRUE` = flagged.
#'
#' @seealso [rb_whack_clean()], the full recipe.
#'
#' @examples
#' # on its own, the forward scan adopts a bad first ping after a gap as its anchor:
#' # it flags 3 pings where rb_whack_clean() flags only the bad one
#' d <- rb_whack_forward(whacky_sim[whacky_sim$vid == 5, ])
#' sum(d$whack2)
#' @export
rb_whack_forward <- function(x, kn_max = 25, max_gap_h = 4, min_dt_s = 10) {
  ms_max      <- kn_max * 0.514444
  max_gap_sec <- max_gap_h * 3600

  # Compiled core (src/whack_forward.cpp) — ~1000x the R for-loop this
  # replaced (2026-09-23), same algorithm, verified byte-identical against it
  # on whacks1 and stress cases (isolated spikes, runs, a >4h gap, an injected
  # NA). `time` is passed as seconds-since-epoch, which is POSIXct's own
  # underlying numeric representation, so no unit conversion is needed here.
  # Speeds use a time step of at least min_dt_s (2026-10-08; was 1e-6 s), in the
  # port and in the test's R reference alike.
  .fwd <- function(lon, lat, time) {
    whack_forward_cpp(lon, lat, as.numeric(time), ms_max, max_gap_sec, min_dt_s)
  }

  grp_vars <- dplyr::group_vars(x)
  out <- x |>
    dplyr::ungroup() |>
    dplyr::arrange(vid, time) |>
    dplyr::group_by(vid) |>
    dplyr::mutate(whack2 = .fwd(lon, lat, time)) |>
    dplyr::ungroup()
  if (length(grp_vars) > 0)
    out <- dplyr::group_by(out, dplyr::across(dplyr::all_of(grp_vars)))
  out
}
