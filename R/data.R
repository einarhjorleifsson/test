#' Seven simulated vessel tracks with known errors
#'
#' One scenario per vessel, made to show what [rb_whack_clean()] does and where it fails.
#' Pings every 5 minutes at 64 N with 10 m of GPS noise; the tracks are drawn in kilometres,
#' so the distances are exact. Built by `data-raw/whacky_sim.R` in the package source.
#'
#' 1. **spike**: one ping moved 13 km off a 10-knot track.
#' 2. **burst aside**: four consecutive pings moved 25 km to the side of the track.
#' 3. **burst ahead**: four consecutive pings moved 25 km along the vessel's heading.
#' 4. **burst near**: four consecutive pings moved 9 km to the side.
#' 5. **bad ping after a gap**: the first ping after 5 hours of silence moved 18 km.
#' 6. **real move after a gap**: after 6 hours of silence the vessel is 47 km away; nothing is wrong.
#' 7. **fast vessel**: fishing at 4 knots, steaming at 30 knots for an hour, fishing again; nothing is wrong.
#'
#' @format A data frame with 412 rows:
#' \describe{
#'   \item{vid}{vessel, 1-7, one per scenario}
#'   \item{scenario}{the scenario's name}
#'   \item{time}{`POSIXct`, UTC}
#'   \item{lon, lat}{decimal degrees}
#'   \item{truth}{`TRUE` for a ping moved off its track on purpose: an error the filter should find}
#' }
"whacky_sim"
