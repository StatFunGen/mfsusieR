# Package load hook + namespace plumbing.
#
# - `@useDynLib` registers the cpp11-compiled symbols (per design.md
#   D14, Phase 3 hot-path kernels in `src/posterior_mixture.cpp`).
# - `.onLoad` will later host S3 method registrations on susieR's
#   internal generics (mirroring `mvsusieR/R/zzz.R`).
#
# Numerical routines do not live here.

#' @useDynLib mfsusieR, .registration = TRUE
#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  invisible(NULL)
}
