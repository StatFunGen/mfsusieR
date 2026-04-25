# Package load hook.
#
# Phase 3 PR group 1 (skeleton): .onLoad is intentionally empty. Subsequent
# PR groups will populate this with:
#   - S3 method registrations on susieR's internal generics
#     (mirroring mvsusieR/R/zzz.R), dispatched on `mf_individual` and
#     `mfsusie` classes,
#   - cached bindings for susieR / fsusieR internal helpers used
#     directly by mfsusieR (per design.md D8 "function-level imports").
#
# Keep this file the canonical place for namespace plumbing; do not put
# numerical routines here.

.onLoad <- function(libname, pkgname) {
  invisible(NULL)
}
