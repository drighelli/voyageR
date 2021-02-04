muted <- function(x) {
  sink(nullfile())
  on.exit(sink())
  x
}
