library(ggplot2)

minimal_lmplot <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point() +
    geom_smooth(method = method, color = "darkred", ...) +
    theme_minimal()
  p
}

minimal_kdeplot <- function(data, mapping, ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_density() +
    theme_minimal()
  p
}
