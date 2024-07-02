.get_gray_copes <- function(avail, at) {
  to_tbl(avail, measure = "cope") |>
    dplyr::left_join(at, by = c("x", "y", "z")) |>
    dplyr::filter(!is.na(.data$label))
}

get_eff_by_roi <- function(test, at) {
  test |>
    tidyr::unnest(avail) |>
    dplyr::mutate(
      data = purrr::map(avail, .get_gray_copes, at = at),
      avail = as.character(avail)
    ) |>
    tidyr::unnest(data)
}

phi <- function(x, y) {
  tp <- mean(x & y)
  fp <- mean(!x & y)
  tn <- mean(!x & !y)
  fn <- mean(x & !y)
  if (((tp & fp) == 0) | ((tp & fn) == 0) | ((fn & tn) == 0) | ((fp & tn) == 0)) {
    return(0)
  }
  cor(x, y)
}

sim.phi <- function(X, ...) {
  n <- nrow(X)
  D <- array(0, dim = c(n, n))
  for (i in 1:(n - 1)) {
    for (j in i:n) {
      D[i, j] <- phi(X[i, ], X[j, ])
    }
  }
  D <- D + t(D)
  return(D)
}

dice <- function(x, y) {
  2 * sum(x & y) / (sum(x) + sum(y))
}

sim.rho <- function(X, ...) {
  n <- nrow(X)
  D <- array(0, dim = c(n, n))
  for (i in 1:(n - 1)) {
    for (j in i:n) {
      D[i, j] <- cor(X[i, ], X[j, ])
    }
  }
  D <- D + t(D)
  return(D)
}
