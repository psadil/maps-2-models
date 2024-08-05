correct_d <- function(N) {
  # https://doi.org/10.1101/865881
  # Han Bossier1∗, Thomas E. Nichols2 & Beatrijs Moerkerke1
  exp((lgamma((N - 1) / 2)) - log(sqrt((N - 1) / 2)) - lgamma((N - 2) / 2))
}

d_var <- function(N, d) {
  h <- correct_d(N)
  ((N - 1) * (1 + N * d^2) / (N * (N - 3)) - d^2 / h^2) * h^2
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


z_to_d <- function(z, n, paired = FALSE) {
  # Will be 1 if TRUE, and 2 if FALSE
  paired <- 2 - paired
  paired * z / sqrt(n)
}

z_to_g <- function(z, n, paired = FALSE) {
  z_to_d(z=z, n=n, paired = paired) * correct_d(n)
}


mask_gray <- function(d, mask = MNITemplate::getMNISegPath(res = "2mm")) {
  m <- to_tbl0(RNifti::readNifti(mask)) |>
    dplyr::filter(value == 2)
  dplyr::semi_join(d, m, by = c("x", "y", "z"))
}

mask <- function(d, mask = MNITemplate::getMNIPath("Brain_Mask", "2mm")) {
  m <- to_tbl0(RNifti::readNifti(mask)) |>
    dplyr::filter(value > 0)
  dplyr::semi_join(d, m, by = c("x", "y", "z"))
}

format_arrow_table <- function() {
  targets::tar_format(
    read = function(path) {
      arrow::read_parquet(path, as_data_frame = FALSE)
    },
    write = function(object, path) {
      arrow::write_parquet(object, path, version = "2.6")
    },
    marshal = function(object) as.data.frame(object),
    unmarshal = function(object) arrow::Table$create(object)
  )
}


to_tbl0 <- function(value, measure = "value") {
  dimnames(value) <- list(
    "x" = seq_len(dim(value)[1]),
    "y" = seq_len(dim(value)[2]),
    "z" = seq_len(dim(value)[3])
  )
  cubelyr::as.tbl_cube(value, met_name = measure) |>
    tibble::as_tibble()
}


to_tbl <- function(file, measure = "value", volumes = NULL) {
  value <- RNifti::readNifti(file, volumes = volumes)
  to_tbl0(value, measure = measure)
}
