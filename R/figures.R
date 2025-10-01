bootstrap_r <- function(.data, times = 2000) {
  boots <- boot::boot(
    data = .data$r,
    statistic = function(x, i) mean(x[i]),
    R = 2000
  )

  ci <- boot::boot.ci(boots, type = "perc")

  tibble::tibble(
    rr = boots$t0,
    lower = ci$percent[4],
    upper = ci$percent[5]
  )
}
