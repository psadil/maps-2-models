

co |>
  ggplot(aes(x=gold, y=cope)) +
  geom_point(alpha = .02) +
  geom_abline(intercept=0, slope=1) +
  geom_smooth()

lm(cope~gold, data=co)
