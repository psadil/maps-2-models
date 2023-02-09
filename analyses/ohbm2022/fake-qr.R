
# fake QR code
p <- crossing(x=seq(0,1,length.out=100), y=seq(0,1,length.out=100)) |> 
  mutate(z=sample(c(0,1),n(),TRUE)) |>
  ggplot(aes(x=x,y=y)) +
  geom_raster(aes(fill=z), show.legend = FALSE) +
  scale_fill_distiller(
    type = "seq",
    direction = -1,
    palette = "Greys") +
  coord_fixed() +
  theme_void()
  
ggsave("analyses/ohbm2022/qr.svg", height = 3, width = 3)
