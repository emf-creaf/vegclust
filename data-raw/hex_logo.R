library(hexSticker)

logo_image <- fs::path("data-raw", "vegclust_logo_image.png")
sticker(
  logo_image,
  package = "vegclust", p_size = 18, p_y = 1.60, p_color = "#000080",
  s_x = 1, s_y = .9, s_width = .52,
  filename = fs::path("data-raw", "vegclust.png"),
  #   url = "emf.creaf.cat", u_size = 6, u_color = "#BFD77A", u_y = .2, u_x = 1.2,
  h_fill = "#FFFFFF", h_color = "#000080"
)
