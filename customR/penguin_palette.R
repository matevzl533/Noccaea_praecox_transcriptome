penguins <- palmerpenguins::penguins
penguin_corp_color <- function(...) {

  penguin_corp_colors <- c(
    `imperial_red`     = "#F94244",
    `orange1` = "#F37430",
    `carrot_orange`     = "#F8961D",
    `orange2`    = "#F98544",
    `saffron`   = "#FAC84F",
    `pistachio`     = "#8FBE6D",
    `jungle_green` = "#45A989",
    `dark_cyan`     = "#508F8A",
    `glaucous`    = "#5B7494",
    `cerulean`   = "#287CA1")

  cols <- c(...)

  if (is.null(cols))
    return (penguin_corp_colors)

  penguin_corp_colors[cols]
}
palette_gen <- function(palette = "main", direction = 1) {
  function(n) {
    if (n > length(penguin_corp_palette(palette)))
      warning("Not enough colors in this palette!")
    else {
      all_colors <- penguin_corp_palette(palette)
      all_colors <- unname(unlist(all_colors))
      all_colors <- if (direction >= 0) all_colors else rev(all_colors)
      color_list <- all_colors[1:n]
    }
  }
}
scale_colour_penguin <- function(palette = "main", direction = 1, ...) {
  ggplot2::discrete_scale(
    "colour", "penguin",
    palette_gen(palette, direction),
    ...
  )
}
penguin_corp_palette <- function(palette = "main", ...) {
  penguin_corp_palettes <- list(
    `main` = penguin_corp_color("imperial_red","orange1","carrot_orange","orange2","saffron","pistachio","jungle_green","dark_cyan","glaucous","cerulean"),
    `reverse-main` = penguin_corp_color("cerulean", "glaucous", "dark_cyan", "jungle_green", "pistachio", "saffron", "orange2", "carrot_orange", "orange1", "imperial_red"),
    `fiver` = penguin_corp_color("orange2","saffron","pistachio","jungle_green","dark_cyan"),
    `middle-main` = penguin_corp_color("pistachio","jungle_green","dark_cyan","glaucous","cerulean","imperial_red","orange1","carrot_orange","orange2","saffron")
    
  )
  penguin_corp_palettes[[palette]]
}
scale_color_penguin <- scale_colour_penguin
scale_fill_penguin <- function(palette = "main", direction = 1, ...) {

  ggplot2::discrete_scale(
    "fill", "penguin",
    palette_gen(palette, direction),
    ...
  )
}