require(tikzDevice)

# define plot saving method
save_tikz_plot <- function(
    plot, filename, width = NA, height = NA, asp = NA
) {
  # automatic scaling
  if (is.na(asp)) asp <- 1.618
  if (is.na(width) && is.na(height)) {
    height <- 3.71
    width <- height * asp
  }
  else if (is.na(width)) {
    width <- height * asp
  }
  else if (is.na(height)) {
    height <- width / asp
  }
  
  # make tex
  tikz(file = filename, width = width, height = height)
  print(plot)
  dev.off()
  
  # patch cropping issues
  lines <- readLines(con = filename)
  lines <- lines[-which(grepl("\\path\\[clip\\]*", lines))]
  lines <- lines[-which(grepl("\\path\\[use as bounding box*", lines))]
  writeLines(lines, con = filename)
}