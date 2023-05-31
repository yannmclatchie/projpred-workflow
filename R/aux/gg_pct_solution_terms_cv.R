# Source:
# <https://github.com/fweber144/modelselection/blob/9e65616c6bf9e234747b2e97a3f602cb5c264c92/bodyfat.Rmd>,
# which in turn is based on
# <https://github.com/avehtari/modelselection/blob/ad0e10ce7d3a47de560609db7abe83777eed28bc/bodyfat.Rmd>.

gg_pct_solution_terms_cv <- function(cvvs) {
  pct_cum <- apply(cvvs[["pct_solution_terms_cv"]][, -1], 2, cumsum)
  pct_cum <- cbind(cvvs[["pct_solution_terms_cv"]][, 1, drop = FALSE], pct_cum)
  rows <- nrow(pct_cum)
  col <- nrow(pct_cum)
  pctch <- round(pct_cum, 2)
  colnames(pctch)[1] <- ".size"
  pct <- get_pct_arr(pctch, 13)
  col_brks <- get_col_brks()
  pct$val_grp <- sapply(pct$val, function(x) sum(x >= col_brks$breaks))
  if (identical(rows, 0)) rows <- pct$var[1]
  pct$sel <- (pct$.size == col) & (pct$var %in% rows)
  brks <- sort(unique(as.numeric(pct$val_grp)))
  pct$val_grp <- factor(pct$val_grp,levels=brks)
  ggobj <- ggplot(pct, aes(x = .size, y = var)) +
    geom_tile(aes(fill = val_grp, color = sel),
              width = 1, height = 1, linewidth = 1) +
    geom_text(aes(label = val, fontface = sel+1), size = 3) +
    coord_cartesian(expand = FALSE) +
    scale_y_discrete(limits = rev(levels(pct$var))) +
    scale_x_discrete(limits = factor(seq(1,col))) +
    scale_color_manual(values = c("white", "black")) +
    labs(x = "Model size", y = "") +
    scale_fill_manual(breaks = brks, values = col_brks$pal[brks+1]) +
    theme(legend.position = "none")
  return(ggobj)
}
