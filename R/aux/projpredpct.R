# Source: <https://github.com/avehtari/modelselection/blob/ad0e10ce7d3a47de560609db7abe83777eed28bc/projpredpct.R>

library(dplyr)
library(tidyr)
library(RColorBrewer)

get_pct_arr <- function(pctch, nv) {
    as_tibble(pctch[1:nv, 1:(nv + 1)]) %>%
        gather("var", "val", colnames(pctch)[1:nv+1], factor_key = T)
}

get_col_brks <- function() {
  list(breaks = seq(5e-3, 1-5e-3, length.out = 7),
       pal = brewer.pal(11, "RdBu")[3:10])
}
