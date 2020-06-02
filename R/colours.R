#colours taken from Michael Bentacour (https://github.com/betanalpha/knitr_case_studies/blob/master/gaussian_processes/gp_part2/gp_utility.R) 
cols_custom <- list(
  light = c("#DCBCBC"),
  light_highlight = c("#C79999"),
  mid = c("#B97C7C"),
  mid_highlight = c("#A25050"),
  dark = c("#8F2727"),
  dark_highlight = c("#7C0000"),
  
  light_trans = c("#DCBCBC80"),
  light_highlight_trans = c("#C7999980"),
  mid_trans = c("#B97C7C80"),
  mid_highlight_trans = c("#A2505080"),
  dark_trans = c("#8F272780"),
  dark_highlight_trans = c("#7C000080"),
  
  light_teal = c("#6B8E8E"),
  mid_teal = c("#487575"),
  dark_teal = c("#1D4F4F")
)

show_cols_custom <- function(x=cols_custom){
  if(!is.list(x)) stop("x must be a list")
  nms <- names(x)
  val <- unlist(x)
  n_col <- length(val)
  
  image(1, seq_len(n_col), matrix(seq_len(n_col), nrow = 1), col = val, axes=FALSE, ylab="", xlab = "")
  text(rep(1, n_col), seq_len(n_col), nms)
}
