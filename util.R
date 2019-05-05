str_partition <- function(s,size) {
  lgt <- str_length(s)
  starts <- seq(1,lgt,size)
  str_sub(s,starts,starts+size-1)
}

str_collapse <- function(s,collapse) s %>%
  str_c(collapse=collapse)

str_cat <- function(num,width) num %>%
  as.character %>%
  str_partition(width) %>%
  str_collapse("\n") %>%
  cat
