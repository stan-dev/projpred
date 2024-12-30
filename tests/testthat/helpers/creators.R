imperfect_alternation <- function(unq_vals, n_tail = 6L, ...) {
  x <- rep_len(unq_vals, ...)
  tail_idxs_x <- tail(seq_along(x), n = n_tail)
  x[tail_idxs_x] <- x[c(tail_idxs_x[-1],
                        tail_idxs_x[1])]
  return(x)
}
