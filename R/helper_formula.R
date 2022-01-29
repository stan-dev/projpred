#' Break up matrix terms
#'
#' Sometimes there can be terms in a formula that refer to a matrix instead of a
#' single predictor. This function breaks up the matrix term into individual
#' predictors to handle separately, as that is probably the intention of the
#' user.
#'
#' @param formula A [`formula`] for a valid model.
#' @param data The original `data.frame` with a matrix as predictor.
#'
#' @return A `list` containing the expanded [`formula`] and the expanded
#'   `data.frame`.
#'
#' @export
break_up_matrix_term <- function(formula, data) {
  formula <- expand_formula(formula, data)
  tt <- terms(formula)
  response <- attr(tt, "response")
  ## when converting the variables to a list the first element is
  ## "list" itself, so we remove it
  variables_list <- as.list(attr(tt, "variables")[-1])

  ## if there is a response, take it out from the variables as
  ## it is located at the first position after "list"
  if (response) {
    variables_list <- variables_list[-1]
  }

  term_labels <- attr(tt, "term.labels")

  mm <- model.matrix(formula, data)
  assign <- attr(mm, "assign")

  new_data <- data
  for (assignee in unique(assign)) {
    if (assignee == 0) { ## intercept
      next
    }
    appearances <- assign[assign == assignee]
    if (length(appearances) > 1) {
      ## check if special term
      current <- term_labels[assignee]
      int <- grepl(":", current)
      mulilevel <- grepl("\\|", current)
      special <- grepl("[a-z]+\\(([a-z]+)\\)", current)
      individual <- !mulilevel & !int

      linear <- individual & !special
      linear_int <- int & !special

      if (linear) {
        ## if linear we can split it
        split_term <- split_linear_term(current, data)

        formula <- update(formula, paste0(
          ". ~ . - ", current, " + ",
          paste(split_term, collapse = " + ")
        ))
        split_matrix <- mm[, assign == assignee]

        new_data_tmp <- as.data.frame(new_data[, colnames(new_data) != current])
        colnames(new_data_tmp) <-
          colnames(new_data)[colnames(new_data) != current]
        new_data <- cbind(new_data_tmp, split_matrix)
      }

      if (linear_int) {
        ## we can also flatten linear interactions
        vars <- strsplit(current, ":")
        split_terms <- lapply(unlist(vars), function(v) {
          split_linear_term(v, data)
        })

        combined_terms <- c()
        for (v1 in split_terms[[1]]) {
          for (v2 in split_terms[[2]]) {
            combined_terms <- c(combined_terms, paste0(v1, ":", v2))
          }
        }

        formula <- update(formula, paste0(
          ". ~ . - ", current,
          " + ",
          paste(combined_terms, collapse = " + ")
        ))
        ## no need to update the data because the interaction terms
        ## do not appear as features
      }
    }
  }

  tryCatch(model.matrix(formula, data = new_data),
           error = function(e) print(e))
  return(list(formula = formula, data = new_data))
}

## Splits a linear term into individual predictors.
## @param term A matrix term.
## @param data The original data frame.
## @return a list of the expanded linear matrix term.
split_linear_term <- function(term, data) {
  appearances <- ncol(data[, term])

  if (appearances > 1) {
    split_term <- sapply(
      1:appearances,
      function(i) paste0(term, i)
    )
  } else {
    split_term <- term
  }
  return(split_term)
}
