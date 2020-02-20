#' This function splits the given formula in response (left) and predictors (right).
#' @param formula Formula object that specifies a model.
#' @return a list containing the plain linear terms, interaction terms, group terms,
#' response and a boolean global intercept indicating whether the intercept is included
#' or not.
extract_terms_response <- function(formula) {
  tt <- terms(formula)
  terms <- attr(tt, "term.labels")
  variables <- as.list(attr(tt, "variables")[-1])
  response <- attr(tt, "response")
  global_intercept <- attr(tt, "intercept") == 1

  hier <- grepl("\\|", terms)
  int <- grepl(":", terms)
  group_terms <- terms[hier]
  interaction_terms <- terms[int & !hier]
  individual_terms <- terms[!hier & !int]

  if (response)
    response <- variables[response]
  else
    response <- NA

  nlist(individual_terms,
        interaction_terms,
        group_terms,
        response=extract_response(response),
        global_intercept)
}

#' At any point inside projpred, the response can be a single object or instead
#' it can represent multiple outputs. In this function we recover the response/s
#' as a character vector so we can index the dataframe.
#' @param response The response as retrieved from the formula object.
#' @return the response as a character vector.
extract_response <- function(response) {
  response_name_ch <- as.character(response[[1]])
  if ("cbind" %in% response_name_ch)
    ## remove cbind
    response_name_ch <- response_name_ch[-which(response_name_ch == "cbind")]
  else
    response_name_ch <- as.character(response)
  return(response_name_ch)
}

#' Because the formula can imply multiple or single response, in this function we
#' make sure that a formula has a single response. Then, in the case of multiple
#' responses we split the formula.
#' @param formula A formula specifying a model.
#' @return a formula or a list of formulas with a single response each.
validate_response_formula <- function(formula) {
  ## if multi-response model, split the formula to have
  ## one formula per model.
  ## only use this for mixed effects models, for glms or gams

  terms <- extract_terms_response(formula)
  response <- terms$response

  if (length(response) > 1)
    return(lapply(response, function(r)
      update(formula, paste0(r, " ~ ."))))
  else
    return(formula)
}

#' By combining different submodels we may arrive to a formula with repeated terms.
#' This function gets rid of duplicated terms.
#' @param formula A formula specifying a model.
#' @return a formula without duplicated structure.
flatten_formula <- function(formula) {
  terms <- extract_terms_response(formula)
  group_terms <- terms$group_terms
  interaction_terms <- terms$interaction_terms
  individual_terms <- terms$individual_terms

  if (length(individual_terms) > 0 ||
      length(interaction_terms) > 0 ||
      length(group_terms) > 0)
    update(formula, paste(c(". ~ ",
                            flatten_individual_terms(individual_terms),
                            flatten_interaction_terms(interaction_terms),
                            flatten_group_terms(group_terms)),
                          collapse=" + "))
  else
    formula
}

#' Remove duplicated linear terms.
#' @param terms A vector of linear terms as strings.
#' @return a vector of unique linear individual terms.
flatten_individual_terms <- function(terms) {
  if (is.null(terms))
    return(terms)
  return(unique(terms))
}

#' Remove duplicated linear interaction terms.
#' @param terms A vector of linear interaction terms as strings.
#' @return a vector of unique linear interaction terms.
flatten_interaction_terms <- function(terms) {
  if (length(terms) < 1)
    return(terms)
  ## TODO: do this right; a:b == b:a.
  return(unique(terms))
}

#' Simplify and flatten group terms by breaking them down in as many terms
#' as varying effects. It also explicitly adds or removes the varying intercept.
#' @param terms A vector of linear group terms as strings.
#' @return a vector of unique group terms.
flatten_group_terms <- function(terms) {
  if (is.null(terms))
    return(terms)
  split_terms <- strsplit(terms, "[ ]*\\|([^\\|]*\\||)[ ]*")
  group_names <- unique(unlist(lapply(split_terms, function(t) t[2])))

  groups <- setNames(lapply(group_names, function(g)
    lapply(split_terms, function(t) if (t[2] == g) t[1] else NA)),
    group_names)
  groups <- lapply(groups, function(g) unlist(g[!is.na(g)]))

  groups <- lapply(1:length(groups), function(i) {
    g <- groups[[i]]
    g_name <- group_names[i]

    partial_form <- as.formula(paste0(". ~ ", paste(g, collapse=" + ")))
    t <- terms(partial_form)
    t.labels <- attr(terms(partial_form), "term.labels")

    if ("1" %in% g)
      attr(t, "intercept") <- 1

    if (length(t.labels) < 1)
      t.labels <- c("1")

    if (!attr(t, "intercept"))
      paste0("(0 + ", paste(t.labels, collapse=" + "), " | ", g_name, ")")
    else
      paste0("(", paste(t.labels, collapse=" + "), " | ", g_name, ")")
  })
  return(unlist(groups))
}

#' Simplify and flatten group terms by breaking them down in as many terms
#' as varying effects. It also explicitly adds or removes the varying intercept.
#' @param formula A formula for a valid model.
#' @param return_group_terms If TRUE, return group terms as well. Default TRUE.
#' @return a vector of all the minimal valid terms that make up for submodels.
break_formula <- function(formula, return_group_terms=TRUE) {
  terms <- extract_terms_response(formula)
  group_terms <- terms$group_terms
  interaction_terms <- terms$interaction_terms
  individual_terms <- terms$individual_terms
  global_intercept <- terms$global_intercept

  ## if there are no group level terms then `terms` is the most basic
  ## decomposition already
  if (length(group_terms) == 0 && length(interaction_terms) == 0)
    return(individual_terms)

  ## if (length(group_terms) > 0 && return_group_terms)
  ##   warning(paste0("The model involves group terms.",
  ##                  " Make sure that the provided mle and proj_predfun",
  ##                  " handle multi-response projections."),
  ##           call. = FALSE)

  if (return_group_terms)
    ## if there are group levels we should split that into basic components
    allterms <- c(individual_terms,
                  unlist(lapply(group_terms, split_group_term)),
                  unlist(lapply(interaction_terms, split_interaction_term)))
  else
    allterms <- c(individual_terms,
                  unlist(lapply(interaction_terms, split_interaction_term)))

  ## exclude the intercept if there is no intercept in the reference model
  if (!global_intercept) {
    allterms_nobias <- unlist(lapply(allterms, function(term) paste0(term, " + 0")))
    return(allterms_nobias)
  } else
    return(allterms)
}

#' Plugs the main effects to the interaction terms to consider jointly for
#' projection.
#' @param term An interaction term as a string.
#' @return a minimally valid submodel for the interaction term including
#' the overall effects.
split_interaction_term <- function(term) {
  ## strong heredity by default
  variables <- unlist(strsplit(term, ":"))
  individual_joint <- paste(variables, collapse = " + ")
  joint_term <- paste(c(individual_joint, term), collapse = " + ")

  joint_term
}

#' Simplify a single group term by breaking it down in as many terms
#' as varying effects. It also explicitly adds or removes the varying intercept.
#' @param term A group term as a string.
#' @return a vector of all the minimally valid submodels for the group term
#' including a single varying effect with and without varying intercept.
split_group_term <- function(term) {
  ## this expands whatever terms() did not expand
  term <- gsub("\\)$", "",
               gsub("^\\(", "",
                    flatten_group_terms(term)))
  if ("\\-" %in% term)
    stop("Use of `-` is not supported, omit variables or use the ",
         "method update on the formula, or write `0 +` to remove ",
         "the intercept.")

  chunks <- strsplit(term, "[ ]*\\|([^\\|]*\\||)[ ]*")[[1]]
  lhs <- as.formula(paste0("~", chunks[1]))
  tt <- extact_terms_response(lhs)
  variables <- c(tt$individual_terms, tt$interaction_terms)

  ## split possible interaction terms
  int_t <- grepl(":", variables)
  ## int_v <- lapply(variables[!is.na(int_t)], split_interaction_term)
  int_v <- variables[int_t]
  lin_v <- variables[!int_t]

  ## variables <- c(int_v, lin_v)

  ## don't add the intercept twice
  variables <- setdiff(variables, "1")
  group <- chunks[2]

  if ("0" %in% variables)
    group_intercept <- FALSE
  else
    group_intercept <- TRUE

  group_terms <- lapply(lin_v, function(v) paste0(v, " + ", "(0 + ", v, " | ", group, ")"))
  group_terms <- c(group_terms, lapply(int_v, function(v)
    paste0(split_interaction_term(v), " + ", "(0 + ", v, " | ", group, ")")))

  if (group_intercept) {
    group_terms <- c(group_terms, list(paste0("(1 | ", group, ")")))
    group_terms <- c(group_terms,
                     lapply(lin_v,
                            function(v)
                              paste0(v, " + ", "(", v, " | ", group, ")")))
    group_terms <- c(group_terms,
                     lapply(int_v,
                            function(v)
                              paste0(split_interaction_term(v), " + ", "(", v, " | ", group, ")")))

    ## add v + ( 1 | group)
    group_terms <- c(group_terms,
                     lapply(lin_v,
                            function(v)
                              paste0(v, " + ", "(1 | ", group, ")")))
    group_terms <- c(group_terms,
                     lapply(int_v,
                            function(v)
                              paste0(split_interaction_term(v), " + ", "(1 | ", group, ")")))
  }

  group_terms
}

#' Checks whether a formula contains group terms or not.
#' @param formula A formula for a valid model.
#' @return TRUE if the formula contains group terms, FALSE otherwise.
formula_contains_group_terms <- function(formula) {
  t <- terms(formula)
  attributes <- attributes(t)
  terms <- attributes$term.labels

  hier <- grepl("\\|", terms)
  group_terms <- terms[hier]

  return(length(group_terms) > 0)
}

#' Subsets a formula by the given terms. It also properly subsets the data.
#' @param formula A formula for a valid model.
#' @param terms A vector of terms to subset.
#' @param data The original data frame for the full formula.
#' @param y The response vector. Default NULL.
#' @param split_formula If TRUE breaks the response down into single response formulas.
#' Default FALSE. It only works if `y` represents a multi-output response.
#' @return a list containing the subset formula and the data.
subset_formula <- function(formula, terms, data, y=NULL, split_formula=FALSE) {
  formula <- update(formula, paste0(". ~ ", paste(terms, collapse=" + ")))

  tt <- extract_terms_response(formula)
  response_name <- tt$response

  if (is.null(response_name))
    stop("No response")

  X <- data[, colnames(data) != response_name]
  cols_X <- colnames(X)

  group_terms <- formula_contains_group_terms(formula)

  if (is.null(y)) {
    y <- data[, response_name]
  }

  if (!is.null(ncol(y)) && ncol(y) > 1) {
    ## remove parens from response
    response_name <- gsub("\\(", "\\.",
                          gsub("\\)", "\\.",
                               response_name))
    cols_y <- unlist(lapply(1:ncol(y), function(n) paste0(response_name, ".", n)))
    if (!split_formula) {
      vector_y <- paste0("cbind(", paste(cols_y, collapse=", "), ")")
      formula <- update(formula, paste0(vector_y, " ~ ."))
    } else {
      formula <- lapply(cols_y, function(response) update(formula, paste0(response, " ~ .")))
    }
  } else
    cols_y <- response_name

  data <- data.frame(y=y, X=X)
  colnames(data) <- c(cols_y, cols_X)
  list(formula=formula, data=data)
}

#' Sometimes there can be terms in a formula that refer to a matrix instead of
#' a single predictor. Because we can handle groups of predictors, this function breaks
#' the matrix term into individual predictors to handle separately, as that is probably
#' the intention of the user.
#' @param formula A formula for a valid model.
#' @param data The original data frame with a matrix as predictor.
#' @return a  list containing the expanded formula and the expanded data frame.
break_up_matrix_term <- function(formula, data) {
  tt <- terms(formula)
  response <- attr(tt, "response")
  variables.list <- as.list(attr(tt, "variables")[-1])

  if (response)
    variables.list <- variables.list[-1]

  term.labels <- attr(tt, "term.labels")
  intercept <- attr(tt, "intercept")

  mm <- model.matrix(formula, data)
  assign <- attr(mm, "assign")

  new_data <- data
  for (assignee in unique(assign)) {
    if (assignee == 0) ## intercept
      next
    appearances <- assign[assign == assignee]
    if (length(appearances) > 1) {
      ## check if special term
      current <- term.labels[assignee]
      int <- grepl(":", current)
      mulilevel <- grepl("\\|", current)
      special <- grepl("[a-z]+\\(([a-z]+)\\)", current)
      individual <- !mulilevel & !int

      linear <- individual & !special
      linear_int <- int & !special

      if (linear) {
        ## if linear we can split it
        split_term <- split_linear_term(current, data)

        formula <- update(formula, paste0(". ~ . - ", current, " + ",
                                          paste(split_term, collapse=" + ")))
        split_matrix <- mm[, assign == assignee]

        new_data_tmp <- as.data.frame(new_data[, colnames(new_data) != current])
        colnames(new_data_tmp) <- colnames(new_data)[colnames(new_data) != current]
        new_data <- cbind(new_data_tmp, split_matrix)
      }

      if (linear_int) {
        ## we can also flatten linear interactions
        vars <- strsplit(current, ":")
        split_terms <- lapply(unlist(vars), function(v)
          split_linear_term(v, data))

        combined_terms <- c()
        for (v1 in split_terms[[1]])
          for (v2 in split_terms[[2]])
            combined_terms <- c(combined_terms, paste0(v1, ":", v2))

        formula <- update(formula, paste0(". ~ . - ", current,
                                          " + ",
                                          paste(combined_terms, collapse=" + ")))
        ## no need to update the data because the interaction terms
        ## do not appear as features
      }
    }
  }

  tryCatch(model.matrix(formula, data=new_data),
           error = function (e) print(e))
  list(formula=formula, data=new_data)
}

#' Splits a linear term into individual predictors.
#' @param term A matrix term.
#' @param data The original data frame.
#' @return a list of the expanded linear matrix term.
split_linear_term <- function(term, data) {
  appearances <- ncol(data[, term])

  if (appearances > 1)
    split_term <- sapply(1:appearances,
                         function(i) paste0(term, i))
  else
    split_term <- term
  split_term
}

#' Utility to count the number of terms in a given formula.
#' @param subformula A formula for a valid model.
#' @return the number of terms in the subformula.
count_terms_in_subformula <- function(subformula) {
  if (!inherits(subformula, "formula"))
    subformula <- as.formula(paste0("~ ", subformula))
  tt <- extract_terms_response(subformula)
  ind_interaction_terms <- length(tt$individual_terms) + length(tt$interaction_terms)
  group_terms <- sum(unlist(lapply(tt$group_terms, count_terms_in_group_term)))
  ind_interaction_terms + group_terms
}

#' Utility to count the number of terms in a given group term.
#' @param term A group term as a string.
#' @return the number of terms in the group term.
count_terms_in_group_term <- function(term) {
  term <- gsub("[\\(\\)]", "", flatten_group_terms(term))

  chunks <- strsplit(term,  "[ ]*\\|([^\\|]*\\||)[ ]*")[[1]]
  lhs <- as.formula(paste0("~", chunks[1]))
  tt <- extact_terms_response(lhs)
  variables <- c(tt$individual_terms, tt$interaction_terms)

  ## split possible interaction terms
  int_t <- grepl(":", variables)
  int_v <- variables[int_t]
  lin_v <- variables[!int_t]

  if ("0" %in% variables) {
    group_intercept <- FALSE
    variables <- setdiff(variables, "0")
  } else {
    group_intercept <- TRUE
    if (!("1" %in% variables))
      variables <- c(variables, "1")
  }

  length(variables)
}

#' Given a list of formulas, sort them by number of terms in them.
#' @param submodels A list of models' formulas.
#' @return a sorted list of submodels by included terms.
sort_submodels_by_size <- function(submodels) {
  size <- lapply(submodels, count_terms_in_subformula)
  df <- as.data.frame(cbind(submodels=submodels, size=unlist(size)))
  df$size <- as.numeric(df$size)
  ordered <- df[order(df$size), ]

  groups <- list()
  for (size in unique(ordered$size))
    groups[[size]] <- as.character(ordered$submodels[ordered$size == size])

  ord_list <- groups
  ord_list_nona <- lapply(ord_list, function(l) l[!is.na(l)])
  ord_list_nona[!is.na(ord_list_nona)]
}

#' Given a refmodel structure, count the number of terms included.
#' @param refmodel A refmodel structure containing at least a formula and a
#' function fetch_data.
#' @param list_of_terms Subset of terms from refmodel.
#' @return number of terms in list_of_terms.
count_variables_chosen <- function(refmodel, list_of_terms) {
  if (length(list_of_terms) == 0)
    return(0)
  form <- subset_formula(refmodel$formula, unique(unlist(list_of_terms)),
                         data=refmodel$fetch_data())$formula
  count_terms_in_subformula(flatten_formula(form))
}

#' Utility that checks if the next submodel is redundant with the current one.
#' @param refmodel A refmodel struct.
#' @param current A list of terms included in the current submodel.
#' @param new The new term to add to the submodel.
#' @return TRUE if the new term results in a redundant model, FALSE otherwise.
is_next_submodel_redundant <- function(refmodel, current, new) {
  old_submodel <- current
  new_submodel <- c(current, new)
  if (count_variables_chosen(refmodel, new_submodel) >
      count_variables_chosen(refmodel, old_submodel))
    FALSE
  else
    TRUE
}

#' Utility to remove redundant models to consider
#' @param refmodel A refmodel struct.
#' @param chosen A list of included terms at a given point of the search.
#' @return a vector of incremental non redundant submodels for all the possible terms
#' included.
reduce_models <- function(refmodel, chosen) {
  Reduce(function(chosen, x)
    if (is_next_submodel_redundant(refmodel, chosen, x)) chosen else c(chosen, x),
    chosen)
}
