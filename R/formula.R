#' This function splits the given formula in response (left) and predictors (right).
#' @param formula Formula object that specifies a model.
#' @return a list containing the plain linear terms, interaction terms, group terms,
#' response and a boolean global intercept indicating whether the intercept is included
#' or not.
extract_terms_response <- function(formula) {
  tt <- terms(formula)
  terms <- attr(tt, "term.labels")
  ## when converting the variables to a list the first element is
  ## "list" itself, so we remove it
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

  response <- extract_response(response)
  nlist(individual_terms,
        interaction_terms,
        group_terms,
        response,
        global_intercept)
}

#' At any point inside projpred, the response can be a single object or instead
#' it can represent multiple outputs. In this function we recover the response/s
#' as a character vector so we can index the dataframe.
#' @param response The response as retrieved from the formula object.
#' @return the response as a character vector.
extract_response <- function(response) {
  if (length(response) > 1)
    stop("Response must have a single element.")
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
  if (length(terms) == 0)
    return(terms)
  return(unique(terms))
}

#' Remove duplicated linear interaction terms.
#' @param terms A vector of linear interaction terms as strings.
#' @return a vector of unique linear interaction terms.
flatten_interaction_terms <- function(terms) {
  if (length(terms) == 0)
    return(terms)
  ## TODO: do this right; a:b == b:a.
  return(unique(terms))
}

#' Unifies group terms removing any duplicates.
#' @param terms A vector of linear group terms as strings.
#' @return a vector of unique group terms.
flatten_group_terms <- function(terms) {
  if (length(terms) == 0)
    return(terms)
  split_terms <- strsplit(terms, "[ ]*\\|([^\\|]*\\||)[ ]*")
  group_names <- unique(unlist(lapply(split_terms, function(t) t[2])))

  group_terms <- setNames(lapply(group_names, function(g)
    lapply(split_terms, function(t) if (t[2] == g) t[1] else NA)),
    group_names)
  group_terms <- lapply(group_terms, function(g) unlist(g[!is.na(g)]))

  group_terms <- lapply(seq_along(group_terms), function(i) {
    g <- group_terms[[i]]
    g_name <- group_names[i]

    partial_form <- as.formula(paste0(". ~ ", paste(g, collapse=" + ")))
    t <- terms(partial_form)
    t.labels <- attr(terms(partial_form), "term.labels")

    if ("1" %in% g)
      attr(t, "intercept") <- 1

    if (length(t.labels) < 1)
      t.labels <- c("1")

    if (!attr(t, "intercept"))
      return(paste0("(0 + ", paste(t.labels, collapse=" + "), " | ", g_name, ")"))
    else
      return(paste0("(", paste(t.labels, collapse=" + "), " | ", g_name, ")"))
  })
  return(unlist(group_terms))
}

#' Simplify and split a formula by breaking it into all possible submodels.
#' @param formula A formula for a valid model.
#' @param return_group_terms If TRUE, return group terms as well. Default TRUE.
#' @return a vector of all the minimal valid terms that make up for submodels.
split_formula <- function(formula, return_group_terms=TRUE) {
  terms <- extract_terms_response(formula)
  group_terms <- terms$group_terms
  interaction_terms <- terms$interaction_terms
  individual_terms <- terms$individual_terms
  global_intercept <- terms$global_intercept

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
    return(unique(allterms_nobias))
  } else
    return(unique(allterms))
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
  ## if ("\\-" %in% term)
  ##   stop("Use of `-` is not supported, omit variables or use the ",
  ##        "method update on the formula, or write `0 +` to remove ",
  ##        "the intercept.")

  chunks <- strsplit(term, "[ ]*\\|([^\\|]*\\||)[ ]*")[[1]]
  lhs <- as.formula(paste0("~", chunks[1]))
  tt <- extract_terms_response(lhs)
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

  group_intercept <- tt$global_intercept

  if (group_intercept) {
    group_terms <- list(paste0("(1 | ", group, ")"))
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
  } else {
    group_terms <- lapply(lin_v, function(v) paste0(v, " + ", "(0 + ", v, " | ", group, ")"))
    group_terms <- c(group_terms, lapply(int_v, function(v)
      paste0(split_interaction_term(v), " + ", "(0 + ", v, " | ", group, ")")))
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

#' Utility to both subset the formula and update the data
#' @param formula A formula for a valid model.
#' @param terms A vector of terms to subset.
#' @param data The original data frame for the full formula.
#' @param y The response vector. Default NULL.
#' @param split_formula If TRUE breaks the response down into single response formulas.
#' Default FALSE. It only works if `y` represents a multi-output response.
#' @return a list including the updated formula and data
subset_formula_and_data <- function(formula, terms, data, y=NULL, split_formula=FALSE) {
  formula <- make_formula(terms, formula=formula)
  tt <- extract_terms_response(formula)
  response_name <- tt$response

  response_cols <- paste0(".", response_name)
  response_ncol <- ncol(y) %||% 1

  if (!is.null(ncol(y)) && ncol(y) > 1) {
    response_cols <- paste0(response_cols, ".", seq_len(ncol(y)))
    if (!split_formula) {
      response_vector <- paste0("cbind(", paste(response_cols, collapse=", "), ")")
      formula <- update(formula, paste0(response_vector, " ~ ."))
    } else {
      formula <- lapply(response_cols, function(response)
        update(formula, paste0(response, " ~ .")))
    }
  }

  data <- data.frame(y, data[, colnames(data) != response_name])
  colnames(data)[1:response_ncol] <- response_cols
  return(nlist(formula, data))
}

#' Subsets a formula by the given terms.
#' @param terms A vector of terms to subset from the right hand side.
#' @return A formula object with the collapsed terms.
make_formula <- function(terms, formula=NULL) {
  if (is.null(formula))
    return(as.formula(paste0(". ~ ", paste(terms, collapse=" + "))))
  return(update(formula, paste0(". ~ ", paste(terms, collapse=" + "))))
}

#' Utility to count the number of terms in a given formula.
#' @param subformula The right hand side of a formula for a valid model
#' either as a formula object or as a string.
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
  tt <- extract_terms_response(lhs)
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
  df <- data.frame(submodels=as.character(submodels), size=unlist(size))
  ordered <- df[order(df$size), ]

  groups <- list()
  for (size in unique(ordered$size))
    groups[[size]] <- as.character(ordered$submodels[ordered$size == size])

  ord_list <- groups
  ## remove NA inside submodels
  ord_list_nona <- lapply(ord_list, function(l) l[!is.na(l)])
  ## remove NA at the submodels level
  ord_list_nona[!is.na(ord_list_nona)]
}

#' Given a refmodel structure, count the number of variables included.
#' @param formula The reference model's formula.
#' @param list_of_terms Subset of terms from formula.
#' @return number of variables in list_of_terms.
count_variables_chosen <- function(list_of_terms) {
  if (length(list_of_terms) == 0)
    return(0)
  formula <- make_formula(list_of_terms)
  count_terms_in_subformula(flatten_formula(formula))
}

#' Utility that checks if the next submodel is redundant with the current one.
#' @param formula The reference models' formula.
#' @param current A list of terms included in the current submodel.
#' @param new The new term to add to the submodel.
#' @return TRUE if the new term results in a redundant model, FALSE otherwise.
is_next_submodel_redundant <- function(current, new) {
  old_submodel <- current
  new_submodel <- c(current, new)
  if (count_variables_chosen(new_submodel) >
      count_variables_chosen(old_submodel))
    FALSE
  else
    TRUE
}

#' Utility to remove redundant models to consider
#' @param refmodel The reference model's formula.
#' @param chosen A list of included terms at a given point of the search.
#' @return a vector of incremental non redundant submodels for all the possible terms
#' included.
reduce_models <- function(chosen) {
  Reduce(function(chosen, x)
    if (is_next_submodel_redundant(chosen, x)) chosen
    else c(chosen, x),
    chosen)
}
