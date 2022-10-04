## This function splits the given formula in response (left) and predictors
## (right).
## @param formula Formula object that specifies a model.
## @return a list containing the plain linear terms, interaction terms, group
##   terms, response and a boolean global intercept indicating whether the
##   intercept is included or not.
extract_terms_response <- function(formula) {
  tt <- terms(formula)
  terms_ <- attr(tt, "term.labels")
  ## when converting the terms_ to a list the first element is
  ## "list" itself, so we remove it
  allterms_ <- as.list(attr(tt, "variables")[-1])
  response <- attr(tt, "response")
  global_intercept <- attr(tt, "intercept") == 1
  offs_attr <- attr(tt, "offset")

  if (response) {
    response <- allterms_[response]
  } else {
    response <- NA
  }

  if (length(offs_attr)) {
    offset_terms <- sapply(allterms_[offs_attr], deparse)
  } else {
    offset_terms <- NULL
  }

  hier <- grepl("\\|", terms_)
  int <- grepl(":", terms_)
  group_terms <- terms_[hier]
  interaction_terms <- terms_[int & !hier]
  individual_terms <- terms_[!hier & !int]
  additive_terms <- parse_additive_terms(individual_terms)
  individual_terms <- setdiff(individual_terms, additive_terms)

  response <- extract_response(response)
  return(nlist(
    individual_terms,
    interaction_terms,
    additive_terms,
    group_terms,
    response,
    global_intercept,
    offset_terms
  ))
}

expand_formula <- function(formula, data) {
  return(formula(terms(formula, data = data)))
}

remove_duplicates <- function(formula) {
  terms <- extract_terms_response(formula)
  linear <- terms$individual_terms
  additive <- unlist(regmatches(
    terms$additive_terms,
    gregexpr("(?<=\\().*?(?=\\))",
             terms$additive_terms,
             perl = TRUE)
  ))
  additive <- trimws(unique(unlist(
    strsplit(paste0(additive, collapse = ","), ",")
  )))
  dups <- linear[!is.na(match(linear, additive))]
  if (length(dups) > 0) {
    update(formula, as.formula(paste0(
      ". ~ . - ",
      paste(dups, collapse = " - ")
    )))
  } else {
    formula
  }
}

## At any point inside projpred, the response can be a single object or instead
## it can represent multiple outputs. In this function we recover the response/s
## as a character vector so we can index the dataframe.
## @param response The response as retrieved from the formula object.
## @return the response as a character vector.
extract_response <- function(response) {
  if (length(response) > 1) {
    stop("Object `response` must not have length greater than 1.")
  }
  response_name_ch <- as.character(response[[1]])
  if ("cbind" %in% response_name_ch) {
    ## remove cbind
    response_name_ch <- response_name_ch[-which(response_name_ch == "cbind")]
  } else {
    response_name_ch <- as.character(response)
  }
  return(response_name_ch)
}

## Parse additive terms (smooth terms) from a list of individual terms. See
## `?init_refmodel` for allowed smooth terms.
## @param terms list of terms to parse
## @return a vector of smooth terms
parse_additive_terms <- function(terms) {
  excluded_terms <- c("te")
  smooth_terms <- c("s", "t2")
  excluded <- unlist(sapply(excluded_terms, function(et) {
    grep(make_function_regexp(et), terms)
  }))
  if (sum(excluded) > 0) {
    stop("te terms are not supported, please use t2 instead.")
  }
  smooth <- sapply(smooth_terms, function(et) {
    terms[grep(make_function_regexp(et), terms)]
  }) %>% unlist() %>% unname()
  return(smooth)
}

make_function_regexp <- function(fname) {
  return(paste0(fname, "\\(.+\\)"))
}

## Because the formula can imply multiple or single response, in this function
## we make sure that a formula has a single response. Then, in the case of
## multiple responses we split the formula.
## @param formula A formula specifying a model.
## @return a formula or a list of formulas with a single response each.
validate_response_formula <- function(formula) {
  ## if multi-response model, split the formula to have
  ## one formula per model.
  ## only use this for mixed effects models, for glms or gams

  terms_ <- extract_terms_response(formula)
  response <- terms_$response

  if (length(response) > 1) {
    return(lapply(response, function(r) {
      update(formula, paste0(r, " ~ ."))
    }))
  } else {
    return(list(formula))
  }
}

## By combining different submodels we may arrive to a formula with repeated
## terms.
## This function gets rid of duplicated terms_.
## @param formula A formula specifying a model.
## @param duplicates if FALSE removes linear terms if their corresponding smooth
## is included. Default TRUE
## @return a formula without duplicated structure.
flatten_formula <- function(formula, duplicates = TRUE) {
  terms_ <- extract_terms_response(formula)
  group_terms <- terms_$group_terms
  interaction_terms <- terms_$interaction_terms
  individual_terms <- terms_$individual_terms
  additive_terms <- terms_$additive_terms

  if (length(individual_terms) > 0 ||
      length(interaction_terms) > 0 ||
      length(group_terms) > 0 ||
      length(additive_terms) > 0) {
    full <- update(
      formula,
      paste(
        c(". ~ ",
          flatten_individual_terms(individual_terms),
          flatten_additive_terms(additive_terms),
          flatten_interaction_terms(interaction_terms),
          flatten_group_terms(group_terms)),
        collapse = " + "
      )
    )
    if (!duplicates)
      remove_duplicates(full)
    else
      full
  } else {
    formula
  }
}

## Remove duplicated linear terms.
## @param terms A vector of linear terms as strings.
## @return a vector of unique linear individual terms.
flatten_individual_terms <- function(terms_) {
  if (length(terms_) == 0) {
    return(terms_)
  }
  return(unique(terms_))
}

## Remove duplicated linear interaction terms.
## @param terms A vector of linear interaction terms as strings.
## @return a vector of unique linear interaction terms.
flatten_interaction_terms <- function(terms_) {
  if (length(terms_) == 0) {
    return(terms_)
  }
  ## TODO: do this right; a:b == b:a.
  return(unique(terms_))
}

## Remove duplicated additive terms.
## @param terms A vector of additive terms as strings.
## @return a vector of unique linear interaction terms.
flatten_additive_terms <- function(terms_) {
  if (length(terms_) == 0) {
    return(terms_)
  }
  return(unique(terms_))
}

## Unifies group terms removing any duplicates.
## @param terms A vector of linear group terms as strings.
## @return a vector of unique group terms.
flatten_group_terms <- function(terms_) {
  if (length(terms_) == 0) {
    return(terms_)
  }
  split_terms_ <- strsplit(terms_, "[ ]*\\|([^\\|]*\\||)[ ]*")
  group_names <- unique(unlist(lapply(split_terms_, function(t) t[2])))

  group_terms <- setNames(
    lapply(group_names, function(g) {
      lapply(split_terms_, function(t) if (t[2] == g) t[1] else NA)
    }),
    group_names
  )
  group_terms <- lapply(group_terms, function(g) unlist(g[!is.na(g)]))

  group_terms <- lapply(seq_along(group_terms), function(i) {
    g <- group_terms[[i]]
    g_name <- group_names[i]

    partial_form <- as.formula(paste0(". ~ ", paste(g, collapse = " + ")))
    t <- terms(partial_form)
    t.labels <- attr(terms(partial_form), "term.labels")

    if ("1" %in% g) {
      attr(t, "intercept") <- 1
    }

    if (length(t.labels) < 1) {
      t.labels <- c("1")
    }

    if (!attr(t, "intercept")) {
      return(paste0(
        "(0 + ", paste(t.labels, collapse = " + "), " | ", g_name, ")"
      ))
    } else {
      return(paste0(
        "(", paste(t.labels, collapse = " + "), " | ", g_name, ")"
      ))
    }
  })
  return(unlist(group_terms))
}

## Simplify and split a formula by breaking it into all possible submodels.
## @param formula A formula for a valid model.
## @param return_group_terms If TRUE, return group terms as well. Default TRUE.
## @param data The reference model data.
## @return a vector of all the minimal valid terms that make up for submodels.
split_formula <- function(formula, return_group_terms = TRUE, data = NULL,
                          add_main_effects = TRUE) {
  terms_ <- extract_terms_response(formula)
  group_terms <- terms_$group_terms
  interaction_terms <- terms_$interaction_terms
  individual_terms <- terms_$individual_terms
  additive_terms <- terms_$additive_terms
  global_intercept <- terms_$global_intercept

  additive <- unlist(regmatches(
    additive_terms,
    gregexpr("(?<=\\().*?(?=\\))",
             terms_$additive_terms,
             perl = TRUE)
  ))
  additive <- trimws(unique(unlist(
    strsplit(paste0(additive, collapse = ","), ",")
  )))
  if (return_group_terms) {
    ## if there are group levels we should split that into basic components
    group_split <- unlist(lapply(group_terms, split_group_term,
                                 add_main_effects = add_main_effects))
    allterms_ <- c(
      unlist(lapply(additive_terms, split_additive_term, data)),
      unlist(lapply(interaction_terms, split_interaction_term,
                    add_main_effects = add_main_effects))
    )
    group_replace <- regmatches(
      group_split,
      gregexpr("\\w+(?![^(]*\\))", group_split, perl = TRUE)
    )
    groups_to_replace <- group_split[unlist(lapply(
      group_replace,
      function(x) length(x) > 0
    ))]
    to_replace <- group_split[match(group_replace, additive) %>%
                                (function(x) !is.na(x))]
    not_replace <- setdiff(group_split, to_replace)

    replacement <- gsub(
      pattern = "(\\w+)(?![^(]*\\))", replacement = "s(\\1)",
      to_replace, perl = TRUE
    )
    group_split <- c(not_replace, replacement)
    nodups <- unique(c(individual_terms, additive))
    allterms_ <- c(allterms_, group_split, nodups)
  } else {
    nodups <- unique(c(individual_terms, additive))
    allterms_ <- c(
      nodups,
      unlist(lapply(additive_terms, split_additive_term, data)),
      unlist(lapply(interaction_terms, split_interaction_term))
    )
  }

  ## exclude the intercept if there is no intercept in the reference model
  if (!global_intercept) {
    allterms_nobias <- unlist(lapply(allterms_, function(term) {
      paste0(term, " + 0")
    }))
    return(unique(allterms_nobias))
  } else {
    return(c("1", unique(allterms_)))
  }
}

## Plugs the main effects to the interaction terms to consider jointly for
## projection.
## @param term An interaction term as a string.
## @return a minimally valid submodel for the interaction term including
## the overall effects.
split_interaction_term <- function(term, add_main_effects = TRUE) {
  ## strong heredity by default
  terms_ <- unlist(strsplit(term, ":"))
  individual_joint <- paste(terms_, collapse = " + ")
  if (add_main_effects) {
    joint_term <- paste(c(individual_joint, term), collapse = " + ")
  } else {
    return(term)
  }
  return(joint_term)
}

## Plugs the main effects to the smooth additive term if `by` argument is
## provided.
## @param term An additive term as a string.
## @param data The reference model data.
## @return a minimally valid submodel for the additive term.
split_additive_term <- function(term, data) {
  out <- mgcv::interpret.gam(as.formula(paste("~", term)))
  if (out$smooth.spec[[1]]$by == "NA") {
    return(term)
  }
  main <- out$smooth.spec[[1]]$by
  fac <- eval_rhs(as.formula(paste("~", main)), data = data)
  if (!is.factor(fac)) {
    return(term)
  }
  joint_term <- paste(c(main, term), collapse = " + ")
  return(c(main, joint_term))
}

## Simplify a single group term by breaking it down in as many terms_
## as varying effects. It also explicitly adds or removes the varying intercept.
## @param term A group term as a string.
## @return a vector of all the minimally valid submodels for the group term
## including a single varying effect with and without varying intercept.
split_group_term <- function(term, add_main_effects = TRUE) {
  ## this expands whatever terms() did not expand
  term <- gsub(
    "\\)$", "",
    gsub(
      "^\\(", "",
      flatten_group_terms(term)
    )
  )
  ## if ("\\-" %in% term) {
  ##   stop("Use of `-` is not supported, omit terms or use the ",
  ##        "method update on the formula, or write `0 +` to remove ",
  ##        "the intercept.")
  ## }

  chunks <- strsplit(term, "[ ]*\\|([^\\|]*\\||)[ ]*")[[1]]
  lhs <- as.formula(paste0("~", chunks[1]))
  tt <- extract_terms_response(lhs)
  terms_ <- c(tt$individual_terms, tt$interaction_terms)

  ## split possible interaction terms
  int_t <- grepl(":", terms_)
  int_v <- terms_[int_t]
  lin_v <- terms_[!int_t]

  ## don't add the intercept twice
  terms_ <- setdiff(terms_, "1")
  group <- chunks[2]

  group_intercept <- tt$global_intercept

  if (group_intercept) {
    group_terms <- list(paste0("(1 | ", group, ")"))
    if (add_main_effects) {
      group_terms <- c(
        group_terms,
        lapply(lin_v, function(v) {
          paste0(v, " + ", "(", v, " | ", group, ")")
        })
      )
      group_terms <- c(
        group_terms,
        lapply(int_v, function(v) {
          paste0(
            split_interaction_term(v, add_main_effects = add_main_effects),
            " + ",
            "(",
            split_interaction_term(v, add_main_effects = add_main_effects),
            " | ",
            group,
            ")"
          )
        })
      )

      ## add v + ( 1 | group)
      group_terms <- c(
        group_terms,
        lapply(lin_v, function(v) {
          paste0(v, " + ", "(1 | ", group, ")")
        })
      )
      group_terms <- c(
        group_terms,
        lapply(int_v, function(v) {
          paste0(
            split_interaction_term(v, add_main_effects = add_main_effects),
            " + ",
            "(1 | ",
            group,
            ")"
          )
        })
      )
    } else {
      group_terms <- c(
        group_terms,
        lapply(lin_v, function(v) {
          paste0("(", v, " | ", group, ")")
        })
      )
      group_terms <- c(
        group_terms,
        lapply(int_v, function(v) {
          paste0(
            "(",
            split_interaction_term(v, add_main_effects = add_main_effects),
            " | ",
            group,
            ")"
          )
        })
      )
    }
  } else {
    group_terms <- lapply(lin_v, function(v) {
      paste0(v, " + ", "(0 + ", v, " | ", group, ")")
    })
    group_terms <- c(group_terms, lapply(int_v, function(v) {
      paste0(
        split_interaction_term(v, add_main_effects = add_main_effects),
        " + ",
        "(0 + ",
        split_interaction_term(v, add_main_effects = add_main_effects),
        " | ",
        group,
        ")"
      )
    }))
  }

  return(group_terms)
}

## Checks whether a formula contains group terms or not.
## @param formula A formula for a valid model.
## @return TRUE if the formula contains group terms, FALSE otherwise.
formula_contains_group_terms <- function(formula) {
  group_terms <- extract_terms_response(formula)$group_terms
  return(length(group_terms) > 0)
}

## Checks whether a formula contains additive terms or not.
## @param formula A formula for a valid model.
## @return TRUE if the formula contains additive terms, FALSE otherwise.
formula_contains_additive_terms <- function(formula) {
  additive_terms <- extract_terms_response(formula)$additive_terms
  return(length(additive_terms) > 0)
}

## Utility to both subset the formula and update the data
## @param formula A formula for a valid model.
## @param terms_ A vector of terms to subset.
## @param data The original data frame for the full formula.
## @param y The response vector. Default NULL.
## @param split_formula If TRUE breaks the response down into single response
##   formulas.
## Default FALSE. It only works if `y` represents a multi-output response.
## @return a list including the updated formula and data
subset_formula_and_data <- function(formula, terms_, data, y = NULL,
                                    split_formula = FALSE) {
  formula <- make_formula(terms_, formula = formula)
  response_name <- extract_terms_response(formula)$response

  response_cols <- paste0(".", response_name)

  if (NCOL(y) > 1) {
    response_cols <- paste0(response_cols, ".", seq_len(ncol(y)))
    if (!split_formula) {
      response_vector <- paste0(
        "cbind(",
        paste(response_cols, collapse = ", "),
        ")"
      )
      formula <- update(formula, paste0(response_vector, " ~ ."))
    } else {
      formula <- lapply(response_cols, function(response) {
        update(formula, paste0(response, " ~ ."))
      })
    }
  } else {
    formula <- update(formula, paste(response_cols, "~ ."))
  }

  ## don't overwrite original y name
  data <- data.frame(.z = y, data)
  colnames(data)[seq_len(NCOL(y))] <- response_cols
  return(nlist(formula, data))
}

## Subsets a formula by the given terms.
## @param terms_ A vector of terms to subset from the right hand side.
## @return A formula object with the collapsed terms.
make_formula <- function(terms_, formula = NULL) {
  if (length(terms_) == 0) {
    terms_ <- c("1")
  }
  if (is.null(formula)) {
    return(as.formula(paste0(". ~ ", paste(terms_, collapse = " + "))))
  }
  return(update(formula, paste0(". ~ ", paste(terms_, collapse = " + "))))
}

## Utility to count the number of terms in a given formula.
## @param formula The right hand side of a formula for a valid model
## either as a formula object or as a string.
## @return the number of terms in the formula.
count_terms_in_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    formula <- as.formula(paste0("~ ", formula))
  }
  tt <- extract_terms_response(formula)
  ind_interaction_terms <- length(tt$individual_terms) +
    length(tt$interaction_terms) +
    length(tt$additive_terms)
  group_terms <- sum(unlist(lapply(tt$group_terms, count_terms_in_group_term)))
  return(ind_interaction_terms + group_terms + tt$global_intercept)
}

## Utility to count the number of terms in a given group term.
## @param term A group term as a string.
## @return the number of terms in the group term.
count_terms_in_group_term <- function(term) {
  term <- gsub("[\\(\\)]", "", flatten_group_terms(term))

  chunks <- strsplit(term, "[ ]*\\|([^\\|]*\\||)[ ]*")[[1]]
  lhs <- as.formula(paste0("~", chunks[1]))
  tt <- extract_terms_response(lhs)
  terms_ <- c(tt$individual_terms, tt$interaction_terms)

  if ("0" %in% terms_) {
    terms_ <- setdiff(terms_, "0")
  } else if (!("1" %in% terms_)) {
    terms_ <- c(terms_, "1")
  }

  return(length(terms_))
}

## Given a list of formulas, sort them by the number of terms in them.
## @param submodels A list of models' formulas.
## @return a sorted list of submodels by included terms.
sort_submodels_by_size <- function(submodels) {
  size <- lapply(submodels, count_terms_in_formula)
  df <- data.frame(submodels = as.character(submodels), size = unlist(size))
  ordered <- df[order(df$size), ]

  search_terms <- list()
  for (size in unique(ordered$size)) {
    search_terms[[size]] <- as.character(
      ordered$submodels[ordered$size == size]
    )
  }

  ord_list <- search_terms
  ## remove NA inside submodels
  ord_list_nona <- lapply(ord_list, function(l) l[!is.na(l)])
  ## remove NA at the submodels level
  return(ord_list_nona[!is.na(ord_list_nona)])
}

## Select next possible terms without surpassing a specific size
## @param chosen A list of currently chosen terms
## @param terms A list of all possible terms
## @param size Maximum allowed size
select_possible_terms_size <- function(chosen, terms, size) {
  if (size < 1) {
    stop("size must be at least 1")
  }

  valid_submodels <- lapply(terms, function(x) {
    ## if we are adding a linear term whose smooth is already
    ## included, we reject it
    terms <- extract_terms_response(make_formula(c(chosen)))
    terms_new <- extract_terms_response(make_formula(x))
    additive <- unlist(regmatches(
      terms$additive_terms,
      gregexpr("(?<=\\().*?(?=\\))",
               terms$additive_terms,
               perl = TRUE)
    ))
    linear <- terms_new$individual_terms
    dups <- setdiff(linear[!is.na(match(linear, additive))], chosen)

    size_crr <- count_terms_chosen(c(chosen, x), duplicates = TRUE) -
      length(dups)
    if (size_crr == size) {
      if (length(dups) > 0) {
        tt <- terms(formula(paste("~", x, "-", paste(dups, collapse = "-"))))
        x <- setdiff(attr(tt, "term.labels"), chosen)
        if (grepl("\\|", x) && !grepl("[()]", x)) {
          x <- paste0("(", x, ")")
        }
      }
      return(x)
    } else {
      return(NA)
    }
  })
  valid_submodels <- unlist(valid_submodels[!is.na(valid_submodels)])
  if (length(chosen) > 0) {
    add_chosen <- paste0(" + ", paste(chosen, collapse = "+"))
    remove_chosen <- paste0(" - ", paste(chosen, collapse = "-"))
  } else {
    add_chosen <- ""
    remove_chosen <- ""
  }
  full_valid_submodels <- unique(unlist(lapply(valid_submodels, function(x) {
    to_character_rhs(flatten_formula(make_formula(
      paste(x, add_chosen, remove_chosen)
    )))
  })))
  return(full_valid_submodels)
}

## Cast a right hand side formula to a character vector.
## @param rhs a right hand side formula of the type . ~ x + z + ...
## @return a character vector containing only the right hand side.
to_character_rhs <- function(rhs) {
  chr <- as.character(rhs)
  return(chr[length(chr)])
}

## Given a refmodel structure, count the number of terms included.
## @param list_of_terms Subset of terms from formula.
## @param duplicates if FALSE removes linear terms if their corresponding smooth
##   is included. Default TRUE
## @param add_icpt Only relevant if `length(list_of_terms) == 0`. A single
##   logical value indicating whether to add the intercept.
## @return number of terms
count_terms_chosen <- function(list_of_terms, duplicates = TRUE,
                               add_icpt = FALSE) {
  if (length(list_of_terms) == 0) {
    if (!add_icpt) {
      return(0)
    } else {
      list_of_terms <- "1"
    }
  }
  formula <- make_formula(list_of_terms)
  return(
    count_terms_in_formula(flatten_formula(formula, duplicates = duplicates))
  )
}

## Helper function to evaluate right hand side formulas in a context
## @param formula Formula to evaluate.
## @param data Data with which to evaluate.
## @return output from the evaluation
eval_rhs <- function(formula, data) {
  eval_el2(formula = formula, data = data)
}

# Helper function to evaluate the second element of a `formula` in a specific
# environment `data`.
#
# @param formula A `formula` whose second element should be evaluated.
# @param data Passed to argument `envir` of eval().
#
# @return The output from eval().
eval_el2 <- function(formula, data) {
  eval(formula[[2]], data, environment(formula))
}

## Extract left hand side of a formula as a formula itself by removing the right
## hand side.
## @param x Formula
## @return updated formula with an intercept in the right hand side.
lhs <- function(x) {
  x <- as.formula(x)
  if (length(x) == 3L) update(x, . ~ 1) else NULL
}

## remove intercept from formula
## @param formula a model formula
## @return the updated formula without intercept
delete.intercept <- function(formula) {
  return(update(formula, . ~ . - 1))
}

## collapse a list of terms including contrasts
## @param formula model's formula
## @param path list of terms possibly including contrasts
## @param data model's data
## @return the updated list of terms replacing the contrasts with the term name
collapse_contrasts_solution_path <- function(formula, path, data) {
  tt <- terms(formula)
  terms_ <- attr(tt, "term.labels")
  for (term in terms_) {
    # TODO: In the following model.matrix() call, allow user-specified contrasts
    # to be passed to argument `contrasts.arg`. The `contrasts.arg` default
    # (`NULL`) uses `options("contrasts")` internally, but it might be more
    # convenient to let users specify contrasts directly. At that occasion,
    # contrasts should also be tested thoroughly (not done until now).
    x <- model.matrix(as.formula(paste("~ 1 +", term)), data = data)
    if (length(attr(x, "contrasts")) == 0) {
      next
    }
    x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
    path <- Reduce(
      function(current, pattern) {
        pattern <- gsub("\\+", "\\\\+", pattern)
        list(current[[1]],
             gsub(pattern, current[[1]], current[[2]]))
      },
      x = colnames(x),
      init = list(term, path)
    )
    path <- unique(path[[length(path)]])
  }
  return(path)
}

split_formula_random_gamm4 <- function(formula) {
  tt <- extract_terms_response(formula)
  if (length(tt$group_terms) == 0) {
    return(nlist(formula, random = NULL))
  }
  parens_group_terms <- unlist(lapply(tt$group_terms, function(t) {
    paste0("(", t, ")")
  }))
  random <- as.formula(paste(
    "~",
    paste(parens_group_terms, collapse = " + ")
  ))
  formula <- update(formula, make_formula(c(
    tt$individual_terms, tt$interaction_terms, tt$additive_terms,
    tt$offset_terms
  )))
  return(nlist(formula, random))
}

# utility to recover the full gam + random formula from a stan_gamm4 model
formula.gamm4 <- function(x) {
  formula <- x$formula
  if (is.null(x$glmod)) {
    return(formula)
  }
  ref <- extract_terms_response(x$glmod$formula)$group_terms
  ref <- unlist(lapply(ref, function(t) {
    paste0("(", t, ")")
  }))
  ref <- paste(ref, collapse = " + ")
  updated <- update(formula, paste(". ~ . + ", ref))
  form <- flatten_formula(update(
    updated,
    paste(
      ". ~",
      paste(split_formula(updated), collapse = " + ")
    )
  ))
  # TODO (GAMMs): Once rstanarm issue #253 has been resolved, we probably need
  # to include offset terms here (in the output of formula.gamm4()) as well.
  return(form)
}
