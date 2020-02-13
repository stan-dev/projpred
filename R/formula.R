library(stringr)

extract_terms_response <- function(formula) {
  tt <- terms(formula)
  terms <- attr(tt, "term.labels")
  variables <- as.list(attr(tt, "variables")[-1])
  response <- attr(tt, "response")
  global_intercept <- attr(tt, "intercept") == 1

  hier <- str_match(terms, "\\|")
  int <- str_match(terms, ":")
  group_terms <- terms[!is.na(hier)]
  int_terms <- terms[!is.na(int) & is.na(hier)]
  individual_terms <- terms[is.na(hier) & is.na(int)]

  if (response)
    response <- variables[response]
  else
    response <- NA

  list(individual_terms=individual_terms,
       int_terms=int_terms,
       group_terms=group_terms,
       response=extract_response(response),
       global_intercept=global_intercept)
}

extract_response <- function(response) {
  response_name_ch <- as.character(response[[1]])
  if ("cbind" %in% response_name_ch)
    ## remove cbind
    response_name_ch <- response_name_ch[-which(response_name_ch == "cbind")]
  else
    response_name_ch <- as.character(response)
  return(response_name_ch)
}

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

flatten_formula <- function(formula) {
  terms <- extract_terms_response(formula)
  group_terms <- terms$group_terms
  int_terms <- terms$int_terms
  individual_terms <- terms$individual_terms

  if (length(individual_terms) > 0 ||
      length(int_terms) > 0 ||
      length(group_terms) > 0)
    update(formula, paste(c(". ~ ",
                            flatten_individual_terms(individual_terms),
                            flatten_interaction_terms(int_terms),
                            flatten_group_terms(group_terms)),
                          collapse=" + "))
  else
    formula
}

flatten_individual_terms <- function(terms) {
  if (length(terms) < 1)
    return(terms)
  return(unique(terms))
}

flatten_interaction_terms <- function(terms) {
  if (length(terms) < 1)
    return(terms)
  ## TODO: do this right
  return(unique(terms))
}

flatten_group_terms <- function(terms) {
  if (length(terms) < 1)
    return(terms)
  split_terms <- strsplit(terms, "[ ]*\\|[ ]*")
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

break_formula <- function(formula, return_group_terms=TRUE) {
  terms <- extract_terms_response(formula)
  group_terms <- terms$group_terms
  int_terms <- terms$int_terms
  individual_terms <- terms$individual_terms
  global_intercept <- terms$global_intercept

  ## if there are no group level terms then `terms` is the most basic
  ## decomposition already
  if (length(group_terms) == 0 && length(int_terms) == 0)
    return(individual_terms)

  if (length(group_terms) > 0 && return_group_terms)
    warning(paste0("The model involves group terms.",
                   " Make sure that the provided mle and proj_predfun",
                   " handle multi-response projections."),
            call. = FALSE)

  if (return_group_terms)
    ## if there are group levels we should split that into basic components
    allterms <- c(individual_terms,
                  unlist(lapply(group_terms, split_group_term)),
                  unlist(lapply(int_terms, split_interaction_term)))
  else
    allterms <- c(individual_terms,
                  unlist(lapply(int_terms, split_interaction_term)))

  ## exclude the intercept if there is no intercept in the reference model
  if (!global_intercept) {
    allterms_nobias <- unlist(lapply(allterms, function(term) paste0(term, " + 0")))
    return(allterms_nobias)
  } else
    return(allterms)
}

split_interaction_term <- function(term) {
  ## strong heredity by default
  variables <- unlist(strsplit(term, ":"))
  individual_joint <- paste(variables, collapse = " + ")
  joint_intercept <- paste(c(individual_joint, term), collapse = " + ")

  joint_intercept
}

split_group_term <- function(term) {
  ## this expands whatever terms() did not expand
  term <- gsub("\\)$", "",
               gsub("^\\(", "",
                    flatten_group_terms(term)))
  if ("\\-" %in% term)
    stop(paste0("Use of `-` is not supported, omit variables or use the ",
                "method update on the formula, or write `0 +` to remove ",
                "the intercept."))

  chunks <- strsplit(term, " \\| ")[[1]]
  variables <- strsplit(chunks[1], " \\+ ")[[1]]

  ## split possible interaction terms
  int_t <- str_match(variables, ":")
  ## int_v <- lapply(variables[!is.na(int_t)], split_interaction_term)
  int_v <- variables[!is.na(int_t)]
  lin_v <- variables[is.na(int_t)]

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

formula_contains_group_terms <- function(formula) {
  t <- terms(formula)
  attributes <- attributes(t)
  terms <- attributes$term.labels

  hier <- str_match(terms, "\\|")
  group_terms <- terms[!is.na(hier)]

  return(length(group_terms) > 0)
}

subset_formula <- function(formula, terms, data, y=NULL, return_list=FALSE) {
  formula <- update(formula, paste0(". ~ ", paste(terms, collapse=" + ")))

  tt <- extract_terms_response(formula)
  response_name <- tt$response

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
    if (!return_list) {
      vector_y <- paste0("cbind(", paste(cols_y, collapse=", "), ")")
      formula <- update(formula, paste0(vector_y, " ~ ."))
    } else {
      formula <- lapply(cols_y, function(response) update(formula, paste0(response, " ~ .")))
    }
  } else
    cols_y <- response_name

  dat <- data.frame(y=y, X=X)
  colnames(dat) <- c(cols_y, cols_X)
  list(formula=formula, dat=dat)
}

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
      int <- str_match(current, ":")
      mulilevel <- str_match(current, "\\|")
      special <- str_match(current, "[a-z]+\\(([a-z]+)\\)")[, 1]
      individual <- is.na(mulilevel) & is.na(int)

      linear <- individual & is.na(special)
      linear_int <- !is.na(int) & is.na(special)

      if (linear) {
        ## if linear we can split it
        term_formula <- update(formula, paste0(". ~ ", current))
        split_term <- split_linear_term(current, term_formula, data)

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
          split_linear_term(v, update(formula, paste0(". ~ ", v)), data))

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

split_linear_term <- function(term, term_formula, data) {
  ## extract individual variables
  ## simple_mm <- model.matrix(term_formula, data=data)
  ## assign <- attr(simple_mm, "assign")
  ## appearances <- assign[assign != 0]

  appearances <- ncol(data[, term])

  if (appearances > 1)
    split_term <- sapply(1:appearances,
                         function(i) paste0(term, i))
  else
    split_term <- term
  split_term
}

count_terms_in_submodel <- function(submodel) {
  if (class(submodel) != "formula")
    submodel <- as.formula(paste0("~ ", submodel))
  tt <- extract_terms_response(submodel)
  ind_int_terms <- length(tt$individual_terms) + length(tt$int_terms)
  group_terms <- sum(unlist(lapply(tt$group_terms, count_terms_in_group_term)))
  if (is.null(group_terms))
    group_terms <- 0
  ind_int_terms + group_terms
}

count_terms_in_group_term <- function(term) {
  term <- gsub("[\\(\\)]", "", flatten_group_terms(term))

  chunks <- strsplit(term, " \\| ")[[1]]
  variables <- strsplit(chunks[1], " \\+ ")[[1]]

  ## split possible interaction terms
  int_t <- str_match(variables, ":")
  int_v <- variables[!is.na(int_t)]
  lin_v <- variables[is.na(int_t)]

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

sort_submodels_by_size <- function(submodels) {
  size <- lapply(submodels, count_terms_in_submodel)
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

count_variables_chosen <- function(refmodel, list_of_terms) {
  if (length(list_of_terms) == 0)
    return(0)
  form <- subset_formula(refmodel$formula, unique(unlist(list_of_terms)),
                         data=refmodel$fetch_data())$formula
  count_terms_in_submodel(flatten_formula(form))
}

is.redundant <- function(refmodel, current, new) {
  old_submodel <- current
  new_submodel <- c(current, new)
  if (count_variables_chosen(refmodel, new_submodel) >
      count_variables_chosen(refmodel, old_submodel))
    FALSE
  else
    TRUE
}

reduce_models <- function(refmodel, chosen) {
  Reduce(function(chosen, x)
    if (is.redundant(refmodel, chosen, x)) chosen else c(chosen, x),
    chosen)
}
