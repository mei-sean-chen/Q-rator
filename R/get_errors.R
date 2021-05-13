#' Get (the set of) Individuals
#'
#' Produces a vector containing the set of all input individuals provided
#' the phenotype data in data frame form.
#' @param phen_data Data frame with phenotype data on all input individuals
#' @keywords
#' input set
#' phenotype data
#' @export
#' @examples
#' get_individuals()
get_individuals <- function(phen_data){
  individuals <- phen_data[,1,drop=F]
}

#' Clean marker data frame
#'
#' Helper function; clears input marker data frame of all cells containing value "NA"
#' @param all_data Marker data set as data frame
#' @keywords
#' input set
#' marker data
#' @export
#' @examples
#' clean_20k8k()
clean_20k8k <- function(all_data){
  return(all_data[, colSums(is.na(all_data)) == 0])
}

#' Get masterlist of indices and parents
#'
#' Helper function; from input marker data frame, obtains subset containing all possible individuals,
#' parent 1 and parent 2 for each.
#' @param all_data Marker data set as data frame
#' @keywords
#' input set
#' marker data
#' @export
#' @examples
#' to_all_indiv()
to_all_indiv <- function(all_data){
  all_data <- clean_20k8k(all_data)
  all_indiv <- all_data[c(5:length(all_data$X.2)),c(3:5), drop=F]
  colnames(all_indiv) <- c("Index", "Parent1", "Parent2")
  return(all_indiv)
}

#' Filter Clerical Errors In Input Set
#'
#' Cross references input individual set (phenotype data) and
#' pool of individuals which have on-file marker data.
#' Produces a 1-column data frame containing all input individuals not found
#' in the master marker data set.
#' @param phen_data Data frame with phenotype data on all input individuals
#' @param all_data Marker data set as data frame
#' @keywords
#' input set
#' marker data
#' phenotype data
#' errors
#' @export
#' @examples
#' makeD_get_errors()
makeD_get_errors <- function(phenotype_data, all_data){
  individuals <- get_individuals(phenotype_data)
  all_indiv <- to_all_indiv(all_data)

  N <- length(individuals$Index)
  error_filter <- rep(F, N)
  for(x in 1:N){
    error_filter[x] <- (individuals$Index[x] %in% all_indiv$Index)
  }
  errors <-       individuals[!error_filter,,drop=F]
  return(errors)
}

#' Get "Clean" Input-Individual Set
#'
#' Using logic of function \code{makeD_get_errors()} obtains the set of
#' all input-individuals whose IDs exist in the master marker data set.
#' @param phen_data Data frame with phenotype data on all input individuals
#' @param all_data Marker data set as data frame
#' @keywords
#' input set
#' marker data
#' @export
#' @examples
#' makeD_cleaned_input_indiv()
makeD_cleaned_input_indiv <- function(phenotype_data, all_data){
  individuals <- get_individuals(phenotype_data)
  all_indiv <- to_all_indiv(all_data)

  N <- length(individuals$Index)
  error_filter <- rep(F, N)
  for(x in 1:N){
    error_filter[x] <- (individuals$Index[x] %in% all_indiv$Index)
  }
  individuals <-  individuals[error_filter,,drop=F]
  return(individuals)
}
