#' Format Out
#'
#' Given phenotype and marker data, creates a FlexQTL-ready data frame.
#' @param phenotype_data Data frame with phenotype data on all input individuals
#' @param all_data Marker data set
#' @keywords out
#' .dat
#' @export
#' @examples
#' format_out(phenotype_data = your_frame, all_data = 20k_8koverlap)
format_out <- function(phenotype_data, all_data){
  
  start <- 2
  end <- ncol(phenotype_data)
  
  total <- fill_pedigree(phenotype_data, all_data)
  N <- length(all_data$X.2)
  filter <- rep(F, N)

  for(x in 1:length(total$Index)){
    filter <- (filter | all_data$X.2 == total$Index[x])
  }
  filter[3:4] <- c(T,T)
  out <- all_data[filter,]
  rm(N, x)

  out$X <- NULL
  out$X.1 <- NULL
  out$Sort <- NULL
  columnA <- append(c("",";"), rep("1", length(out$X.2)-2))
  out <- out %>% add_column(X.1 = columnA, .before = "X.2")
  for(j in start:end){
    if(j == end){
      phen.insert <- c(";", colnames(phenotype_data)[j], rep("-",length(out$X.2)-2))
    } else {
      phen.insert <- c("", colnames(phenotype_data)[j], rep("-",length(out$X.2)-2))
    }
    for(i in 3:length(out$X.2)){
      if (out$X.2[i] %in% phenotype_data$Index){
        phen.insert[i] <- as.character(phenotype_data[as.character(phenotype_data$Index)
                                                 == as.character(out$X.2[i]), j])
      }
    }
    out <- out %>% add_column(phen.insert, .before = "X1")
  }
  out <- out[, colSums(is.na(out)) == 0]
  return(out)
}

#' Subset Data File
#'
#' Helper function; obtains subset of 'out' based on a provided set of marker IDs by set subtraction
#' @param out Data frame with phenotype data on all input individuals
#' @param markers Vector, set of marker IDs to subtract from .dat
#' @keywords out
#' .dat
#' @export
#' @examples
#' subset_out_markers(out = my_dat, markers = exclude_these)
subset_out_markers <- function(out, markers){
  filter <- rep(T, length(out))
  for(x in 1:length(out)){
    if(as.character(out[1, x]) %in% markers){
      filter[x] <- F
    }
  }
  return(out[,filter])
}

#' Update Data File By Summary
#'
#' Produces a subset of all marker data based on individuals provided by the summary set.
#' Output may be a subset or supserset of the existing .dat data frame depending on the summary.
#' @param phenotype_data Data frame with phenotype data on all input individuals
#' @param all_data Marker data set as data frame
#' @param summary SNP or marker summary as data frame
#' @keywords out
#' .dat
#' @export
#' @examples
#' update_out_by_summary()
update_out_by_summary <- function(phenotype_data, all_data, summary){
  
  negative_set <- excluded_markers(summary)$Markers
  out <- format_out(phenotype_data, all_data)
  filter <- rep(T, length(out))
  for(x in 1:length(out)){
    if(as.character(out[1, x]) %in% negative_set){
      filter[x] <- F
    }
  }
  return(out[,filter])
}

#' Get Parents
#'
#' Given a properly-formatted data frame, with individuals in column 1, parents in columns
#' 2 and 3, this function returns vector containing the set of all parents
#' @param indiv_and_parents Data frame, input individuals set and parents of each individual
#' @keywords
#' pedigree
#' .dat
#' @export
#' @examples
#' get_parents()
get_parents <- function(indiv_and_parents){
  parents <- c(levels(factor(indiv_and_parents$Parent1)),
               levels(factor(indiv_and_parents$Parent2)))
  return(parents)
}


#' Fill Pedigree
#'
#' Given a phenotype and marker data, assembles a data frame with the pedigree of each input
#' individual represented in column 1. Each individual has their parents represented in columns
#' 2 and 3.
#' @param phenotype_data Data frame with phenotype data on all input individuals
#' @param all_data Marker data set as data frame
#' @keywords
#' pedigree
#' .dat
#' @export
#' @examples
#' fill_pedigree()
fill_pedigree <- function(phenotype_data, all_data){

  total <- input_align_parents(phenotype_data, all_data)
  parents <- get_parents(total)
  all.ind <- to_all_indiv(all_data)

  while(length(parents) > 1){
    filter <- rep(F, length(all.ind$Index))
    for(x in 1:length(parents)){
      if(!(parents[x] %in% total$Index) & !(parents[x] == '-')){
        filter <- (filter | all.ind$Index == parents[x])
      }
    }
    temp <- all.ind[filter,]
    parents <- c(levels(factor(temp$Parent1)), levels(factor(temp$Parent2)))
    total <- rbind(total, temp)
  }
  return(total)
}

#' Align Input Indices and Parents
#'
#' Given a phenotype and marker data, assembles a data frame with every input individual
#' in column 1, and parents 1 and 2 in columns 2 and 3. Erroneous individual IDs/indices are
#' not included.
#' @param phenotype_data Data frame with phenotype data on all input individuals
#' @param all_data Marker data set as data frame
#' @keywords
#' phenotype data
#' marker data
#' formatting
#' .dat
#' @export
#' @examples
#' input_align_parents()
input_align_parents <- function(phenotype_data, all_data){
  all_indiv <- to_all_indiv(all_data)
  input_indiv <- makeD_cleaned_input_indiv(phenotype_data, all_data)

  N <- length(all_indiv$Index)
  filter <- rep(F, N)
  for(x in 1:N){
    filter[x] <- all_indiv$Index[x] %in% input_indiv$Index
  }
  return(all_indiv[filter,])
}
