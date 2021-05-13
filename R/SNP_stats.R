#' Assemble marker set from map data
#'
#' Helper function; returns vector of all SNPs/markers provided a map file as a data frame
#' @param map Data frame with map data (group number and loci)
#' @keywords
#' map
#' marker data
#' @export
#' @examples
#' map_to_SNPs()
map_to_SNPs <- function(map){
  SNPs <- map[,1]
  SNPs <- SNPs[!is.na(SNPs)]
}

#' Assemble set of parent-homozygous markers
#'
#' Helper function; returns set (as vector) of all parent-homozygous SNPs/markers for
#' removal from curated output frame. Assumes input summary is correct.
#' @param summary SNP or marker summary as data frame
#' @keywords
#' summary
#' marker data
#' @export
#' @examples
#' to_extract_parent_homozygous()
to_extract_parent_homozygous <- function(summary){

  extract_these <- c()
  for(x in 1:length(summary$`A A`)){
    thisRow <- summary[x,]
    if(
      (thisRow[2] + thisRow[8] == 1) |
      (thisRow[4] + thisRow[8] == 1)
    ){
      extract_these <- append(extract_these,
                                 as.character(summary$Markers[x]))
    }
  }
  return(extract_these)
}


#' Assemble parent marker data
#'
#' Creates a data frame subset of the curated output .dat frame, which
#' contains only parents of the input set, and corresponding marker data
#' @param phenotype_data Data frame with phenotype data on all input individuals
#' @param all_data Marker data set as data frame
#' @param start First column of data in phenotype data frame
#' @param end Last column of data in phenotype data frame
#' @keywords
#' summary
#' marker data
#' @export
#' @examples
#' first_parent_SNP_data()
first_parent_SNP_data <- function(phenotype_data, all_data, start, end){
  out <- format_out(phenotype_data, all_data, start, end)

  filtered.all.ind <- input_align_parents(phenotype_data, all_data)
  FirstParents <- get_parents(filtered.all.ind)

  outIndices <- out$X.2
  filter <- rep(F, length(outIndices))
  for(g in 1:length(FirstParents)){
    if(!(FirstParents[g] == "-") & !(FirstParents[g] == "")
       & FirstParents[g] %in% outIndices){
      filter <- (filter | out$X.2 == FirstParents[g])
    }
  }
  filter[1:2] <- c(T,T)
  firstParentSNPData <- out[filter,]
  firstParentSNPData <- firstParentSNPData[, colSums(is.na(firstParentSNPData)) == 0]
  return(firstParentSNPData)
}
#' Assemble parent-homozygous marker data
#'
#' Creates a data frame subset of the curated output .dat frame, which
#' contains only parents of the input set, and only homozygous markers
#' @param phenotype_data Data frame with phenotype data on all input individuals
#' @param all_data Marker data set as data frame
#' @param start First column of data in phenotype data frame
#' @param end Last column of data in phenotype data frame
#' @keywords
#' marker data
#' pedigree
#' @export
#' @examples
#' parent_SNPs_culled()
parent_SNPs_culled <- function(phenotype_data, all_data, map, start, end){
  firstParentSNPData <- first_parent_SNP_data(phenotype_data, all_data, start, end)
  summary <- parent_SNP_summary(phenotype_data, all_data, map, start, end)
  extract_these <- to_extract_parent_homozygous(summary)
  filter <- c(F, rep(T, 3), rep(F, length(firstParentSNPData[2,]) - 4))
  for(x in 12:length(filter)){
    filter[x] <- (firstParentSNPData[1,x] %in% extract_these)
  }
  firstParentSNPData <- firstParentSNPData[, filter]
  return(firstParentSNPData)
}

#' Parent Marker Summary
#'
#' Creates a data frame that summarizes SNP frequencies for the set of all
#' parents of the input individual set provided by the phenotype data.
#' @param phenotype_data Data frame with phenotype data on all input individuals
#' @param all_data Marker data set as data frame
#' @param start First column of data in phenotype data frame
#' @param end Last column of data in phenotype data frame
#' @keywords
#' marker data
#' allele frequencies
#' summary
#' pedigree
#' @export
#' @examples
#' parent_SNPs_summary()
parent_SNP_summary <- function(phenotype_data, all_data, map, start, end){

  firstParentSNPData <- first_parent_SNP_data(phenotype_data, all_data, start, end)

  Markers <- map_to_SNPs(map)

  toMatrix <- rep(0, 7*length(Markers))
  summary <- matrix(toMatrix, nrow = length(Markers), ncol = 7)
  summary <- data.frame(summary)
  colnames(summary) <- c("A A", "A B", "B B", "A C", "B C", "C C", "- -")

  for(s in 1:length(Markers)){
    thisRow <- rep(0, 7)
    thisColumn <- firstParentSNPData[, s+11]
    for(t in 1:length(thisColumn)){
      switch(as.character(thisColumn[t]),
             "A A" = {
               thisRow[1] <- thisRow[1] + 1
             },
             "A B" = {
               thisRow[2] <- thisRow[2] + 1
             },
             "B B" = {
               thisRow[3] <- thisRow[3] + 1
             },
             "A C" = {
               thisRow[4] <- thisRow[4] + 1
             },
             "B C" = {
               thisRow[5] <- thisRow[5] + 1
             },
             "C C" = {
               thisRow[6] <- thisRow[6] + 1
             },
             "- -" = {
               thisRow[7] <- thisRow[7] + 1
             },
             {
               thisRow[7] <- thisRow[7] + 1
             }
      )
    }
    for(r in 1:7){
      thisRow[r] <- thisRow[r] / as.double(length(thisColumn))
    }
    summary[s,] <- thisRow
  }
  summary <- summary %>% add_column(Markers, .before = "A A")

  return(summary)
}

#' Input Set Marker Summary
#'
#' Creates a data frame that summarizes SNP frequencies for the set of all
#' input individuals provided by the phenotype data.
#' @param phenotype_data Data frame with phenotype data on all input individuals
#' @param all_data Marker data set as data frame
#' @param start First column of data in phenotype data frame
#' @param end Last column of data in phenotype data frame
#' @keywords
#' marker data
#' allele frequencies
#' summary
#' @export
#' @examples
#' input_SNP_summary()
input_SNP_summary <- function(phenotype_data, all_data, map, start, end){
  out <- format_out(phenotype_data, all_data, start, end)
  Markers <- map_to_SNPs(map)

  map <- map_format(map)

  inputIndices <- get_individuals(phenotype_data)

  filter <- rep(F, length(out$X.1))
  filter[1:2] <- c(T,T)
  for(i in 3:length(filter)){
    filter[i] <- out$X.2[i] %in% inputIndices$Index
  }
  inputIndicesData <- droplevels(out[filter,])

  toMatrix <- rep(0, 7*length(Markers))
  input.SNPs.frame <- matrix(toMatrix, nrow = length(Markers), ncol = 7)
  input.SNPs.frame <- data.frame(input.SNPs.frame)
  colnames(input.SNPs.frame) <- c("A A", "A B", "B B", "A C", "B C", "C C", "- -")
  for(s in 1:length(Markers)){
    thisRow <- rep(0, 7)
    thisColumn <- inputIndicesData[, s+11]
    for(t in 1:length(thisColumn)){
      switch(as.character(thisColumn[t]),
             "A A" = {
               thisRow[1] <- thisRow[1] + 1
             },
             "A B" = {
               thisRow[2] <- thisRow[2] + 1
             },
             "B B" = {
               thisRow[3] <- thisRow[3] + 1
             },
             "A C" = {
               thisRow[4] <- thisRow[4] + 1
             },
             "B C" = {
               thisRow[5] <- thisRow[5] + 1
             },
             "C C" = {
               thisRow[6] <- thisRow[6] + 1
             },
             "- -" = {
               thisRow[7] <- thisRow[7] + 1
             },
             {
               thisRow[7] <- thisRow[7] + 1
             }
      )
    }
    for(r in 1:7){
      thisRow[r] <- thisRow[r] / as.double(length(thisColumn))
    }
    input.SNPs.frame[s,] <- thisRow
  }
  input.SNPs.frame <- input.SNPs.frame %>% add_column(Markers, .before = "A A")
  input.SNPs.frame <- input.SNPs.frame %>% add_column(map$Locus, .before = "A A")
  input.SNPs.frame <- input.SNPs.frame %>% add_column(map$Group, .before = "A A")
  colnames(input.SNPs.frame)[2:3] <- c("Locus", "Group")
  filter <- rep(T, length(input.SNPs.frame$Markers))

  Included <- rep(1, length(input.SNPs.frame$Markers))
  input.SNPs.frame <- input.SNPs.frame %>% add_column(Included)

  return(input.SNPs.frame)
}

#' Cull by locus
#'
#' Creates a soft subset of the marker summary frame based on marker location.
#' Removes redundant markers; defaults to keeping first marker at location. Returns
#' summary frame with column "Included" modified.
#' @param summary SNP or marker summary as data frame
#' @keywords
#' marker data
#' loci
#' summary
#' @export
#' @examples
#' cull_by_locus()
cull_by_locus <- function(summary){
  filter <- rep(T, length(summary$Markers))
  for(i in 2:length(summary$Markers)){
    if(summary$Group[i] == summary$Group[i - 1]){
      if(summary$Locus[i] == summary$Locus[i - 1]){
        filter[i] <- F
      }
    }
  }

  for(i in 1:length(filter)){
    if(filter[i]){
      summary$Included[i] <- 1
    }
    else {
      summary$Included[i] <- 0
    }
  }
  return(summary)
}

#' Assemble locus-culled set
#'
#' Helper function; creates a set of markers as a vector excluded by the logic of \code{cull_by_locus()}
#' @param map Marker map as data frame
#' @keywords
#' marker data
#' locus
#' loci
#' map
#' @export
#' @examples
#' to_extract_by_locus()
to_extract_by_locus <- function(map){
  map <- map_format(map)
  filter <- rep(FALSE, length(map$Marker))
  for(i in 2:length(map$Marker)){
    if(map$Group[i] == map$Group[i - 1]){
      if(map$Locus[i] == map$Locus[i - 1]){
        filter[i] <- TRUE
      }
    }
  }

  culled <- map$Marker[filter]
  return(culled)
}

#' Summary-included Markers
#'
#' Helper function; produces a true subset of the input summary based on Included
#' column, values of TRUE or 1 make the subset.
#' @param summary SNP or marker summary as data frame
#' @keywords
#' map
#' summary
#' @export
#' @examples
#' included_markers()
included_markers <- function(summary){
  filter <- rep(F, length(summary$Markers))
  for(i in 1:length(summary$Markers)){
    if(summary$Included[i]){
      filter[i] <- TRUE
    }
  }
  included <- summary[filter,]
  return(included)
}

#' Summary-excluded Markers
#'
#' Helper function; produces a true subset of the input summary based on Included
#' column, values of FALSE or 0 make the subset.
#' @param summary SNP or marker summary as data frame
#' @keywords
#' map
#' summary
#' @export
#' @examples
#' excluded_markers()
excluded_markers <- function(summary){
  filter <- rep(F, length(summary$Markers))
  for(i in 1:length(summary$Markers)){
    if(!summary$Included[i]){
      filter[i] <- TRUE
    }
  }
  excluded <- summary[filter,]
  return(excluded)
}

#' Summary Set-Subtraction
#'
#' Soft subset of SNPs/markers included in the input summary frame by
#' set subtraction, summary - markers. Markers must be a vector.
#' @param summary SNP or marker summary as data frame
#' @param markers Vector of marker IDs to exclude from returned output
#' @keywords
#' map
#' summary
#' marker data
#' @export
#' @examples
#' subtraction_input_summary(summary = your_summary_table, markers = set_to_subtract)
subtraction_input_summary <- function(summary, markers){
  for(i in 1:length(summary$Markers)){
    if(summary$Markers[i] %in% markers){
      summary$Included[i] <- 0
    }
  }
  return(summary)
}


#' Marker Count By Group
#'
#' Obtains a vector of marker counts by Group number, useful to .par file.
#' Group number corresponds to vector index.
#' @param summary SNP or marker summary as data frame
#' @keywords
#' map
#' summary
#' .par constituent
#' @export
#' @examples
#' marker_count()
marker_count <- function(summary){
  groups <- summary$Group[length(summary$Markers)]
  groups <- as.numeric(as.character(groups))
  SNP.count <- rep(0, groups)
  for(i in 1:length(summary$Group)){
    thisIndex <- as.numeric(as.character(summary$Group[i]))
    if(summary$Included[i] == 1){
      SNP.count[thisIndex] <- 1 + SNP.count[thisIndex]
    }
  }
  return(SNP.count)
}

#' Subset Summary By Input
#'
#' Updates input summary based on user-input string "add" or "remove"
#' and vector containing markers.
#' Best paired with functions which update dependents of the
#' input summary.
#' Requires another function to update the .dat data frame by input summary.
#' @param summary SNP or marker summary as data frame
#' @param map Marker map as data frame
#' @param add_remove String, "add" or "remove"
#' @param markers Vector of marker IDs to add or remove to input summary, Markers must be present in map set
#' @keywords
#' map
#' summary
#' .par constituent
#' markers
#' @export
#' @examples
#' subset_marker_data()
subset_marker_data <- function(summary, map, add_remove, markers){

  SNPs <- map_to_SNPs(map)


  for(SNP in markers){
    if(add_remove == "add"){
      if (SNP %in% SNPs){
        x <- 1
        while(SNP != SNPs[x]){
          x <- x + 1
        }
        if (summary$Included[x] == 1){
          print(paste("SNP is already present:", SNP, sep=" "))
        }
        else {
          #Add SNP back in, must update relevantSNPs
          summary$Included[x] <- 1
        }
      }
      else {
        print(paste("Invalid SNP:", SNP, sep=" "))
      }
    } else if(add_remove == "remove"){
      if (SNP %in% SNPs){
        x <- 1
        while(SNP != SNPs[x]){
          x <- x + 1
        }
        if (summary$Included[x] == 1){
          #Remove the SNP
          summary$Included[x] <- 0
        }
        else {
          print(paste("SNP is not present:", SNP, sep=" "))
        }
      }
      else {
        print(paste("Invalid SNP:", SNP, sep=" "))
      }
    } else {
      print("Invalid Input")
    }
    return(summary)
  }
}





















