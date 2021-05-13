#' Format Map Helper Frame
#'
#' Produces a helper frame with SNP/marker location and group number paired to
#' each marker ID in a 3-column table
#' @param map Marker map as data frame
#' @keywords
#' map
#' summary
#' markers
#' @export
#' @examples
#' map_format()
map_format <- function(map){
  x <- 1
  to_map <- rep("0", 3*length(map$X))
  to_map <- matrix(to_map, nrow = length(to_map) / 3, ncol = 3)
  colnames(to_map) <- c("Marker", "Locus", "Group")
  for(i in 1:length(map$X)){
    if(!is.na(map$X[i])){
      to_map[i,] <- c(map$X[i],
                   as.character(map$Group.1[i]),
                   x)
    }
    else {
      x <- x + 1
    }
  }
  rm(x, i)
  to_map <- data.frame(to_map)
  filter <- rep(T, length(to_map$Marker))
  for(i in 1:length(to_map$Marker)){
    filter[i] <- (to_map[i,] != c(0, 0, 0))[1]
  }
  to_map <- to_map[filter,]
  return(to_map)
}

#' Update .map
#'
#' Returns data frame version of .map for use in
#' in FQTL, obtains image of current SNP set with group and locus data provided existing,
#' complete .map and input set summary. Assumes correctness of summary. Can produce default map is summary
#' parameter is not subset.
#' @param map Marker map as data frame
#' @param summary SNP or marker summary as data frame
#' @keywords
#' map
#' summary
#' markers
#' loci
#' @export
#' @examples
#' map_format()
map_final <- function(summary, map){
  exclusions_as_vector <- excluded_markers(summary)$Markers

  filter <- rep(T,length(map$X))
  for(i in 1:length(filter)){
    filter[i] <- !((as.character(map$X[i]) %in% exclusions_as_vector))
  }
  map_subset_update <- map[filter,]
  colnames(map_subset_update) <- c("", "Group 1")
  return(map_subset_update)
}


