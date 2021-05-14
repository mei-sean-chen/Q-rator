#' Adjust map resolution
#'
#' Given an integer, adjusts map resolution and subsets input summary accordingly
#' @param interval Integer sets cM interval for map resolution
#' @param map Data frame with map data (group number and loci)
#' @param summary SNP or marker summary as data frame
#' @keywords
#' map
#' marker data
#' @export
#' @examples
#' map_to_SNPs()
summary_resolution_adjust <- function(interval, map, summary){
  
  map_frame <- map_format(map)
  
  filter <- rep(FALSE, length(map_frame$Marker))
  previous_loc <- as.numeric(0)
  counter <- 0.0
  
  for(i in 1:length(filter)){
    this_loc <- as.numeric(map_frame$Locus[i])
    if(isTRUE(previous_loc > this_loc)){
      previous_loc <- 0
      counter <- 0
    }
    if(previous_loc == 0){
      filter[i] <- T
    } else if(isTRUE(abs(this_loc - counter) > 1)){
      filter[i] <- F
    }
    else{
      filter[i] <- T
    }
    if(isTRUE(this_loc > counter)){
      counter <- counter + interval
      
    }
    previous_loc <- this_loc
  }
  
  for(i in 1:length(filter)){
    if(!filter[i] & (summary$Included[i] == 1)){
      summary$Included[i] <- 0
    }
  }
  return(summary)
}