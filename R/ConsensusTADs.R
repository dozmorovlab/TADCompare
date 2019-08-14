#' Consensus boundary identification
#'
#' @import dplyr
#' @import magrittr
#' @import PRIMME
#' @param cont_mats List of contact matrices in either sparse 3 column, n x n or n x (n+3)
#' form where the first three columns are coordinates in BED format.
#' If an x n matrix is used, the column names must correspond to the start
#' point of the corresponding bin. Required.
#' @param resolution Resolution of the data. Used to assign TAD boundaries
#' to genomic regions. If not provided, resolution will be estimated from
#' column names of the first matrix. Default is "auto"
#' @param z_thresh Threshold for boundary score. Higher values result in a
#' higher threshold for differential TADs. Default is 2.
#' @window_size Size of sliding window for TAD detection, measured in bins.
#' Results should be consistent Default is 15.
#' @param gap_thresh Required \% of 0s before a region will be considered a gap
#' and excluded. Default is .8
#' @return A list containing consensus TAD boundaries and overall scores
#'  \itemize{
#'  \item Consensus - Data frame containing location of all consensus
#'  boundaries
#'  \item All_Regions - Data frame containing consensus scores for all regions
#' }
#' @export
#' @details Given two sparse 3 column, n x n , or n x (n+3) contact matrices,
#' TADCompare identifies differential TAD boundaries. Using a novel boundary
#' score metric, TADCompare simultaneously identifies TAD boundaries and
#' tests for differential enrichment between datasets. Strength of interactions
#' are provided using raw boundary scores and p-values.
#' @examples
#' #Read in data
#' data("rao_chr22_rep")
#' data("rao_chr22_prim")
#' cont_mats = list(rao_chr22_rep, rao_chr22_prim)
#' #Find consensus TAD boundaries
#' diff_list <- TADCompare(cont_mats, resolution = 50000)

ConsensusTADs = function(cont_mats, resolution,
                      z_thresh = 2, window_size = 15, gap_thresh = .8) {

  #Calculate boundary scores
  bound_scores = lapply(seq_len(length(cont_mats)), function(x) {

    dist_sub = Single_Dist(cont_mats[[x]], resolution,
                           window_size = window_size,
                           gap_thresh = gap_thresh)
    dist_sub = data.frame(Sample = paste("Sample", x), dist_sub[,c(2:3)])
    dist_sub
  })

  coord_sum = lapply(bound_scores, function(x) x[,2])
  shared_cols = Reduce(intersect, coord_sum)
  bound_scores = lapply(bound_scores, function(x) x %>%
                          dplyr::filter(as.numeric(Coordinate) %in% as.numeric(shared_cols)))

  #Bind boundary scores
  score_frame = dplyr::bind_rows(bound_scores)

  colnames(score_frame)[1] = "Sample"
  base_sample = score_frame %>% filter(Sample == "Sample 1")

  #Get differences in boundary scores

  score_frame = score_frame %>% dplyr::group_by(Sample)  %>%
    dplyr::mutate(Diff_Score = (base_sample$Boundary-Boundary)) %>%
    dplyr::ungroup() %>% dplyr::mutate(Diff_Score = (Diff_Score-
                           mean(Diff_Score, na.rm = TRUE))/
                           sd(Diff_Score, na.rm =TRUE)) %>% dplyr::ungroup() %>% dplyr::mutate(
                         TAD_Score = (Boundary-mean(Boundary, na.rm = TRUE))/
                           sd(Boundary, na.rm = TRUE))

  score_frame = score_frame %>% dplyr::mutate(Differential = ifelse(abs(Diff_Score)>z_thresh,
                                       "Differential", "Non-Differential"),
                                       Differential = ifelse(is.na(Diff_Score),
                                                             "Non-Differential",
                                                             Differential))

  #Getting a frame summarizing boundaries

  TAD_Frame = score_frame %>% dplyr::select(Sample, Coordinate, TAD_Score) %>%
    dplyr::arrange(as.numeric(gsub("Sample", "", Sample))) %>%
    dplyr::mutate(Sample = factor(Sample, levels = unique(Sample)))


  #Spread into wide format
  TAD_Frame = tidyr::spread(as.data.frame(TAD_Frame),
                            key = Sample,
                            value = TAD_Score)

  #Get median of boundary scores for each row

  TAD_Frame = TAD_Frame %>%
    dplyr::mutate(Consensus_Score =
    apply(TAD_Frame %>% dplyr::select(-Coordinate) %>% as.matrix(.)
          ,1, median))

  return(list(Consensus = TAD_Frame,
              All_Regions = score_frame))
}


