#' Time-varying TAD boundary analysis
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
#' @param groupings Vector identifying which group each contact matrix belongs
#' to. Each group will be combined using consensus boundary scores. Default is
#' NULL
#' @return A list containing consensus TAD boundaries and overall scores
#'  \itemize{
#'  \item Consensus - Data frame containing location of all consensus
#'  boundaries
#'  \item All_Regions - Data frame containing consensus scores for all regions
#' }
#' @export
#' @details Given a list of sparse 3 column, n x n, or n x (n+3) contact
#' matrices, TimeCompare identifies TADs and classifies time points based on
#' boundary change. A data frame of time points with at least one boundary
#' and the corresponding classification is returned.
#' @examples
#' #Read in data
#' data("rao_chr22_rep")
#' data("rao_chr22_prim")
#' cont_mats = list(rao_chr22_rep, rao_chr22_prim)
#' #Find time varying TAD boundaries
#' diff_list <- TimeCompare(cont_mats, resolution = 50000)

TimeCompare = function(cont_mats, resolution,
                      z_thresh = 2, window_size = 15, gap_thresh = .8,
                      groupings = NULL) {

  #Calculate boundary scores
  bound_scores = lapply(seq_len(length(cont_mats)), function(x) {

    dist_sub = Single_Dist(cont_mats[[x]], resolution, window_size = window_size)
    dist_sub = data.frame(Sample = paste("Sample", x), dist_sub[,c(2:3)])
    dist_sub
  })

  coord_sum = lapply(bound_scores, function(x) x[,2])
  shared_cols = Reduce(intersect, coord_sum)
  bound_scores = lapply(bound_scores, function(x) x %>%
                          filter(as.numeric(Coordinate) %in% as.numeric(shared_cols)))

  #Bind boundary scores
  score_frame = bind_rows(bound_scores)

  colnames(score_frame)[1] = "Sample"
  base_sample = score_frame %>% filter(Sample == "Sample 1")

  if (!is.null(groupings)) {
    #Map groupings to samples
    Group_Frame = data.frame(Groups = groupings,
                             Sample = unique(score_frame$Sample))

    #Join to replace sample with group
    score_frame = left_join(score_frame, Group_Frame) %>%
      dplyr::select(Sample = Groups, Coordinate, Boundary)

    score_frame = score_frame %>% group_by(Sample, Coordinate) %>%
      mutate(Boundary = median(Boundary)) %>% distinct()
  }

  #Get differences

  score_frame = score_frame %>% group_by(Sample)  %>%
    mutate(Diff_Score = (base_sample$Boundary-Boundary)) %>%
    ungroup() %>% mutate(Diff_Score = (Diff_Score-
                           mean(Diff_Score, na.rm = TRUE))/
                           sd(Diff_Score, na.rm =TRUE)) %>% ungroup() %>% mutate(
                         TAD_Score = (Boundary-mean(Boundary, na.rm = TRUE))/
                           sd(Boundary, na.rm = TRUE))

  score_frame = score_frame %>% mutate(Differential = ifelse(abs(Diff_Score)>z_thresh,
                                       "Differential", "Non-Differential"),
                                       Differential = ifelse(is.na(Diff_Score),
                                                             "Non-Differential",
                                                             Differential))

  #Getting a frame summarizing boundaries

  TAD_Frame = score_frame %>% dplyr::select(Sample, Coordinate, TAD_Score) %>%
    arrange(as.numeric(gsub("Sample", "", Sample))) %>%
    mutate(Sample = factor(Sample, levels = unique(Sample)))


  #Spread into wide format
  TAD_Frame = tidyr::spread(as.data.frame(TAD_Frame),
                            key = Sample,
                            value = TAD_Score)

  #Get median of boundary scores for each row

  TAD_Frame = TAD_Frame %>%
    mutate(Consensus_Score =
    apply(TAD_Frame %>% dplyr::select(-Coordinate) %>% as.matrix(.)
          ,1, median))

  Differential_Points = score_frame %>% filter(Differential == "Differential")

  #Pulling out consensus for classification

  TAD_Iden = TAD_Frame[,c(-1, -ncol(TAD_Frame))]>3

  #Classify time trends

  #All_TADs = apply(TAD_Iden, 1, function(x) all(x == TRUE))
  All_Non_TADs = apply(TAD_Iden, 1, function(x) all(x == FALSE))

  #Split into 4 groups for summarization
  Num_Points = seq_len(ncol(TAD_Iden))

  four_groups = split(Num_Points, ggplot2::cut_number(Num_Points,4))

  #Get summary for each group and put back together

  Full_Sum = do.call(cbind.data.frame,
          lapply(four_groups, function(x) (rowSums(as.matrix(TAD_Iden[,x]))/
    ncol(as.matrix(TAD_Iden[,x])))>=.5))

  #Highly Common TADs
  Common_TADs = (Full_Sum[,1] == Full_Sum[,2]) &
     (Full_Sum[,1] == Full_Sum[,ncol(Full_Sum)])

  Common_TADs = ifelse(Common_TADs, "Highly Common TAD", NA)

  #Early Appearing TADs
  Early_Appearing = (Full_Sum[,1] != Full_Sum[,2]) &
    (Full_Sum[,1] ==FALSE) &
    (Full_Sum[,2] == Full_Sum[,ncol(Full_Sum)])

  Early_Appearing = ifelse(Early_Appearing, "Early Appearing TAD", NA)

  #Late Appearing TADs
  Late_Appearing = (Full_Sum[,1] == Full_Sum[,2]) &
    (Full_Sum[,1] ==FALSE) &
    (Full_Sum[,1] != Full_Sum[,ncol(Full_Sum)])

  Late_Appearing = ifelse(Late_Appearing, "Late Appearing TAD", NA)

  #Early disappearing TADs
  Early_Disappearing = (Full_Sum[,1] != Full_Sum[,2]) &
    (Full_Sum[,1] == TRUE) &
    (Full_Sum[,2] == Full_Sum[,ncol(Full_Sum)])

  Early_Disappearing = ifelse(Early_Disappearing, "Early Disappearing TAD", NA)

  #Late disapearing TADs
  Late_Disappearing = (Full_Sum[,1] == Full_Sum[,2]) &
    (Full_Sum[,1] == TRUE) &
    (Full_Sum[,1] != Full_Sum[,ncol(Full_Sum)])

  Late_Disappearing = ifelse(Late_Disappearing, "Late Disappearing TAD", NA)

  #Dynamic TAD
  Dynamic = (Full_Sum[,1] != Full_Sum[,2]) &
    (Full_Sum[,1] == Full_Sum[,ncol(Full_Sum)])

  Dynamic = ifelse(Dynamic, "Dynamic TAD", NA)

  TAD_Cat = dplyr::coalesce(as.character(Common_TADs),
                     as.character(Early_Appearing),
                     as.character(Late_Appearing),
                     as.character(Early_Disappearing),
                     as.character(Late_Disappearing),
                     as.character(Dynamic))

  TAD_Frame = TAD_Frame %>% dplyr::mutate(Category = TAD_Cat)

  return(list(Differential_Points = Differential_Points,
              Full = score_frame,
              Consensus = TAD_Frame))
}


