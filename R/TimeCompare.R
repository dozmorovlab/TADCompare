#' Time-varying TAD boundary analysis
#'
#' @import dplyr
#' @import magrittr
#' @import PRIMME
#' @importFrom HiCcompare sparse2full
#' @param cont_mats List of contact matrices in either sparse 3 column, n x n or n x (n+3)
#' form where the first three columns are coordinates in BED format.
#' If an x n matrix is used, the column names must correspond to the start
#' point of the corresponding bin. Required.
#' @param resolution Resolution of the data. Used to assign TAD boundaries
#' to genomic regions. If not provided, resolution will be estimated from
#' column names of the first matrix. Default is "auto"
#' @param z_thresh Threshold for boundary score. Higher values result in a
#' higher threshold for differential TADs. Default is 3.
#' @window_size Size of sliding window for TAD detection, measured in bins.
#' Results should be consistent Default is 15.
#' @param gap_thresh Required \% of 0s before a region will be considered a gap
#' and excluded. Default is .8
#' @param groupings Vector identifying which group each contact matrix belongs
#' to. Each group will be combined using consensus boundary scores. Default is
#' NULL
#' @return A list containing consensus TAD boundaries and overall scores
#'  \itemize{
#'  \item Time_TADs - Data frame containing all regions with a TAD boundary
#'  at one or more time point
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
                      z_thresh = 3, window_size = 15, gap_thresh = .8,
                      groupings = NULL) {

  #Get dimensions of first contact matrix
  row_test = dim(cont_mats[[1]])[1]
  col_test = dim(cont_mats[[1]])[2]

  if (row_test == col_test) {
    if (all(is.finite(cont_mats[[1]])) == FALSE) {
      stop("Contact matrix 1 contains non-numeric entries")
    }

  }

  if (col_test == 3) {


    #Convert sparse matrix to n x n matrix

    message("Converting to n x n matrix")

    cont_mats = lapply(cont_mats, HiCcompare::sparse2full)

    if (all(is.finite(cont_mats[[1]])) == FALSE) {
      stop("Contact matrix 1 contains non-numeric entries")
    }

    if (resolution == "auto") {
      message("Estimating resolution")
      resolution = as.numeric(names(table(as.numeric(colnames(cont_mats[[1]]))-
                                            lag(
                                              as.numeric(
                                                colnames(
                                                  cont_mats[[1]])))))[1])
    }

  } else if (col_test-row_test == 3) {

    message("Converting to n x n matrix")

    cont_mats = lapply(cont_mats, function(x) {
      #Find the start coordinates based on the second column of the bed file portion of matrix

      start_coords = x[,2]

      #Calculate resolution based on given bin size in bed file

      resolution = as.numeric(x[1,3])-as.numeric(x[1,2])

      #Remove bed file portion

      x = as.matrix(x[,-c(seq_len(3))])

      if (all(is.finite(x)) == FALSE) {
        stop("Contact matrix contains non-numeric entries")
      }


      #Make column names correspond to bin start

      colnames(x) = start_coords
      return(x)
    })

  } else if (col_test!=3 & (row_test != col_test) & (col_test-row_test != 3)) {

    #Throw error if matrix does not correspond to known matrix type

    stop("Contact matrix must be sparse or n x n or n x (n+3)!")

  } else if ( (resolution == "auto") & (col_test-row_test == 0) ) {

    message("Estimating resolution")

    #Estimating resolution based on most common distance between loci

    resolution = as.numeric(names(table(as.numeric(colnames(cont_mats[[1]]))-
                                          lag(
                                            as.numeric(
                                              colnames(
                                                cont_mats[[1]])))))[1])
  }

  #Calculate boundary scores
  bound_scores = lapply(seq_len(length(cont_mats)), function(x) {

    dist_sub = Single_Dist(cont_mats[[x]], resolution, window_size = window_size)
    dist_sub = data.frame(Sample = paste("Sample", x), dist_sub[,c(2:3)])
    dist_sub
  })


  #Reduce matrices to only include shared regions
  coord_sum = lapply(bound_scores, function(x) x[,2])
  shared_cols = Reduce(intersect, coord_sum)
  bound_scores = lapply(bound_scores, function(x) x %>%
                          filter(as.numeric(Coordinate) %in% as.numeric(shared_cols)))

  #Bind boundary scores
  score_frame = bind_rows(bound_scores)

  #Set column names for base sample
  colnames(score_frame)[1] = "Sample"
  base_sample = score_frame %>% filter(Sample == "Sample 1")

  #Check if user specified groups
  if (groupings == TRUE) {
    #Map groupings to samples
    Group_Frame = data.frame(Groups = groupings,
                             Sample = unique(score_frame$Sample))

    #Join to replace sample with group
    score_frame = left_join(score_frame, Group_Frame) %>%
      dplyr::select(Sample = Groups, Coordinate, Boundary)

    score_frame = score_frame %>% group_by(Sample, Coordinate) %>%
      mutate(Boundary = median(Boundary)) %>% distinct()
  }

  #Get differences and convert to boundary scores

  score_frame = score_frame %>% group_by(Sample)  %>%
    mutate(Diff_Score = (base_sample$Boundary-Boundary)) %>%
    ungroup() %>% mutate(Diff_Score = (Diff_Score-
                           mean(Diff_Score, na.rm = TRUE))/
                           sd(Diff_Score, na.rm =TRUE)) %>% ungroup() %>% mutate(
                         TAD_Score = (Boundary-mean(Boundary, na.rm = TRUE))/
                           sd(Boundary, na.rm = TRUE))

  #Determine if boundaries are differential or non-differential
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

  #Subset differential frame to only include differential points

  Differential_Points = score_frame %>% filter(Differential == "Differential")

  #Pulling out consensus for classification

  TAD_Iden = TAD_Frame[,c(-1, -ncol(TAD_Frame))]>3

  #Classify time trends
  All_Non_TADs = apply(TAD_Iden, 1, function(x) all(x == FALSE))

  #Split into 4 groups for summarization
  Num_Points = seq_len(ncol(TAD_Iden))

  four_groups = split(Num_Points, ggplot2::cut_number(Num_Points,4))

  #Get summary for each group and put back together

  Full_Sum = do.call(cbind.data.frame,
          lapply(four_groups, function(x) (rowSums(as.matrix(TAD_Iden[,x]))/
    ncol(as.matrix(TAD_Iden[,x])))>=.5))

  #Define Highly Common TADs
  Common_TADs = (Full_Sum[,1] == Full_Sum[,2]) &
     (Full_Sum[,1] == Full_Sum[,ncol(Full_Sum)])

  Common_TADs = ifelse(Common_TADs, "Highly Common TAD", NA)

  #Define Early Appearing TADs
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

  TAD_Frame_Sub = TAD_Frame %>%
    dplyr::filter_at(dplyr::vars(`Sample 1`:Consensus_Score), dplyr::any_vars(.>3))

  return(list(TAD_Bounds = TAD_Frame_Sub,
              All_Bounds = TAD_Frame))
}


