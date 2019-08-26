#' Differential TAD boundary detection
#'
#' @import dplyr
#' @import magrittr
#' @import PRIMME
#' @import ggplot2
#' @param cont_mat1 Contact matrix in either sparse 3 column, n x n or n x (n+3)
#' form where the first three columns are coordinates in BED format.
#' If an x n matrix is used, the column names must correspond to the start
#' point of the corresponding bin. If large mode is selected, then
#' this matrix must be a tab-seperated n x n or n x (n+3) and it should be the
#' path to a contact matrix. Required.
#' @param cont_mat2 Second contact matrix, used for differential comparison,
#' must be in same format as cont_mat1. Required
#' @param resolution Resolution of the data. Used to assign TAD boundaries
#' to genomic regions. If not provided, resolution will be estimated from
#' column names of matrix. Default is "auto"
#' @param z_thresh Threshold for boundary score. Higher values result in a
#' higher threshold for differential TADs. Default is 2.
#' @param gap_thresh Required \% of 0s before a region will be considered a gap
#' and excluded. Default is .8
#' @param pre_TADs Set TRUE if a pre-defined set of TAD boundaries will be
#' tested (Defined by TADs parameter). If FALSE, TADCompare will test all
#' possible boundaries. Default is FALSE
#' @param TADs Data frame containing TAD boundaries to be tested. Column with
#' boundaries must be named "end". Only works if pre_TADs is TRUE. Default is
#' NULL
#' @return A list containing differential TAD characteristics
#'  \itemize{
#'  \item TAD_Frame - Data frame containing any region where a TAD boundary
#'  was detected.
#'  \item Diff_Loci - Data frame containing any regions with differential
#'  loci
#'  \item Boundary_Scores - Boundary scores for the entire genome
#' }
#' @export
#' @details Given two sparse 3 column, n x n , or n x (n+3) contact matrices,
#' TADCompare identifies differential TAD boundaries. Using a novel boundary
#' score metric, TADCompare simultaneously identifies TAD boundaries and
#' tests for differential enrichment between datasets. Strength of interactions
#' are provided using raw boundary scores and p-values.
#' @examples
#' #Read in data
#' data("rao_chr20_25_rep")
#' data("rao_chr20_25_prim")
#' #Find differential TADs
#' diff_list <- TADCompare(rao_chr20_25_rep, rao_chr20_25_prim,
#' resolution = 25000)

TADCompare = function(cont_mat1, cont_mat2, resolution = "auto",
                      z_thresh = 2, window_size = 25,
                      gap_thresh = .8, pre_TADs = FALSE,
                      TADs = NULL) {

  #Pulling out dimensions to test for matrix type
  row_test = dim(cont_mat1)[1]
  col_test = dim(cont_mat1)[2]

  if (row_test == col_test) {
    if (all(is.finite(cont_mat1)) == FALSE) {
      stop("Contact matrix 1 contains non-numeric entries")
    }

    if (all(is.finite(cont_mat2)) == FALSE) {
      stop("Contact matrix 2 contains non-numeric entries")
    }
  }

  if (col_test == 3) {


    #Convert sparse matrix to n x n matrix

    message("Converting to n x n matrix")

    if (nrow(cont_mat1) == 1) {
      stop("Matrix 1 is too small to convert to full")
    }

    if (nrow(cont_mat2) == 1) {
      stop("Matrix 2 is too small to convert to full")
    }

    cont_mat1 = HiCcompare::sparse2full(cont_mat1)
    cont_mat2 = HiCcompare::sparse2full(cont_mat2)

    if (all(is.finite(cont_mat1)) == FALSE) {
      stop("Contact matrix 1 contains non-numeric entries")
    }

    if (all(is.finite(cont_mat2)) == FALSE) {
      stop("Contact matrix 2 contains non-numeric entries")
    }
    if (resolution == "auto") {
      message("Estimating resolution")
      resolution = as.numeric(names(table(as.numeric(colnames(cont_mat1))-lag(as.numeric(colnames(cont_mat1)))))[1])
    }

  } else if (col_test-row_test == 3) {

    message("Converting to n x n matrix")

    #Find the start coordinates based on the second column of the bed file portion of matrix

    start_coords = cont_mat1[,2]

    #Calculate resolution based on given bin size in bed file

    resolution = as.numeric(cont_mat1[1,3])-as.numeric(cont_mat1[1,2])

    #Remove bed file portion

    cont_mat1 = as.matrix(cont_mat1[,-c(seq_len(3))])
    cont_mat2 = as.matrix(cont_mat2[,-c(seq_len(3))])

    if (all(is.finite(cont_mat1)) == FALSE) {
      stop("Contact matrix 1 contains non-numeric entries")
    }

    if (all(is.finite(cont_mat2)) == FALSE) {
      stop("Contact matrix 2 contains non-numeric entries")
    }

    #Make column names correspond to bin start

    colnames(cont_mat1) = start_coords
    colnames(cont_mat2) = start_coords


  } else if (col_test!=3 & (row_test != col_test) & (col_test-row_test != 3)) {

    #Throw error if matrix does not correspond to known matrix type

    stop("Contact matrix must be sparse or n x n or n x (n+3)!")

  } else if ( (resolution == "auto") & (col_test-row_test == 0) ) {
    message("Estimating resolution")

    #Estimating resolution based on most common distance between loci

    resolution = as.numeric(names(table(as.numeric(colnames(cont_mat1))-lag(as.numeric(colnames(cont_mat1)))))[1])
  }

  #Make sure contact matrices only include shared columns
  coord_sum = list(colnames(cont_mat1), colnames(cont_mat2))

  #Only include shared columns in analysis
  shared_cols = Reduce(intersect, coord_sum)

  cont_mat1 = cont_mat1[colnames(cont_mat1) %in% shared_cols,
                        colnames(cont_mat1) %in% shared_cols]

  cont_mat2 = cont_mat2[colnames(cont_mat2) %in% shared_cols,
                        colnames(cont_mat2) %in% shared_cols]

  #Set maximize size of sliding window

  window_size = window_size

  #Remove full gaps from matrices

  non_gaps = which(colSums(cont_mat1) !=0 & (colSums(cont_mat2) !=0))

  #Remove gaps
  cont_mat1 = cont_mat1[non_gaps,non_gaps]
  cont_mat2 = cont_mat2[non_gaps,non_gaps]

  #Defining window size
  max_end = window_size
  max_size = window_size/ceiling(200000/resolution)
  min_size = ceiling(200000/resolution)
  Group_over = bind_rows()

  start = 1
  end = max_end

  end_loop = 0

  #If window is larger than matrix make it equal to matrix size
  if (end+window_size>nrow(cont_mat1)) {
    end = nrow(cont_mat1)
  }

  #Pre-allocate vectors
  point_dists1 = c()
  point_dists2 = c()
  Regions = c()

  while (end_loop == 0) {
    #Subsetting
    sub_filt1 = cont_mat1[start:end, start:end]
    sub_filt2 = cont_mat2[start:end, start:end]

    #Removing gap regions from sub_matrices

    Per_Zero1 = colSums(sub_filt1 ==0)/nrow(sub_filt1)
    Per_Zero2 = colSums(sub_filt2 ==0)/nrow(sub_filt2)

    #Remove columns with more zeros than threshold
    sub_gaps1 = Per_Zero1<gap_thresh
    sub_gaps2 = Per_Zero2<gap_thresh

    comp_rows = sub_gaps1 & sub_gaps2
    sub_filt1 = sub_filt1[ comp_rows, comp_rows]
    sub_filt2 = sub_filt2[ comp_rows, comp_rows]

    #Slide window to end if window size is less than 2
    if ( (length(sub_filt1) == 0) | (length(sub_filt1) == 1) ) {
      start = start+max_end
      end = end+max_end

    } else {

    #Getting degree matrices

    dr1 = rowSums(abs(sub_filt1))
    dr2 = rowSums(abs(sub_filt2))

    #Creating the normalized laplacian

    Dinvsqrt1 = diag((1/sqrt(dr1+2e-16)))
    Dinvsqrt2 = diag((1/sqrt(dr2+2e-16)))

    P_Part1 = crossprod(as.matrix(sub_filt1), Dinvsqrt1)
    sub_mat1 = crossprod(Dinvsqrt1, P_Part1)

    P_Part2 = crossprod(as.matrix(sub_filt2), Dinvsqrt2)
    sub_mat2 = crossprod(Dinvsqrt2, P_Part2)

    #Reading names

    colnames(sub_mat1) = colnames(sub_mat2) = colnames(sub_filt1)

    #Find gaps at 2mb and remove


    #Get first two eigenvectors

    Eigen1 = PRIMME::eigs_sym(sub_mat1, NEig = 2)

    eig_vals1 = Eigen1$values
    eig_vecs1 = Eigen1$vectors

    #Get order of eigenvalues from largest to smallest

    large_small1 = order(-eig_vals1)

    eig_vals1 = eig_vals1[large_small1]
    eig_vecs1 = eig_vecs1[,large_small1]

    #Repeat for matrix 2

    Eigen2 = eigs_sym(sub_mat2, NEig = 2)

    eig_vals2 = Eigen2$values
    eig_vecs2 = Eigen2$vectors

    #Get order of eigenvalues from largest to smallest

    large_small2 = order(-eig_vals2)

    eig_vals2 = eig_vals2[large_small2]
    eig_vecs2 = eig_vecs2[,large_small2]


    #Normalize the eigenvectors

    norm_ones = sqrt(dim(sub_mat1)[2])

    for (i in 1:dim(eig_vecs1)[2]) {
      eig_vecs1[,i] = (eig_vecs1[,i]/sqrt(sum(eig_vecs1[,i]^2)))  * norm_ones
      if (eig_vecs1[1,i] !=0) {
        eig_vecs1[,i] = -1*eig_vecs1[,i] * sign(eig_vecs1[1,i])
      }
    }


    for (i in 1:dim(eig_vecs2)[2]) {
      eig_vecs2[,i] = (eig_vecs2[,i]/sqrt(sum(eig_vecs2[,i]^2)))  * norm_ones
      if (eig_vecs2[1,i] !=0) {
        eig_vecs2[,i] = -1*eig_vecs2[,i] * sign(eig_vecs2[1,i])
      }
    }

    eps = 2.2204e-16

    n = dim(eig_vecs1)[1]
    k = dim(eig_vecs1)[2]

    #Project eigenvectors onto a unit circle

    vm1 = matrix(kronecker(rep(1,k), as.matrix(sqrt(rowSums(eig_vecs1^2)))),n,k)
    eig_vecs1 = eig_vecs1/vm1

    vm2 = matrix(kronecker(rep(1,k), as.matrix(sqrt(rowSums(eig_vecs2^2)))),n,k)
    eig_vecs2 = eig_vecs2/vm2

    #Get distance between points on circle

    point_dist1 = sqrt(rowSums( (eig_vecs1-rbind(NA,eig_vecs1[-nrow(eig_vecs1),]))^2  ))
    point_dist2 = sqrt(rowSums( (eig_vecs2-rbind(NA,eig_vecs2[-nrow(eig_vecs2),]))^2  ))

    #Remove NA entry at start of windows

    point_dists1 = c(point_dists1, point_dist1[-1])
    point_dists2 = c(point_dists2, point_dist2[-1])

    #Assign to regions based on column names
    Regions = c(Regions, colnames(sub_filt1)[-1])

    }

    #Test if we've reached end of matrix
    if (end == nrow(cont_mat1)) {
      end_loop = 1
    }

    #Set new start and end for window
    start = end
    end = end+max_end

    if ( (end + max_end) >nrow(cont_mat1)) {
      end = nrow(cont_mat1)
    }

    if (start == end | start>nrow(cont_mat1)) {
      end_loop = 1
    }
  }
  #Calculating the difference between log gaps

  dist_diff = (point_dists1)-(point_dists2)
  #Getting the z-scores
  sd_diff = (dist_diff-mean(dist_diff, na.rm = TRUE))/(sd(dist_diff, na.rm = TRUE))

  TAD_Score1 = (point_dists1-mean(point_dists1, na.rm = TRUE))/
    (sd(point_dists1, na.rm = TRUE))

  TAD_Score2 = (point_dists2-mean(point_dists2, na.rm = TRUE))/
    (sd(point_dists2, na.rm = TRUE))

  #Get areas with high z-scores
  gaps = which(abs(sd_diff)>z_thresh)

  #Put differential regions into a data frame
  diff_loci = data.frame(Region = as.numeric(Regions)[gaps],Gap_Score = sd_diff[gaps])

  #Return differential TAD boundaries
  diff_loci = diff_loci %>% mutate(P_value = 2*(1-pnorm(abs(Gap_Score))), TAD_Region = ifelse(Gap_Score>0, "Matrix 1", "Matrix 2"))
  Gap_Scores = data.frame(Region = as.numeric(Regions), TAD_Score1 = TAD_Score1,
                          TAD_Score2 =TAD_Score2,Gap_Score = sd_diff)
  TAD_Frame = data.frame(Boundary = as.numeric(Regions),
                         Gap_Score = sd_diff,
                         TAD_Score1, TAD_Score2)

  #Assign labels to boundary type and identify which matrix has the boundary

  TAD_Frame = TAD_Frame %>%
    filter( (TAD_Score1>3) | TAD_Score2>3) %>%
    mutate(Differential = ifelse(abs(Gap_Score)>z_thresh, "Differential",
                                 "Non-Differential"),
           Enriched_In = ifelse(Gap_Score>0, "Matrix 1", "Matrix 2")) %>%
    arrange(Boundary) %>%
    mutate(Bound_Dist = pmin(abs(Boundary-lag(Boundary))/resolution,
           abs((Boundary-lead(Boundary)))/resolution)) %>%

    mutate(Differential = ifelse( (Differential == "Differential") &
                                    (Bound_Dist<=5) & !is.na(Bound_Dist),
                                  "Shifted", Differential)) %>%
    dplyr::select(-Bound_Dist)

  #Classifying merged-split
  TAD_Frame = TAD_Frame %>%
    mutate(Type = ifelse( (Differential == "Differential") &
                            (lag(Differential) == "Non-Differential") &
                            (lead(Differential) == "Non-Differential"),
                          ifelse(Enriched_In == "Matrix 1", "Split", "Merge"),
                          Differential))

  #See if TADs are pre-specified

  if (pre_TADs == TRUE) {

    #Subset to only include preset TADs
    Full_TADs = unique(c(TADs$end))
    TAD_Locs = match(Full_TADs, Regions)
    TAD_Frame = data.frame(Boundary = Full_TADs, Gap_Score = sd_diff[TAD_Locs])
    TAD_Frame = TAD_Frame %>% arrange(Boundary)

    #Get gap scores and subset based on pre-defined TADs
    Ind_Score1 = ((point_dists1-mean(point_dists1, na.rm = TRUE))/
                    sd(point_dists1, na.rm = TRUE))[sort(TAD_Locs)]

    Ind_Score2 = ((point_dists2-mean(point_dists2, na.rm = TRUE))/
                    sd(point_dists2, na.rm = TRUE))[sort(TAD_Locs)]

    #Define boundayr types
    TAD_Frame = TAD_Frame %>%
      mutate(Differential = ifelse(abs(Gap_Score)>z_thresh, "Differential",
                                "Non-Differential"),
             Enriched_In = ifelse(Gap_Score>0, "Matrix 1", "Matrix 2")) %>%
      arrange(Boundary)

    #Testing if boundaries are within 5 bins in front (Shifted)

    TAD_Frame = TAD_Frame %>%
      mutate(Bound_Dist = (Boundary-lag(Boundary))/resolution) %>%
      mutate(Differential = ifelse( (Differential == "Differential") &
                                    (Bound_Dist<=5) & !is.na(Bound_Dist),
                                    "Shifted", Differential)) %>%
                                    dplyr::select(-Bound_Dist)

    #Testing if boundaries are within 5 bins behind

    TAD_Frame = TAD_Frame %>%
      mutate(Bound_Dist = abs((Boundary-lead(Boundary))/resolution)) %>%
      mutate(Differential = ifelse( (Differential == "Differential") &
                                      (Bound_Dist<=5) & !is.na(Bound_Dist),
                                    "Shifted", Differential)) %>%
      dplyr::select(-Bound_Dist) %>% filter(!is.na(Differential))

    #Defining merged and split boundaries

    TAD_Frame = TAD_Frame %>%
      mutate(Type = ifelse( (Differential == "Differential") &
                              (lag(Differential) == "Non-Differential") &
                              (lead(Differential) == "Non-Differential"),
                            ifelse(Enriched_In == "Matrix 1", "Split", "Merge"),
                            Differential) )

    #Appending TAD Scores

    TAD_Frame = cbind.data.frame(TAD_Frame, Ind_Score1, Ind_Score2)
    colnames(TAD_Frame)[c(ncol(TAD_Frame)-1, ncol(TAD_Frame))] =
      c("TAD_Score1", "TAD_Score2")


    return(list(TAD_Frame = TAD_Frame,
                Gap_Scores = Gap_Scores))
  }

  #Add up-down enrichment of TAD boundaries
  TAD_Frame = TAD_Frame %>%
    mutate(Type = ifelse( (TAD_Score1>3) &
                            (TAD_Score2>3) &
                            (Differential == "Differential"),
                          "Strength Change", Type))

  TAD_Size = TAD_Frame %>% group_by(Enriched_In) %>%
    transmute(Boundary_Distance = Boundary-lag(Boundary))

  Size_Plot = ggplot(TAD_Size, aes(x = Enriched_In, y = Boundary_Distance,
                                   fill = Enriched_In)) + geom_boxplot() +
    theme_bw(base_size = 24)

  TAD_Sum = TAD_Frame %>% group_by(Enriched_In, Type) %>% summarise(Count = n())

  Count_Plot = ggplot(TAD_Sum, aes(x = Enriched_In, y = Count, fill = Type)) +
    geom_bar(stat="identity") + theme_bw(base_size = 24)

  return(list(TAD_Frame =TAD_Frame,
              Diff_Loci = diff_loci,
              Gap_Scores = Gap_Scores,
              Size_Plot = Size_Plot,
              Count_Plot = Count_Plot ))
}





