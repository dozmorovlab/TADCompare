#' Consensus boundary identification
#'
#' @import dplyr
#' @import magrittr
#' @import PRIMME
#' @param cont_mats List of contact matrices in either sparse 3 column,
#' n x n or n x (n+3) form where the first three columns are coordinates in
#' BED format. See "Input_Data" vignette for more information. 
#' If an x n matrix is used, the column names must correspond to
#' the start point of the corresponding bin. Required.
#' @param resolution Resolution of the data. Used to assign TAD boundaries
#' to genomic regions. If not provided, resolution will be estimated from
#' column names of the first matrix. Default is "auto"
#' @param z_thresh Threshold for boundary score. Higher values result in a
#' higher threshold for differential TADs. Default is 3.
#' @param window_size Size of sliding window for TAD detection, measured in bins.
#' Results should be consistent Default is 15.
#' @param gap_thresh Required \% of non-zero entries before a region will
#' be considered non-informative and excluded. Default is .2
#' @return A list containing consensus TAD boundaries and overall scores
#'  \itemize{
#'  \item Consensus - Data frame containing location of all consensus
#'  boundaries. Coordinate is the region of the genome, Sample columns
#'  correspond to individual boundary scores. Consensus_Score is consensus
#'  boundary score.
#'  \item All_Regions - Data frame containing consensus scores for all regions.
#'  All columns are identiical to the Consensus object.
#' }
#' @export
#' @details Given a list of sparse 3 column, n x n , or n x (n+3) contact matrices,
#' ConsensusTADs provides the set of consensus TAD boundaries across
#' them. Consensus TADs are defined by the consensus boundary score,
#' a score measuring TAD boundary likelihood across all matrices.
#' @examples
#' # Read in data
#' data("time_mats")
#' # Find consensus TAD boundaries
#' diff_list <- ConsensusTADs(time_mats, resolution = 50000)


ConsensusTADs = function(cont_mats, resolution,
                      z_thresh = 3, window_size = 15, gap_thresh = .2) {

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
      #Find the start coordinates based on the second column of the
      #bed file portion of matrix

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

    } else if (col_test!=3 & (row_test != col_test) &
               (col_test-row_test != 3)) {

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

    dist_sub = .single_dist(cont_mats[[x]], resolution,
                           window_size = window_size,
                           gap_thresh = gap_thresh)
    dist_sub = data.frame(Sample = paste("Sample", x), dist_sub[,c(2,3)])
    dist_sub
  })

  #Reduce matrices to only include shared regions
  coord_sum = lapply(bound_scores, function(x) x[,2])
  shared_cols = Reduce(intersect, coord_sum)
  bound_scores = lapply(bound_scores, function(x) x %>%
                          dplyr::filter(as.numeric(Coordinate) %in%
                                          as.numeric(shared_cols)))

  #Bind boundary scores
  score_frame = dplyr::bind_rows(bound_scores)

  colnames(score_frame)[1] = "Sample"
  base_sample = score_frame %>% filter(Sample == "Sample 1")

  #Get differences in boundary scores and convert to z-scores

  score_frame = score_frame %>% dplyr::group_by(Sample)  %>%
    dplyr::mutate(Diff_Score = (base_sample$Boundary-Boundary)) %>%
    dplyr::ungroup() %>% dplyr::mutate(Diff_Score = (Diff_Score-
                           mean(Diff_Score, na.rm = TRUE))/
                           sd(Diff_Score, na.rm =TRUE)) %>%
                          dplyr::ungroup() %>% dplyr::mutate(
                         TAD_Score = (Boundary-mean(Boundary, na.rm = TRUE))/
                           sd(Boundary, na.rm = TRUE))

  #Identify differential boundaries
  score_frame = score_frame %>% dplyr::mutate(Differential =
                                       ifelse(abs(Diff_Score)>z_thresh,
                                       "Differential", "Non-Differential"),
                                       Differential = ifelse(is.na(Diff_Score),
                                                            "Non-Differential",
                                                            Differential))

  #Getting a frame summarizing boundaries and subset to only have significant

  score_frame = score_frame %>%
    dplyr::select(Sample, Coordinate, TAD_Score) %>%
    dplyr::arrange(as.numeric(gsub("Sample", "", Sample))) %>%
    dplyr::mutate(Sample = factor(Sample, levels = unique(Sample)))

  #Filtering for TADs specifically
  TAD_Frame = score_frame %>% group_by(Coordinate) %>% filter(any(TAD_Score>3))


  #Spread the score frame and TAD frame into wide format
  TAD_Frame = tidyr::spread(as.data.frame(TAD_Frame),
                            key = Sample,
                            value = TAD_Score)

  score_frame = tidyr::spread(as.data.frame(score_frame),
                              key = Sample,
                              value = TAD_Score)

  #Get median of boundary scores for each row

  TAD_Frame = TAD_Frame %>%
    dplyr::mutate(Consensus_Score =
    apply(TAD_Frame %>% dplyr::select(-Coordinate) %>% as.matrix(.)
          ,1, median))

  score_frame = score_frame %>% ungroup() %>%
    dplyr::mutate(Consensus_Score =
                    apply(score_frame %>% dplyr::select(-Coordinate) %>%
                            as.matrix(.)
                          ,1, median))

  #Get a subset of TAD_Frame that contains only regions with significant TADs


  #Return consensus TAD boundaries and scores for all regions
  return(list(Consensus = TAD_Frame,
              All_Regions = score_frame))
}

.single_dist = function(cont_mat1,
                        resolution,
                        window_size = 15,
                        gap_thresh = .2) {

  #Remove full gaps from matrices

  non_gaps = which(colSums(cont_mat1) !=0)

  cont_mat_filt = cont_mat1[non_gaps,non_gaps]

  #Setting window size parameters
  max_end = window_size

  start = 1
  end = max_end

  end_loop = 0

  if ( (end+window_size)>nrow(cont_mat_filt)) {
    end = nrow(cont_mat_filt)
  }

  point_dists1 = rbind()

  while (end_loop == 0) {

    #Get windowed portion of matrix
    sub_filt1 = cont_mat_filt[seq(start, end, 1), seq(start,end,1)]

    #Identify columns with more than the gap threshold of zeros
    Per_Zero1 = (colSums(sub_filt1 != 0)/nrow(sub_filt1))<gap_thresh

    #Remove rows and dcolumns with zeros above gap threshold
    sub_filt1 = sub_filt1[!Per_Zero1, !Per_Zero1]

    #Test if matrix is too small to analyze
    if ((length(sub_filt1) == 0) | (length(sub_filt1) == 1)) {
      if (end == nrow(cont_mat1)) {
        end_loop = 1
      }

      #Move window to next point
      start = end
      end = end+max_end

      #Check if window is overlapping end of matrix and shorten if so
      if ( (end + max_end) >nrow(cont_mat_filt)) {
        end = nrow(cont_mat_filt)
      }

      #Check if matrix starts at end of matrix and kill loop
      if (start == end | start>nrow(cont_mat_filt)) {
        end_loop = 1
      }
      next
    }

    #Creating the normalized laplacian

    #Calculate row sums (degrees)
    dr1 = rowSums(abs(sub_filt1))

    #Perturbation factor for stability
    Dinvsqrt1 = diag((1/sqrt(dr1+2e-16)))

    #Form degree matrix
    P_Part1 = crossprod(as.matrix(sub_filt1), Dinvsqrt1)
    sub_mat1 = crossprod(Dinvsqrt1, P_Part1)

    #sub_mat = crossprod(diag(dr^-(1/2)), as.matrix(sub_filt))

    #Set column names to match original matrix
    colnames(sub_mat1) = colnames(sub_filt1)

    #Get first two eigenvectors

    Eigen1 = eigs_sym(sub_mat1, NEig = 2)

    #Pull out eigenvalues and eigenvectors
    eig_vals1 = Eigen1$values
    eig_vecs1 = Eigen1$vectors

    #Get order of eigenvalues from largest to smallest

    large_small1 = order(-eig_vals1)

    eig_vals1 = eig_vals1[large_small1]
    eig_vecs1 = eig_vecs1[,large_small1]

    #Normalize the eigenvectors

    norm_ones = sqrt(dim(sub_mat1)[2])

    for (i in seq_len(dim(eig_vecs1)[2])) {
      eig_vecs1[,i] = (eig_vecs1[,i]/sqrt(sum(eig_vecs1[,i]^2)))  * norm_ones
      if (eig_vecs1[1,i] !=0) {
        eig_vecs1[,i] = -1*eig_vecs1[,i] * sign(eig_vecs1[1,i])
      }
    }

    #Get rows and columns of contact matrix

    n = dim(eig_vecs1)[1]
    k = dim(eig_vecs1)[2]

    #Project eigenvectors onto a unit circle

    vm1 = matrix(kronecker(rep(1,k),
                           as.matrix(sqrt(rowSums(eig_vecs1^2)))),n,k)
    eig_vecs1 = eig_vecs1/vm1

    #Get distance between points on circle

    point_dist1 = sqrt(
      rowSums( (eig_vecs1-rbind(NA,eig_vecs1[-nrow(eig_vecs1),]))^2)
      )

    #Match column names (Coordinates) with eigenvector distances
    point_dist1 = cbind( match(colnames(sub_mat1),colnames(cont_mat1)),
                         as.numeric(colnames(sub_mat1)), point_dist1)

    #Combine current distances with old distances
    point_dists1 = rbind(point_dists1, point_dist1[-1,])

    #Check if window is at the end of matrix and kill loop
    if (end == nrow(cont_mat1)) {
      end_loop = 1
    }

    #Reset window
    start = end
    end = end+max_end

    #Check if window is near end of matrix and expand to end if true
    if ( (end + max_end) >nrow(cont_mat_filt)) {
      end = nrow(cont_mat_filt)
    }

    #Check if start of window overlaps end of matrix and kill if true
    if (start == end | start>nrow(cont_mat_filt)) {
      end_loop = 1
    }
  }

  #Convert to data frame
  point_dists1 = as.data.frame(point_dists1)
  colnames(point_dists1) = c("Index", "Coordinate","Boundary")
  return(point_dists1 = point_dists1)
}




