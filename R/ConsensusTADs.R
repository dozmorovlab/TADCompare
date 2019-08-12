#' Differential TAD boundary detection
#'
#' @import dplyr
#' @import magrittr
#' @import PRIMME
#' @importFrom GenomicRanges GRanges GRangesList start
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
#'

ConsensusTADs = function(cont_mats, resolution,
                      z_thresh = 2, window_size = 25, gap_thresh = .8,
                       Pre_Def = TRUE) {

  require(PRIMME)
  require(Matrix)
  require(dplyr)
  require(cluster)
  require(softImpute)

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

  #Get differences

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
    mutate(Consensus_Score =
    apply(TAD_Frame %>% dplyr::select(-Coordinate) %>% as.matrix(.)
          ,1, median))

  return(list(Full = score_frame,
              Consensus = TAD_Frame))
}


Single_Dist = function(cont_mat1, resolution, window_size = 25, gap_thresh = .8) {

  require(PRIMME)
  require(Matrix)
  require(dplyr)
  require(cluster)

  #Set maximize size of sliding window

  #resolution = as.numeric(colnames(cont_mat1)[2])-as.numeric(colnames(cont_mat1)[1])
  window_size = window_size

  #Remove full gaps from matrices

  non_gaps = which(colSums(cont_mat1) !=0)

  cont_mat_filt = cont_mat1[non_gaps,non_gaps]

  max_end = window_size

  start = 1
  end = max_end

  end_loop = 0

  if ( (end+window_size)>nrow(cont_mat_filt)) {
    end = nrow(cont_mat_filt)
  }

  point_dists1 = rbind()

  while (end_loop == 0) {
    sub_filt1 = cont_mat_filt[start:end, start:end]

    Per_Zero1 = (colSums(sub_filt1 ==0)/nrow(sub_filt1))>gap_thresh

    sub_filt1 = sub_filt1[!Per_Zero1, !Per_Zero1]

    #sub_filt1[sub_filt1==0] =NA

    if ((length(sub_filt1) == 0) | (length(sub_filt1) == 1)) {
      if (end == nrow(cont_mat1)) {
        end_loop = 1
      }

      start = end
      end = end+max_end

      if ( (end + max_end) >nrow(cont_mat_filt)) {
        end = nrow(cont_mat_filt)
      }

      if (start == end | start>nrow(cont_mat_filt)) {
        end_loop = 1
      }
      next
    }

    #imputed = softImpute(sub_filt1)

    #sub_filt1 = complete(sub_filt1, imputed)

    #Creating the normalized laplacian

    dr1 = rowSums(abs(sub_filt1))

    Dinvsqrt1 = diag((1/sqrt(dr1+2e-16)))

    P_Part1 = crossprod(as.matrix(sub_filt1), Dinvsqrt1)
    sub_mat1 = crossprod(Dinvsqrt1, P_Part1)

    #sub_mat = crossprod(diag(dr^-(1/2)), as.matrix(sub_filt))

    colnames(sub_mat1) = colnames(sub_filt1)

    #Find gaps at 2mb and remove


    #Remove gaps if true

    # if (length(sub_gaps) != 0) {
    #   sub_mat = sub_mat[!sub_gaps, !sub_gaps]
    # }


    #Get first two eigenvectors

    Eigen1 = eigs_sym(sub_mat1, NEig = 2)

    eig_vals1 = Eigen1$values
    eig_vecs1 = Eigen1$vectors

    #Get order of eigenvalues from largest to smallest

    large_small1 = order(-eig_vals1)

    eig_vals1 = eig_vals1[large_small1]
    eig_vecs1 = eig_vecs1[,large_small1]

    #Normalize the eigenvectors

    norm_ones = sqrt(dim(sub_mat1)[2])

    for (i in 1:dim(eig_vecs1)[2]) {
      eig_vecs1[,i] = (eig_vecs1[,i]/sqrt(sum(eig_vecs1[,i]^2)))  * norm_ones
      if (eig_vecs1[1,i] !=0) {
        eig_vecs1[,i] = -1*eig_vecs1[,i] * sign(eig_vecs1[1,i])
      }
    }


    eps = 2.2204e-16

    n = dim(eig_vecs1)[1]
    k = dim(eig_vecs1)[2]

    #Project eigenvectors onto a unit circle

    vm1 = matrix(kronecker(rep(1,k), as.matrix(sqrt(rowSums(eig_vecs1^2)))),n,k)
    eig_vecs1 = eig_vecs1/vm1


    #Get distance between points on circle

    point_dist1 = sqrt(rowSums( (eig_vecs1-rbind(NA,eig_vecs1[-nrow(eig_vecs1),]))^2  ))

    #Get number of rows without gap threshold of contacts and make NA

    # point_dist1[Per_Zero1] = NA
    #
    # point_dist1[which(Per_Zero1)+1] = NA


    point_dist1 = cbind( match(colnames(sub_mat1),colnames(cont_mat1) ) , as.numeric(colnames(sub_mat1)), point_dist1)

    point_dists1 = rbind(point_dists1, point_dist1[-1,])

    if (end == nrow(cont_mat1)) {
      end_loop = 1
    }

    start = end
    end = end+max_end

    if ( (end + max_end) >nrow(cont_mat_filt)) {
      end = nrow(cont_mat_filt)
    }

    if (start == end | start>nrow(cont_mat_filt)) {
      end_loop = 1
    }
  }

  point_dists1 = as.data.frame(point_dists1)
  colnames(point_dists1) = c("Index", "Coordinate","Boundary")
  return(point_dists1 = point_dists1)
}
