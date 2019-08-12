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

  TAD_Cat = coalesce(as.character(Common_TADs),
                     as.character(Early_Appearing),
                     as.character(Late_Appearing),
                     as.character(Early_Disappearing),
                     as.character(Late_Disappearing),
                     as.character(Dynamic))

  TAD_Frame = TAD_Frame %>% dplyr::mutate(Category = TAD_Cat)

  return(list(Differential_Points = Differential_Points,
              Full = score_frame,
              Consensus = TAD_Frame,
              )
}


Single_Dist = function(cont_mat1, resolution, window_size = 15, gap_thresh = .8) {

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

Pairwise_z = function(x, z_thresh) {
  poss_combs = combn(1:ncol(x),2)
  lab_list = c()
  z_frame = cbind()
  for (i in 1:ncol(poss_combs)) {
    index_1 = poss_combs[1,i]
    index_2 = poss_combs[2,i]
    diff = x[,index_1]-x[,index_2]
    z_diff = ifelse(((diff-mean(diff))/sd(diff))>z_thresh,
                    "Differential",
                    "Non-Differential")
    lab_list = c(lab_list,paste0(index_1,"-", index_2))
    z_frame = cbind(z_frame, z_diff)
  }
  colnames(z_frame) = lab_list
  return(z_frame)
}

