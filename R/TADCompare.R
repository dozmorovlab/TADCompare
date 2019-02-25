#' Algorithm for differential TAD detection
#'
#' @param cont_mat1 The first contact matrix to compare in either sparse 3 column, n x n or n x (n+3)
#' form where the first three columns are coordinates in BED format.
#' If an x n matrix is used, the column names must correspond to the start
#' point of the corresponding bin. Required.
#' @param cont_mat2 The second contact matrix to compare in either sparse 3 column, n x n or n x (n+3)
#' form where the first three columns are coordinates in BED format.
#' If an x n matrix is used, the column names must correspond to the start
#' point of the corresponding bin. Required.
#' @param resolution Resolution of the contact matrix. Required
#' @param z_thresh Threshold for boundary score. Corresponds to the z-score of
#' differences between eigenvector gaps. Default is 2
#' @param window_size Size of the sliding window. Defaults to 2000000/resolution
#' @param gap_threshold Corresponds to the percentage of zeros allowed before
#' a column/row is removed from the analysis. 1=100\%, .7 = 70\%, etc. Default is .8.
#' @return \item{diff_loci}{Data frame containing differential regions and their
#' corresponding gap score. Additionally, the matrix in which the TAD boundary
#' exists is listed.}
#' @return \item{sd_diff}{Vector containing all boundary scores}
#' @export
#' @details Given two contact matrices, TADCompare returns a dataframe listing
#' regions which contain differential TADs. It also returns a continuous score
#' of boundary differences that can be used to compare regions of the genome.
#' @examples
#' #Read in data
#' data("rao_chr20_25_prim")
#' data("rao_chr20_25_rep")
#' #Find differential TADs
#' diff_tads <- TADCompare(rao_chr20_25_prim, rao_chr20_25_rep, resolution = 25000)
#'
TADCompare = function(cont_mat1, cont_mat2, resolution, z_thresh = 2, window_size = 2000000/resolution, gap_thresh = .8) {

  require(PRIMME)
  require(Matrix)
  require(dplyr)
  require(cluster)

  #Set maximize size of sliding window

  #resolution = as.numeric(colnames(cont_mat1)[2])-as.numeric(colnames(cont_mat1)[1])
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

  if (end+window_size>nrow(cont_mat1)) {
    end = nrow(cont_mat1)
  }

  point_dists1 = c()
  point_dists2 = c()
  Regions = c()

  while (end_loop == 0) {
    #Subsetting
    sub_filt1 = cont_mat1[start:end, start:end]
    sub_filt2 = cont_mat2[start:end, start:end]

    #Removing gap regions

    Per_Zero1 = colSums(sub_filt1 ==0)/nrow(sub_filt1)
    Per_Zero2 = colSums(sub_filt2 ==0)/nrow(sub_filt2)

    sub_gaps1 = Per_Zero1<gap_thresh
    sub_gaps2 = Per_Zero2<gap_thresh

    comp_rows = sub_gaps1 & sub_gaps2
    sub_filt1 = sub_filt1[ comp_rows, comp_rows]
    sub_filt2 = sub_filt2[ comp_rows, comp_rows]

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


    #sub_mat = crossprod(diag(dr^-(1/2)), as.matrix(sub_filt))

    #Reading names

    colnames(sub_mat1) = colnames(sub_mat2) = colnames(sub_filt1)

    #Find gaps at 2mb and remove

    if (length(sub_mat1) == 0) {
      start = start+max_end
      end = end+max_end

      } else {


    #Remove gaps if true

    # if (length(sub_gaps) != 0) {
    #   sub_mat = sub_mat[!sub_gaps, !sub_gaps]
    # }


    # Per_Zero1 = colSums(sub_mat1 ==0)/nrow(sub_mat1)
    # Per_Zero2 = colSums(sub_mat1 ==0)/nrow(sub_mat2)

    #Get first two eigenvectors

    Eigen1 = eigs_sym(sub_mat1, NEig = 2)

    eig_vals1 = Eigen1$values
    eig_vecs1 = Eigen1$vectors

    #Get order of eigenvalues from largest to smallest

    large_small1 = order(-eig_vals1)

    eig_vals1 = eig_vals1[large_small1]
    eig_vecs1 = eig_vecs1[,large_small1]

    #Matrix 2

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

    #Get number of rows without 80% of contacts and make NA to control sparsity

    point_dists1 = c(point_dists1, point_dist1)
    point_dists2 = c(point_dists2, point_dist2)

    Regions = c(Regions, colnames(sub_filt1))

      }
    #Test if we've reached end of matrix
    if (end == nrow(cont_mat1)) {
      end_loop = 1
    }

    #Set new start and end for window
    start = start+max_end
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
  sd_diff = dist_diff/(sd(dist_diff, na.rm = TRUE))

  #Get areas with high z-scores
  gaps = which(abs(sd_diff)>z_thresh)

  #Put differential regions into a data frame
  diff_loci = data.frame(Region = as.numeric(Regions)[gaps],Gap_Score = sd_diff[gaps])

  #Return differential TAD boundaries
  diff_loci = diff_loci %>% mutate(TAD_Region = ifelse(Gap_Score>0, "Matrix 1", "Matrix 2"))

  return(list(diff_loci = diff_loci, dist_diff = sd_diff))
}


