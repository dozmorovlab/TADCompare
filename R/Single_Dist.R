Single_Dist = function(cont_mat1, resolution, window_size = 15, gap_thresh = .8) {

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
    sub_filt1 = cont_mat_filt[start:end, start:end]

    #Identify columns with more than the gap threshold of zeros
    Per_Zero1 = (colSums(sub_filt1 ==0)/nrow(sub_filt1))>gap_thresh

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

    for (i in 1:dim(eig_vecs1)[2]) {
      eig_vecs1[,i] = (eig_vecs1[,i]/sqrt(sum(eig_vecs1[,i]^2)))  * norm_ones
      if (eig_vecs1[1,i] !=0) {
        eig_vecs1[,i] = -1*eig_vecs1[,i] * sign(eig_vecs1[1,i])
      }
    }

    #Get rows and columns of contact matrix

    n = dim(eig_vecs1)[1]
    k = dim(eig_vecs1)[2]

    #Project eigenvectors onto a unit circle

    vm1 = matrix(kronecker(rep(1,k), as.matrix(sqrt(rowSums(eig_vecs1^2)))),n,k)
    eig_vecs1 = eig_vecs1/vm1

    #Get distance between points on circle

    point_dist1 = sqrt(rowSums( (eig_vecs1-rbind(NA,eig_vecs1[-nrow(eig_vecs1),]))^2  ))

    #Match column names (Coordinates) with eigenvector distances
    point_dist1 = cbind( match(colnames(sub_mat1),colnames(cont_mat1) ) , as.numeric(colnames(sub_mat1)), point_dist1)

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
