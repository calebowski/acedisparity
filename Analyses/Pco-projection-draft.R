X1 <- matrix(sample(0:2, size = 25, replace = TRUE), nrow = 5, ncol = 5)
## 2. change second matrix
X2 <- X1
X2[1,] <- sample(0:2, size = 5, replace = TRUE)
## 3. Compute distance matrices
D1 <- dist(X1)
D2 <- dist(X2)



## 5. Gram matrix and eigen decomposition (for D1)
dim <- nrow(X1)

similarities <- -0.5 * as.matrix(D1)^2

center_mat <- diag(dim) - (1/dim) * matrix(1, dim, dim)

gram_matrix <- center_mat %*% similarities %*% center_mat


## check if this is mathematically equal to eigen algorithm
eigen_result <- svd(gram_matrix)
eigenvector <- eigen_result$u

ord_mds_eigen <- gram_matrix %*% eigenvector


###### projection

similarities_two <- -0.5 * as.matrix(D2)^2
gram_two <- center_mat %*% similarities_two %*% center_mat
X2_proj <- gram_two %*% eigen_result$vector