# 1st Exercise

# --------------------------------
# Calc the nth power of a matrix
potM <- function(M, n=1){
	ret=M
	i=1
	while (i < n) {
		ret = ret%*%M
		i=1+i
	}
	return (ret)
}

# Question matrix
m <- matrix(
	c(0.3, 0.3, 0.4,
	  0.2, 0.6, 0.2,
	  0.0, 0.0, 1.0), 
	3, 3, byrow=T)


# --------------------------------
# 1. Solving by potM
m3 <- potM(m, 3)

# P(x_3=2|x_0=0) = m^3[0,2]
m3[1, 3]


# --------------------------------
# 2. Solving by eigen-stuff

# Eigen-stuff of a matrix
egstuff <- eigen(t(m))
egval <- egstuff$values
egvec <- egstuff$vectors

eigenValDiag <- diag(egval**3)

# m^n <- Q * V^n * Q^{-1}, where
#	Q is the eigenvector matrix (in the columns) of m
#	V is a diagonal matrix with the eigenvalues of m

# Output (must be equal to m3=potM(m, 3))
stc <- t(egvec %*% eigenValDiag %*% solve(egvec))

# Print output (must be equal to part 1)
stc[1, 3]

