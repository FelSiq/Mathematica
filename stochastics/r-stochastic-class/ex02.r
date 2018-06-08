# 2nd Exercise

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

greatest.estate <- 4 

matriz <- function ( greatest.estate ) {
	M <- matrix (0 , nrow = greatest.estate +1 , ncol = greatest.estate +1)
	M [1 , 2] <- 1
	M [ greatest.estate +1 , greatest.estate ] <- 1
	for ( j in 2: greatest.estate ) {



	M [j , j -1] <- (j -1) / greatest.estate
	M [j , j +1] <- ( greatest.estate - (j -1) ) / greatest.estate} 

	return ( M ) 
}

# Question matrix
m <- matriz(greatest.estate)

# --------------------------------
# 1.
for (i in 37:40) {
	t <- potM(m, i)
	print(t)
}

#       [,1] [,2] [,3] [,4]  [,5]
# [1,] 0.000  0.5 0.00  0.5 0.000
# [2,] 0.125  0.0 0.75  0.0 0.125
# [3,] 0.000  0.5 0.00  0.5 0.000
# [4,] 0.125  0.0 0.75  0.0 0.125
# [5,] 0.000  0.5 0.00  0.5 0.000
#       [,1] [,2] [,3] [,4]  [,5]
# [1,] 0.125  0.0 0.75  0.0 0.125
# [2,] 0.000  0.5 0.00  0.5 0.000
# [3,] 0.125  0.0 0.75  0.0 0.125
# [4,] 0.000  0.5 0.00  0.5 0.000
# [5,] 0.125  0.0 0.75  0.0 0.125
#       [,1] [,2] [,3] [,4]  [,5]
# [1,] 0.000  0.5 0.00  0.5 0.000
# [2,] 0.125  0.0 0.75  0.0 0.125
# [3,] 0.000  0.5 0.00  0.5 0.000
# [4,] 0.125  0.0 0.75  0.0 0.125
# [5,] 0.000  0.5 0.00  0.5 0.000
#       [,1] [,2] [,3] [,4]  [,5]
# [1,] 0.125  0.0 0.75  0.0 0.125
# [2,] 0.000  0.5 0.00  0.5 0.000
# [3,] 0.125  0.0 0.75  0.0 0.125
# [4,] 0.000  0.5 0.00  0.5 0.000
# [5,] 0.125  0.0 0.75  0.0 0.125

# 2.
# The matrix has a stationary state
# of period 2, as it goes in and back
# to the exatcly same two states shown
# above, depending if the state index
# is even or odd.

# 3.
e <- eigen(t(m))
e$values
# [1]  1.000000e+00 -1.000000e+00 -5.000000e-01  5.000000e-01  1.343484e-16

# The process has a recurrent stationary state,
# with periodicity two because it has a two 
# unitary eigenvalues:
# 
# let P be the transition matrix and V the
# stationary vector:
# V = [v[0] ... v[n]]^{T} and, by definition,
# V[j] = sum(V[i]*P(i,j)), for all j
#
# So, V = P^{t} * V => P^{t} * V = 1.0 * V, 
# which means that if V exists it is a P^{t}
# eigenvector and 1.0 is it's correspondent
# eigenvalue.
#
# In conclusion, m^{t} has two unitary eigenvalue,
# so there is two correspondent eigenvectors which
# represents m recurrent states.

mt<-t(m)
sys<-mt-diag(1, ncol(m)) # V * m = 1.0 * V => (m - I) * V = 0
sys[1,]<-rep(1, ncol(m)) # Change any equation to sum(V[i])=1
b<-c(1, rep(0, ncol(m)-1)) # Null vector, except to the y[i]=1 above

# Solving V=m^{t} * V <=> (m^{t} - I) * V = 0
v<-solve(sys, b)
v
# [1] 0.0625 0.2500 0.3750 0.2500 0.0625

mt%*%v == v
# 3.
# You can't get the same result with
# matrix power because the periodicity
# of the state matrix is not 1.

# 4.
# The values obtained with the stationary system
# is the average value of any row of the two 
# recurrent states.
(potM(m, 39)+potM(m, 40))/2
