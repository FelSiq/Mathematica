from matrix import matrix, rref

"""
Input: Symmetric (A = A^T) and Positive Definite (x^T * A * x > 0 for all x in R^n) square matrix A[n x n].

This matrix A can be decomposed in a lower triangular matrix L: A = L * L^T
Or a upper triangular matrix U: A = U^T * U
Remembering that L = U^T and, of course, U = L^T.

Both L and U have positive diagonal (l[i,i], u[i,i] > 0 for all i = 1, ..., n).

Output: L matrix.

Note: Cholesky Method is a good way to check if a symmetric matrix M 
is Positive Definite (algorithm complexity is O(n^3)).
"""

def cholesky(M):
	n, m = matrix.mshape(M)

	if n != m:
		print('Error: M matrix must be square for Cholesky Decomposition.')
		return None

	if not matrix.missymmetric(M):
		print('Error: M matrix must be symmetric for Cholesky Decompostion.')
		return None

	L = matrix.minit(n, n)

	for k in range(n):
		s = 0.0
		for j in range(k):
			s = L[k][j]**2.0

		if M[k][k] >= s:
			L[k][k] = (M[k][k] - s)**0.5
		else:
			print('Error: M matrix is not positive definite (can\'t perform Cholesky Decomposition)!')
			return None

		for j in range(k+1, n):
			s = 0.0
			for i in range(k):
				s = L[j][k] * L[k][i]
			L[j][k] = 1.0/L[k][k] * (M[j][k] - s) 

	return L

if __name__ == '__main__':
	M1 = [
		[2, -1, 0],
		[-1, 2, -1],
		[0, -1, 2]
	]
	M2 = [
		[6, 4, -3],
		[4, -2, 0],
		[-3, 0, 1]
	]
	M3 = [
		[6, 4, -3],
		[4, -2, 0],
		[-3, 1, 1]
	]

	matrix.mprint(cholesky(M1))
	matrix.mprint(cholesky(M2))
	matrix.mprint(cholesky(M3))
