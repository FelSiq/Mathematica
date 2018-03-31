from matrix import matrix, rref, determinant
import copy
import numpy as np

"""
This Algorithm implements the Jacobi Method for linear system solving.
It's a iterative method.
"""

def relativeError(a, b, delta=1e-9):
	return max([abs(a[i] - b[i]) for i in range(min(len(a), len(b)))])/(delta + max(a))

def jacobi(A, b, x0, itMax=100, epsilon=1e-4, showError=False):
	n = len(x0)
	nrow, ncol = matrix.mshape(A)
	if not (nrow == ncol == n):
		print('Error: A matrix must be square and compatible with the dimension of x vector.')
		return None

	if not determinant.reducedDet(A):
		print('Error: A must be non-singular (det(A) != 0.0).')
		return None

	xCur = copy.copy(x0)
	M = copy.deepcopy(A)
	Dinv = matrix.minit(n, n)
	for i in range(n):
		Dinv[i][i] = 1.0/A[i][i]
		M[i][i] = 0.0 # M is (A - D), D is diag(A).

	err = epsilon * 2.0
	i = 0
	if showError:
		print('# :\t Relative error:')
	while i < itMax and err > epsilon:
		i += 1
		xNew = matrix.mmult(Dinv, (matrix.msub(b, matrix.mmult(M, xCur))))
		err = relativeError(matrix.mtransp(xNew), matrix.mtransp(xCur))
		xCur = xNew
		if showError:
			print(i, ':\t', err)
	return xCur	

if __name__ == '__main__':
	M = [
		[2, 1],
		[3, 4]
	]

	b = matrix.mtransp([1, -1])

	x0 = matrix.mtransp([0, 0])

	matrix.mprint(jacobi(M, b, x0, showError=True))
