from matrix import matrix, rref, determinant
import sys
import copy
import numpy as np

"""
This Algorithm implements the Jacobi Method for linear system solving.
It's a iterative method.
"""

def calcError(a, b, delta=1e-9, relative=True):
	err = max([abs(a[i] - b[i]) for i in range(min(len(a), len(b)))])
	return err/((delta + max(a)) if relative else 1.0)

def jacobi(A, b, x0, itMax=100, epsilon=1e-4, showError=False, relativeError=True, delta=1.0e-10):
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
	
	# C is the Jacobi Method iteration matrix (C = -D^-1 * (A - D))
	C = -np.matrix(Dinv) * np.matrix(M)

	# Infinite Norm of the C matrix (highest sum of absolute values of each row)
	CinfNorm = abs(C).sum(1).max()

	# This is a error estimation auxiliary constant value
	errAux = CinfNorm / (1.0 - CinfNorm + delta)	

	# Multiple value of which the error should be printed.
	itToPrint = max(int(itMax/20), 1)

	# Init Jacobi method iterations
	err = epsilon + 1.0
	i = 0
	if showError:
		print('# :\t', 'Relative' if relativeError else 'Infinite Norm', 'error:')
	while i < itMax and err > epsilon:
		i += 1
		xNew = matrix.mmult(Dinv, (matrix.msub(b, matrix.mmult(M, xCur))))
		err = errAux * calcError(matrix.mtransp(xNew), matrix.mtransp(xCur), relative=relativeError)
		xCur = xNew
		if showError and not (i % itToPrint):
			print(i, ':\t', err)
	print('Process finished @ ieration', i, 'with', err, 'error value.')	

	return xCur

if __name__ == '__main__':
	if len(sys.argv) <= 5:
		print('usage:', sys.argv[0], '<A filepath> <b filepath> <x0 filepath> <itMax # (integer)> <max error (float)>')
		exit(1)	

	try:
		M = matrix.mread(sys.argv[1])
		b = matrix.mread(sys.argv[2])
		x0 = matrix.mread(sys.argv[3])
	except:
		print('Error reading input matrices.')
		exit(4)

	try:
		itMax = int(sys.argv[4])
		if (itMax <= 0):
			raise Exception
	except:
		print('Error: \'itMax\' parameter must be a natural number!')
		exit(2)

	try:
		errMax = float(sys.argv[5])
		if (errMax < 0.0):
			raise Exception
	except:
		print('Error: \'maxError\' parameter must be a non-negative floating point!')
		exit(3)

	matrix.mprint(jacobi(M, b, x0, itMax=itMax,
		epsilon=errMax, showError=True, relativeError=False))
