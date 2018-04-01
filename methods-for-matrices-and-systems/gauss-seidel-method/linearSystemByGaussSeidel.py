from enum import IntEnum
import numpy as np
import sys
import copy

"""
INFORMATION:

This program solve a Ax = b system using Gauss-Seidel algorithm.

"""

def error(a, b, delta=1e-13, relative=False):
	denominator = (max(a) + delta) if relative else 1.0
	return max(a - b)/denominator

def minit(rownum, colnum, val=0.0):
	return np.matrix([[val] * colnum for i in range(rownum)])

def gaussSeidel(A, b, x0, itMax=100, epsilon=1e-5, showError=False, relativeError=False):
	if type(A) != np.matrix:
		print('Error: expecting np.matrix as input matrix \'A\'.')
		return None
	
	ArowNum, AcolNum = A.shape
	browNum, bcolNum = b.shape
	x0rowNum, x0colNum = x0.shape

	if x0colNum != 1 or bcolNum != 1:
		print('Error: \'b\' and \'x0\' must be a vector',
			'(n rows and 1 single column).')
		return None

	if not (ArowNum == AcolNum == browNum == x0rowNum):
		print('Error: matrix \'A\' must be square with n dimensions,',
			'and vectors \'b\' and \'x0\' must be size compatible (in R^n).')
		return None

	n = ArowNum
	L = minit(n, n)
	U = minit(n, n)

	for i in range(n):
		for j in range(i, n):
			L[j, i] = A[j, i]
			U[i, j] = A[i, j] if i !=j else 0.0

	try:
		Linv = L.I
	except:
		print('Error: \'L\' (lower triangular matrix) is singular!')
		return None

	xCur = copy.deepcopy(x0)
	i = 0
	err = epsilon * 2.0

	errorType = 'Relative' if relativeError else 'Infinite Norm'

	if showError:
		print('#        :', errorType, 'error:')

	while i < itMax and err > epsilon:
		i += 1
		xNew = Linv * (b - U * xCur)
		err = error(xNew, xCur, relative=relativeError)
		xCur = xNew
		if showError:
			print('{val:<{fill}}'.format(val=i, fill=8), ':', float(err))
	
	print('Process finished. ', end='')
	if err <= epsilon:
		print('(algorithm converged to a', errorType, 'error <= ' + str(epsilon) + ')')
	else:
		print('(max number of iterations reached: ' + str(itMax) + ')')

	return xCur

def __readfile__(filepath):
	values = []
	with open(filepath, 'r') as f:
		for lines in f:
			values.append(lines.split())
	return np.matrix(values, dtype=float)
			

def initProblemMatrices(n=10, Afp=None, bfp=None, x0fp=None):
	if Afp and bfp:
		A = __readfile__(Afp)
		b = __readfile__(bfp)
	else:
		# Default A: pentadiagonal matrix defined on problem specification
		# Default b: 1.0/i for i = 1,..., n 
		A = [[0.0] * n for i in range(n)]
		for i in range(n-3):
			A[i][i] = 4.0
			A[i][i+3] = -1.0
			A[i+3][i] = -1.0
			A[i][i+1] = -1.0
			A[i+1][i] = -1.0
		for i in range(n-3, n-1):
			A[i][i] = 4.0
			A[i+1][i] = -1.0
			A[i][i+1] = -1.0
		A[n-1][n-1] = 4.0
		A = np.matrix(A)
		b = np.matrix([[1.0/(1.0 + i)] for i in range(n)]) 
	if x0fp:
		x0 = __readfile__(x0fp) 
	else:
		# Default x0: zero vector in R^n
		x0 = np.matrix([[0.0] for i in range(n)]) 
	return A, b, x0

def __showHelp__():
	print('Gauss-Seidel Method for Linear System Solving by Felipe Alves Siqueira')
	print('Program usage:', sys.argv[0], 
		'[-a matrixFilepath] [-b vectorFilepath] [-x x0Filepath]',
		'[-n matrixSize] [-i itMax] [-e epsilon]')
	print('Options:',
		'-a:\tfull filepath of the coefficient matrix \'A\'. Default is the Pentadiagonal Matrix defined in problem specification.',
		'-b:\tfull filepath of the solution vector \'b\'. Default is 1.0/i for i= {1, ..., n}',
		'-x:\tinitial approximation of Gauss-Seidel solution. Default is the zero vector in R^n.',
		'-n:\tdimension of \'A\', \'b\' and \'x0\' for default values when filepaths was not specified. Default is 10.',
		'-i:\tmaximum iteration of Gauss-Seidel method. Default is 100.',
		'-e:\tmaximum error on a single Gauss-Seidel iteration. Default is 1e^-5 as \'Infinite Norm error\' (||a - b||[inf]).', 
		'-s:\tshow error each iteration of Gauss-Seidel method.',
		'-r:\tset error to be \'Relative Error\' (||a - b||[inf]/||a||[inf]) instead of \'Infinite Norm Error\' (||a - b||[inf]).',
		sep='\n')

def translateArgs(argv):
	param = {
		'-n': 10, 
		'-i': 100, 
		'-e': 1e-5,  
		'-a': None, 
		'-x': None, 
		'-b': None, 
		'-h': False, 
		'--help': False, 
		'-s': False,
		'-r': False
	}

	i = 1
	while i < len(argv):
		if argv[i] in param:
			if type(param[argv[i]]) != bool:
				param[argv[i]] = argv[i + 1]
				i += 1
			else:
				param[argv[i]] = True
		else:
			print('Warning: unknown option ', argv[i], 
				', ignoring it.', sep='\'')
			i += 1
		i += 1

	if param['-h'] or param['--help']:
		__showHelp__()
		return [None] * 8
	else:
		Afp = param['-a']
		bfp = param['-b']
		x0fp = param['-x']
		n = param['-n']
		epsilon = param['-e']
		itMax = param['-i']
		showError = param['-s']
		relativeError = param['-r']
		try:
			n = int(n)
			if n <= 0:
				raise Exception
		except:
			print('Warning: -n option must be a natural number.')

		try:
			epsilon = float(epsilon)
			if epsilon < 0:
				raise Exception
		except:
			print('Warning: -f option must be a non-negative real number.')

		try:
			itMax = int(itMax)
			if itMax <= 0:
				raise Exception
		except:
			print('Warning: -i option must be a natural number.' )

		try:
			if (Afp and not bfp) or (not Afp and bfp):
				raise Exception
		except:
			print('Warning: you must specify options -a (coefficients matrix)',
				'and -b (solution vector) together.')

		return n, epsilon, itMax, Afp, bfp, x0fp, showError, relativeError

if __name__ == '__main__':
	n, epsilon, itMax, Afp, bfp, x0fp, showError, relativeError = translateArgs(sys.argv)
	if n or Afp:
		A, b, x0 = initProblemMatrices(n, Afp, bfp, x0fp)
		if len(A) and len(b) and len(x0):
			r = gaussSeidel(A, b, x0, itMax, epsilon, showError, relativeError)
			if type(r) is np.matrix:
				print(r)

