from enum import IntEnum
import numpy as np
import sys

"""
INFORMATION:
This program solve a Ax = b system using Gauss-Seidel numeric method for
Linear System solving. It's very similar to the Jacobi Method, but it's
faster on a single-core run (Jacobi Method's better with parallel programming).

Gauss-Seidel numeric method is a anytime algorithm, wich means that it
always have a answer, independently of how much time was wasted running
at til a certain execution point. Howhever, the more iterations the al-
gorithm perform, the more precise will be the output.
"""

# Main function of this program. Effectively runs the Gauss-Seidel iterative method
def gaussSeidel(A, b, x0, itMax=100, epsilon=1e-5, showError=False, relativeError=False):
	if type(A) != np.matrix:
		print('Error: expecting np.matrix as input matrix \'A\'.')
		return None
	
	ArowNum, AcolNum = A.shape
	browNum, bcolNum = b.shape
	x0rowNum, x0colNum = x0.shape

	# Check if x0 and b is a column matrix (or a vector)
	if x0colNum != 1 or bcolNum != 1:
		print('Error: \'b\' and \'x0\' must be a vector',
			'(n rows and 1 single column).')
		return None

	# Check if matrix A, b and x0 dimensions are compatible 
	# in between each other for a equation system. In this case,
	# A must be square and have the exactly same dimension of
	# b vector and x0 vector.
	if not (ArowNum == AcolNum == browNum == x0rowNum):
		print('Error: matrix \'A\' must be square with n dimensions,',
			'and vectors \'b\' and \'x0\' must be size compatible (in R^n).')
		return None

	n = ArowNum

	# Init L (Lower-triangular) and U (Upper-triangular matrices)
	L = minit(n, n)
	U = minit(n, n)

	# L is defined as a lower-triangular matrix with the same elements
	# inside and below the main diagonal of A matrix, i.e:
	#	[a11 	0	0	...	0  ]
	#	[a21	a22	0	...	0  ]
	# L =	[a31	a32	a33	...	0  ]
	#	[...	...	...	...	0  ]
	#	[an1	an2	an3	...	ann]
	#
	# U is defined as a upper-triangular matrix that have all elements
	# strictly above the main diagonal of A, which means that U = A - L.
	for i in range(n):
		for j in range(i, n):
			L[j, i] = A[j, i]
			U[i, j] = A[i, j] if i !=j else 0.0

	# Calculate the Inverse of L (it's needed to Gauss-Seidel method).
	try:
		Linv = L.I
	except:
		print('Error: \'L\' (lower triangular matrix) is singular!')
		return None

	# Init some necessary variables
	xCur = x0
	i = 0
	err = epsilon * 2.0 + 1.0

	# Stuff for user interface
	errorType = 'Relative' if relativeError else 'Infinite Norm'
	if showError:
		print('#        :', errorType, 'error:')

	# Gauss-Seidel numeric method implementation
	while i < itMax and err > epsilon:
		i += 1
		# Calculate a (theoretically, supposing the process will converge)
		# better x approximation based on the old x.
		xNew = Linv * (b - U * xCur)

		# Calculate the iteration error (defaults to Infinite Norm error,
		# but can be used the Relative Error instead. Check error function
		# for more detail.)
		err = error(xNew, xCur, relative=relativeError)

		# Update current x approximation.
		xCur = xNew

		# User interface stuff (print each iteration error if '-s' 
		# options is selected).
		if showError:
			print('{val:<{fill}}'.format(val=i, fill=8), ':', float(err))
		# Repeat till one of the stop criteria is reached.	

	# User interface statistics
	print('Process finished. ', end='')
	if err <= epsilon:
		print('(algorithm converged to a', errorType, 'error', 
			float(err), '<=', epsilon, 'at iteration #' + str(i) + ')')
	else:
		print('(max number of iterations reached: ' + str(itMax) + ')')

	# Return the last x approximation
	return xCur

# Calculates the Infinite Norm or Relative Error of a Gauss-Seidel iteration
def error(a, b, delta=1e-13, relative=False):
	denominator = (max(a) + delta) if relative else 1.0
	return max(a - b)/denominator

# Init a rownum x colnum matrix with all values equals to 'val'
def minit(rownum, colnum, val=0.0):
	return np.matrix([[val] * colnum for i in range(rownum)])

# Read matrix values from a given filepath
def __readfile__(filepath):
	values = []
	with open(filepath, 'r') as f:
		for lines in f:
			values.append(lines.split())
	return np.matrix(values, dtype=float)

# Gather all data needed to perform the official problem specifications experiments
def initProblemMatrices(n=10, Afp=None, bfp=None, x0fp=None):
	if Afp and bfp:
		print('Note: Reading user A (\''+Afp+'\') and b (\''+bfp+'\') matrices...')
		# In case the user specify a particular A and b matrices
		A = __readfile__(Afp)
		b = __readfile__(bfp)
	else:
		print('Note: Using predefined pentadiagonal A (with n = '+str(n)+') and standard b matrices.',
			'(Tip: change with \'-a\' and \'b\' options).')
		# Default A: pentadiagonal matrix defined on problem specification
		# Default b: 1.0/i for i = 1,..., n 
		A = [[0.0] * n for i in range(n)]
		if n > 1:
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
		print('Note: reading user initial approximation x0 (\''+fx0+'\')...')
		# In case user specify a particular starting approximation of x
		x0 = __readfile__(x0fp) 
	else:
		print('Note: using null vector as initial approximation x0 (Tip: use \'-x\' option to change).')
		bnRow, bnCol = b.shape
		# Default x0: zero vector in R^n
		x0 = np.matrix([[0.0] for i in range(bnRow)]) 
	return A, b, x0

# Show problem help (used when user call the program using '-h' or '--help' options)
def __showHelp__():
	print('Gauss-Seidel Method for Linear System Solving by Felipe Alves Siqueira')
	print('Program usage:', sys.argv[0], 
		'[-a matrixFilepath] [-b vectorFilepath] [-x x0Filepath]',
		'[-n matrixSize] [-i itMax] [-e epsilon] [-s (show error each iteration)] [-r (use Relative Error)]')
	print('Options:',
		'-a:\tfull filepath of the coefficient matrix \'A\'. Default is the Pentadiagonal Matrix defined in problem specification.',
		'-b:\tfull filepath of the solution vector \'b\'. Default is 1.0/i for i= {1, ..., n}',
		'-x:\tinitial approximation of Gauss-Seidel solution. Default is the zero vector in R^n.',
		'-n:\tdimension of \'A\', \'b\' and \'x0\' for default values when filepaths was not specified. Default is 10.',
		'-i:\tmaximum iteration of Gauss-Seidel method. Default is 100.',
		'-e:\tmaximum error on a single Gauss-Seidel iteration. Default is 1e^-10 as \'Infinite Norm error\' (||a - b||[inf]).', 
		'-s:\tshow error each iteration of Gauss-Seidel method.',
		'-r:\tset error to be \'Relative Error\' (||a - b||[inf]/||a||[inf]) instead of \'Infinite Norm Error\' (||a - b||[inf]).',
		sep='\n')

# Translate user command line arguments into real program parameters
def translateArgs(argv):
	param = {
		'-n': 10, 
		'-i': 100, 
		'-e': 1e-10,  
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

# Program driver
if __name__ == '__main__':
	# Get needed parameters based on user command line arguments
	n, epsilon, itMax, Afp, bfp, x0fp, showError, relativeError = translateArgs(sys.argv)
	if n or Afp:
		# User interface messages
		print('Parameters choosen:\n',
			'\tEpsilon =', epsilon, '\n',
			'\tError =', 'Relative' if relativeError else 'Infinite Norm', '\n',
			'\tMax Iterations =', itMax)

		# Fill A, b and x0 matrices with properly data
		A, b, x0 = initProblemMatrices(n, Afp, bfp, x0fp)
		if len(A) and len(b) and len(x0):
			# Runs Gauss-Seidel numeric method, trying to solve Ax=b, if possible
			print('\nStarted Gauss-Seidel method...')
			r = gaussSeidel(A, b, x0, itMax, epsilon, showError, relativeError)
			if type(r) is np.matrix:
				# Print problem answer.
				print(r)

