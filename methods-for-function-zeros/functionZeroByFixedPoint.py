"""
This script receives a f and g fuctions, where f(x) = 0 <=> g(x) = x, alongside a starting point x0,
and tries to use the Iterative Fixed Point Theorem to find a root z of f which is, by definition, a Fixed Point for g.
"""

from rpy2 import robjects as ro
import re
import sys

def funcFixedPoint(g, a, itMax=100, maxError=1.0e-10, var='x'):

	print('Started fixed point method...\n\tMax iteration #:', 
		itMax, '\n\tMax error value:', maxError, '\n')

	xVal = a
	findVar = re.compile(r'\b' + var + r'\b')

	itToPrint = max(1, int(itMax/10))
	i = 0
	err = maxError + 1.0
	while i < itMax and err > maxError:
		i += 1
		xPrev = xVal
		xVal = ro.r(findVar.sub('('+str(xVal)+')', g))[0]
		err = abs(xPrev - xVal)
		if i % itToPrint == 0:
			print(i, 'it error:', err)
	print('Process finished @ iteration', i, 'with error', err, '\n') 
	return xVal

if __name__ == '__main__':
	if not (3 < len(sys.argv) < 7):
		print('usage:', sys.argv[0], '<f function> <g function> <x0>' 
			+ ' <max Error (optional)> <iteration number (optional)>')
		exit(1)

	f = sys.argv[1]
	g = sys.argv[2]

	if not (f and  g):
		print('Error: invalid input functions.')
		exit(6)

	try:
		a = float(sys.argv[3])
	except:
		print('Error: \'x0\' argument should be a real number!')
		exit(2)

	itNum = 100
	if len(sys.argv) >= 6:
		try:
			itNum = int(sys.argv[5])
			if itNum <= 0:
				raise Exception
		except:
			print('Error: maximum iteration number should be a positive integer!')
			exit(3)

	maxError =1.0e-10
	if len(sys.argv) >= 5:
		try:
			maxError = float(sys.argv[4])
			if maxError < 0.0:
				raise Exception
		except:
			print('Error: maximum error should be a non-negative real number!')
			exit(4)

	z = funcFixedPoint(g, a, itMax=itNum, maxError=maxError)
	print('f Root/g Fixed Point = z =', z)
	print('f(z) =', ro.r(re.sub(r'\bx\b', '('+str(z)+')', f))[0])
