from rpy2 import robjects as ro
import numpy as np
import sys
import re

def genPhiVectors(x, f, var='x'):
	n = len(x)
	m = len(f)

	phiVectors = [[0.0] * n for i in range(m)]

	for i in range(m):
		for j in range(n):
			phiVectors[i][j] = ro.r(re.sub(r'\b' + var + r'\b',
				str(x[j]), f[i]))[0]

	return phiVectors

def genNormalEquationSystem(phiVectors, y):
	n = len(phiVectors)
	nes = np.array([[0.0] * n for i in range(n)])
	b = [0.0] * n

	for i in range(n):
		b[i] = np.dot(y, phiVectors[i])
		for j in range(n):
			nes[i, j] = np.dot(phiVectors[i], 
				phiVectors[j])
	return nes, np.array(b)

def mmq(x, y, f):
	pv = genPhiVectors(x, f)
	nes, b = genNormalEquationSystem(pv, y)
	coeffs = np.linalg.solve(nes, b)
	approxFun = ''
	for i in range(len(coeffs)):
		approxFun += '(' + str(coeffs[i]) + ')*' + f[i] + ' + '
	return approxFun[:-2]

if __name__ == '__main__':
	if len(sys.argv) < 4:
		print('usage:', sys.argv[0], 
			'<x values separated by spaces>',
			'<y values separated by spaces>', 
			'<Base functions separated by spaces>',
			'[(optional) input values for approximated function]')
		exit(1)
	# k input funcions
	# m+1 points
	x = list(map(float, sys.argv[1].split(' ')))
	y = list(map(float, sys.argv[2].split(' ')))
	f = sys.argv[3].split(' ')

	approxFunc = mmq(x, y, f)
	print('Approximated function: ', approxFunc)

	if len(sys.argv) >= 5:
		vlist = list(map(float, sys.argv[4].split(' ')))
		print('Evaluating given values on the approximated function:')
		for val in vlist:
			print('f(' + str(val) + ')=', 
				ro.r(re.sub(r'\bx\b', str(val), approxFunc))[0])
