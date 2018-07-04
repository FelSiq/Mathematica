import matplotlib.pyplot as plt
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
		approxFun += '(' + str(coeffs[i]) + ')*(' + f[i] + ') + '

	return approxFun[:-3]

def plot(fun, xlim=(-10.0, 10.0), num=1000, points=(None, None)):
    x=[]
    y=[]
    inc=(xlim[1]-xlim[0])/num
    for i in range(num):
        val=i*inc+xlim[0]
        x.append(val)
        y.append(float(ro.r(re.sub(r'\bx\b', str(val), fun))[0]))

    plt.plot(x, y)
    if points is not None and points[0] is not None and points[1] is not None:
        plt.plot(points[0], points[1], 'o')
    plt.show()

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
	f = [re.sub(r'\bx\b', '(x)', f) for f in sys.argv[3].split(' ')]

	approxFunc = mmq(x, y, f)
	print('Approximated function: ', approxFunc)

	if len(sys.argv) >= 5:
		vlist = list(map(float, sys.argv[4].split(' ')))
		print('Evaluating given values on the approximated function:')
		for val in vlist:
			print('f(' + str(val) + ')=', 
				ro.r(re.sub(r'\bx\b', str(val), approxFunc))[0])

	# Result plot
	intInc=(max(x)-min(x))*0.15
	plot(approxFunc, xlim=(min(x)-intInc, max(x)+intInc), points=(x, y))	
