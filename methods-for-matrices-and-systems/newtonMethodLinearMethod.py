from matrix import rref
from rpy2 import robjects as ro
import regex
import sys

"""
This script uses the Newton's method for Linear System Solving.
"""

def solveEquations(m, vals):
	a=[]
	for i in range(len(m)):
		a.append([])
		for j in range(len(m[i])):
			curEq = m[i][j]
			counter = 1
			for v in vals:
				curEq = regex.sub(r'\bx'+str(counter)+r'\b', str(v), curEq)
				counter += 1
			a[i].append(ro.r(curEq)[0])
	return a

def __jacobianMatrix__(eq):
	ro.r('library(Deriv)')
	jacobian = []
	n = len(eq)
	for i in range(n):
		jacobian.append([])
		for m in range(n):
			jacobian[i].append(ro.r('Deriv(\''+str(eq[i][0])+'\', x=\'x'+str(m+1)+'\')')[0])
	return jacobian

def newtonMethod(x0, eq, itNum=100, showIt=True):
	jacobi = __jacobianMatrix__(eq)
	xk = list(map(float, x0.split()))
	n = len(xk)
	if showIt:
		print('Started Newton\'s method...')
	for i in range(1,1+itNum):
		jacobiXk = solveEquations(jacobi, xk)
		fXk = solveEquations(eq, xk)
		fXk = [[-l[0]] for l in fXk]

		system = [jacobiXk[l] + fXk[l] for l in range(n)]

		dk = [q[n] for q in rref.rref(system, returnReduced=True)[0]]

		xk = [xk[l] + dk[l] for l in range(n)]

		if not i % 20:
			print('iteration:\t', i)
	if showIt:
		print('Process finished.')
	return xk, fXk

if __name__ == '__main__':
	if len(sys.argv) < 4:
		print('usage:', sys.argv[0], '<x0> <#Iteration> <equations...>')
		exit(1)

	try:
		n = int(sys.argv[2])
		if n < 0:
			raise Exception
	except:
		print('Error: # iterations must be a positive integer value.')
		exit(2)

	eqs = [[e] for e in sys.argv[3:]]
	x, y = newtonMethod(sys.argv[1], eqs, n)

	print('\nResults:\n-\tx^(k+1)=', x,'^ T\n-\tf(x^(k+1))=', [r[0] for r in y],'^ T')

