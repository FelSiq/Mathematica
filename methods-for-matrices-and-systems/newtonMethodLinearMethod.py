from rpy2 import robjects as ro
import regex
import sys

"""
This script uses the Newton's method for Linear System Solving.
"""

def solveEquations(eq, vals):
	a=[]
	for e in eq:
		curEq = e
		counter = 1
		for v in vals:
			curEq = regex.sub(r'\bx'+str(counter)+r'\b', str(v), curEq)
			counter += 1
		a.append(ro.r(curEq)[0])
	return a

def __jacobianMatrix__(eq):
	ro.r('library(Deriv)')
	jacobian = []
	n = len(eq)
	for i in range(n):
		jacobian.append([])
		for m in range(n):
			jacobian[i].append(ro.r('Deriv(\''+str(eq[i])+'\', x=\'x'+str(m+1)+'\')')[0])
	return jacobian

def newtonMethod(x0, eq, itNum=100):
	jacobi = __jacobianMatrix__(eq)
	xk = x0
	for i in range(itNum):
		jacobiXk = solveEquations(jacobi, xk)
		fXk = solveEquations(eq, xk)
		# Calculate dk
		xk = xk + dk
	return xk

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

	r = newtonMethod(sys.argv[1], n, sys.argv[3:])
	print('result:', r)

