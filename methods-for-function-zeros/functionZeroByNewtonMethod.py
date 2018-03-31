from rpy2 import robjects as ro
import re
import sys
"""
INFO:
This script uses Newton's Method to find a approximation of a function f zero (root).
It does not check any sequence convergence criteria.
"""


# a.k.a. "Tangent Method"
def newtonMethod(fun, x0, it=10, var=r'x'):
	xCur = x0
	ro.r('library(Deriv)')
	funDev = ro.r('Deriv(\'' + fun + '\')')[0]
	
	subVarRegex = re.compile(r'\b' + var + r'\b')

	for i in range(it):
		funVal = ro.r(subVarRegex.sub(str(xCur), fun))[0]
		derivVal = ro.r(subVarRegex.sub(str(xCur), funDev))[0]
		xCur -= funVal/derivVal 
	return xCur

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print('usage:', sys.argv[0], '<function> <x0> <itNum (optional)>')
		exit(1)

	itNum = 10
	if len(sys.argv) == 4:
		itNum = int(sys.argv[3])

	approxZero = newtonMethod(sys.argv[1], float(sys.argv[2]), itNum)

	print('~z:', approxZero)
	print('f(~z):', ro.r(re.sub(r'\bx\b', str(approxZero), sys.argv[1]))[0])
