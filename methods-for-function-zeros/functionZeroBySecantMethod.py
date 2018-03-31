from rpy2 import robjects as ro
import re
import sys
"""
INFO:
This script uses Secant Method to find a approximation of a function f zero (root).
It does not check any sequence convergence criteria.
"""

def secantMethod(fun, a, b, it=10, var=r'x'):
	xA = a
	xB = b

	subVarRegex = re.compile(r'\b' + var + r'\b')

	i = 0
	while i < it:
		i += 1
		funValA = ro.r(subVarRegex.sub('('+str(xA)+')', fun))[0]
		funValB = ro.r(subVarRegex.sub('('+str(xB)+')', fun))[0]
		if funValA - funValB:
			aux = xA
			xA -= funValA * (xA - xB)/(funValA - funValB)
			xB = aux
		else:
			i = it
	return xA

if __name__ == '__main__':
	if len(sys.argv) < 4:
		print('usage:', sys.argv[0], '<function> <a> <b> <itNum (optional)>')
		exit(1)

	itNum = 10
	if len(sys.argv) == 5:
		itNum = int(sys.argv[4])

	approxZero = secantMethod(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), itNum)

	print('~z:', approxZero)
	print('f(~z):', ro.r(re.sub(r'\bx\b', '('+str(approxZero)+')', sys.argv[1]))[0])
