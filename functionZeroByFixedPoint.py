"""
This script receives a f and g fuctions, where f(x) = 0 <=> g(x) = x, alongside a starting point x0,
and tries to use the Iterative Fixed Point Theorem to find a root z of f which is, by definition, a Fixed Point for g.
"""

from rpy2 import robjects as ro
import sys

def funcFixedPoint(g, a, it=100):
	xVal = a
	for i in range(it):
		xVal = ro.r(g.replace('x', str(xVal)))[0]
	return xVal

if __name__ == '__main__':
	if not (4 < len(sys.argv) < 7):
		print('usage:', sys.argv[0], '<f function> <g function> <x0>' 
			+ ' <iteration number (optional)>')
		exit(1)

	f = sys.argv[1]
	g = sys.argv[2]
	a = float(sys.argv[3])

	itNum = 100
	if len(sys.argv) == 5:
		itNum = int(sys.argv[4])

	z = funcFixedPoint(g, a, itNum)
	print('f Root/g Fixed Point = z =', z)
	print('f(z) =', ro.r(f.replace('x', str(z)))[0])
