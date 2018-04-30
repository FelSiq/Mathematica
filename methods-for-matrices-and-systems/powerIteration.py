import sys
import numpy as np

def readMatrix(filepath, sep=' '):
	m = []
	with open(filepath, 'r') as f:
		for line in f:
			m.append(list(map(float, line.split(sep))))
	return np.matrix(m)

def powerIteration(A, x0, itMax=1e3):
	x=x0
	r=0.0

	pm=x0.argmax()
	if x[pm, 0] != 1.0:
		print('Error: ||x0||[infinity] should be 1.0.')
		return None

	for i in range(itMax):
		y=A*x
		if i < itMax - 1:
			pm=y.argmax()
			x=y/y[pm, 0]
		else:
			r = y[pm, 0]
	return r

if __name__ == '__main__':
	if len(sys.argv) <= 3:
		print('usage:', sys.argv[0], '<A filepath> <x0 filepath> <# of iterations>')
		exit(1)	
	
	try:
		A = readMatrix(sys.argv[1])
		x0 = readMatrix(sys.argv[2])
	except:
		print('Error at reading input matrices.')
		exit(2)

	try:
		itMax = int(sys.argv[3])
		if (itMax <= 0):
			raise Exception
	except:
		print('Error: parameter \'itMax\' must be a positive integer!')
		exit(3)
	
	r = powerIteration(A, x0, itMax)
	print('Highest A eigenvalue absolute value approximation:', r)
	
