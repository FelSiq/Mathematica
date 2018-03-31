# This is the linear algebra approach for linear regression with Least Squares method
from matrix import rref, matrix

"""
Least Squase implementation based on the Projection idea:

if Ax = b don't have a solution (more equations than variables), then
(A^T)Ax = (A^T)b <=> x = (A^T * A)^-1 * A^T * b is a good approximation
for x.

Sample input (7 data points at R^5):
1 2 3 4 5
1 6 7 8 9
2 3 4 5 6
1 2 3 4 1
1 9 0 8 7
1 2 3 4 5
6 4 3 1 2

Sample output:
best data descriptor hyperplane coefficients are:
a: 1.0
b: -0.561
c: 0.045
d: 1.848
e: -3.742

"""

def leastSquares(x, y):
	xTransp = matrix.mtransp(x)
	newSystem = matrix.mmult(xTransp, x)
	newLabels = matrix.mmult(xTransp, y)

	# Augmented matrix for the system
	for i in range(len(newSystem)):
		newSystem[i] += newLabels[i]

	reduced, redutor = rref.rref(newSystem)
	
	coeffs = list()
	# Create solutions
	for s in reduced:
		coeffs.append(s[-1])
	return coeffs

def printSolution(coeffs):
	dataDimension = len(coeffs) - 1
	if dataDimension <= 0:
		print('No geometric data descriptor found (no dimensions available)')
	else:
		separator = ''
		if dataDimension == 1:
			separator = 'line'
		elif dataDimension == 2:
			separator = 'plane'
		else:
			separator = 'hyperplane'

		print('best data descriptor', separator, 'coefficients are:')
		letter = ord('a')
		for i in coeffs:
			print (chr(letter) + ':', round(i, 3))
			letter += 1

# Program driver
if __name__ == '__main__':
	x = list()
	y = list()
	flag = True
	while flag:
		try:
			newData = list(map(float, input().split(' ')))
			y.append([newData[-1]])
			newData[-1] = 1.0
			x.append(newData)
		except:
			flag = False

	coeffs = leastSquares(x, y)

	printSolution(coeffs)
