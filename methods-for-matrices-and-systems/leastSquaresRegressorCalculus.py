"""

Expected output: 
n data points (x0, ..., xm) with m dimension, one at each input line.

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

# This is the calculus approach for linear regression with Least Squares method
from matrix import rref

# This solution uses partial derivatives concept in order do minimize the
# squared error sum.
def leastSquares(x, y):
	system = list()

	varNum = len(x[0])
	# For eeach variable
	for i in range(varNum):
		equation = [0] * (varNum + 1)
		# for each equation
		for j in range(len(x)):
			derivativeCoef = x[j][i]
			for k in range(varNum):
				equation[k] += derivativeCoef * x[j][k] 
			equation[-1] += x[j][i] * y[j]
		system.append(equation)

	reduced, redutor = rref.rref(system)
	
	coeffs = list()
	# Create solutions
	for s in reduced:
		coeffs.append(s[-1])
	return coeffs

def printSolution(coeffs):
	dataDimension = len(coeffs) - 1
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
			y.append(newData[-1])
			newData[-1] = 1.0
			x.append(newData)
		except:
			flag = False

	coeffs = leastSquares(x, y)

	printSolution(coeffs)
