# This is the calculus approach for linear regression with Least Squares method
with open('rref.py', 'r') as f:
	exec(f.read())

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

	reduced, redutor = rref(system)
	
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

	print('best data separator', separator, 'coefficients are:')
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