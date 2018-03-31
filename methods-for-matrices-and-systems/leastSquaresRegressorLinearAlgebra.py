# This is the linear algebra approach for linear regression with Least Squares method
with open('rref.py', 'r') as f:
	exec(f.read())

def leastSquares(x, y):
	xTransp = mtransp(x)
	newSystem = mmult(xTransp, x)
	newLabels = mmult(xTransp, y)

	# Augmented matrix for the system
	for i in range(len(newSystem)):
		newSystem[i] += newLabels[i]

	reduced, redutor = rref(newSystem)
	
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
			y.append([newData[-1]])
			newData[-1] = 1.0
			x.append(newData)
		except:
			flag = False

	coeffs = leastSquares(x, y)

	printSolution(coeffs)