"""
Expected input: n vectors coordinates, each one in a column, separated by whitespaces
Example:
1 1 5
2 1 4
3 3 9

where v1 = [1 2 3]^T, v2 = [1 1 3]^T and v3 = [5 4 9]^T

General case (m vectors with n dimensions:
x1	y1	z1	..	k1
x2	y2	z2	..	k2
x3	y3	z3	..	k3
..	..	..	..	..
xn	yn	zn	..	kn

where v1 = [x1 x2 x3 .. xn]^T, ..., vm = [k1 k2 k3 .. kn]^T

"""
from matrix import matrix, rref

def checkIndenpendency(vecs):
	vCopy = matrix.mcopy(vecs)
	reducted, reductor = rref.rref(vCopy)
	
	dep = 0
	for r in reducted:
		dep += sum(r) == 0
	
	return dep == 0

def vecnorm(vec):
	norm = 0
	for i in vec:
		norm += i[0]**2
	return norm**0.5

def gso(vecs):
	# 1. Check indenpendency two-by-two
	colVecs = matrix.mtransp(vecs)

	if not checkIndenpendency(colVecs):
		print('E: These vectors aren\'t indenpendent.')
		return None

	for i in range(len(colVecs)):
		for j in range(len(vecs)):
			colVecs[i][j] = [colVecs[i][j]]

	newVectors = list()
	# 2. Apply Gram process with Schmidt operation
	for i in range(len(colVecs)):
		newBaseVec = [[k[0]] for k in colVecs[i]]
		for j in range(i):
			coeff1 = matrix.mmult(matrix.mtransp(newVectors[j]), newBaseVec)
			coeff2 = matrix.mmult(matrix.mtransp(newVectors[j]), newVectors[j])
			curProj = matrix.mmconst(newVectors[j], coeff1/coeff2)
			newBaseVec = matrix.msub(newBaseVec, curProj)

		newVectors.append(matrix.mmconst(newBaseVec, 1.0/vecnorm(newBaseVec)))

	return newVectors

if __name__ == '__main__':
	vecs = list()
	flag = True
	while flag:
		try:
			vecs.append(list(map(float, input().split(' '))))
		except:
			flag = False
	newvecs = gso(vecs)

	if newvecs:
		letter = ord('a')
		for v in newvecs:
			print(chr(letter) + ':', v)
			letter += 1

		# Testing if they're orthogonal two by two
		# test: if i == j then != 0, 0 otherwise
		# In other words, this must print a diagonal matrix.
		# Because the vectors should be unitary, then this diag.
		# matrix must be the Identity Matrix.
		print('\nTruth proof (should always be the identity matrix):')
		nVecs = len(newvecs)
		resultMat = matrix.minit(nVecs, nVecs)
		for i in range(nVecs):
			for j in range(nVecs):
				resultMat[i][j] = matrix.mmult(matrix.mtransp(newvecs[i]), newvecs[j])

		matrix.mprint(resultMat)
