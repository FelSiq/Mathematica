with open('matrix.py', 'r') as f:
	exec(f.read())

def rref(m, reduced = True, returnReduced = True):
	# if pivot, in module, is smaller than DELTA, assume it's zero
	DELTA = 1.0e-10

	mCol = len(m[0])
	mRow = len(m)

	mReduced = mcopy(m)
	matReverse = midentity(mRow)

	shift = 0
	i = 0
	while i < mRow and i + shift < mCol:
		pivot = mReduced[i][i + shift]

		# Check for row exchanges
		p = i + 1
		while abs(pivot) < DELTA and p < mRow:
			mReduced = mrowswap(mReduced, i, p)
			matReverse = mrowswap(matReverse, i, p)
			pivot = mReduced[i][i + shift]
			p += 1

		# If pivot is still zero, after all possible row exchanges, 
		# step to the next column and repeat the algorithm, if possible
		if abs(pivot) < DELTA:
			shift += 1
		else:
			# Just some row operations with the augmented matrix
			for j in range(i * (not reduced), mRow):
				if j != i:
					multiplier = mReduced[j][i + shift]/pivot
					for k in range(mCol):
						mReduced[j][k] -= mReduced[i][k] * multiplier 
					for k in range(mRow):
						matReverse[j][k] -= matReverse[i][k] * multiplier
			# Set the pivots to 1
			if reduced:
				for k in range(mCol):
					mReduced[i][k] /= pivot
				for k in range(mRow):
					matReverse[i][k] /= pivot
			i += 1

	if (returnReduced):
		return mReduced, matReverse
	return matReverse

# Program drive
"""
if __name__ == '__main__':
	mat = [
		[1, 4, 5],
		[2, -1, 1],
		[3, 0, 7],
	]


	# If reduced = False, the the output should be the Upper Triangular Echelon Form Matrix.
	# Otherwise, it will be the Row Reduced Echelon Form diagonal (if square) matrix.
	reducted, elimination = rref(mat, reduced = True)
	print('Matrix:')
	mprint(mat)
	print('\n\nReducted (should be identity for reversible square matrix):')
	mprint(reducted)
	print('\n\nElimination (should be inverse for reversible square matrix):')
	mprint(elimination)
	print('\n\nElimination x Mat == Reducted?:')
	mprint(mmult(elimination, mat))

	identity, redInverse = rref(elimination)

	print('\n\nShould be identity:')
	mprint(identity)
	print('\n\nElimination reversed:')
	mprint(redInverse)

	print('\n\nElimination^-1 x Reducted == Matrix? :')
	mprint(mmult(redInverse, reducted))
"""
