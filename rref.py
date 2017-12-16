with open('matrix.py', 'r') as f:
	exec(f.read())

def rref(m):
	# if pivot, in module, is smaller than DELTA, assume it's zero
	DELTA = 1.0e-10

	mCol = len(m[0])
	mRow = len(m)

	mCopy = mcopy(m)
	mAug = midentity(mRow)

	shift = 0
	i = 0
	while i < mRow and i + shift < mCol:
		pivot = mCopy[i][i + shift]

		# Check for row exchanges
		p = i + 1
		while abs(pivot) < DELTA and p < mRow:
			mCopy = mrowswap(mCopy, i, p)
			mAug = mrowswap(mAug, i, p)
			pivot = mCopy[i][i + shift]
			p += 1

		# If pivot is still zero, after all possible row exchanges, 
		# step to the next column and repeat the algorithm, if possible
		if abs(pivot) < DELTA:
			shift += 1
		else:
			# Just some row operations with the augmented matrix
			for j in range(mRow):
				if j != i:
					multiplier = mCopy[j][i + shift]/pivot
					for k in range(mCol):
						mCopy[j][k] -= mCopy[i][k] * multiplier 
					for k in range(mRow):
						mAug[j][k] -= mAug[i][k] * multiplier
			# Set the pivots to 1
			for k in range(mCol):
				mCopy[i][k] /= pivot
			for k in range(mRow):
				mAug[i][k] /= pivot
			i += 1

	return mCopy, mAug

# Program drive
"""if __name__ == '__main__':
	mat = [
		[1, 0, 2, 7],
		[2, -1, -2, 7],
		[3, 0, 6, 6],
	]

	reducted, elimination = rref(mat)
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
	mprint(mmult(redInverse, reducted))"""
