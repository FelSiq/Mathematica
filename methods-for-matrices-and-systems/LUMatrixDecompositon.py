with open('matrix.py') as f:
	exec(f.read())

# Using DoLittle technique l[i,i] = 1, for all i in [0, n], l in L 
def LU(M):
	n = len(M)
	m = len(M[0])

	if n != m:
		print('Error: expecting M to be a square matrix.')
		return None, None

	L = minit(n, n)
	U = minit(n, n)

	#DoLittle Technique
	for i in range(n):
		L[i][i] = 1.0

	for k in range(n):
		for j in range(k, n):
			c = 0
			for r in range(k):
				c += L[k][r] * U[r][j]
			U[k][j] = M[k][j] - c
		for i in range(k+1, n):
			c = 0
			for r in range(k):
				c += L[i][r] * U[r][k]
			L[i][k] = (M[i][k] - c)/U[k][k]
	return L, U
	

if __name__ == '__main__':
	M = [
		[ 2.0,  3.0, -1.0],
		[ 0.0,  3.0,  4.0],
		[ 2.0, -3.0,  8.0]
	]

	L, U = LU(M)
	
	print('Original Matrix (M):')
	mprint(M)
	print('L:')
	mprint(L)
	print('U:')
	mprint(U)
	print('LU - M (must be all equal to 0):')
	mprint(msub(mmult(L, U), M))
