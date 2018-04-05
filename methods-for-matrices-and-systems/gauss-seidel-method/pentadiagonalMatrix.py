import sys
"""
This subscript generates a pentadiagonal matrix, used for test case input
of Gauss-Seidel numeric method.
"""

if __name__ == '__main__':
	n = 10
	if len(sys.argv) > 1:
		try:
			n = int(sys.argv[1])
			if n <= 0:
				n = 10
				raise Exception
		except:
			print('Warning: matrix dimension \'n\' must be natural value.')

	A = [[0.0] * n for i in range(n)]
	if n > 1:
		for i in range(n-3):
			A[i][i] = 4.0
			A[i][i+3] = -1.0
			A[i+3][i] = -1.0
			A[i][i+1] = -1.0
			A[i+1][i] = -1.0
		for i in range(n-3, n-1):
			A[i][i] = 4.0
			A[i+1][i] = -1.0
			A[i][i+1] = -1.0
	A[n-1][n-1] = 4.0
	
	for n in A:
		for m in n:
			print(m, end=' ')
		print()

