"""
This program uses the Laplace Algorithm to calculate a square Matrix Determinant.
It is a recursive algorithm. It's not very efficient, but still a interesting idea.
"""
def laplacedet(mat):
	if not len(mat):
		return 1

	result = 0.0
	for i in range(len(mat)):
		matcopy = list()
		for j in range(len(mat)):
			if j != i:
				matcopy.append(mat[j][1:])
		result += (-1) ** i  * mat[i][0] * laplacedet(matcopy)
	return result

if __name__ == '__main__':
	mat = [
		[2, 0, 15, -19, 0],
		[3, -3, -68, -1, 11],
		[0, -2.3, -1.234, 19, -5],
		[9, -1, -99, 0, -2],
		[1, 2, 1, -3, -1]
	]
	print(laplacedet(mat))
