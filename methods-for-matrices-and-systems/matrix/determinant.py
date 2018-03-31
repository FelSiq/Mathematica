from matrix import rref

def reducedDet(m):
	if len(m) != len(m[0]):
		print('Matrix should be square to have a properly determinant.')
		return None

	reduced, elim = rref.rref(m, reduced = False)

	det = 1.0
	for i in range(len(reduced)):
		det *= reduced[i][i]
	return det

if __name__ == '__main__':
	mat = [
		[1, 4, 5, 9, -29, 27],
		[2, -1, 1, 1, 23, -1],
		[3, 0, 7, -2, 17, 0],
		[1, 2, 1, -1, 0, 0],
		[-30, 20, -1, -1, 17, 1],
		[9, 9, -1, -2, -3, 0],
	]

	print(reducedDet(mat))
