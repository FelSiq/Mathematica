def mtransp(m):
	mRow = len(m)
	mCol = len(m[0])

	t = minit(mCol, mRow)
	for i in range(mRow):
		for j in range(mCol):
			t[j][i] = m[i][j]

	return t

def minit(rows, cols, n = 0):
	return [[n] * cols for i in range(rows)]

def midentity(order = 1):
	if order <= 0:
		print('E: Identity matrix order should be at least 1.')
		return None
	m = minit(order, order)
	for i in range(order):
		m[i][i] = 1
	return m

def mmult(m1, m2):
	m1Col = len(m1[0])
	m2Col = len(m2[0])
	m1Row = len(m1)
	m2Row = len(m2)
	
	if m1Col != m2Row:
		print('E: can\'t multiply a matrix with', m1Col, 'columns with another one with', m2Row, 'rows.')
	m3 = minit(m1Row, m2Col)

	for i in range(m1Row):
		for j in range(m2Col):
			for k in range(m1Col):
				m3[i][j] += m1[i][k] * m2[k][j]
	return m3

def mcopy(m):
	mcopy = list()
	for i in m:
		mcopy.append([k for k in i])
	return mcopy

def mrowswap(m, row1, row2):
	if not (0 <= row1 < len(m) and 0 <= row2 < len(m)):
		print('E: index for row exchange invalid.')
		return None

	mSwap = midentity(len(m))
	mCopy = mcopy(m)
	mSwap[row1], mSwap[row2] = mSwap[row2], mSwap[row1]
	
	return mmult(mSwap, mCopy)

def mprint(m):
	for i in m:
		for j in i:
			print('{value: <{space}}'.format(value = round(j, 6), space = 15), end = '')
		print()