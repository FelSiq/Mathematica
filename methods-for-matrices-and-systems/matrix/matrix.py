def missymmetric(m):
	nrow, ncol = mshape(m)

	if nrow != ncol:
		return False

	for i in range(1, nrow):
		for j in range(0, i):
			if m[i][j] != m[j][i]:
				return False
	return True

def mshape(m):
	nrow = len(m)
	ncol = 0
	if nrow:
		if type(m[0]) == list:
			ncol = len(m[0])
		else:
			# Then it's just a (transposed) vector
			ncol = nrow
			nrow = 1
	return (nrow, ncol)

def mtransp(m):
	mRow, mCol = mshape(m)

	t = minit(mCol, mRow)
	if mRow != 1 and mCol != 1:
		for i in range(mRow):
			for j in range(mCol):
				t[j][i] = m[i][j]
	elif mRow != 1 and mCol == 1:
		for i in range(mRow):
			t[i] = m[i][0]
	elif mRow == 1 and mCol != 1:
		for i in range(mCol):
			t[i][0] = m[i]
	else:
		t[0][0] = m[0][0]
	return t

def minit(rows, cols, n=0):
	if rows > 1:
		return [[n] * cols for i in range(rows)]
	return [n] * cols

def mread(filepath, sep=None):
	m = []
	with open(filepath, 'r') as f:
		for lines in f:
			data = list(map(float, lines.split(sep=sep)))
			m.append(data)
	return m		

def midentity(order=1):
	if order <= 0:
		print('E: Identity matrix order should be at least 1.')
		return None
	m = minit(order, order)
	for i in range(order):
		m[i][i] = 1
	return m

def mmult(m1, m2):
	m1Row, m1Col = mshape(m1)
	m2Row, m2Col = mshape(m2)
	
	if m1Col != m2Row:
		print('E: can\'t multiply a matrix with', m1Col, 'columns with another one with', m2Row, 'rows.')
	m3 = minit(m1Row, m2Col)

	for i in range(m1Row):
		for j in range(m2Col):
			for k in range(m1Col):
				m3[i][j] += m1[i][k] * m2[k][j]
	# If matrix resultant is a 1x1, then it was just a dot product and, then, return just the value
	return m3[0][0] if (m1Row == 1 and m2Col == 1) else m3

def msub(m1, m2):
	m1Row, m1Col = mshape(m1)
	m2Row, m2Col = mshape(m2)
	
	if m1Col != m1Col or m2Row != m2Row:
		print('E: can\'t subtract two matrices with different dimensions.')
	
	m3 = minit(m1Row, m2Col)
	for i in range(m1Row):
		for j in range(m2Col):
			m3[i][j] = m1[i][j] - m2[i][j]
	
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

def mmconst(m, const):
	mRow, mCol = mshape(m)

	n = minit(mRow, mCol)
	for i in range(mRow):
		for j in range(mCol):
			n[i][j] = m[i][j] * const

	return n

def mprint(m):
	if m:
		for i in m:
			if i:
				for j in i:
					print('{value: <{space}}'.format(value = round(j, 6) if type(j) != str else j, 
						space = 15), end = '')
				print()

