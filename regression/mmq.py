import matplotlib.pyplot as plt
from rpy2 import robjects as ro
import numpy as np
import sys
import re

def genPhiVectors(x, f, var="x"):
	n = len(x)
	m = len(f)

	phiVectors = [[0.0] * n for i in range(m)]

	for i in range(m):
		for j in range(n):
			phiVectors[i][j] = ro.r(re.sub(r"\b" + var + r"\b",
				str(x[j]), f[i]))[0]

	return phiVectors

def genNormalEquationSystem(phiVectors, y, F):
	n = len(phiVectors)
	nes = np.array([[0.0] * n for i in range(n)])
	b = [0.0] * n

	Fvalues=list(ro.r(re.sub(r"\bx\b", \
		"c(" + ",".join(map(str, y)) + ")", F)))

	for i in range(n):
		b[i] = np.dot(Fvalues, phiVectors[i])
		for j in range(n):
			nes[i, j] = np.dot(phiVectors[i], 
				phiVectors[j])

	return nes, np.array(b)

def mmq(x, y, f, F, printSys=True):
	pv = genPhiVectors(x, f)
	nes, b = genNormalEquationSystem(pv, y, F)

	if printSys:
		print("Normal Equation System:\n", nes, sep="", end="\n\n")

	coeffs = np.linalg.solve(nes, b)
	approxFun = ""

	for i in range(len(coeffs)):
		approxFun += "(" + str(coeffs[i]) + ")*(" + f[i] + ") + "

	return approxFun[:-3]

def plot(fun, xlim=(-10.0, 10.0), num=1000, points=(None, None)):
    x=[]
    y=[]
    inc=(xlim[1]-xlim[0])/num
    for i in range(num):
        val=i*inc+xlim[0]
        x.append(val)
        y.append(float(ro.r(re.sub(r"\bx\b", str(val), fun))[0]))

    plt.plot(x, y)
    if points is not None and points[0] is not None and points[1] is not None:
        plt.plot(points[0], points[1], "o")
    plt.show()

if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("usage:", sys.argv[0], 
			"<x values separated by spaces>",
			"<y values separated by spaces>", 
			"<Base functions separated by spaces>",
			"<F function>",
			"[(optional) input values for approximated function]")
		exit(1)
	# k input funcions
	# m+1 points
	x = list(map(float, sys.argv[1].split(" ")))
	y = list(map(float, sys.argv[2].split(" ")))
	f = [re.sub(r"\bx\b", "(x)", f) for f in sys.argv[3].split(" ")]
	F = re.sub(r"\bx\b", "(x)", sys.argv[4])

	approxFunc = mmq(x, y, f, F)
	print("Approximated function: ", re.sub(r"\(x\)", "f(x)", F), "=", approxFunc)

	if len(sys.argv) >= 6:
		vlist = list(map(float, sys.argv[5].split(" ")))
		print("Evaluating given values on the approximated function:")
		for val in vlist:
			print(re.sub(r"\bx\b", str(val), F), "=", 
				ro.r(re.sub(r"\bx\b", str(val), approxFunc))[0])

	# Result plot
	intInc=(max(x)-min(x))*0.15


	y_t=list(ro.r(re.sub(r"\bx\b", "c("+",".join(map(str, y))+")", F)))
	plot(approxFunc, xlim=(min(x)-intInc, max(x)+intInc), points=(x, y_t))	
