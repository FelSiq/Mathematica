import sympy
import sys
import copy

def readMatrix(filepath):
	m = []
	with open(filepath) as f:
	    for l in f:
		    m.append(list(map(float, l.split())))
	return m

def solveTransitionMatrix(m, startProbs):
	system = copy.deepcopy(m)

	# Move the mi to the 'right' and the 1.0 step to the 'left'
	for i in range(len(system)):
		if system[i][i] != 1.0:
			system[i][i] -= 1.0
			system[i] += [-1.0]
		else:
			system[i] += [0.0]
	steps = (sympy.Matrix(system).rref()[0]).col(-1)

	E = sum([steps[i] * startProbs[i] for i in range(len(m))])

	return E

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print('usage:', sys.argv[0], '<T matrix> <(optional) start state prob (p1 ... pn)>')
		exit(1)
	m = readMatrix(sys.argv[1])
	startProbs = []
	if len(sys.argv) >= 3:
		startProbs = list(map(float, sys.argv[2:]))
		if abs(sum(startProbs) - 1.0) > 1.0e-5:
			print('Sum of starting probabilities of states must be 1.0.')
			exit(2)
	else:
		aux = 1.0/len(m)
		print('Initing all starting probabilities \'P(X0=i)\' with equal values ('+str(aux)+')')
		startProbs = [aux for p in range(len(m))]
	E = solveTransitionMatrix(m, startProbs)
	print('Avg # of steps til absortion:', E)
