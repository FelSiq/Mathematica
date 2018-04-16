import sympy
import sys
import numpy as np
import copy

def readMatrix(filepath):
    m = []
    with open(filepath) as f:
            for l in f:
                    m.append(list(map(float, l.split())))
    return m

def solveTransitionMatrix(m, start, stop):
    system = copy.deepcopy(m)
    E = 0.0
    P = 0.0

    # Move the mi to the 'right' and the 1.0 step to the 'left'
    for i in range(len(system)):
        system[i][i] -= 1.0
        system[i] += [-1.0]
    steps = (sympy.Matrix(system).rref()[0]).col(-1)

    print(steps)

    E = steps[start]

    return {'Expected # steps': E, 'Probability P(T|X0='+str(start)+')': P}

if __name__ == '__main__':
    if len(sys.argv) < 3:
            print('usage:', sys.argv[0], '<T matrix> <init state> <stop states (x1 ... xn)>')
            exit(1)
    m = readMatrix(sys.argv[1])
    s = int(sys.argv[2])
    r = solveTransitionMatrix(m, s, sys.argv[3:])
    print(r)
