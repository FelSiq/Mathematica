from rpy2 import robjects as ro
import matplotlib.pyplot as plt
import numpy as np
import sys
import re

def calcLk(x, k):
    lk = ''
    for j in range(len(x)):
        if j != k:
            lk += '(x - ' + str(x[j]) + ')/(' + (str(x[k] - x[j])) + ')*'
    return lk[:-1]

def lagrangeInterpol(x, y):
    pol = ''
    for k in range(len(x)):
        pol += calcLk(x, k) + '*' + str(y[k]) + '+'
    return pol[:-1]

def plot(pol, xlim=(-10.0, 10.0), num=1000, points=(None, None)):
    x=[]
    y=[]
    inc=(xlim[1]-xlim[0])/num
    for i in range(num):
        val=i*inc+xlim[0]
        x.append(val)
        y.append(float(ro.r(re.sub(r'\bx\b', str(val), pol))[0]))

    plt.plot(x, y)
    if points is not None and points[0] is not None and points[1] is not None:
        plt.plot(points[0], points[1], 'o')
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('usage:', sys.argv[0], '<x values separated by',
		'spaces> <y values separated by spaces>',
		'[(optional) lagrange polynomial input]')
        exit(1)

    x = list(map(float, sys.argv[1].split(' ')))
    y = list(map(float, sys.argv[2].split(' ')))

    if len(x) != len(y):
        print('Lengths of x and y values must be the same.')
        exit(2)

    p = lagrangeInterpol(x, y)

    print('Lagrange polynomial:', p)

    if len(sys.argv) >= 4:
        res = ro.r(re.sub(r'\bx\b', sys.argv[3], p))[0]
        print('\nResult:', res)

    # Result plot
    intInc=(max(x)-min(x))*0.15
    plot(p, xlim=(min(x)-intInc, max(x)+intInc), points=(x, y))

