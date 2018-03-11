import random
import math
import sys

def coprime(a, b):
	for i in range(2, min(math.floor(a/2), math.floor(b/2))):
		if a % i == 0 and b % i == 0:
			return False
	return True

def pi(it=5000, maxval=1e+5):
	cpFreq = 0
	for i in range(it):
		a = random.randint(0, maxval)
		b = random.randint(0, maxval)
		cpFreq += coprime(a, b)

	pi = (6.0 * it/cpFreq)**0.5
	return pi

if __name__ == '__main__':
	print(pi())
