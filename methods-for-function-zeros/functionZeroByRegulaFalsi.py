from rpy2 import robjects as ro
import regex
import sys

"""
This method is a variation of Bissection Method. Differently from the original
Bissection Method (which only considers signal(f(a)*f(b) each iteration), 
the Regula Falsi (or False Position)  does consider the values of f(a) and f(b) 
when calculating the f((a+b)/2), i.e, it calculates f((a*f(a) + b*f(b))/(f(a) + (f(b)) 
which is, of course, an weighted mean between a and b.

The motivation behind this is the assumption that, in f function, continuous in [a,b] interval,
if exists a z real that f(z) = 0, then z is closer to |f(m)| than |f(n)| if and only if
|f(m)| < |f(n)|, for all m, n in [a, b].
"""

def bissection(function, a, b, n=100, var='x'):
	subVarRe = regex.compile(r'\b' + var + r'\b')
	for i in range(n):
		fa = ro.r(subVarRe.sub('('+str(a)+')', function))[0] 
		fb = ro.r(subVarRe.sub('('+str(b)+')', function))[0]
		if fb != fa: 
			curXMean = (a*fb - b*fa)/(fb - fa) 
			curImgVal = ro.r(subVarRe.sub('('+str(curXMean)+')', function))[0]
			if curImgVal * fa < 0:
				b = curXMean
			else:
				a = curXMean
		else:
			break;

	return curXMean

if __name__ == '__main__':
	if len(sys.argv) < 4:
		print('usage:', sys.argv[0], '<function> <a> <b> (interval is [a, b]) <# of iterations? (optional)>')
		exit(1)
	function = sys.argv[1]
	try:
		a = float(sys.argv[2])
		b = float(sys.argv[3])
	except:
		print('Error: \'a\' and \'b\' args must be a real number.')
		exit(2)
	try:
		n = 100
		if len(sys.argv) == 5:
			n = int(sys.argv[4])
		if n <= 0:
			raise Exception()
	except:
		print('Error: number of iterations must be a positive natural value.')
		exit(3)

	root = bissection(function, a, b, n)
	print('root (z):', root, '\tAbsError:', abs(ro.r(regex.sub(r'\bx\b', '('+str(root)+')', function))[0])) 
	
