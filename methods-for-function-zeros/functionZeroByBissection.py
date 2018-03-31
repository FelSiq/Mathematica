from rpy2 import robjects as ro
import regex
import sys

def bissection(function, a, b, n=100, retList=False, var='x'):
	start = a
	end = b
	xMeans = list()
	subVarRe = regex.compile(r'\b' + var + r'\b')
	for i in range(n):
		curXMean = (start + end) * 0.5		
		xMeans.append(curXMean)
		
		val1 = ro.r(subVarRe.sub('('+str(curXMean)+')', function))[0] 
		val2 = ro.r(subVarRe.sub('('+str(start)+')', function))[0] 

		if val1 * val2 > 0:
			# The sign does not change between f(curXMean) and f(start) i.e. the 
			# exact solution z is on the 'other' interval.
			# Then z e [curXMean, end]
			start = curXMean
		else:
			# The sign between f(curXMean) and f(start) changes, so the
			# exact solution z is on 'this' interval.
			# Then z e [start, curXMean]
			end = curXMean
	if retList:
		return xMeans
	return xMeans.pop()

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
	
