from rpy2 import robjects as ro
import regex

def bissection(function, a, b, n=100, retList=False, var='x'):
	start = a
	end = b
	xMeans = list()
	for i in range(n):
		curXMean = (start + end) * 0.5		
		xMeans.append(curXMean)
		
		val1 = ro.r(regex.sub(var, str(curXMean), function))[0] 
		val2 = ro.r(regex.sub(var, str(start), function))[0] 

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
	function = '3*x^5 + 2*x^3 - 2*x^2 + 7*x + 4.5'	
	root = bissection(function, -600, 600)
	print('root (z):', root, '\tAbsError:', abs(ro.r(regex.sub('x', str(root), function))[0])) 
	
