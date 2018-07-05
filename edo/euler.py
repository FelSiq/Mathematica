from rpy2 import robjects as ro
import sys
import re

"""
	This program solves a IVP EDO using Euler's method.
"""

def euler_method(f, a, b, y_0, h, y_val):
	it_num=1+int(y_val/h)
	y_i=y_0

	re_sub_x=re.compile(r"\bx\b")
	re_sub_y=re.compile(r"\by\b")

	print("y_0\t:", y_0)
	for i in range(it_num):
		x_val=a+h*i
		f_aux=re_sub_x.sub(str(x_val), f)
		f_aux=re_sub_y.sub(str(y_i), f_aux)

		y_i=y_i + h * ro.r(f_aux)[0]
		print("y_"+str(i+1)+"\t:", y_i)

	return y_i

if __name__ == "__main__":
	if len(sys.argv) < 7:
		print("usage:", sys.argv[0], "<function f> <x_a> <x_b> <y_0> <h> <Eval y at>")
		exit(1)

	try:
		a=float(sys.argv[2])
		b=float(sys.argv[3])
		y_0=float(sys.argv[4])
		h=float(sys.argv[5])
		y_val=float(sys.argv[6])

	except:
		print("All parameters should be real numbers.")
		exit(2)

	
	f=re.sub(r"\b([xy])\b", r"(\1)", sys.argv[1])
	print("Given function: y' = ", f)

	if not a <= y_val <= b:
		print("Wanted y value not in [x_a, x_b] interval.")
		exit(4)

	if (y_val/h) % 1 > 1.0e-6:
		print("Wanted y not a h's multiple.")
		exit(5)

	res=euler_method(f, a, b, y_0, h, y_val)

	print("Result\t:", res)
