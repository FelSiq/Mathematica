import math
import regex

class expressionSolver:
	def __init__(self):
		self.reg_Sin = regex.compile(r'sin\{([^}]+)\}')
		self.reg_Cos = regex.compile(r'cos\{([^}]+)\}')
		self.reg_Log = regex.compile(r'log\{([^,]+),\s*([^}]+)\}')
		self.reg_Mon = regex.compile(r'([^(\s]+)\s*^\s*([^)]+)')
		self.reg_Fac = regex.compile(r'([^(\s]+)\s*!')
		self.reg_Low = regex.compile(r'([^(\s]+)\s*+\s*([^\s)]+)')
		self.reg_High = regex.compile(r'([^(\s]+)\s*([*\\])\s*([^\s)]+)')
		self.reg_X = regex.compile('x')
		self.reg_Pi = regex.compile('pi')
		self.reg_Euler = regex.compile('e')
		self.reg_Par = regex.compile(r'\(([^)]+)\)')
		self.reg_Minus = regex.compile(r'([^\s\(\+]+)\s*\-\s*([^\s)]+)')
		self.reg_Pre = regex.compile(r'(cos|sin|log)\((.+?)\)')
	
	def __solveExpr__(self, function):
		try:
			return float(regex.sub('\s', '', function))
		except:
			function = self.reg_Minus.sub(r'\1+-\2', function)
			# Parenthesis solving
			compFunc = self.reg_Par.findall(function)
			for p in compFunc:
				function = self.reg_Par.sub(str(self.__solveExpr__(p)), function, count=1)

			# Post-parenthesis stuff
			compFunc = self.reg_Log.findall(function)
			for e in compFunc:
				argument = self.__solveExpr__(e[0])
				base = self.__solveExpr__(e[1])
				function = self.reg_Log.sub(str(math.log(argument, base)), function, count=1) 

			compFunc = self.reg_Sin.findall(function)
			for e in compFunc:
				function = self.reg_Sin.sub(str(math.sin(self.__solveExpr__(e))), function, count=1) 

			compFunc = self.reg_Cos.findall(function)
			for e in compFunc:
				function = self.reg_Cos.sub(str(math.cos(self.__solveExpr__(e))), function, count=1) 

			compFunc = self.reg_Fac.findall(function)
			for e in compFunc:
				function = self.reg_Fac.sub(str(math.factorial(self.__solveExpr__(e))), function, count=1) 

			compFunc = self.reg_Mon.search(function)
			while compFunc:
				curVal = self.__solveExpr__(compFunc.group(1)) ** self.__solveExpr__(compFunc.group(2))
				function = self.reg_Mon.sub(str(curVal), function, count=1)
				compFunc = self.reg_Mon.search(function)
			
			compFunc = self.reg_High.search(function)
			while compFunc:
				operator = compFunc.group(2)
				if operator == '*':
					curVal = self.__solveExpr__(compFunc.group(1)) * self.__solveExpr__(compFunc.group(3))
				else:
					curVal = self.__solveExpr__(compFunc.group(1)) / self.__solveExpr__(compFunc.group(3))
				function = self.reg_High.sub(str(curVal), function, count=1)
				compFunc = self.reg_High.search(function)

			compFunc = self.reg_Low.search(function)
			while compFunc:
				operandA = compFunc.group(1)
				operandB = compFunc.group(2)
				print(operandA, '\t|\t', operandB, '\t|\t', function)
				curVal = self.__solveExpr__(operandA) + self.__solveExpr__(operandB)
				function = self.reg_Low.sub(str(curVal), function, count=1)
				compFunc = self.reg_Low.search(function)
			return self.__solveExpr__(function)

	def __subVariable__(self, function, x):
		number = self.reg_X.sub(str(x), function)
		number = self.reg_Pi.sub(str(math.pi), number)
		number = self.reg_Euler.sub(str(math.e), number)
		return number

	def __preprocessing__(self, function):
		preprocessedFunc = function
		while self.reg_Pre.search(preprocessedFunc):
			preprocessedFunc = self.reg_Pre.sub(r'\1{\2}', preprocessedFunc)
		return preprocessedFunc

	def solve(self, function, varValue=0.0):
		function = self.__subVariable__(function, varValue)
		function = self.__preprocessing__(function)
		return self.__solveExpr__(function)

if __name__ == '__main__':
	function = 'log(e, e)'
	arithSolver = expressionSolver()
	print(arithSolver.solve(function, 2.0))
