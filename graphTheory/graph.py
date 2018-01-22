import math

class graph:
	def __init__(self, filepath=None, separator=' ', directed=True, dimension=0, label=None):
		self.edgeList = {}
		self.dimension = max(0, dimension)
		self.label = label

		with open(filepath, mode = 'r') as file:
			for line in file:
				nodeA, nodeB, weight = line.split(separator)
				weight = float(weight)

				if not nodeA in self.edgeList:
					self.edgeList[nodeA] = {
						'coords': [0.0] * self.dimension, 
						'adj': {}}
				if not nodeB in self.edgeList:
					self.edgeList[nodeB] = {
						'coords': [0.0] * self.dimension,
						'adj': {}}

				self.edgeList[nodeA]['adj'][nodeB] = weight
				if not directed:
					self.edgeList[nodeB]['adj'][nodeA] = weight

	def print(self, coordinates=True, decorate=True):
		nodeNames = self.edgeList.keys()
		if self.dimension and coordinates:
			if decorate:
				for n in nodeNames:
					print('\u256DNODE: \'', n, 
						'\':\n\u251C\u2500\u2500\u2500\u2500\u257C Coordinates: ', self.edgeList[n]['coords'], 
						'\n\u2514\u2500\u2500\u2500\u2500\u257C Adjacents:   ',self.edgeList[n]['adj'], sep ='')
			else:
				for n in nodeNames:
					print(n, ':', self.edgeList[n])
		else:
			for n in nodeNames:
				print(n, ':', self.edgeList[n]['adj'])

	def setCoord(self, coordinates):
		for c in coordinates.keys():
			coordDim = len(coordinates[c])
			if coordDim < self.dimension:
				coordinates[c] = coordinates[c] + [0.0] * (self.dimension - coordDim)
			if coordDim > self.dimension:
				print('warning:', c, 'coordinates reside in a larger dimension than graph space.',
					'Projection on dim(' + str(self.dimension)+ ') will be used.')
			self.edgeList[c]['coords'] = coordinates[c][:self.dimension]


	def _genOutput(self):
		return {'distance': 0.0, 'path': [], 
			'iterations' : 0, 'totalChildrens' : 0}

	def _buildPathAndDistance(self, predVec, start, end):
		totalDistance = 0.0
		curNode = end
		path = []

		while curNode != start and curNode != '':
			path.insert(0, curNode)
			prevNode = predVec[curNode]
			totalDistance += self.edgeList[prevNode]['adj'][curNode]
			curNode = prevNode

		if curNode == start:
			path.insert(0, start)
			return path, totalDistance
		return [], math.inf

	def _updateOutput(self, output, statisticOutput):
		if not statisticOutput:
			output.pop('iterations')
			output.pop('totalChildrens')
		return output

	def _blindSearch(self, start, end, depthFirst, prune=True, statisticOutput=False):
		output = self._genOutput()

		predVec = {key : '' for key in self.edgeList}
		visitedVec = {key : False for key in self.edgeList}

		deque = [start]
		curNode = ''
		while len(deque):
			output['iterations'] += statisticOutput
			curNode = deque.pop()
			visitedVec[curNode] = True

			if curNode != end:
				for a in self.edgeList[curNode]['adj']:
					if not visitedVec[a]:
						output['totalChildrens'] += statisticOutput
						predVec[a] = curNode
						if depthFirst:
							deque.append(a)
						else:
							deque.insert(0, a) 
			else:
				deque.clear()

		output['path'], output['distance'] = self._buildPathAndDistance(predVec, start, end)

		return self._updateOutput(output, statisticOutput)

	def dfs(self, start, end, prune=True, statisticOutput=False):
		return self._blindSearch(start, end, True, prune, statisticOutput)

	def bfs(self, start, end, prune=True, statisticOutput=False):
		return self._blindSearch(start, end, False, prune, statisticOutput)

	def hillClimbing(self, start, end):
		if not self.dimension:
			print('error: hillClimbing works only with',
				'coordinate systems (graph is in space with null dimension).')
			return None
		return None
	def beamSearch(self, start, end):
		if not self.dimension:
			print('error: beamSearch works only with',
				'coordinate systems (graph is in space with null dimension).')
			return None
		return None
	def branchAndBound(self, start, end, prune=False, admissibleHeuristic=False):
		if admissibleHeuristic and not self.dimension:
			print('error: a admissible heuristic works only with',
				'coordinate systems (graph is in space with null dimension).')
			return None
		return None

	def Astar(self, start, end):
		return self.branchAndBound(start, end, prune=True, admissibleHeuristic=True)

if __name__ == '__main__':
	G = graph('1.in', dimension = 3)
	G.setCoord({'B': [1,2,3], 'C':[5,5,5], 'D':[10,15,10]})
	G.print(decorate=True)
	print(G.dfs('A', 'D', statisticOutput=True))
	print(G.bfs('A', 'D', statisticOutput=True))