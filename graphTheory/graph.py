import matplotlib.pyplot as plt
import operator
import math
import sys

class graph:
	def __init__(self, filepath=None, separator=' ', directed=True, 
		geometrical=True, dimension=0, label=None):
		self.edgeList = {}
		self.label = label
		self.geometrical = geometrical

		if dimension < 0:
			print('w: graph space dimension set to 0 (can\'t use negative values).')
		self.dimension = max(0, dimension)

		with open(filepath, mode = 'r') as file:
			for line in file:
				tokens = line.strip().split(separator)
				if tokens[0] == 'edge':
					nodeA = tokens[1]
					nodeB = tokens[2]
					if not nodeA in self.edgeList:
						self.edgeList[nodeA] = {
							'coord': [0.0] * self.dimension, 
							'adj': {}}
					if not nodeB in self.edgeList:
						self.edgeList[nodeB] = {
							'coord': [0.0] * self.dimension,
							'adj': {}}

					if geometrical:
						weight = self._calcDist(nodeA, nodeB)
					else:
						weight = float(tokens[3])

					self.edgeList[nodeA]['adj'][nodeB] = weight
					if not directed:
						self.edgeList[nodeB]['adj'][nodeA] = weight

				elif tokens[0] == 'node':
					nodeLabel = tokens[1]
					if geometrical:
						nodeCoord = list(map(float, tokens[2:]))
						newCoordDim = len(nodeCoord)
						if self.dimension > 0:
							if newCoordDim > self.dimension:
								print('w: \'' + nodeLabel + 
									'\' node coordinates has more dimensions than graph. I\'ll use the projection.')
								nodeCoord = nodeCoord[:self.dimension]
							elif newCoordDim < self.dimension:
								print('w: \'' + nodeLabel + 
									'\' node coordinates has less dimensions than graph. Filling with \'0\'.')
								nodeCoord += [0.0] * (self.dimension - newCoordDim)
						else:
							self.dimension = len(nodeCoord)
					else:
						nodeCoord = []

					if not nodeLabel in self.edgeList:
						self.edgeList[nodeLabel] = {
							'coord': nodeCoord, 
							'adj': {}}
					else:
						self.edgeList[nodeLabel]['coord'] = nodeCoord
				else:
					print('e: unknown token \'' + tokens[0] + '\'. Ignoring it.')

	def print(self, coordinates=True, decorate=True):
		nodeNames = self.edgeList.keys()
		if self.geometrical and coordinates:
			if decorate:
				for n in nodeNames:
					print('\u256DNODE: \'', n, 
						'\':\n\u251C\u2500\u2500\u2500\u2500\u257C Coordinates: ', self.edgeList[n]['coord'], 
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
			self.edgeList[c]['coord'] = coordinates[c][:self.dimension]


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
			if prevNode != '':
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

	def _blindSearch(self, start, end, depthFirst, prune=True, 
		statisticOutput=False, lexicographical=True, plot=False, plotSpeed=0.2, 
		title = 'Generic Blind Search Algorithm'):
		output = self._genOutput()

		if plot:
			self.plot(show=False)
			plt.title(title)
			plt.ion()

		predVec = {key : '' for key in self.edgeList}
		visitedVec = {key : False for key in self.edgeList}

		deque = [start]
		curNode = start
		while len(deque):
			if plot:
				self._plotColorLocation(curNode, color='black')

			output['iterations'] += statisticOutput
			curNode = deque.pop()
			visitedVec[curNode] = True

			if plot:
				self._plotCurrentPath(predVec, curNode)
				self._plotColorLocation(curNode, color='Red', time=plotSpeed)

			if curNode != end:
				adjList = list(self.edgeList[curNode]['adj'].keys())
				if lexicographical:
					adjList.sort(reverse=depthFirst)

				for a in adjList:
					if (prune and not visitedVec[a]) or not (prune or self._checkPreds(predVec, curNode, a)):
						output['totalChildrens'] += statisticOutput
						predVec[a] = curNode
						if depthFirst:
							deque.append(a)
						else:
							deque.insert(0, a) 
			else:
				deque.clear()

		output['path'], output['distance'] = self._buildPathAndDistance(predVec, start, end)

		if plot:
			if visitedVec[end]:
				self._plotCurrentPath(predVec, end)
			plt.pause(3.0)

		return self._updateOutput(output, statisticOutput)

	def _calcDist(self, a, b):
		aCoord = self.edgeList[a]['coord']
		bCoord = self.edgeList[b]['coord']
		return sum([(aCoord[i] - bCoord[i])**2.0 for i in range(min(len(aCoord), len(bCoord)))])**0.5

	def _checkPreds(self, predVec, start, query):
		found = False
		curNode = start
		while not found and curNode != '':
			found = curNode == query
			curNode = predVec[curNode]
		return found

	def _informedSearch(self, start, end, depthFirst, prune=True, 
		statisticOutput=False, keptChildrens=0, plot=False, plotSpeed=0.2, 
		title='Generic Informed Search algorithm'):
		output = self._genOutput()

		if plot:
			self.plot(show=False)
			plt.title(title)
			plt.ion()

		predVec = {key : '' for key in self.edgeList}
		visitedVec = {key : False for key in self.edgeList}

		deque = [(start, 0.0)]
		curNode = start
		while len(deque):

			if plot:
				self._plotColorLocation(curNode, color='black')

			output['iterations'] += statisticOutput
			curNode = deque.pop()[0]
			visitedVec[curNode] = True

			if plot:
				self._plotCurrentPath(predVec, curNode)
				self._plotColorLocation(curNode, color='Red', time=plotSpeed)

			if curNode != end:
				auxList = []
				adjList = list(self.edgeList[curNode]['adj'].keys())

				for a in adjList:
					if (prune and not visitedVec[a]) or not (prune or self._checkPreds(predVec, curNode, a)):
						output['totalChildrens'] += statisticOutput
						predVec[a] = curNode
						newChildren = (a, self._calcDist(a, end))
						auxList.append(newChildren)
				
				# Sort based on which children node is heuristically closer to the objective,
				# then use lexicographical order in case of draw.

				if depthFirst:
					auxList.sort(key=operator.itemgetter(1, 0), reverse=True)
					deque = deque + auxList
				else:
					deque = auxList + deque
					deque.sort(key=operator.itemgetter(1))
					deque[:keptChildrens]

			else:
				deque.clear()

		output['path'], output['distance'] = self._buildPathAndDistance(predVec, start, end)

		if plot:
			if visitedVec[end]:
				self._plotCurrentPath(predVec, end)
			plt.pause(3.0)

		return self._updateOutput(output, statisticOutput)

	def dfs(self, start, end, prune=True, lexicographical=True, 
		statisticOutput=False, plot=False, plotSpeed=0.2):
		return self._blindSearch(start, end, True, prune, 
			statisticOutput, lexicographical, plot, plotSpeed, 'Depth First Search')

	def bfs(self, start, end, prune=True, lexicographical=True, 
		statisticOutput=False, plot=False, plotSpeed=0.2):
		return self._blindSearch(start, end, False, prune, 
			statisticOutput, lexicographical, plot, plotSpeed, 'Breadth First Search')

	"""
	Hill climbing is a DFS that proceed first into children node that is closest to the goal.
	"""
	def hillClimbing(self, start, end, prune=True, 
		statisticOutput=False, plot=False, plotSpeed=0.2):
		if not self.geometrical:
			print('error: hillClimbing works only with',
				'coordinate systems (graph is in space with null dimension).')
			return None
		return self._informedSearch(start, end, True, prune, 
			statisticOutput, 0, plot, plotSpeed, 'Hill Climbing algorithm')

	"""
	Beam Search is like BFS, but it does only keep w childrens that is closer to the goal per search tree level.
	"""
	def beamSearch(self, start, end, keptChildren=2, prune=True, 
		statisticOutput=False, plot=False, plotSpeed=0.2):
		if not self.geometrical:
			print('error: beamSearch works only with',
				'coordinate systems (graph is in space with null dimension).')
			return None
		return self._informedSearch(start, end, False, prune, 
			statisticOutput, keptChildren, plot, plotSpeed, 'Beam Search algorithm')

	def _plotColorLocation(self, curNode, color='Blue', time=0.0):
		curPosition = self.edgeList[curNode]['coord']
		plt.scatter(curPosition[0], curPosition[1], color=color)
		if time > 0.0:
			plt.pause(time)

	def _plotCurrentPath(self, predVec, start):
		curNode = start
		curPosition = self.edgeList[curNode]['coord']

		for k in self.edgeList.keys():
			for a in self.edgeList[k]['adj']:
				kCoord = self.edgeList[k]['coord']
				aCoord = self.edgeList[a]['coord']
				plt.plot([kCoord[0], aCoord[0]], [kCoord[1], aCoord[1]], color='black')

		while curNode != '':
			prevPosition = curPosition
			curNode = predVec[curNode]
			if curNode != '':
				curPosition = self.edgeList[curNode]['coord']
				plt.plot([prevPosition[0], curPosition[0]], [prevPosition[1], curPosition[1]], color='Red')


	"""
	This is Branch And Bound algorithm. To be fair, it's exactly the same as Dijkstra's Algorithm.
	"""
	def branchAndBound(self, start, end, prune=False, admissibleHeuristic=False, 
		statisticOutput=False, plot=False, plotSpeed=0.2, title='\'Branch and Bound\' algorithm'):
		if admissibleHeuristic and not self.geometrical:
			print('error: a admissible heuristic works only with',
				'coordinate systems (graph is in space with null dimension).')
			return None
		
		if plot:
			self.plot(show=False)
			plt.title(title)
			plt.ion()

		output = self._genOutput()

		predVec = {key : '' for key in self.edgeList}
		visitedVec = {key : False for key in self.edgeList}

		minHeap = [(start, 0.0, 0.0)]
		minPathLen = math.inf
		curNode = start

		while len(minHeap):
			output['iterations'] += statisticOutput

			if plot:
				self._plotColorLocation(curNode, color='black')

			curItem = minHeap.pop()
			curNode = curItem[0]

			if plot:
				self._plotCurrentPath(predVec, curNode)
				self._plotColorLocation(curNode, color='Red', time=plotSpeed)

			accumulatedDist = curItem[1]
			visitedVec[curNode] = True

			if curNode != end:
				adjList = list(self.edgeList[curNode]['adj'].keys())

				for a in adjList:
					if (prune and not visitedVec[a]) or not (prune or self._checkPreds(predVec, curNode, a)):
						totalDist = accumulatedDist + self.edgeList[curNode]['adj'][a]
						if totalDist < minPathLen:
							output['totalChildrens'] += statisticOutput
							predVec[a] = curNode
							newChildren = (a, totalDist, self._calcDist(a, end))
							minHeap.append(newChildren)
					
				# Sort based on which children node is heuristically closer to the objective
				if admissibleHeuristic:
					minHeap.sort(key=lambda item : item[1] + item[2], reverse=True)
				else:
					minHeap.sort(key=lambda item : item[1], reverse=True)

			else:
				minPathLen = min(minPathLen, accumulatedDist)

		output['path'], output['distance'] = self._buildPathAndDistance(predVec, start, end)

		if plot:
			if visitedVec[end]:
				self._plotCurrentPath(predVec, end)
			plt.pause(3.0)

		return self._updateOutput(output, statisticOutput)

	"""
	A* is just Branch and Bound (or Dijkstra's algorhtm) without revisiting nodes and with a admissible heuristic.

	A admissible heuristic is a value that is always less or equal than the true value, i.e, it never overestimate
	the total value/cost/distance/etc to the goal. For example, if you take a straight line between two points, 
	desconsidering if there is any obstacle or other constraint between then, it is a admissible heuristic, because 
	the straight line is the smallest distance between two points, no matter if this straight path is or isn't possible
	to be done.
	"""
	def Astar(self, start, end, statisticOutput=False, plot=False, plotSpeed=0.2):
		return self.branchAndBound(start, end, True, True, 
			statisticOutput, plot, plotSpeed, 'A* algorithm')

	def plot(self, show=True, time=2.0):
		x =[]
		y =[]

		if show:
			plt.ion()

		nodes = self.edgeList.keys()
		for k in nodes:
			node = self.edgeList[k]['coord']
			x.append(node[0])
			y.append(node[1])
		for k in nodes:
			for a in self.edgeList[k]['adj']:
				kCoord = self.edgeList[k]['coord']
				aCoord = self.edgeList[a]['coord']
				plt.plot([kCoord[0], aCoord[0]], [kCoord[1], aCoord[1]], color='black')
		plt.scatter(x=x, y=y)

		if show:
			plt.pause(min(0.0, time))

if __name__ == '__main__':

	if len(sys.argv) <= 4:
		print('usage: <graph filepath> <start> <end> <plot? 0/1>')
		exit(1)

	start = sys.argv[2]
	end = sys.argv[3]
	plot = bool(int(sys.argv[4]))

	G = graph(sys.argv[1], geometrical = True, directed=False)
	G.print()

	print('\n')
	print('DFS:', G.dfs(start, end, prune=False, lexicographical=True, statisticOutput=True, plot=plot))
	print('BFS:', G.bfs(start, end, prune=False, lexicographical=True, statisticOutput=True, plot=plot))
	print('HC:', G.hillClimbing(start, end, statisticOutput=True, plot=plot))
	print('BS:', G.beamSearch(start, end, statisticOutput=True, keptChildren=2, plot=plot))
	print('B&B:', G.branchAndBound(start, end, statisticOutput=True, plot=plot))
	print('A*:', G.Astar(start, end, statisticOutput=True, plot=plot))
