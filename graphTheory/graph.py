import matplotlib.pyplot as plt
import operator
import math
import sys

class graph:
	def _addNode(self, label, coord=None, adj=None, hcost={'enabled': False, 'cost': 0.0}):
		if not label in self.edgeList:
			self.edgeList[label] = {
				'coord': coord if coord else [0.0] * self.dimension,
				'adj': adj if adj else {},
				'hcost': hcost}

	def __init__(self, filepath=None, separator=' ', directed=True, 
		geometrical=True, dimension=0, label=None):
		self.edgeList = {}
		self.label = label
		self.geometrical = geometrical
		self.plotPathList = []

		if dimension < 0:
			print('w: graph space dimension set to 0 (can\'t use negative values).')
		self.dimension = max(0, dimension)

		with open(filepath, mode = 'r') as file:
			for line in file:
				tokens = line.strip().split(separator)
				if tokens[0] == 'edge':
					nodeA = tokens[1]
					nodeB = tokens[2]
					self._addNode(nodeA)
					self._addNode(nodeB)

					if geometrical and len(tokens) < 4:
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

					self._addNode(nodeLabel)
					self.edgeList[nodeLabel]['coord'] = nodeCoord

				elif tokens[0] == 'hcost':
					curNode = tokens[1]
					self._addNode(curNode)
					self.edgeList[curNode]['hcost'] = {'enabled': True, 'cost': float(tokens[2])}

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

	def _searchPlotSetup(self, title, start, end):
		xSize, ySize = self.plot(show=False, returnDimensions = True)
		circleRadius = min(xSize, ySize) * 0.1
		plt.gcf().gca().add_artist(plt.Circle(tuple(self.edgeList[end]['coord']), circleRadius, color='blue', fill=False))
		plt.gcf().gca().add_artist(plt.Circle(tuple(self.edgeList[start]['coord']), circleRadius, color='green', fill=False))
		plt.title(title)
		plt.ion()


	def _genOutput(self):
		return {'distance': 0.0, 'path': [], 'iterations' : 0, 
		'totalChildrens' : 0, 'visitOrder' : []}

	def _calculateTotalDist(self, curPath='', goal=''):
		totalDistance = 0.0
		if curPath != '':
			prevNode = curPath[0]
			for p in curPath[1:]:
				totalDistance += self.edgeList[prevNode]['adj'][p]
				prevNode = p

		if prevNode == goal:
			return curPath, totalDistance
		return '', math.inf

	def _updateOutput(self, output, statisticOutput):
		if not statisticOutput:
			output.pop('iterations')
			output.pop('totalChildrens')
		return output

	def _calcDist(self, a, b):
		aCoord = self.edgeList[a]['coord']
		bCoord = self.edgeList[b]['coord']
		return sum([(aCoord[i] - bCoord[i])**2.0 for i in range(min(len(aCoord), len(bCoord)))])**0.5

	def dfs(self, start, end, prune=True, lexicographical=True, 
		statisticOutput=False, plot=False, plotSpeed=0.2):
		return self._blindSearch(start, end, True, prune, 
			statisticOutput, lexicographical, plot, plotSpeed, 'Depth First Search')

	def bfs(self, start, end, prune=True, lexicographical=True, 
		statisticOutput=False, plot=False, plotSpeed=0.2):
		return self._blindSearch(start, end, False, prune, 
			statisticOutput, lexicographical, plot, plotSpeed, 'Breadth First Search')

	"""
	Hill climbing is a DFS that proceed only into children node that is closest to the goal.
	It's a Beam Search with keptChildren = 1.
	"""
	def hillClimbing(self, start, end, prune=True, 
		statisticOutput=False, plot=False, plotSpeed=0.2, lexicographical=True):
		if not self.geometrical:
			print('error: hillClimbing works only with',
				'coordinate systems (graph is in space with null dimension).')
			return None
		return self._informedSearch(start, end, prune, 
			statisticOutput, 1, plot, plotSpeed, 'Hill Climbing algorithm',lexicographical)

	"""
	Beam Search is like BFS, but it does only keep w childrens that is 
	closer to the goal per search tree level.
	"""
	def beamSearch(self, start, end, keptChildren=3, prune=True, 
		statisticOutput=False, plot=False, plotSpeed=0.2, lexicographical=True):
		if not self.geometrical:
			print('error: beamSearch works only with',
				'coordinate systems (graph is in space with null dimension).')
			return None
		return self._informedSearch(start, end, prune, 
			statisticOutput, keptChildren, plot, plotSpeed, 'Beam Search algorithm',lexicographical)

	def _plotColorLocation(self, curNode, color='Blue', time=0.0):
		curPosition = self.edgeList[curNode]['coord']
		plt.scatter(curPosition[0], curPosition[1], color=color)
		if time > 0.0:
			plt.pause(time)

	def _plotCurrentPath(self, start='', predVec=None, curPath=''):
		while len(self.plotPathList):
			path = self.plotPathList.pop()
			plt.plot(path[0], path[1], color='black')

		if predVec:
			curNode = start
			curPosition = self.edgeList[curNode]['coord']
		else:
			curNode = curPath[0]
			curPosition = self.edgeList[curPath[0]]['coord']

		i = 1
		while (predVec and curNode != '') or (i < len(curPath)):
			prevPosition = curPosition
			curNode = predVec[curNode] if predVec else curPath[i]
			i += 1
			if curNode != '':
				curPosition = self.edgeList[curNode]['coord']
				plotXs = [prevPosition[0], curPosition[0]]
				plotYs = [prevPosition[1], curPosition[1]]
				self.plotPathList.append((plotXs, plotYs))
				plt.plot(plotXs, plotYs, color='red')

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

	def plot(self, show=True, time=2.0, returnDimensions=False):
		if not self.geometrical:
			print('E: can\'t plot a non-geometrical graph.')
			return None

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
		plt.scatter(x=x, y=y, color='blue')

		if show:
			plt.pause(min(0.0, time))

		if returnDimensions:
			xLims = plt.gca().get_xlim()
			yLims = plt.gca().get_ylim() 
			return xLims[1] - xLims[0], yLims[1] - yLims[0]

	def _blindSearch(self, start, end, depthFirst, prune=True, 
		statisticOutput=False, lexicographical=True, plot=False, plotSpeed=0.2, 
		title = 'Generic Blind Search Algorithm'):
		output = self._genOutput()

		if plot:
			self._searchPlotSetup(title + (' (prunned)' if prune else ''), start, end)

		predVec = {key : '' for key in self.edgeList}
		visitedVec = {key : False for key in self.edgeList}

		deque = [(start, start)]
		curNode = start
		prevNode = start
		while len(deque):
			curItem = deque.pop()
			aux = curItem[0]

			if not prune or not visitedVec[aux]:
				if plot:
					self._plotColorLocation(prevNode, color='black')
				prevNode = curNode
				curNode = curItem[0]
				curPath = curItem[1]

				visitedVec[curNode] = True
				output['iterations'] += statisticOutput
				output['visitOrder'].append(curNode)

				if plot:
					self._plotCurrentPath(curPath=curPath)
					self._plotColorLocation(curNode, color='Red', time=plotSpeed)

				done = False
				if curNode != end and not done:
					adjList = list(self.edgeList[curNode]['adj'].keys())
					if lexicographical:
						adjList.sort(reverse=depthFirst)

					for a in adjList:
						if curPath.find(a) == -1:
							output['totalChildrens'] += statisticOutput
							
							newItem = (a, curPath + a)
							predVec[a] = curNode

							if depthFirst:
								deque.append(newItem)
							else:
								deque.insert(0, newItem)
							if a == end:
								done = True
				else:
					deque.clear()

		output['path'], output['distance'] = self._calculateTotalDist(curPath=curPath, goal=end)

		if plot:
			self._plotCurrentPath(curPath=curPath)
			plt.pause(3.0)

		return self._updateOutput(output, statisticOutput)

	def _informedSearch(self, start, end, prune=True, 
		statisticOutput=False, keptChildrens=1, plot=False, plotSpeed=0.2, 
		title='Generic Informed Search algorithm', lexicographical=True):
		output = self._genOutput()

		if plot:
			self._searchPlotSetup(title + (' (prunned)' if prune else ''), start, end)

		visitedVec = {key : False for key in self.edgeList}

		curLevelNodes = [(start, 0.0, 0.0, start)]
		childrenNodes = []
		curNode = start
		curPath = start
		prevNode = start

		while len(curLevelNodes) or len(childrenNodes):
			if not len(curLevelNodes):
				# Sort based on which children node is heuristically closer to the objective,
				# then use lexicographical order in case of draw.
				childrenNodes.sort(key=lambda item : item[1] + item[2], reverse=True)

				size = len(childrenNodes)
				if size > keptChildrens:
					childrenNodes = childrenNodes[size - keptChildrens:size]
				if lexicographical:
					curLevelNodes = sorted(childrenNodes, key=operator.itemgetter(0), reverse=True)
				else:
					curLevelNodes = [i for i in a]
				childrenNodes.clear()

			curItem = curLevelNodes.pop()
			aux = curItem[0]
			if not prune or not visitedVec[aux]:
				if plot:
					self._plotColorLocation(prevNode, color='black')

				prevNode = curNode
				curNode = curItem[0]
				curCost = curItem[1]
				curPath = curItem[3]
				visitedVec[curNode] = True
				output['iterations'] += statisticOutput
				output['visitOrder'].append(curNode)

				if plot:
					self._plotCurrentPath(curPath=curPath)
					self._plotColorLocation(curNode, color='Red', time=plotSpeed)

				if curNode != end:
					adjList = list(self.edgeList[curNode]['adj'].keys())

					for a in adjList:
						if curPath.find(a) == -1:
							output['totalChildrens'] += statisticOutput
							if self.edgeList[a]['hcost']['enabled']:
								heuristicCost = self.edgeList[a]['hcost']['cost']
							else:
								heuristicCost = self._calcDist(a, end)
							newChildren = (a, curCost + self.edgeList[curNode]['adj'][a], heuristicCost, curPath + a)
							childrenNodes.append(newChildren)

				else:
					curLevelNodes.clear()
					childrenNodes.clear()

		output['path'], output['distance'] = self._calculateTotalDist(curPath=curPath, goal=end)

		if plot:
			self._plotCurrentPath(curPath=curPath)
			plt.pause(3.0)

		return self._updateOutput(output, statisticOutput)

	"""
	This is Branch And Bound algorithm. To be fair, it's exactly the same as Dijkstra's Algorithm.
	"""
	def branchAndBound(self, start, end, prune=False, admissibleHeuristic=False, 
		statisticOutput=False, plot=False, plotSpeed=0.2, title='Branch and Bound algorithm', lexicographical=True):
		if admissibleHeuristic and not self.geometrical:
			print('error: a admissible heuristic works only with',
				'coordinate systems (graph is in space with null dimension).')
			return None
		
		if plot:
			self._searchPlotSetup(title + (' (prunned)' if prune else ''), start, end)

		output = self._genOutput()

		visitedVec = {key : False for key in self.edgeList}

		minHeap = [(start, 0.0, 0.0, start)]
		minPathLen = math.inf
		curNode = start
		prevNode = start
		curPath = start
		winnerPath = start

		while len(minHeap):
			curItem = minHeap.pop()
			aux = curItem[0]

			if not prune or not visitedVec[aux]:
				if plot:
					self._plotColorLocation(prevNode, color='black')

				prevNode = curNode
				curNode = curItem[0]
				curPath = curItem[3]
				output['iterations'] += statisticOutput
				output['visitOrder'].append(curNode)

				if plot:
					self._plotCurrentPath(curPath=curPath)
					self._plotColorLocation(curNode, color='Red', time=plotSpeed)

				accumulatedDist = curItem[1]
				visitedVec[curNode] = True

				if accumulatedDist < minPathLen:
					if curNode != end:
						adjList = list(self.edgeList[curNode]['adj'].keys())

						for a in adjList:
							if curPath.find(a) == -1:
								totalDist = accumulatedDist + self.edgeList[curNode]['adj'][a]
								if totalDist < minPathLen:
									output['totalChildrens'] += statisticOutput
									newChildren = (a, totalDist, self._calcDist(a, end), curPath + a)
									minHeap.append(newChildren)
									if a == end:
										winnerPath = curPath + a
										minPathLen = totalDist
							
						# Sort based on which children node is heuristically closer to the objective
						if admissibleHeuristic:
							minHeap.sort(key=lambda item : item[1] + item[2], reverse=True)
						else:
							minHeap.sort(key=lambda item : item[1], reverse=True)
					else:
						winnerPath = curPath
						minPathLen = accumulatedDist

		output['path'], output['distance'] = self._calculateTotalDist(curPath=winnerPath, goal=end)

		if plot:
			self._plotCurrentPath(curPath=winnerPath)
			plt.pause(3.0)

		return self._updateOutput(output, statisticOutput)


def printSearchOutput(algorithm, outputDic):
	keys = sorted(outputDic.keys())
	print('\u256D', algorithm)
	for k in range(len(keys)):
		curKey = keys[k]
		print(('\u2514' if k == len(keys) - 1 else '\u251C') 
			+ '\u2500\u2500\u2500\u2500\u257C', curKey, ':', outputDic[curKey])

if __name__ == '__main__':

	if len(sys.argv) <= 4:
		print('usage: <graph filepath> <start> <end> <plot? 0/1> <plotDelay (optional)>')
		exit(1)

	start = sys.argv[2]
	end = sys.argv[3]
	plot = bool(int(sys.argv[4]))

	plotDelay = 0.025
	if len(sys.argv) >= 6:
		plotDelay = float(sys.argv[5])

	G = graph(sys.argv[1], geometrical = True, directed=False)
	G.print()

	if not (start in G.edgeList and end in G.edgeList):
		print('E: invalid start/end indexes for search.') 
		exit(2)

	print('\n')

	searches = {'BFS':False, 'DFS':False, 'HC':False, 'BS':True, 'BB':False, 'AStar':False}

	if searches['BFS']:
		printSearchOutput('BFS', G.bfs(start, end, prune=False, lexicographical=True, 
			statisticOutput=True, plot=plot, plotSpeed=plotDelay))
	if searches['DFS']:
		printSearchOutput('DFS', G.dfs(start, end, prune=True, lexicographical=True, 
			statisticOutput=True, plot=plot, plotSpeed=plotDelay))
	if searches['HC']:
		printSearchOutput('HC', G.hillClimbing(start, end, 
			statisticOutput=True, plot=plot, plotSpeed=plotDelay))
	if searches['BS']:
		printSearchOutput('BS', G.beamSearch(start, end, 
			statisticOutput=True, keptChildren=2, plot=plot, plotSpeed=plotDelay))
	if searches['BB']:
		printSearchOutput('B&B', G.branchAndBound(start, end, 
			statisticOutput=True, plot=plot, plotSpeed=plotDelay, prune=True))
	if searches['AStar']:
		printSearchOutput('A*', G.Astar(start, end, statisticOutput=True, plot=plot))
