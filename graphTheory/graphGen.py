import random
import sys

def dist(a, b):
	return sum([(a[i] - b[i])**2.0 for i in range(min(len(a), len(b)))])**0.5

if __name__ == '__main__':
	if len(sys.argv) <= 4:
		print('usage: <# nodes> <# edges> <max dist> <max radius>')
		exit(1)

	graph = {}
	edges = []

	nodeNum = int(sys.argv[1])
	edgeNum = int(sys.argv[2])
	maxDist = float(sys.argv[3])
	maxRadius = float(sys.argv[4])
	proximityDelta = maxDist * 0.05

	for nodeID in range(nodeNum):
		nodeX = round(random.random() * maxDist, 2)
		nodeY = round(random.random() * maxDist, 2)
		graph[nodeID] = (nodeX, nodeY)
		print('node', nodeID, nodeX, nodeY)

	for _ in range(edgeNum):

		it = 0
		flag = True
		while flag and it < 500:
			it += 1
			A = random.randint(0, nodeNum-1)
			B = A
			while B == A:
				B = random.randint(0, nodeNum-1)

			if dist(graph[A], graph[B]) <= maxRadius and not (A,B) in edges and not (B,A) in edges:
				
				slope = (graph[A][1] - graph[B][1])/(graph[A][0] - graph[B][1])
				intersect = graph[A][1] - slope * graph[A][0]
				
				noIntersection = 2
				for n in graph:
					if n != A and n != B:
						curNode = graph[n]
						noIntersection += (slope * curNode[0] + intersect - curNode[1] > proximityDelta)

				if noIntersection == len(graph):
					edges.append((A,B))
					flag = False
			if it == 500:
				maxRadius *= 1.1

		print('edge', A, B)