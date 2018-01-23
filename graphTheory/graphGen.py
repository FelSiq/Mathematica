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

	for nodeID in range(nodeNum):
		nodeX = round(random.random() * maxDist, 2)
		nodeY = round(random.random() * maxDist, 2)
		graph[nodeID] = (nodeX, nodeY)
		print('node', nodeID, nodeX, nodeY)

	for _ in range(edgeNum):

		it = 0
		while True and it < 500:
			it += 1
			A = random.randint(0, nodeNum-1)
			B = A
			while B == A:
				B = random.randint(0, nodeNum-1)

			if dist(graph[A], graph[B]) <= maxRadius and not (A,B) in edges and not (B,A) in edges:
				edges.append((A,B))
				break
			if it == 500:
				maxRadius *= 1.1

		print('edge', A, B)