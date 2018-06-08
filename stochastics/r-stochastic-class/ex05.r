# 5th and last exercise
im.resol=800

# Load data files
x1 <- read.csv('coords/base01coord_x.csv')$x
y1 <- read.csv('coords/base01coord_y.csv')$x
x2 <- read.csv('coords/base02coord_x.csv')$x
y2 <- read.csv('coords/base02coord_y.csv')$x

# Trying a simple 2D visualization
jpeg(
	'images/5q-2d-points.jpeg', 
	width=2*im.resol, 
	height=im.resol)
par(mfrow=c(1,2))
plot(x1, y1, main='Plot for Coords Set #1')
plot(x2, y2, main='Plot for Coords Set #2')
dev.off()
# Not a good vizualization, ai first

# Trying a more sofisticated approach
library(plot3D) # Self-explanatory
library(arules) # To discretize data
intervalNum=10

jpeg(
	'images/5q-3d-density.jpeg', 
	width=2*im.resol, 
	height=im.resol)
par(mfrow=c(1,2))

x1d<-discretize(
	x1, 
	method='interval', 
	breaks=intervalNum, 
	labels=F)
y1d<-discretize(
	y1, 
	method='interval', 
	breaks=intervalNum,
	labels=F)

hist3D(
	z=table(x1d, y1d), 
	border='black',
	main='Density for Coords Set #1',
	mai = c(0.1, 0.0, 0., 0.0)
)

x2d<-discretize(
	x2,
	method='interval', 
	breaks=intervalNum, 
	labels=F)
y2d<-discretize(
	y2,
	method='interval', 
	breaks=intervalNum,
	labels=F)

hist3D(
	z=table(x2d, y2d), 
	border='black',
	main='Density for Coords Set #2',
	mai = c(1, 0.1, 0.1, 0.1)
)

dev.off()
