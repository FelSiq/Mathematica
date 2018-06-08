# 3rd exercise
im.resol<-800

circleRandomWalk <- function(n=100, p=0.5, s=6, s0=0) {
	ret<-s0
	
	curState=s0
	for (i in 1:n) {
		curState <- curState + (if (runif(1) > p) 1 else -1)
		curState <- curState %% s
		ret<-c(ret, curState)
	}

	return (ret)
}

cat('Generating line plots...\n')
jpeg(
	filename='images/line.jpeg', 
	width=im.resol,
	height=im.resol)
par(mfrow=c(2,2))
for (i in c(10, 25, 50, 100)) {
	r<-circleRandomWalk(n=i)
	plot(
		r, 
		type='l', 
		main=paste(i, '-Steps', sep=''), 
		xlab='Steps', 
		ylab='State')
}
dev.off()

cat('Generating histograms...\n')
jpeg(
	filename='images/hist.jpeg',
	width=im.resol,
	height=im.resol)
par(mfrow=c(2,2))
for (i in c(10, 25, 50, 100)) {
	r<-circleRandomWalk(n=i)
	hist(
		r, 
		main=paste(i, '-Steps', sep=''), 
		breaks=24,
		xlab='State')
}
dev.off()


# A more extreme case
cat('Running a bigger case, this may take a while...\n')
jpeg(
	filename='images/inf-hist.jpeg',
	width=im.resol,
	height=im.resol)
par(mfrow=c(1,1))
hist(
	circleRandomWalk(1e+5), 
	breaks=24,
	main='1e+5-Steps',
	xlab='State')
dev.off()
