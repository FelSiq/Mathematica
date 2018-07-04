data=read.table("results.dat", sep="&")
jpeg("plot.jpeg", width=640, height=640)
plot(cbind(seq(1, 10), data$V2), 
	col='blue', 
	type='l', 
	main='Gráfico de x_{n+1} e Erro E_{n+1} por Iteração n',
	ylab='x_{n+1}', 
	ylim=c(min(data$V1), max(data$V2)),
	lwd=4,
	xlab='Iteração') 
lines(cbind(seq(1,10), data$V1), lwd=4, lty=2, col='red')

# Subtitles
middle=sum(c(min(data$V1), max(data$V2)))/2
text(rep(9, 2), 
	c(middle, middle+0.10), 
	labels=c('x_{n+1}', 'E_{n+1}'))
lines(x=c(7.2, 8.2), y=c(middle, middle), lwd=4, col='blue')
lines(x=c(7.2, 8.2), y=c(middle+0.10, middle+0.10), lwd=4, lty=2, col='red')
dev.off()	
