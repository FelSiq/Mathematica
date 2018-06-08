# 4th exercise

# Generate the time of each bus arrival
genBusTimes <- function(start, n=1000, lambda=5.5) {
	t<-start
	times<-NULL
	for (i in 1:n) {
		t<-t+rexp(1, lambda)
		times<-c(times, t)
	}
	return (times)
}


days<-1000
days<-1

# ------------------------
# First item

catarinaArrival<-5+45/60
waitTimes<-NULL
for (i in 1:days) {
	times <- genBusTimes(start=5)

	# Of course, Catarina will take the 
	# earliest bus after her arrival, always.
	earliestBusTime<-min(times[times >= catarinaArrival])

	curWaitTime <- earliestBusTime-catarinaArrival

	# Wait time, in hours
	waitTimes<-c(waitTimes, curWaitTime)
}

# Catarina's Average wait time in one thousand days
# (in hours)
mean(waitTimes)
# Result: ~0.1871888 h



days<-1000

# ------------------------
# 2nd question Item
catarinaArrival<-8+30/60
intervalTimes<-NULL
waitTimes<-NULL
for (i in 1:days) {
	times <- genBusTimes(start=5)

	# Of course, Catarina will take the 
	# earliest bus after her arrival, always.
	earliestBusIndex<-which(min(
		times[times >= catarinaArrival]) == times)

	# First order difference between the time of the
	# Catarina's bus and the immediately previous bus.
	intTime <- times[earliestBusIndex] - 
		times[earliestBusIndex-1]

	curWaitTime <- times[earliestBusIndex] - 
		catarinaArrival

	# Interval cardinality, in hours
	intervalTimes<-c(intervalTimes, intTime)

	# Wait time, in hours
	waitTimes<-c(waitTimes, curWaitTime)
}

# Average value of first order difference
# between Catarina's bus and the previous one
intAvg<-mean(intervalTimes)
# Result: ~0.3638214 h

# ------------------------
# 3.
2.0*mean(waitTimes)	 # ~ 0.3607519 h
intAvg			 # ~ 0.3638214 h

# The two values are approximately the same.
