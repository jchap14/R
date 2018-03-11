empty <- c()

for (i in 1:nrow(log2.rlog.counts)) {row <- c(mean(log2.rlog.counts[i,1:2])) 
empty <-rbind(empty, row)}

log2.rlog.counts <- cbind(log2.rlog.counts,empty)