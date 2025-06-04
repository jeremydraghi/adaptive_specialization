setwd("/Users/user/Desktop/outputs")
expName = "E"
dt = read.table(paste0("table_", expName, ".txt"), header=TRUE)
maxTime = 101980
nTs = length(seq(from = 20, to = maxTime, by = 20))
nIds = dim(dt)[1]
ids = 1:nIds

nPrefs = length(unique(dt$pref))
for(pr in 1:nPrefs)
{
  if(nPrefs > 1)
  {
    pref = FALSE
    if(pr == 1) pref = TRUE
    probs = sort(unique(dt$prob[dt$pref == pref]))
  } else {
    probs = sort(unique(dt$prob))
    pref = dt$pref[1]
  }
  nSurvived = rep(0, length(probs))
  diversity = matrix(0, ncol=length(probs), nrow=nTs)
  diversityReps = rep(0, length(probs))
  picky = rep(0, length(probs))
  reps = rep(0, length(probs))

  pIDs = which(dt$pref == pref)
  for(i in pIDs)
  {
    d = try(read.table(paste0("types_", expName, "_", formatC(i, width=3, flag="0"), ".txt"), header=FALSE), silent=TRUE)
    if(is.data.frame(d))
    {
    treat = which(probs == dt$prob[i])
    fits = tapply(d$V3 * d$V4, d$V1, sum) / tapply(d$V3, d$V1, sum)
    if(max(d$V1) == maxTime)
    {
      nSurvived[treat] = nSurvived[treat] + 1
      #diversity[,treat] = diversity[,treat] + tapply(d$V1, d$V1, length)
      #diversityReps[treat] = diversityReps[treat] + 1
      # ds = d[which(d$V1 == 100000),]
      # ranges = rep(0, dim(ds)[1])
      # coverage = rep(0, 4)
      # pickiness = 0
      # for(j in 1:dim(ds)[1])
      # {
      #   pattern = as.integer(strsplit(as.character(ds$V2[j]), "")[[1]])
      #   if(length(pattern) < 4) pattern = c(rep(0, 4 - length(pattern)), pattern)
      #   coverage = coverage + (pattern+1) %% 2
      #   ranges[j] = 4 - sum(pattern)
      #   pickiness = sum(pattern) * ds$V3[j]
      # }
      # picky[treat] = picky[treat] + pickiness / sum(ds$V3)
      # if(nPrefs > 1) 
      # {
      #   write.table(cbind(rep(dt$prob[i], length(ranges)), ranges, ds$V4), paste0("endpoint_fitnesses_",expName,"_", pr, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
      #   write.table(cbind(dt$prob[i], sum(coverage > 0)), paste0("endpoint_coverage_",expName,"_", pr, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
      # } else {
      #   write.table(cbind(rep(dt$prob[i], length(ranges)), ranges, ds$V4), paste0("endpoint_fitnesses_",expName,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
      #   write.table(cbind(dt$prob[i], sum(coverage > 0)), paste0("endpoint_coverage_",expName,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
      # }
      
    } else {
      if(tail(d, 1)[3] > 300)
      {
        print(i)
        print(tail(d, 3))
      }
    }
    reps[treat] = reps[treat] + 1
    #d2 = read.table(paste0("fits_", expName, "_", formatC(i, width=3, flag="0"), ".txt"), header=FALSE)
    #points(d2$V1, rowMeans(d2[,3:6]), type="l", lty="dotted", col=cols[treat]) 
    }
  }
  
  if(nPrefs > 1) 
  {
    write.table(cbind(probs, nSurvived, reps), paste0("survival_", expName,"_", pr, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
  } else {
    write.table(cbind(probs, nSurvived, reps), paste0("survival_", expName, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
  }
}

