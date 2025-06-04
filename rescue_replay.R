setwd("/Users/draghi/Desktop/outputs/")
expName = "K"
dt = read.table(paste0("table_", expName, ".txt"), header=TRUE)
sample = read.table(paste0("types_", expName, "_", formatC(1, width=3, flag="0"), ".txt"), header=FALSE)
maxTime = 101980
nTs = length(seq(from = 20, to = maxTime, by = 20))
nIds = dim(dt)[1]
ids = 1:nIds
probs = sort(unique(dt$prob))
target = 15
foundList = rep(0, length(probs))

outfile = paste0("replay_", expName, ".txt")

for(i in ids)
{
  if(dt$pref[i] == TRUE)
  {
    d = read.table(paste0("types_", expName, "_", formatC(i, width=3, flag="0"), ".txt"), header=FALSE)
    if(max(d$V1) == maxTime)
    {
      types = d$V2[d$V1 == maxTime]
      counts = d$V3[d$V1 == maxTime]
      if((any(types == 0) == FALSE) || counts[which(types == 0)] < 10000)
      {
        firstType = d$V2[which(d$V2 > 0 & d$V3 >= 500)[1]]
        dx = d[which(d$V2 == firstType),]
        firstX = which(dx$V3 >= 500)[1]
        initialCount = dx$V3[firstX]
        firstTime = dx$V1[firstX]
        
        obs = firstX
        cursor = firstX
        maxLength = dim(dx)[1]
        chain = TRUE
        while(chain)
        {
          if(cursor < maxLength & dx$V1[cursor + 1] - dx$V1[cursor] == 20)
          {
            cursor = cursor + 1
            obs = c(obs, cursor)
          } else {
            chain = FALSE
          }
        }
        
        summedCounts = sum(dx$V3[obs])
        plot(dx$V1[obs], dx$V3[obs], main= paste0(i, ", ", summedCounts))
        whichProb = which(probs == dt$prob[i])
        if(summedCounts > 500000 & foundList[whichProb] < target & whichProb > 2)
        {
          foundList[whichProb] = foundList[whichProb] + 1
          write.table(cbind(i, dt[i,], firstTime, formatC(firstType, width = 4, format = "d", flag = "0"), initialCount), outfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
        }
      }
    }
  }
}
