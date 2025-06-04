library(viridis)
setwd("/Users/user/Desktop/outputs")
expName = "A"
dt = read.table(paste0("replay_", expName, ".txt"), header=FALSE, col.names = c("id", "selection", "prob", "mu", "pref", "seed", "time", "pheno", "initial"))
nIds = dim(dt)[1]
ids = 1:nIds
probs = sort(unique(dt$prob))
cols = magma(length(probs)+1)[-(length(probs)+1)]
par(mfrow=c(1,2))
plot(c(0.5, 35), c(0.5, 35), type="l", asp=1, xlab="Eco", ylab="Eco + evo", log="xy")
points(rep(1,2), c(0.5, 35), type="l", lty="dotted")
points( c(0.5, 35),rep(1,2), type="l", lty="dotted")

# timecodes = rep(0, 162)
# countcodes = rep(0, 162)
# for(i in ids)
# {
#   v = read.table(paste0("validation_", expName, "_", formatC(i, width=3, flag="0"), ".txt"), header=FALSE)
#   timecodes[i] = v$V1[1]
#   countcodes[i] = v$V3[which(v$V3 >= 500 & v$V3 < 2000)]
# }

for(i in 1:162)
{
  d = read.table(paste0("trials_", expName, "_", formatC(i, width=3, flag="0"), ".txt"), header=FALSE)
  base = dt$initial[i]
  evoRatio = mean(d$V3[d$V2 == 0]) / base
  ecoRatio = mean(d$V3[d$V2 == 1]) / base
  treat = which(probs == dt$prob[i])
  col = "black"
  if(mean(d$V3[d$V2 == 1]) > 14000) col="gray"
  points(ecoRatio, evoRatio, pch=16, col=col, cex=0.8)
  
  sdX = sd(d$V3[d$V2 == 1]) / base
  sdY = sd(d$V3[d$V2 == 0]) / base
  seX = 1.96 * sdX / sqrt(500)
  seY = 1.96 * sdY / sqrt(500)
  points(rep(ecoRatio, 2), evoRatio + c(-1, 1) * seY, type="l", col=col)
  points(ecoRatio + c(-1,1) * seX, rep(evoRatio,2), type="l", col=col)
}
text(0.58, 30, "A", cex=1.6)

plot(range(probs), c(0.75, 3.25), type="n", xlab="Prob. env. change", ylab="Response ratio")
text(2.2e-5, 3.15, "B", cex=1.6)
points(range(probs), rep(1, 2), type="l")
xs = NULL
ys = NULL
reps = rep(0, 162)
for(i in 1:162)
{
  d = read.table(paste0("trials_", expName, "_", formatC(i, width=3, flag="0"), ".txt"), header=FALSE)
  base = dt$initial[i]
  evoRatio = mean(d$V3[d$V2 == 0] / base)
  ecoRatio = mean(d$V3[d$V2 == 1] / base)
  reps[i] = length(d$V3[d$V2 == 0])
  finalCount =  mean(d$V3[d$V2 == 1])
  treat = which(probs == dt$prob[i])
  if(finalCount <= 14000)
  {
    points(dt$prob[i] + rnorm(1, 0, 2.5e-7), evoRatio / ecoRatio, pch=16, cex=0.85)
    xs = c(xs, dt$prob[i])
    ys = c(ys, evoRatio / ecoRatio)
  }
}
summary(lm(ys ~ xs))
print(reps)
