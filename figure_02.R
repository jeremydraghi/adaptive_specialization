library(prevalence)

trials = 100

logit = function(p)
{
  return(log(p / (1 - p)))
}

findInflection = function(xs, ys)
{
  mod = glm(ys ~ xs, family=binomial(link = "logit"))
  return(-mod$coefficients[1]/mod$coefficients[2])
}

d_np1 = read.table("/Users/draghi/Desktop/outputs/survival_A_2.txt", header=FALSE, col.names=c("prob", "survived", "reps"))
d_np2 = read.table("/Users/draghi/Desktop/outputs/survival_K_2.txt", header=FALSE, col.names=c("prob", "survived", "reps"))

d_np = data.frame(prob = sort(unique(d_np1$prob, d_np2$prob)))
d_np$survived = rep(0, dim(d_np)[1])
d_np$reps = rep(0, dim(d_np)[1])
for(i in 1:dim(d_np)[1])
{
  reps = 0
  survived = 0
  x1 = which(d_np1$prob == d_np$prob[i])
  x2 = which(d_np2$prob == d_np$prob[i])
  if(length(x1) == 1)
  {
    reps = reps + d_np1$reps[x1]
    survived = survived + d_np1$survived[x1]
  }
  if(length(x2) == 1)
  {
    reps = reps + d_np2$reps[x2]
    survived = survived + d_np2$survived[x2]
  }
  d_np$reps[i] = reps
  d_np$survived[i] = survived
}

d_p1 = read.table("/Users/draghi/Desktop/outputs/survival_A_1.txt", header=FALSE, col.names=c("prob", "survived", "reps"))
d_p2 = read.table("/Users/draghi/Desktop/outputs/survival_K_1.txt", header=FALSE, col.names=c("prob", "survived", "reps"))

d_p = data.frame(prob = sort(unique(d_p1$prob, d_p2$prob)))
d_p$survived = rep(0, dim(d_p)[1])
d_p$reps = rep(0, dim(d_p)[1])
for(i in 1:dim(d_p)[1])
{
  reps = 0
  survived = 0
  x1 = which(d_p1$prob == d_p$prob[i])
  x2 = which(d_p2$prob == d_p$prob[i])
  if(length(x1) == 1)
  {
    reps = reps + d_p1$reps[x1]
    survived = survived + d_p1$survived[x1]
  }
  if(length(x2) == 1)
  {
    reps = reps + d_p2$reps[x2]
    survived = survived + d_p2$survived[x2]
  }
  d_p$reps[i] = reps
  d_p$survived[i] = survived
}

d_s = read.table("/Users/draghi/Desktop/outputs/survival_C.txt", header=FALSE, col.names=c("prob", "survived", "reps"))
# Remove excessive treatments with low replication
d_s = d_s[1:18,]
probs = sort(unique(c(d_np$prob, d_p$prob, d_s$prob)))



par(mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(1.5, 0.5, 0), cex.lab=1.3, cex.axis=1)
plot(c(0, 0.00013), c(0,1), type="n", xlab="Prob. of environmental change", ylab="Proportion survived")
for(i in 1:length(d_np$prob))
{
  x = propCI(d_np$survived[i], d_np$reps[i], method="jeffreys")
  points(rep(d_np$prob[i],2)-2.5e-7, c(x$lower, x$upper), type="l", col="orange")
}
#points(d_np$prob-2.5e-7, d_np$survived / d_np$reps, pch=16, col="white", cex=0.8)
points(d_np$prob-2.5e-7, d_np$survived / d_np$reps, type="b", pch=16, col="orange")
xs = rep(d_np$prob, times = d_np$reps)
ys = NULL
for(i in 1:dim(d_np)[1])
{
  ys = c(ys, rep(1, d_np$survived[i]), rep(0, d_np$reps[i] - d_np$survived[i]))
}
p50Gen = findInflection(xs, ys)
#points(rep(p50Gen,2), c(0, 1), type="l", lty="dotted", col="orange")
ps = rep(0, trials)
for(i in 1:trials)
{
  resample = sample(length(xs), length(xs), replace=TRUE)
  ps[i] = findInflection(xs[resample], ys[resample])
}
print("Generalists: ") 
print(p50Gen)
print(quantile(ps, probs = c(0.025, 0.975)))

for(i in 1:length(d_p$prob))
{
  x = propCI(d_p$survived[i], d_p$reps[i], method="jeffreys")
  points(rep(d_p$prob[i],2)+2.5e-7, c(x$lower, x$upper), type="l", col="dodgerblue")
}
points(d_p$prob+2.5e-7, d_p$survived / d_p$reps, type="b", col="dodgerblue", pch=16)
xs = rep(d_p$prob, times = d_p$reps)
ys = NULL
for(i in 1:dim(d_p)[1])
{
  ys = c(ys, rep(1, d_p$survived[i]), rep(0, d_p$reps[i] - d_p$survived[i]))
}
p50Evo = findInflection(xs, ys)
#points(rep(p50Evo,2), c(0, 1), type="l", lty="dotted", col="dodgerblue")
ps = rep(0, trials)
for(i in 1:trials)
{
  resample = sample(length(xs), length(xs), replace=TRUE)
  ps[i] = findInflection(xs[resample], ys[resample])
}
print("Evolving: ") 
print(p50Evo)
print(quantile(ps, probs = c(0.025, 0.975)))

for(i in 1:length(d_s$prob))
{
  x = propCI(d_s$survived[i], d_s$reps[i], method="jeffreys")
  points(rep(d_s$prob[i],2)+5e-7, c(x$lower, x$upper), type="l", col="chartreuse3")
}
points(d_s$prob+5e-7, d_s$survived / d_s$reps, type="b", col="chartreuse3", pch=16)
xs = rep(d_s$prob, times = d_s$reps)
ys = NULL
for(i in 1:dim(d_s)[1])
{
  ys = c(ys, rep(1, d_s$survived[i]), rep(0, d_s$reps[i] - d_s$survived[i]))
}
p50Sp = findInflection(xs, ys)
#points(rep(p50Sp,2), c(0, 1), type="l", lty="dotted", col="chartreuse3")
ps = rep(0, trials)
for(i in 1:trials)
{
  resample = sample(length(xs), length(xs), replace=TRUE)
  ps[i] = findInflection(xs[resample], ys[resample])
}
print("Specialists: ") 
print(p50Sp)
print(quantile(ps, probs = c(0.025, 0.975)))

legend("topright", c("Generalists", "Niche evolves", "Specialists"), col=c("orange", "dodgerblue", "chartreuse3"), lwd=2, pch=16)

dh = read.table("/Users/draghi/Desktop/outputs/survival_X.txt", header=FALSE, col.names=c("prob", "survived", "reps"))
points(dh$prob, dh$survived / dh$reps, col="darkorchid")
for(i in 1:length(dh$prob))
{
  x = propCI(dh$survived[i], dh$reps[i], method="jeffreys")
  points(rep(dh$prob[i],2), c(x$lower, x$upper), type="l", col="darkorchid")
}
