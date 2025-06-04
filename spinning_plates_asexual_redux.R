library(doParallel)
library("flock")
library("fastmatch")

calcFits = function(geno, opt)
{
  sames = sum(geno == opt)
  if(selection == "birth") return(c(fitMap[1 + sames], baseDeath))
  if(selection == "death") return(c(1, baseDeath * fitMap[1 + sames]))
}

# Total number of sites for adults (controls density dependence)
K = 20000
# K sites divided evenly into demes
demes = 4
# Number of sites in each deme
k = K / demes
# Maximum number of time units, relative to birth rate when rare
burnIn = 50000
treatment = 50000
coda = 2000
maxTime = burnIn + treatment + coda

baseProb = 0.000005

# Abundance threshold for a genotype to be recorded
threshold = 100

# Number of loci per deme
L = 25

maxBatch = 200
maxProportion = 0.05

# indices of preference loci
prefLoci = 1:demes
totalLoci = demes * (L + 1)
totalFitLoci = demes * L
fitLoci = (demes + 1):totalLoci

# map of loci assigning them to each environment
subMap = matrix(fitLoci, ncol=demes, byrow=FALSE)

# map of optima assigning them to each environment
optimaMap = matrix(1:totalFitLoci, ncol=demes, byrow=FALSE)

# Number of alleles per locus
A = 20
mutantVector = 1:(A-1)
# Performance increment per match
s = 0.1
fitMap = (1 - s)^(L:0)
# Search cost
c = 0.01

setwd("/Users/user/Desktop/outputs/")

# Base rate of death; base rate of birth is 1.
baseDeath = 0.2

expName = "B"
replayFile = paste0("replay_", expName, ".txt")
dt = read.table(replayFile, header=FALSE, col.names = c("id", "selection", "prob", "mu", "pref", "seed", "time", "pheno", "initial"))
nTreatments = dim(dt)[1]

hubC = makeForkCluster(16)
registerDoParallel(hubC)

temp = foreach(r=1:nTreatments) %dopar%
  {
    # Load parameters from treats and tableFile
    set.seed(dt$seed[r])
    mu = dt$mu[r]
    PREF = dt$pref[r]
    treatmentProb = dt$prob[r]
    selection = dt$selection[r]
    stopTime = dt$time[r]
    U = totalFitLoci * mu
    envProb = baseProb
    
    popfile = paste0("pop_", expName, "_", formatC(r, width=3, flag="0"), ".txt")
    envfile = paste0("env_", expName, "_", formatC(r, width=3, flag="0"), ".txt")
    valfile = paste0("validation_", expName, "_", formatC(r, width=3, flag="0"), ".txt")
    vacantfile = paste0("vacant_", expName, "_", formatC(r, width=3, flag="0"), ".txt")
    envChangesFile = paste0("envChanges_", expName, "_", formatC(r, width=3, flag="0"), ".txt")
    pop = matrix(0, nrow=K, ncol = totalLoci)
    envs = rep(1:demes, each = k)
    birthRates = rep(1, K)
    deathRates = rep(1, K)
    totalBirth = K
    totalDeath = K
    optima = rep(0, totalFitLoci)
    vacant = NULL
    N = K
    
    t = 0
    nextT = 1
    phase = "burnIn"
    transition = burnIn
    
    nBuffer = 1e6
    idBuffer = sample(1:K, nBuffer, replace=TRUE)
    pBuffer = runif(nBuffer)
    cursor = 1
    preChanges = rep(0, demes)
    postChanges = rep(0, demes)
    
    while(N > 1 & t < stopTime)
    {
      if(t >= transition)
      {
        if(phase == "treatment")
        {
          phase = "coda"
          transition = maxTime + 1
          envProb = baseProb
        }
        if(phase == "burnIn")
        {
          phase = "treatment"
          transition = burnIn + treatment
          envProb = treatmentProb
        }
      }
      if(t >= nextT)
      {
        periods = ceiling(t - nextT)
        stopifnot(N == (K - length(vacant)))
        # Check if environments are changing
        toChange = rbinom(demes, L, periods * envProb)
        if(any(toChange > 0))
        {
          if(phase == "burnIn") preChanges = preChanges + toChange
          if(phase == "treatment") postChanges = postChanges + toChange
          for(j in 1:demes)
          {
            if(toChange[j] > 0)
            {
              changes = sample(1:L, toChange[j], replace=FALSE)
              for(i in changes)
              {
                optima[optimaMap[i,j]] = (optima[optimaMap[i,j]] + sample(mutantVector,1)) %% A
              }
              ids = (1:K)[-vacant]
              ids = ids[which(envs[ids] == j)]
              for(i in ids)
              {
                newFits = calcFits(pop[i,subMap[,envs[i]]], optima[optimaMap[,envs[i]]])
                birthRates[i] = newFits[1]
                deathRates[i] = newFits[2]
              }
            }
          }
        }
        totalBirth = sum(birthRates)
        totalDeath = sum(deathRates)
        nextT = nextT + 1
      }
      batch = min(c(maxBatch, ceiling((K - length(vacant)) * maxProportion)))
      # Compute waiting time and identity of next event
      # Birth rates are reduced to exclude offspring that die before ever seeing a site; this reduction
      # is compensated below by presenting each offspring with a site before they have a chance to die.
      p = length(vacant) / K
      modifiedBirth = totalBirth * p / (c + p)
      gTotal = modifiedBirth + totalDeath
      t = t + rexp(1, gTotal / batch)
      
      modBirthRates = birthRates *  p / (c + p)
      rates = deathRates + modBirthRates
      rates = rates / max(rates)
      ids = rep(0, batch)
      i = 1
      while(i <= batch)
      {
        candidate = idBuffer[cursor]
        if(pBuffer[cursor] < rates[candidate])
        {
          ids[i] = candidate
          i = i + 1
        }
        cursor = cursor + 1
        if(cursor > nBuffer)
        {
          cursor = 1
          idBuffer = sample(1:K, nBuffer, replace=TRUE)
          pBuffer = runif(nBuffer)
        }
      }
      events = rbinom(batch, 1, modBirthRates[ids] / (modBirthRates[ids] + deathRates[ids]))
      dups = duplicated(ids)
      if(any(dups))
      {
        events = events[-which(dups)]
        ids = ids[-which(dups)]
      }
      
      for(i in 1:length(events))
      {
        event = events[i]
        # Death event
        if(event == 0)  
        {
          dead = ids[i]
          totalBirth = totalBirth - birthRates[dead]
          totalDeath = totalDeath - deathRates[dead]
          birthRates[dead] = 0
          deathRates[dead] = 0
          vacant = c(vacant, dead)
          p = length(vacant) / K
          N = N - 1
        } else {
          # Birth event
          birth = ids[i]
          # Assign fate to numbered spot, or -1 if organism dies during search
          fate = 0
          
          if(length(vacant) > 0)
          {
            # Mutate preferences before search begins
            preference = pop[birth, 1:demes]
            if(PREF)
            {
              prefMutants = rbinom(demes, 1, mu)
              if(any(prefMutants == 1)) preference[which(prefMutants == 1)] = (preference[which(prefMutants == 1)] + 1) %% 2
            }
          } else {
            fate = -1
          }
          
          # Examine sites until organism chooses one (fate == [1..K]) or dies (fate == -1). Don't bother             # mutating performance traits until the organism finds an empty spot
          while(fate == 0)
          {
            if(length(vacant) > 1)
            {
              spot = sample(vacant, 1)
            } else {
              spot = vacant
            }
            if(preference[envs[spot]] == 0)
            {
              fate = spot
            } else {
              if(rbinom(1, 1, c / (c + p)) == 1) fate = -1
            }
          }
          if(fate > 0)
          {
            vacant = vacant[-match(fate, vacant)]
            p = length(vacant) / K
            pop[fate,] = pop[birth,]
            pop[fate,1:demes] = preference
            # Mutation in performance traits
            mutants = rpois(1, U)
            if(mutants > 0)
            {
              hits = sample(fitLoci, mutants, replace=TRUE)
              for(i in hits)
              {
                pop[fate, i] = (pop[fate, i] + sample(mutantVector,1)) %% A
              }
            }
            # Calculate fitness components here
            newFits = calcFits(pop[fate,subMap[,envs[fate]]], optima[optimaMap[,envs[fate]]])
            birthRates[fate] = newFits[1]
            deathRates[fate] = newFits[2]
            totalBirth = totalBirth + birthRates[fate]
            totalDeath = totalDeath + deathRates[fate]
            N = N + 1
          }
        }
      }
    }
    write.table(optima, envfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
    write.table(pop, popfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
    write.table(vacant, vacantfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
    write.table(cbind(preChanges, postChanges), envChangesFile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
    if(N < K)
    {
      testPop = pop[-vacant,]
      testBR = birthRates[-vacant]
      testEnvs = envs[-vacant]
    } else {
      testPop = pop
      testBR = birthRates
      testEnvs = envs
    }
    meanB = rep(0, demes)
    tmp = tapply(testBR, testEnvs, mean)
    meanB[as.integer(names(tmp))] = tmp
    
    codes = apply(testPop[,1:demes], 1, paste0, collapse="")
    codeTable = table(codes)
    sigs = names(codeTable)[which(codeTable >= threshold)]
    toPrint = NULL
    for(ty in sigs)
    {
      toPrint = rbind(toPrint, c(nextT, ty, codeTable[which(names(codeTable) == ty)], mean(testBR[which(codes == ty)])))
    }
    write.table(toPrint, valfile, row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)
  }  
stopCluster(hubC)




