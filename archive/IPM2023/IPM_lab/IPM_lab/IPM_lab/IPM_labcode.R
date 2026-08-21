
################ Integrated Population Model ##########################




### Load Data ----------------- 


library(IPMbook)
library(jagsUI)

data(woodchat5)
str(woodchat5)



### Bundle data and produce data overview -------------------


marr <- marrayAge(woodchat5$ch, woodchat5$age)

jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
                  rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=woodchat5$repro[,1],
                  year=woodchat5$repro[,2], age=woodchat5$repro[,3], C=woodchat5$count, pNinit=dUnif(1, 300))
str(jags.data)



### Write JAGS model file -----------------------



cat(file="model4.txt", "
model {
# Priors and linear models
for (t in 1:(n.occasions-1)){
logit.sj[t] ~ dnorm(mu.sj, tau.sj)
sj[t] <- ilogit(logit.sj[t]) # Back-transformation from logit scale
logit.sa[t] ~ dnorm(mu.sa, tau.sa)
sa[t] <- ilogit(logit.sa[t]) # Back-transformation from logit scale
p[t] <- mean.p
}

for (t in 1:n.occasions){
log.f[1,t] ~ dnorm(mu.f[1], tau.f[1])
f[1,t] <- exp(log.f[1,t]) # Back-transformation from log scale
log.f[2,t] ~ dnorm(mu.f[2], tau.f[2])
f[2,t] <- exp(log.f[2,t]) # Back-transformation from log scale
}
mean.sj ~ dunif(0, 1)
mu.sj <- logit(mean.sj) # Logit transformation
mean.sa ~ dunif(0, 1)
mu.sa <- logit(mean.sa) # Logit transformation
sigma.sj ~ dunif(0, 3)
tau.sj <- pow(sigma.sj, -2)
sigma.sa ~ dunif(0, 3)
tau.sa <- pow(sigma.sa, -2)
for (j in 1:2){
mean.f[j] ~ dunif(0, 10)
mu.f[j] <- log(mean.f[j]) # Log transformation
sigma.f[j] ~ dunif(0, 3)
tau.f[j] <- pow(sigma.f[j], -2)
}
mean.p ~ dunif(0, 1)
sigma ~ dunif(0.5, 100)
tau <- pow(sigma, -2)
# Population count data (state-space model)
# Model for initial stage-spec. population sizes: discrete uniform priors
N[1,1] ~ dcat(pNinit)
N[2,1] ~ dcat(pNinit)
# Process model over time: our model of population dynamics
for (t in 1:(n.occasions-1)){
N[1,t+1] ~ dpois(N[1,t] * f[1,t] / 2 * sj[t] + N[2,t] * f[2,t] / 2 * sj[t])
N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
}

# Observation model
for (t in 1:n.occasions){
C[t] ~ dnorm(N[1,t] + N[2,t], tau)
}

# Productivity data (Poisson regression model)
for (i in 1:length(J)){
J[i] ~ dpois(f[age[i],year[i]])
}

# Capture-recapture data (CJS model with multinomial likelihood)
# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
}

# Define the cell probabilities of the m-arrays
for (t in 1:(n.occasions-1)){
# Main diagonal
q[t] <- 1 - p[t] # Probability of non-recapture
pr.j[t,t] <- sj[t] * p[t]
pr.a[t,t] <- sa[t] * p[t]
# Above main diagonal
for (j in (t+1):(n.occasions-1)){
pr.j[t,j] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
pr.a[t,j] <- prod(sa[t:j]) * prod(q[t:(j-1)]) * p[j]
} #j
# Below main diagonal
for (j in 1:(t-1)){
pr.j[t,j] <- 0
pr.a[t,j] <- 0
} #j
} #t

# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
}

# Derived parameters
# Annual population growth rate (added 0.001 to avoid possible division by 0)
for (t in 1:(n.occasions-1)){
ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t] + 0.001)
}
# Total population size
for (t in 1:n.occasions){
Ntot[t] <- N[1,t] + N[2,t]
}
}
")



# Run the Model -----------------


# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}
# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.f", "mean.p", "sigma.sj", "sigma.sa", "sigma.f",
                "sigma", "sj", "sa", "f", "N", "ann.growth.rate", "Ntot")
# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 2; na <- 1000
# Call JAGS (ART 1 min), check convergence and summarize posteriors
out4 <- jags(jags.data, inits, parameters, "model4.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
             n.thin=nt, n.adapt=na, parallel=TRUE)



# Results ---------------------


#traceplot(out4) # Warning: there are a lot of estimated parameters
print(out4, 3)



mag <- 1.25
cex.tif <- mag * 1.25
lwd.tif <- 3 * mag
op <- par(mar=c(4, 4, 3, 0), las=1, cex=cex.tif, lwd=lwd.tif)
u <- col2rgb("grey82")
T <- length(woodchat5$count)
col.pol <- rgb(u[1], u[2], u[3], alpha=100, maxColorValue=255)
plot(out4$mean$Ntot, type="n",
     ylim=range(c(out4$q2.5$Ntot, out4$q97.5$Ntot, woodchat5$count)),
     ylab="Population size", xlab="Year", las=1, cex=1.5, axes=FALSE)
axis(2, las=1, lwd=lwd.tif)
axis(2, at=c(90, 110, 130, 150), labels=NA, tcl=-0.25, lwd=lwd.tif)
axis(1, at=1:T, labels=NA, tcl=-0.25, lwd=lwd.tif)
axis(1, at=c(5, 10, 15, 20), labels=c(5, 10, 15, 20), tcl=-0.5, lwd=lwd.tif)
polygon(c(1:T, T:1), c(out4$q2.5$Ntot, out4$q97.5$Ntot[T:1]), border=NA, col=col.pol)
points(out4$mean$Ntot, type="b", col="black", pch=16, lty=1, lwd=lwd.tif)
points(woodchat5$count, type="b", col="blue", pch=1, lty=2, lwd=lwd.tif)
legend("topleft", legend=c("Observed population counts", "Estimated population size"),
       pch=c(1, 16), lwd=c(lwd.tif, lwd.tif), col=c("blue", "black"), lty=c(2, 1), bty="n")
par(op)
