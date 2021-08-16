sink("cjs-tsm.jags")
cat("
model {
# Priors and constraints
for (t in 1:(n.occasions-1)){
logit(phi.trans[t]) <- mu.t + epsilon.t[t]
epsilon.t[t] ~ dnorm(0, tau.t)T(-15,15) # Range restriction
phi.t[t] <- 1/(1+exp(-mu.t-epsilon.t[t]))

logit(phi.app[t]) <- mu + iv[1]*beta[1]*dep.mass[t] + iv[2]*beta[2]*snow.cov[t] + epsilon[t]
epsilon[t] ~ dnorm(0, tau)T(-15,15) # Range restriction
phi[t] <- 1/(1+exp(-mu-iv[1]*beta[1]*dep.mass[t]-iv[2]*beta[2]*snow.cov[t]-epsilon[t]))

logit(psight[t]) <- mu.p + epsilon.p[t]
epsilon.p[t] ~ dnorm(0, tau.p)T(-15,15) # Range restriction
p[t] <- 1/(1+exp(-mu.p-epsilon.p[t]))
}

mu.t <- log(mean.phi.t / (1-mean.phi.t))
mean.phi.t ~ dunif(0, 1) # Prior for mean transience
sigma.t ~ dunif(0, 5) # Prior on sd of temp. var
tau.t <- pow(sigma.t, -2)
sigma.t2 <- pow(sigma.t, 2)

mu <- log(mean.phi / (1-mean.phi))
mean.phi ~ dunif(0, 1) # Prior for mean survival
sigma ~ dunif(0, 5) # Prior on sd of temp. var
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)

for (i in 1:2){
beta[i] ~ dnorm(0, 0.001)T(-10, 10) # Prior for slope parameter
iv[i] ~ dbern(0.5) # Prior for indicator variable
}

mu.p <- log(mean.p / (1-mean.p))
mean.p ~ dunif(0, 1) # Prior for mean survival
sigma.p ~ dunif(0, 5) # Prior on sd of temp. var
tau.p <- pow(sigma.p, -2)
sigma.p2 <- pow(sigma.p, 2)

# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
marr.t[t,1:n.occasions] ~ dmulti(pr.t[t,], rel.t[t])
marr[t,1:n.occasions] ~ dmulti(pr[t,], rel[t])
}

# Define the cell probabilities of the m-arrays
# Main diagonal
for (t in 1:(n.occasions-1)){
q[t] <- 1-p[t] # Probability of non-recapture
pr.t[t,t] <- phi.trans[t]*p[t]
pr[t,t] <- phi[t]*p[t]

# Above main diagonal
for (j in (t+1):(n.occasions-1)){
pr.t[t,j] <- phi.trans[t]*prod(phi[(t+1):j])*prod(q[t:(j-1)])*p[j]
pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
} # j

# Below main diagonal
for (j in 1:(t-1)){
pr.t[t,j] <- 0
pr[t,j] <- 0
} # j
} # t

# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
pr.t[t,n.occasions] <- 1-sum(pr.t[t,1:(n.occasions-1)])
pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
} # t

# Assess model fit using Freeman-Tukey statistic
# Compute fit statistics for observed data
for (t in 1:(n.occasions-1)){
   for (j in 1:n.occasions){
      expmarr[t,j] <- rel[t]*pr[t,j]
      expmarr.t[t,j] <- rel.t[t]*pr.t[t,j]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      E.org.t[t,j] <- pow((pow(marr.t[t,j], 0.5)-pow(expmarr.t[t,j], 0.5)), 2)
      } #j
   } #t
   
# Generate replicate data and compute fit stats from them
for (t in 1:(n.occasions-1)){
   marr.new[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   marr.t.new[t,1:n.occasions] ~ dmulti(pr.t[t, ], rel.t[t])
   for (j in 1:n.occasions){
      E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      E.t.new[t,j] <- pow((pow(marr.t.new[t,j], 0.5)-pow(expmarr.t[t,j], 0.5)), 2)
      } #j
   } #t
   
fit <- sum(E.org[,])
fit.t <- sum(E.org.t[,])
fit.new <- sum(E.new[,])
fit.t.new <- sum(E.t.new[,])
}
",fill = TRUE)
sink()