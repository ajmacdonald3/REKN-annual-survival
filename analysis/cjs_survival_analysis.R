################################################################################
# REKN SURVIVAL ANALYSIS IN JAGS
#
################################################################################

# load packages

#library(RMark)
library(tidyverse)
library(lubridate)
library(reshape)
library(bb2enchist)
library(R2ucare)
library(cowplot)
library(corrr)
library(R2jags)
library(jagsUI)

################################################################################

# GOF testing with R2ucare
# read in and format data
rekn <- read_inp("./data/rekn_enchist.inp")

rekn_hist <- rekn$encounter_histories
rekn_freq <- rekn$sample_size

# run GOF tests

test3sr(rekn_hist, rekn_freq) # transience - significant
test3sm(rekn_hist, rekn_freq) # resighting rates - not significant
test2ct(rekn_hist, rekn_freq) # equal resightability - significant
test2cl(rekn_hist, rekn_freq) # equal resightability (before and after) - not significant
overall_CJS(rekn_hist, rekn_freq) # significant

# results suggest fitting CJS model incorporating transience
# some suggestion for trap-dependence effect
# calculate new GOF test for this model

# subtract the components 3SR and 2CT to the CJS test statistic
stat_new <- overall_CJS(rekn_hist, rekn_freq)$chi2 - (test3sr(rekn_hist, rekn_freq)$test3sr[[1]])
#+  test2ct(rekn_hist, rekn_freq)$test2ct[[1]])

# calculate degrees of freedom associated with the new test statistic
df_new <- overall_CJS(rekn_hist, rekn_freq)$degree_of_freedom - (test3sr(rekn_hist, rekn_freq)$test3sr[[2]])
#+  test2ct(rekn_hist, rekn_freq)$test2ct[[2]])

# compute p-value
1 - pchisq(stat_new, df_new)

# chat = chi-square/df = 52.11/31 = 1.681
# adjust for estimated overdispersion

################################################################################

# load data
rekn_enchist <- readRDS("data/rekn_enchist.rds")
std_covs <- readRDS("data/survival_covariates_standardized.rds")

# convert BirdID to row names so a numeric matrix can be generated
rownames(rekn_enchist) <- NULL

enchist <- rekn_enchist %>% 
  column_to_rownames(var = "FlagID")

enchist <- sapply(enchist, as.numeric)

# function to create a m-array based on capture-histories (CH)
m_array <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

# Create separate m-arrays for first capture and subsequent recaptures
CH <- enchist

cap <- apply(CH, 1, sum)
ind <- which(cap >= 2)
CH.R <- CH[ind,] # First period CH recaptured at least once
CH.N <- CH[-ind,] # First period CH never recaptured
# Remove first capture
first <- numeric()
for (i in 1:dim(CH.R)[1]){
  first[i] <- min(which(CH.R[i,]==1))
}
CH.R1 <- CH.R
for (i in 1:dim(CH.R)[1]){
  CH.R1[i,first[i]] <- 0
}
# Create m-array of those recaptured at least once
CH.marray <- m_array(CH.R1)
# Create CH matrix for first period, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.R1)[1]){
  second[i] <- min(which(CH.R1[i,]==1))
}
CH.R2 <- matrix(0, nrow = dim(CH.R)[1], ncol = dim(CH.R)[2])
for (i in 1:dim(CH.R)[1]){
  CH.R2[i,first[i]] <- 1
  CH.R2[i,second[i]] <- 1
}
# Create m-array for these
CH.R.marray <- m_array(CH.R2)
# The last column ought to show the number of transients not recaptured
# again and should all be zeros, since all of them are released as "residents"
CH.R.marray[,dim(CH)[2]] <- 0
# Create the m-array for transients never recaptured and add it to the
# previous m-array
CH.N.marray <- m_array(CH.N)
CH.T.marray <- CH.R.marray + CH.N.marray

# pull out covariates
dep_mass <- std_covs %>% 
  dplyr::select(dep_mass_std) %>% pull()

snow_cov <- std_covs %>% 
  dplyr::select(mean_snowc_std) %>% pull()

################################################################################

# source JAGS model code
source("./analysis/JAGS_cjs_tsm.R")

# Bundle data
jags.data <- list(marr.t = CH.T.marray, marr = CH.marray, n.occasions = dim(CH)[2],
                  rel.t = rowSums(CH.T.marray), rel = rowSums(CH.marray),
                  dep.mass = dep_mass, snow.cov = snow_cov)

# Initial values
inits <- function(){list(mean.phi.t = runif(1, 0, 1), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1),
                         sigma.t = runif(1, 0, 5), sigma = runif(1, 0, 5), sigma.p = runif(1, 0, 5),
                         iv = rep(1, 2), beta = runif(2, -5, 5))}  

# Parameters monitored
parameters <- c("mean.phi.t", "phi.t", "sigma.t2", "mean.phi", "phi", "sigma2",
                "mean.p", "p", "sigma.p2",
                "beta", "iv",
                "fit", "fit.new", "fit.t", "fit.t.new")

# MCMC settings
ni <- 100000
nt <- 10
nb <- 50000
nc <- 3

# Call JAGS from R (BRT 1 min)
cjs <- jagsUI::jags(jags.data, inits, parameters, "cjs-tsm.jags",
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(cjs, digits = 3)

saveRDS(cjs$summary, file = paste0("./analysis-output/cjs-tsm-summary", Sys.Date(), ".rds"))
write.csv(cjs$summary, file = paste0("./analysis-output/cjs-tsm-summary", Sys.Date(), ".csv"))

saveRDS(cjs$sims.list, file = paste0("./analysis-output/cjs-tsm-simslist", Sys.Date(), ".rds"))

# traceplots
pdf(file = paste0("./analysis-output/traceplots/cjs-tsm-traceplots", Sys.Date(), ".pdf"))
par(mfrow = c(3, 1))
traceplot(cjs)
dev.off()

# parameter identifiability checks
sims.list <- cjs$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/",
                        i,"-","check.png"),
      width=4, height=3, units="in", res=600)
  
  print(ggplot(sims.list, aes(sims.list[,i])) +
          geom_density() +
          #geom_hline(yintercept = 1, linetype = "dashed") +
          xlab(i))
  
  dev.off()
  
}

# evaluation of fit
mean(sims.list$fit.new > sims.list$fit)
mean(sims.list$fit.t.new > sims.list$fit.t)

ppcheck <- ggplot(sims.list, aes(x = fit, y = fit.new)) +
  geom_point(alpha = 0.3, pch = 16) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Discrepancy actual data") +
  ylab("Discrepancy ideal data") +
  theme(legend.position = "none")

png(filename = "figures/cjs-tsm-ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

ppcheck.t <- ggplot(sims.list, aes(x = fit.t, y = fit.t.new)) +
  geom_point(alpha = 0.3, pch = 16) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Discrepancy actual data") +
  ylab("Discrepancy ideal data") +
  theme(legend.position = "none")

png(filename = "figures/cjs-tsm-ppcheck-t.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck.t)

dev.off()

################################################################################

# plot results

windowsFonts(Times=windowsFont("TT Times New Roman"))

# set custom theme for all plots
theme_cust <- function() {
  theme_classic(base_family = "Times") %+replace%
    theme(axis.title.x = element_text(size=10),
          axis.text.x  = element_text(size=8, colour = "black"),
          axis.title.y = element_text(size=10, angle = 90, margin = margin(t = 0, r = 5, b = 0, l = 0)),
          axis.text.y = element_text(size=8, colour = "black"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          strip.text.x = element_text(size=12, face = "bold"),
          legend.text = element_text(size=8),
          legend.key.height = unit(1, "mm"),
          plot.title = element_text(size = 12, hjust = 0, vjust = 1.5),
          #panel.border = element_rect(size =0.5, fill = "transparent"),
          plot.margin = margin(10, 10, 10, 15))
}

# format data
cjs_res <- readRDS("./analysis-output/cjs-tsm-summary2021-08-14.rds")

cjs_res <- cjs_res %>% 
  as.data.frame %>% 
  rownames_to_column() %>% 
  dplyr::rename(parameter = rowname) %>% 
  dplyr::rename(lcl = `2.5%`) %>% 
  dplyr::rename(median = `50%`) %>% 
  dplyr::rename(ucl = `97.5%`) %>% 
  select(parameter, mean, lcl, median, ucl)

phi_res <- cjs_res %>% 
  filter(str_detect(parameter, "phi\\[")) %>% 
  mutate(year = 2009:2017) %>% 
  mutate(var = "phi")
  
# apparent annual survival
phi_plot <- ggplot() +
  geom_rect(data = filter(cjs_res, parameter == "mean.phi"),
                    aes(xmin=-Inf, xmax=Inf, ymin=lcl, ymax=ucl), fill = "grey50", alpha=0.3) +
  geom_hline(data = filter(cjs_res, parameter == "mean.phi"),
             aes(yintercept=mean)) +
  geom_errorbar(data = phi_res, aes(x=as.factor(year), ymin=lcl, ymax=ucl),
                width=0, size=0.5, colour="black", linetype=1) +
  geom_line(data = phi_res, aes(x=as.factor(year), y=mean, group = var),
            linetype="dashed", size=0.5) +
  geom_point(data = phi_res, aes(x=as.factor(year), y=mean),
             size=3) +
  #scale_fill_manual(values = c("black", "white")) +
  #scale_x_continuous(labels = function(x) format(as.Date(x, origin = "1970-01-01"), "%d %b")) +
  ylim(0, 1) +
  ylab("Annual apparent survival") +
  theme_cust() +
  theme(#legend.position = c(0.4,0.15),
        #legend.title = element_blank(),
        #legend.background = element_blank(),
        axis.title.x = element_blank())

png(filename = paste0("figures/phi-plot.png"),
    width=6, height=4, units="in", res=600)

plot(phi_plot)

dev.off()

