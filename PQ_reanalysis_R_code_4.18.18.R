 #Set seed for random num generator to allow reproduction of results
set.seed(8675309)

#POST HOC DESIGN ANALYSIS

#install.packages("Exact")
library(Exact)

#Table 1, Panel A
#Column 1 (maximum possible)
power.exact.test(.14, .00, 96, 60, alpha=0.05, method="fisher")
# 0.99
#Column 2 (strong)
power.exact.test(.11, .03, 96, 60, alpha=0.05, method="fisher")
# 0.37
#Column 3 (moderate)
power.exact.test(.10, .05, 96, 60, alpha=0.05, method="fisher")
# 0.14

#Table 1, Panel B
#Column 1 (maximum possible)
power.exact.test(.23, .00, 22, 134, alpha=0.05, method="fisher")
# 0.98
#Column 2 (strong)
power.exact.test(.14, .01, 22, 134, alpha=0.05, method="fisher")
# 0.66
#Column 3 (moderate)
power.exact.test(.09, .02, 22, 134, alpha=0.05, method="fisher")
# 0.28

#Table 1, Panel C
#Column 1 (maximum possible)
power.exact.test(.13, .00, 124, 32, alpha=0.05, method="fisher")
# 0.66
#Column 2 (strong)
power.exact.test(.12, .03, 124, 32, alpha=0.05, method="fisher")
# 0.22
#Column 3 (moderate)
power.exact.test(.11, .06, 124, 32, alpha=0.05, method="fisher")
# 0.06

#Table 1, Panel D
#Column 1 (maximum possible)
power.exact.test(.23, .00, 22, 32, alpha=0.05, method="fisher")
# 0.78
#Column 2 (strong)
power.exact.test(.18, .03, 22, 32, alpha=0.05, method="fisher")
# 0.37
#Column 3 (moderate)
power.exact.test(.14, .06, 22, 32, alpha=0.05, method="fisher")
# 0.11


#BAYESIAN REANALYSIS 

#first, install program JAGS on computer from http://mcmc-jags.sourceforge.net/
  #then do in R:
#install.packages("devtools")
#install.packages("rjags")
#devtools::install_github("rasmusab/bayesian_first_aid")
#library(devtools)

library(rjags)
library(BayesianFirstAid)

#TABLE 2, COLUMN 3 (Proportion Difference & Exact Tests)
#Fisher exact test of equality for two binomial proportions (n1-success, n1-failure, n2-success, n2-failure)

  #Panel 2.A, Col3 (Observed)
  fisher.test(matrix(c(7, 89, 4, 56), ncol=2))
    #OR=1.10, p=1 (fail to reject null of no difference in callbacks by vignette, given observed data)
  
  #Panel 2.B, Col3 (Observed)
  fisher.test(matrix(c(2, 20, 9, 125), ncol=2))
    #OR=1.39, p=.655 (fail to reject null of no difference in callbacks by vignette, given observed data)
  
  #Panel 2.C, Col3 (Observed)
  fisher.test(matrix(c(10, 114, 1, 31), ncol=2))
    #OR=2.71, p=.463 (fail to reject null of no difference in callbacks by vignette, given observed data)

  #Panel 2.D, Col3 (Observed)
  fisher.test(matrix(c(2, 20, 1, 31), ncol=2))
    #OR=3.03, p=.560 (fail to reject null of no difference in callbacks by vignette, given observed data)
  
  
#TABLE 2, COLUMNS 4-6 (Proportion Difference & Exact Tests)
  #Col. 4: Bayes.prop.test estimates  rel freq of callbacks by employer willingness  
  #Col. 5: Also estimates diff in rel freq & HDIs, aka, callback diffs by "more willing" vs "less willing" employers  
  #Col. 6: Also estimates probability that group 1 has higher relative frequency of callbacks than group 2
    #where group 1 = employers reportedly "more willing" to hire in vignette
    #where group 2 = employers reportedly "less willing" to hire in vignette  

  #Panel 2.A, Cols.4-6 (Observed)
  set.seed(8675309)
  n_callbacks2A <- c(7, 4)
  n_employers2A <- c(96, 60)
  bayes.prop.test(n_callbacks2A, n_employers2A)
  fit2A <- bayes.prop.test(n_callbacks2A, n_employers2A)
  plot(fit2A)
  summary(fit2A)
  fit2A$stats

  model.code(fit2A)
  #Code for Estimating Probability at/above counterfactual critical values (Figure 2):
    ### Model code for the Bayesian First Aid  ###
    ### alternative to the test of proportions ###
    require(rjags)
    
    # Setting up the data
    x <- c(7, 4) 
    n <- c(96, 60) 
    
    # The model string written in the JAGS language
    model_string <- "model {
    for(i in 1:length(x)) {
      x[i] ~ dbinom(theta[i], n[i])
      theta[i] ~ dbeta(1, 1)
      x_pred[i] ~ dbinom(theta[i], n[i])
    }
  }"
    
    # Running the model
    model2A <- jags.model(textConnection(model_string), data = list(x = x, n = n), 
                        n.chains = 3, n.adapt=1000)
    samples2A <- coda.samples(model2A, c("theta", "x_pred"), n.iter=5000)
    
    # Inspecting the posterior
    plot(samples2A)
    summary(samples2A)
    
    # Extract the mcmc samples as a matrix and compare the thetas 
    # of the groups. The following shows the median and 95%
    # credible interval for the difference between Group 1 and Group 2.
    samp_mat2A <- as.matrix(samples2A)
    quantile(samp_mat2A[, "theta[1]"] - samp_mat2A[, "theta[2]"], c(0.025, 0.5, 0.975))
    
    #calculate difference in thetas btw grp1 and grp2 from mcmc samples
    diff2A <- (samp_mat2A[, "theta[1]"] - samp_mat2A[, "theta[2]"])
    #plot distribution of diff12
    hist(diff2A)
    #number of diff12 estimates from mcmc samples > 0
    min2A <- length(diff2A[diff2A > 0.00])
    mod2A <- length(diff2A[diff2A >= 0.05])
    str2A <- length(diff2A[diff2A >= 0.08])
    max2A <- length(diff2A[diff2A >= 0.14])
    #total number of diff12 estimates from mcmc samples
    dend2A <- length(diff2A)
    #calc proportion of diff12 > 0
    Pdmin2A <- min2A/dend2A
    Pdmod2A <- mod2A/dend2A
    Pdstr2A <- str2A/dend2A
    Pdmax2A <- max2A/dend2A
    #Probability group 1 > group 2 (i.e., more likely callbacks > less likely callbacks)
    Pdmin2A
    Pdmod2A
    Pdstr2A
    Pdmax2A
  
  
  #Panel 2.B, Cols.4-6 (Observed)
  set.seed(8675309)
  n_callbacks2B <- c(2, 9)
  n_employers2B <- c(22, 134)
  bayes.prop.test(n_callbacks2B, n_employers2B)
  fit2B <- bayes.prop.test(n_callbacks2B, n_employers2B)
  plot(fit2B)
  fit2B$stats
  
  model.code(fit2B)
  #Code for Estimating Probability at/above counterfactual critical values (Figure 2):
  ### Model code for the Bayesian First Aid  ###
    ### alternative to the test of proportions ###
    require(rjags)
    
    # Setting up the data
    x <- c(2, 9) 
    n <- c(22, 134) 
    
    # The model string written in the JAGS language
    model_string <- "model {
    for(i in 1:length(x)) {
      x[i] ~ dbinom(theta[i], n[i])
      theta[i] ~ dbeta(1, 1)
      x_pred[i] ~ dbinom(theta[i], n[i])
    }
  }"
    
    # Running the model
    model2B <- jags.model(textConnection(model_string), data = list(x = x, n = n), 
                          n.chains = 3, n.adapt=1000)
    samples2B <- coda.samples(model2B, c("theta", "x_pred"), n.iter=5000)
    
    # Inspecting the posterior
    plot(samples2B)
    summary(samples2B)
    
    # Extract the mcmc samples as a matrix and compare the thetas 
    # of the groups. The following shows the median and 95%
    # credible interval for the difference between Group 1 and Group 2.
    samp_mat2B <- as.matrix(samples2B)
    quantile(samp_mat2B[, "theta[1]"] - samp_mat2B[, "theta[2]"], c(0.025, 0.5, 0.975))
    
    #calculate difference in thetas btw grp1 and grp2 from mcmc samples
    diff2B <- (samp_mat2B[, "theta[1]"] - samp_mat2B[, "theta[2]"])
    #plot distribution of diff12
    hist(diff2B)
    #number of diff12 estimates from mcmc samples > 0
    min2B <- length(diff2B[diff2B > 0.00])
    mod2B <- length(diff2B[diff2B >= 0.07])
    str2B <- length(diff2B[diff2B >= 0.12])
    max2B <- length(diff2B[diff2B >= 0.23])
    #total number of diff12 estimates from mcmc samples
    dend2B <- length(diff2B)
    #calc proportion of diff12 > 0
    Pdmin2B <- min2B/dend2B
    Pdmod2B <- mod2B/dend2B
    Pdstr2B <- str2B/dend2B
    Pdmax2B <- max2B/dend2B
    #Probability group 1 > group 2 (i.e., more likely callbacks > less likely callbacks)
    Pdmin2B
    Pdmod2B
    Pdstr2B
    Pdmax2B
  
  
  #Panel 2.C, Cols.4-6 (Observed)
  set.seed(8675309)
  n_callbacks2C <- c(10, 1)
  n_employers2C <- c(124, 32)
  bayes.prop.test(n_callbacks2C, n_employers2C)
  fit2C <- bayes.prop.test(n_callbacks2C, n_employers2C)
  plot(fit2C)
  fit2C$stats
  
  model.code(fit2C)
  #Code for Estimating Probability at/above counterfactual critical values (Figure 2):
    ### Model code for the Bayesian First Aid  ###
    ### alternative to the test of proportions ###
    require(rjags)
    
    # Setting up the data
    x <- c(10, 1) 
    n <- c(124, 32) 
    
    # The model string written in the JAGS language
    model_string <- "model {
    for(i in 1:length(x)) {
      x[i] ~ dbinom(theta[i], n[i])
      theta[i] ~ dbeta(1, 1)
      x_pred[i] ~ dbinom(theta[i], n[i])
    }
  }"
    
    # Running the model
    model2C <- jags.model(textConnection(model_string), data = list(x = x, n = n), 
                          n.chains = 3, n.adapt=1000)
    samples2C <- coda.samples(model2C, c("theta", "x_pred"), n.iter=5000)
    
    # Inspecting the posterior
    plot(samples2C)
    summary(samples2C)
    
    # Extract the mcmc samples as a matrix and compare the thetas 
    # of the groups. The following shows the median and 95%
    # credible interval for the difference between Group 1 and Group 2.
    samp_mat2C <- as.matrix(samples2C)
    quantile(samp_mat2C[, "theta[1]"] - samp_mat2C[, "theta[2]"], c(0.025, 0.5, 0.975))
    
    #calculate difference in thetas btw grp1 and grp2 from mcmc samples
    diff2C <- (samp_mat2C[, "theta[1]"] - samp_mat2C[, "theta[2]"])
    #plot distribution of diff12
    hist(diff2C)
    #number of diff12 estimates from mcmc samples > 0
    min2C <- length(diff2C[diff2C > 0.00])
    mod2C <- length(diff2C[diff2C >= 0.05])
    str2C <- length(diff2C[diff2C >= 0.09])
    max2C <- length(diff2C[diff2C >= 0.13])
    #total number of diff12 estimates from mcmc samples
    dend2C <- length(diff2C)
    #calc proportion of diff12 > 0
    Pdmin2C <- min2C/dend2C
    Pdmod2C <- mod2C/dend2C
    Pdstr2C <- str2C/dend2C
    Pdmax2C <- max2C/dend2C
    #Probability group 1 > group 2 (i.e., more likely callbacks > less likely callbacks)
    Pdmin2C
    Pdmod2C
    Pdstr2C
    Pdmax2C
  
  #Panel 2.D, Cols.4-6 (Observed)
  set.seed(8675309)
  n_callbacks2D <- c(2, 1)
  n_employers2D <- c(22, 32)
  bayes.prop.test(n_callbacks2D, n_employers2D)
  fit2D <- bayes.prop.test(n_callbacks2D, n_employers2D)
  plot(fit2D)
  fit2D$stats
  
  model.code(fit2D)
  #Code for Estimating Probability at/above counterfactual critical values (Figure 2):
    ### Model code for the Bayesian First Aid  ###
    ### alternative to the test of proportions ###
    require(rjags)
    
    # Setting up the data
    x <- c(2, 1) 
    n <- c(22, 32) 
    
    # The model string written in the JAGS language
    model_string <- "model {
    for(i in 1:length(x)) {
      x[i] ~ dbinom(theta[i], n[i])
      theta[i] ~ dbeta(1, 1)
      x_pred[i] ~ dbinom(theta[i], n[i])
    }
  }"
    
    # Running the model
    model2D <- jags.model(textConnection(model_string), data = list(x = x, n = n), 
                          n.chains = 3, n.adapt=1000)
    samples2D <- coda.samples(model2D, c("theta", "x_pred"), n.iter=5000)
    
    # Inspecting the posterior
    plot(samples2D)
    summary(samples2D)
    
    # Extract the mcmc samples as a matrix and compare the thetas 
    # of the groups. The following shows the median and 95%
    # credible interval for the difference between Group 1 and Group 2.
    samp_mat2D <- as.matrix(samples2D)
    quantile(samp_mat2D[, "theta[1]"] - samp_mat2D[, "theta[2]"], c(0.025, 0.5, 0.975))
    
    #calculate difference in thetas btw grp1 and grp2 from mcmc samples
    diff2D <- (samp_mat2D[, "theta[1]"] - samp_mat2D[, "theta[2]"])
    #plot distribution of diff12
    hist(diff2D)
    #number of diff12 estimates from mcmc samples > 0
    min2D <- length(diff2D[diff2D > 0.00])
    mod2D <- length(diff2D[diff2D >= 0.07])
    str2D <- length(diff2D[diff2D >= 0.15])
    max2D <- length(diff2D[diff2D >= 0.23])
    #total number of diff12 estimates from mcmc samples
    dend2D <- length(diff2D)
    #calc proportion of diff12 > 0
    Pdmin2D <- min2D/dend2D
    Pdmod2D <- mod2D/dend2D
    Pdstr2D <- str2D/dend2D
    Pdmax2D <- max2D/dend2D
    #Probability group 1 > group 2 (i.e., more likely callbacks > less likely callbacks)
    Pdmin2D
    Pdmod2D
    Pdstr2D
    Pdmax2D
    
#Reproduce Figure 2: 
  fit2A.data = as.data.frame(fit2A)
  fit2A.data$tdif=fit2A.data$theta1-fit2A.data$theta2
  fit2A.data$cat4=1
  fit2A.data$cat4alpha="2A"
  fit2A.data$cat4range="0 (-.09, .09)"

  fit2B.data = as.data.frame(fit2B)
  fit2B.data$tdif=fit2B.data$theta1-fit2B.data$theta2
  fit2B.data$cat4=2
  fit2B.data$cat4alpha="2B"
  fit2B.data$cat4range=".04 (-.08, .20)"
  
  fit2C.data = as.data.frame(fit2C)
  fit2C.data$tdif=fit2C.data$theta1-fit2C.data$theta2
  fit2C.data$cat4=3 
  fit2C.data$cat4alpha="2C"
  fit2C.data$cat4range=".03 (-.07, .12)"
  
  fit2D.data = as.data.frame(fit2D)
  fit2D.data$tdif=fit2D.data$theta1-fit2D.data$theta2
  fit2D.data$cat4=4
  fit2D.data$cat4alpha="2D"
  fit2D.data$cat4range=".06 (-.09, .22)"
  
  fit.data=rbind(fit2A.data, fit2B.data, fit2C.data, fit2D.data)  

  #Density Plots
  library(ggplot2)

  #Create Figure 2:
  tdif.facet= 
  ggplot()+
  facet_wrap(~cat4alpha, ncol=1, strip.position="left") + 
  geom_density(data=fit.data, aes(x=tdif), fill="lightblue", color="lightblue") + 
    theme(axis.line.x=element_line(linetype="solid"),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),
          legend.key=element_blank(),
          panel.background=element_blank(),
          panel.spacing=unit(0, "points"),
          panel.grid=element_blank(),
          strip.text.y=element_text(angle=180),
          strip.background=element_blank(),
          strip.placement=("outside"),
          panel.border=element_rect(color="black", fill=NA, size=0.5)) +
    #Add line at median of distribution
    geom_segment(data=fit2A.data, aes(x=median(tdif), xend=median(tdif), y=0, yend=9.1), color="black", linetype="solid", size=1) +
    geom_segment(data=fit2B.data, aes(x=median(tdif), xend=median(tdif), y=0, yend=5.8), color="black", linetype="solid", size=1) +
    geom_segment(data=fit2C.data, aes(x=median(tdif), xend=median(tdif), y=0, yend=9), color="black", linetype="solid", size=1) +
    geom_segment(data=fit2D.data, aes(x=median(tdif), xend=median(tdif), y=0, yend=5.2), color="black", linetype="solid", size=1) +
    #Add lines representing 95% HDI:
    geom_segment(data=fit2A.data, aes(x=-.09, xend=.09, y=0, yend=0), color="black", linetype="solid", size=1) +
    geom_segment(data=fit2B.data, aes(x=-.08, xend=.20, y=0, yend=0), color="black", linetype="solid", size=1) +
    geom_segment(data=fit2C.data, aes(x=-.07, xend=.12, y=0, yend=0), color="black", linetype="solid", size=1) +
    geom_segment(data=fit2D.data, aes(x=-.09, xend=.22, y=0, yend=0), color="black", linetype="solid", size=1) +
    #Add lines at "moderate" and "strong" critical values:
    geom_segment(data=fit2A.data, mapping=aes(x=.05, xend=.05, y=0, yend=Inf), color="black", linetype="dashed", size=.5) +
    geom_segment(data=fit2A.data, mapping=aes(x=.08, xend=.08, y=0, yend=Inf), color="black", linetype="dashed", size=.5) +
    geom_segment(data=fit2B.data, mapping=aes(x=.07, xend=.07, y=0, yend=Inf), color="black", linetype="dashed", size=.5) +
    geom_segment(data=fit2B.data, mapping=aes(x=.12, xend=.12, y=0, yend=Inf), color="black", linetype="dashed", size=.5) +
    geom_segment(data=fit2C.data, mapping=aes(x=.05, xend=.05, y=0, yend=Inf), color="black", linetype="dashed", size=.5) +
    geom_segment(data=fit2C.data, mapping=aes(x=.09, xend=.09, y=0, yend=Inf), color="black", linetype="dashed", size=.5) +
    geom_segment(data=fit2D.data, mapping=aes(x=.07, xend=.07, y=0, yend=Inf), color="black", linetype="dashed", size=.5) +
    geom_segment(data=fit2D.data, mapping=aes(x=.15, xend=.15, y=0, yend=Inf), color="black", linetype="dashed", size=.5) +
    #Add labels for "Maximum," "Strong," and "Moderate":
    geom_text(data=fit2A.data, aes(x=.05, y=9, label="a. Moderate"), size=2.5, hjust=0, nudge_x=.01) +
    geom_text(data=fit2A.data, aes(x=.08, y=7.5, label="b. Strong"), size=2.5, hjust=0, nudge_x=.01) +
    geom_text(data=fit2B.data, aes(x=.07, y=9, label="a"), size=2.5, hjust=0, nudge_x=.01) +
    geom_text(data=fit2B.data, aes(x=.12, y=9, label="b"), size=2.5, hjust=0, nudge_x=.01) +
    geom_text(data=fit2C.data, aes(x=.05, y=9, label="a"), size=2.5, hjust=0, nudge_x=.01) +
    geom_text(data=fit2C.data, aes(x=.09, y=9, label="b"), size=2.5, hjust=0, nudge_x=.01) +
    geom_text(data=fit2D.data, aes(x=.07, y=9, label="a"), size=2.5, hjust=0, nudge_x=.01) +
    geom_text(data=fit2D.data, aes(x=.15, y=9, label="b"), size=2.5, hjust=0, nudge_x=.01) +
    #Add dotted line at zero across all graphs
    geom_vline(xintercept=c(0), linetype="dotted") +
    scale_x_continuous(breaks=seq(-.2,.5,.1))
  tdif.facet


  
  
#SUPPLEMENTAL BAYESIAN REANALYSIS   

  #Specify Jeffreys prior rather than uniform prior 
    #After estimating bayes.prop.test, use model.code (e.g., model.code(fit2A)) command to output code 
    #Then change "theta[i] ~ dbeta(1, 1)" to "theta[i] ~ dbeta(.5, .5)"

  #Results of robustness check (code for models below)    
    #NOTE: random draws may make results  vary slightly (code below does not recognize set.seed command)

    #2A, UNIFORM PRIOR, theta[i] ~ dbeta(1, 1), P(g1>g2) = .52
    #2A, JEFFREYS PRIOR, CHANGE TO: theta[i] ~ dbeta(.5, .5), P(g1>g2) = .55
    
    #2B, UNIFORM PRIOR, theta[i] ~ dbeta(1, 1), P(g1>g2) = .76
    #2B, JEFFREYS PRIOR, CHANGE TO: theta[i] ~ dbeta(.5, .5), P(g1>g2) = .68
    
    #2C, UNIFORM PRIOR, theta[i] ~ dbeta(1, 1), P(g1>g2) = .76
    #2C, JEFFREYS PRIOR, CHANGE TO: theta[i] ~ dbeta(.5, .5), P(g1>g2) = .83
    
    #2D, UNIFORM PRIOR, theta[i] ~ dbeta(1, 1), P(g1>g2) = .82
    #2D, JEFFREYS PRIOR, CHANGE TO: theta[i] ~ dbeta(.5, .5), P(g1>g2) = .83
    
    
#Panel 2A (Observed) - Supplemental (Jeffreys prior)

    ### Model code for the Bayesian First Aid  ###
    ### alternative to the test of proportions ###
    require(rjags)
    
    # Setting up the data
    x <- c(7, 4) 
    n <- c(96, 60) 
    
    # The model string written in the JAGS language
    model_string <- "model {
    for(i in 1:length(x)) {
    x[i] ~ dbinom(theta[i], n[i])
    theta[i] ~ dbeta(.5, .5)
    x_pred[i] ~ dbinom(theta[i], n[i])
    }
    }"

    # Running the model
    model <- jags.model(textConnection(model_string), data = list(x = x, n = n), 
                        n.chains = 3, n.adapt=1000)
    samples <- coda.samples(model, c("theta", "x_pred"), n.iter=5000)
    
    # Inspecting the posterior
    plot(samples)
    summary(samples)
    
    # You can extract the mcmc samples as a matrix and compare the thetas 
    # of the groups. For example, the following shows the median and 95%
    # credible interval for the difference between Group 1 and Group 2.
    samp_mat <- as.matrix(samples)
    quantile(samp_mat[, "theta[1]"] - samp_mat[, "theta[2]"], c(0.025, 0.5, 0.975))
    
    ######
    #calculate difference in thetas btw grp1 and grp2 from mcmc samples
    diff12a <- (samp_mat[, "theta[1]"] - samp_mat[, "theta[2]"])
    #plot distribution of diff12
    hist(diff12a)
    #number of diff12 estimates from mcmc samples > 0
    numa <- length(diff12a[diff12a > 0])
    #total number of diff12 estimates from mcmc samples
    dena <- length(diff12a)
    #calc proportion of diff12 > 0
    Pa <- numa/dena
    #Probability group 1 > group 2 (i.e., more likely callbacks > less likely callbacks)
    Pa
    
    
#Panel 2B (Observed) - Supplemental (Jeffreys prior)

  ### Model code for the Bayesian First Aid  ###
  ### alternative to the test of proportions ###
  require(rjags)
  
  # Setting up the data
  x <- c(2, 9) 
  n <- c(22, 134) 
  
  # The model string written in the JAGS language
  model_string <- "model {
  for(i in 1:length(x)) {
  x[i] ~ dbinom(theta[i], n[i])
  theta[i] ~ dbeta(.5, .5)
  x_pred[i] ~ dbinom(theta[i], n[i])
  }
  }"

  # Running the model
  model <- jags.model(textConnection(model_string), data = list(x = x, n = n), 
                      n.chains = 3, n.adapt=1000)
  samples <- coda.samples(model, c("theta", "x_pred"), n.iter=5000)
  
  # Inspecting the posterior
  plot(samples)
  summary(samples)
  
  # You can extract the mcmc samples as a matrix and compare the thetas 
  # of the groups. For example, the following shows the median and 95%
  # credible interval for the difference between Group 1 and Group 2.
  samp_mat <- as.matrix(samples)
  quantile(samp_mat[, "theta[1]"] - samp_mat[, "theta[2]"], c(0.025, 0.5, 0.975))    
  
  ######
  #calculate difference in thetas btw grp1 and grp2 from mcmc samples
  diff12b <- (samp_mat[, "theta[1]"] - samp_mat[, "theta[2]"])
  #plot distribution of diff12
  hist(diff12b)
  #number of diff12 estimates from mcmc samples > 0
  numb <- length(diff12b[diff12b > 0])
  #total number of diff12 estimates from mcmc samples
  denb <- length(diff12b)
  #calc proportion of diff12 > 0
  Pb <- numb/denb
  #Probability group 1 > group 2 (i.e., more likely callbacks > less likely callbacks)
  Pb
  

#Panel 2.C (Observed) - Supplemental (Jeffreys prior)

  ### Model code for the Bayesian First Aid  ###
  ### alternative to the test of proportions ###
  require(rjags)
  
  # Setting up the data
  x <- c(10, 1) 
  n <- c(124, 32) 
  
  # The model string written in the JAGS language
  model_string <- "model {
  for(i in 1:length(x)) {
  x[i] ~ dbinom(theta[i], n[i])
  theta[i] ~ dbeta(.5, .5)
  x_pred[i] ~ dbinom(theta[i], n[i])
  }
  }"

  # Running the model
  model <- jags.model(textConnection(model_string), data = list(x = x, n = n), 
                      n.chains = 3, n.adapt=1000)
  samples <- coda.samples(model, c("theta", "x_pred"), n.iter=5000)
  
  # Inspecting the posterior
  plot(samples)
  summary(samples)
  
  # You can extract the mcmc samples as a matrix and compare the thetas 
  # of the groups. For example, the following shows the median and 95%
  # credible interval for the difference between Group 1 and Group 2.
  samp_mat <- as.matrix(samples)
  quantile(samp_mat[, "theta[1]"] - samp_mat[, "theta[2]"], c(0.025, 0.5, 0.975))
  
  ######
  #calculate difference in thetas btw grp1 and grp2 from mcmc samples
  diff12c <- (samp_mat[, "theta[1]"] - samp_mat[, "theta[2]"])
  #plot distribution of diff12
  hist(diff12c)
  #number of diff12 estimates from mcmc samples > 0
  numc <- length(diff12c[diff12c > 0])
  #total number of diff12 estimates from mcmc samples
  denc <- length(diff12c)
  #calc proportion of diff12 > 0
  Pc <- numc/denc
  #Probability group 1 > group 2 (i.e., more likely callbacks > less likely callbacks)
  Pc
  
  
#Panel 2D (Observed) - Supplemental (Jeffreys prior)

  ### Model code for the Bayesian First Aid  ###
  ### alternative to the test of proportions ###
  require(rjags)
  
  # Setting up the data
  x <- c(2, 1) 
  n <- c(22, 32) 
  
  # The model string written in the JAGS language
  model_string <- "model {
  for(i in 1:length(x)) {
  x[i] ~ dbinom(theta[i], n[i])
  theta[i] ~ dbeta(.5, .5)
  x_pred[i] ~ dbinom(theta[i], n[i])
  }
  }"

  # Running the model
  model <- jags.model(textConnection(model_string), data = list(x = x, n = n), 
                      n.chains = 3, n.adapt=1000)
  samples <- coda.samples(model, c("theta", "x_pred"), n.iter=5000)
  
  # Inspecting the posterior
  plot(samples)
  summary(samples)
  
  # You can extract the mcmc samples as a matrix and compare the thetas 
  # of the groups. For example, the following shows the median and 95%
  # credible interval for the difference between Group 1 and Group 2.
  samp_mat <- as.matrix(samples)
  quantile(samp_mat[, "theta[1]"] - samp_mat[, "theta[2]"], c(0.025, 0.5, 0.975))
  
  
  ######
  #calculate difference in thetas btw grp1 and grp2 from mcmc samples
  diff12d <- (samp_mat[, "theta[1]"] - samp_mat[, "theta[2]"])
  #plot distribution of diff12
  hist(diff12d)
  #number of diff12 estimates from mcmc samples > 0
  numd <- length(diff12d[diff12d > 0])
  #total number of diff12 estimates from mcmc samples
  dend <- length(diff12d)
  #calc proportion of diff12 > 0
  Pd <- numd/dend
  #Probability group 1 > group 2 (i.e., more likely callbacks > less likely callbacks)
  Pd
  