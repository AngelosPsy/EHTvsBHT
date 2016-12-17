#' ---
#' title: "Bayesian Hypothesis Testing versus Equivalence Hypothesis Testing"
#' author: " "
#' ---
#' Please note that although this version has been partially tested, a most thorough
#' look needs to be taken.
#' 
#' # Bayesian Hypothesis Testing versus Equivalent Hypothesis testing
#' This document includes a few functions that could serve as the start for the
#' further development of a head-to-head test between Bayesian Hypothesis Testing versus Equivalent Hypothesis testing.
#' In this document we compare the results between the default Bayes factor t-test (http://link.springer.com/article/10.3758/PBR.16.2.225)
#' and equivalence hypothesis testing for two sample independent t-test.
#' 
#' Clear workplace
rm(list = ls(all = T))
#' Load R packages
suppressPackageStartupMessages(library(effsize))
suppressPackageStartupMessages(library(BayesFactor))
suppressPackageStartupMessages(library(dplyr))
#' Load data here. To see how the data were generated, please see below.
load("/Users/angelos-miltiadiskrypotos/Dropbox/eqBayes/tmp_11_27_2016.RData")
#' Function of equivalence testing. This is the modified function
#' of TOSTER::TOSTtwo (https://github.com/Lakens/TOSTER/blob/master/R/TOSTtwo.R).
#' I modified it a bit so that I avoid the printing of the messages and the
#' generation of plots.
TOSTtwo <-function(m1, m2, sd1, sd2, n1, n2, eqbound_d, alpha = 0.05, var.equal = TRUE){
  low_eqbound_d = -eqbound_d
  high_eqbound_d = eqbound_d
  if(missing(alpha)) {
    alpha<-0.05
  }
  if(missing(var.equal)) {
    var.equal<-FALSE
  }
  if(var.equal==TRUE) {
    sdpooled<-sqrt((((n1 - 1)*(sd1^2)) + (n2 - 1)*(sd2^2))/((n1+n2)-2)) #calculate sd pooled
    low_eqbound<-low_eqbound_d*sdpooled
    high_eqbound<-high_eqbound_d*sdpooled
    degree_f<-n1+n2-2
    t1<-(abs(m1-m2)-low_eqbound)/(sdpooled*sqrt(1/n1 + 1/n2))  #students t-test lower bound
    p1<-pt(t1, degree_f, lower=FALSE) 
    t2<-(abs(m1-m2)-high_eqbound)/(sdpooled*sqrt(1/n1 + 1/n2)) #students t-test upper bound
    p2<-pt(t2, degree_f, lower=TRUE) 
    t<-(m1-m2)/(sdpooled*sqrt(1/n1 + 1/n2))
    pttest<-2*pt(-abs(t), df=degree_f)
    LL90<-(m1-m2)-qt(1-alpha, n1+n2-2)*(sdpooled*sqrt(1/n1 + 1/n2))
    UL90<-(m1-m2)+qt(1-alpha, n1+n2-2)*(sdpooled*sqrt(1/n1 + 1/n2))
    LL95<-(m1-m2)-qt(1-(alpha/2), n1+n2-2)*(sdpooled*sqrt(1/n1 + 1/n2))
    UL95<-(m1-m2)+qt(1-(alpha/2), n1+n2-2)*(sdpooled*sqrt(1/n1 + 1/n2))
  } else {
    sdpooled<-sqrt((sd1^2 + sd2^2)/2) #calculate sd root mean squared for Welch's t-test
    low_eqbound<-low_eqbound_d*sdpooled
    high_eqbound<-high_eqbound_d*sdpooled
    degree_f<-(sd1^2/n1+sd2^2/n2)^2/(((sd1^2/n1)^2/(n1-1))+((sd2^2/n2)^2/(n2-1))) #degrees of freedom for Welch's t-test
    t1<-(abs(m1-m2)-low_eqbound)/sqrt(sd1^2/n1 + sd2^2/n2) #welch's t-test upper bound
    p1<-pt(t1, degree_f, lower=FALSE) #p-value for Welch's TOST t-test
    t2<-(abs(m1-m2)-high_eqbound)/sqrt(sd1^2/n1 + sd2^2/n2) #welch's t-test lower bound
    p2<-pt(t2, degree_f, lower=TRUE) #p-value for Welch's TOST t-test
    t<-(m1-m2)/sqrt(sd1^2/n1 + sd2^2/n2) #welch's t-test NHST
    pttest<-2*pt(-abs(t), df=degree_f) #p-value for Welch's t-test
    LL90<-(m1-m2)-qt(1-alpha, degree_f)*sqrt(sd1^2/n1 + sd2^2/n2) #Lower limit for CI Welch's t-test
    UL90<-(m1-m2)+qt(1-alpha, degree_f)*sqrt(sd1^2/n1 + sd2^2/n2) #Upper limit for CI Welch's t-test
    LL95<-(m1-m2)-qt(1-(alpha/2), degree_f)*sqrt(sd1^2/n1 + sd2^2/n2) #Lower limit for CI Welch's t-test
    UL95<-(m1-m2)+qt(1-(alpha/2), degree_f)*sqrt(sd1^2/n1 + sd2^2/n2) #Upper limit for CI Welch's t-test
  }
  ptost<-max(p1,p2) #Get highest p-value for summary TOST result
  ttost<-ifelse(abs(t1) < abs(t2), t1, t2) #Get lowest t-value for summary TOST result
  results<-data.frame(t1,p1,t2,p2,degree_f,LL90,UL90, t, pttest)
  return(results)
}
#' Wrapper function for Bayes factor and TOST in which only Mean, SD, and N are provided as input

TOSTone<-function(m,mu,sd,n, eqbound_d, alpha){
  low_eqbound_d = -eqbound_d
  high_eqbound_d = eqbound_d
  if(missing(alpha)) {
    alpha<-0.05
  }
  low_eqbound<-low_eqbound_d*sd
  high_eqbound<-high_eqbound_d*sd
  degree_f<-n-1
  t1<-(m-mu-low_eqbound)/(sd/sqrt(n))# t-test
  p1<-pt(t1, degree_f, lower=FALSE) 
  t2<-(m-mu-high_eqbound)/(sd/sqrt(n)) #t-test
  p2<-pt(t2, degree_f, lower=TRUE) 
  t<-(m-mu)/(sd/sqrt(n))
  pttest<-2*pt(-abs(t), df=degree_f)
  LL90<-(m-mu)-qt(1-alpha, degree_f)*(sd/sqrt(n))
  UL90<-(m-mu)+qt(1-alpha, degree_f)*(sd/sqrt(n))
  LL95<-(m-mu)-qt(1-(alpha/2), degree_f)*(sd/sqrt(n))
  UL95<-(m-mu)+qt(1-(alpha/2), degree_f)*(sd/sqrt(n))
  ptost<-max(p1,p2) #Get highest p-value for summary TOST result
  ttost<-ifelse(abs(t1) < abs(t2), t1, t2) #Get lowest t-value for summary TOST result
  results<-data.frame(t1,p1,t2,p2,degree_f,LL90,UL90, t, pttest)
  colnames(results) <- c("t-value 1","p-value 1","t-value 2","p-value 2","df", paste("Lower Limit ",100*(1-alpha*2),"% CI",sep=""),
                         paste("Upper Limit ",100*(1-alpha*2),"% CI",sep=""),
                         "t",  "p")
  return(results)
}

bayesEqTwo <- function(m1, m2, s1, s2, n1, n2,
                     eqVar = FALSE, eqbound_d = .5,
                     alpha = 0.05, 
                     rscale = sqrt(2)/2){
  c.d <- (m2 - m1)/(sqrt((s1^2 + s2^2) /2))
  tost <- TOSTtwo(m1 = m1, m2 = m2, sd1 = s1, sd2 = s2,
               n1 = n1, n2 = n2, eqbound_d = eqbound_d, alpha = alpha,
               var.equal = eqVar)
  b.f <- BayesFactor::ttest.tstat(as.numeric(tost[8]), n1 = n1, n2 = n2, rscale = rscale, simple = TRUE)
  res <- data.frame("T-test stat" = tost[8], "P-value" = tost[9],
                    "Cohen' d" = c.d,
                    "BF10" = b.f,
                    "BF01" = 1/b.f, rscale = rscale, 
                    "t1" = tost[1], "p1" = tost[2],
                    "t2" = tost[3], "p2" = tost[4],
                    "df" = tost[5], "LL90" = tost[6],
                    "UL90" = tost[7], 
                    eqbound_d = eqbound_d,
                    m1 = m1, m2 = m2, s1 = s1, s2 = s2, n1 = n1, n2 = n2)
  rownames(res) <- NULL
  return(res)
}


bayesEqOne <- function(m, mu, sd, n,eqbound_d = .5, 
                       alpha = 0.05, rscale = sqrt(2)/2){
  tost <- TOSTone(m = m, mu = mu, n = n, sd = sd, eqbound_d = eqbound_d, alpha = alpha)
  b.f <- BayesFactor::ttest.tstat(as.numeric(tost[8]), n1 = n, n2 = 0, rscale = rscale, simple = TRUE)
  res <- data.frame("T-test stat" = tost[8], "P-value" = tost[9],
                    "BF10" = b.f,
                    "BF01" = 1/b.f, rscale = rscale, 
                    "t1" = tost[1], "p1" = tost[2],
                    "t2" = tost[3], "p2" = tost[4],
                    "df" = tost[5], "LL90" = tost[6],
                    "UL90" = tost[7], 
                    eqbound_d = eqbound_d,
                    m = m, mu = mu, sd = sd, n = n)
  rownames(res) <- NULL
  return(res)
}
#' Function for plotting BFs and p values for equivalence testing
bfEqPlot <- function(bf01, p1, p2, bf01E1 = NULL, bf01E2 = NULL, ...){
  layout(t(1:2))
  plot(unlist(bf01), ylab = "BF01", xlab = "Sample Size", type = 'b', xaxt = "n", ...)
  axis(1, at = seq(1, 96, 5), seq(5, 100, 5))
  if(!is.null(bf01E1)){
    lines(unlist(bf01E1), col = "red", type = "b")
  }
  if (!is.null(bf01E2)){
    lines(unlist(bf01E2), col = "blue", type = "b")
  }
  plot(unlist(p1), ylab = "p-Values (EHT)", type = 'b', xlab = "Sample Size", xaxt = "n",
       ylim = c(0, 1))
  axis(1, at = seq(1, 96, 5), seq(5, 100, 5))
  lines(unlist(p2), col = "red", type = "b")
  abline(h = 0.05, lty = 2, lwd = 2, col = "blue")
  legend("topleft", legend = c("p Value Lower", "p Value Higher", "0.05") , lty = c(1, 1, 2), col = c("black", "red", "blue"), bty = "n")
}

#' We are going to create the tmp object that has our results for different
#' sample size, r scale parameters, and choice of EHT for bounds. It takes long
#' to create this object. I have commented out the lines for now and load
#' the data later on.
#tmp <- matrix(NA, nrow = 115200, ncol = 20)
#colnames(tmp) <- c("T-test stat", "P-value", "Cohen' d", "BF10", "BF01", 
#                   "rscale", "t1", "p1", "t2", "p2", "df", "LL90", "UL90",
#                   "eqbound_d", "m1", "m2", "sd1", "sd2", "n1", "n2")
#int = 1
#for (i in seq(0.1, 2, 0.1)){ # Mean difference
#  for (j in seq(5, 100, 1)){ # Sample size
#    for (k in seq(0.1, 2, 0.1)){ # bounds
#      for (rs in c(.707, 1, 1.4)){ #' rscale
#        tmp[int, ] = as.numeric(bayesEqTwo(m1=i,m2=0.0,s1=1,s2=1,n1=j,n2=j, rscale = rs, eqbound_d=k, eqVar = TRUE))
#      int = int + 1
#    }
#  }
#}
#}
#tmp <- data.frame(tmp)

#setwd("/Users/angelos-miltiadiskrypotos/Dropbox/eqBayes")
#save.image(paste0("tmp", Sys.Date(), ".RData"))

#' Below is a series of plots where we see the BF01
#' and the p values for the EHT, while we change the sample
#' size and the eqbound_d. The last parameter goes from
#' .1 to 2 in steps of .1. Regarding the Bayesian plots (on the left side),
#' we have the BF01 when we have an rscale of .707 (black line), 1 (red line),
#' or 1.4 (black). These plots are the same all the time but we just repeat them
#' to allow the easy comparison.
for (i in seq(.1, 2, .1)){
ptmp <- tmp %>% filter(eqbound_d == i & m1 == .1)
  bfEqPlot(ptmp$BF01[ptmp$rscale == .707], ptmp$p1[ptmp$rscale == .707], 
           ptmp$p2[ptmp$rscale == .707],
           ptmp$BF01[ptmp$rscale == 1],
           ptmp$BF01[ptmp$rscale == 1.4], ylim = c(0, 10))
  title(paste0("eqbound = ", i))
}

#' An easy conclusion is that someone has to be careful
#' when choosing the eqbound limits as this will affect
#' the conclusions -- logical :). However, we see that the BF is able
#' to gather quicker evidence for the null hypothesis quicker, and quicker as the
#' distribution over the H1 becomes wider (i.e., the rscale increases).