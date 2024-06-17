#########################################
# R-Code and illustration for 
# "Agglomerative Hierarchical Clustering for Selecting Valid Instrumental Variables"
# by Nicolas Apfel and Xiaoran Liang
# Version: 1.0
# Date: 17/06/2024
# Author: Nicolas Apfel
#########################################

rm(list=ls())
# Packages ----
# Install packages if necessary
library(haven)
library(AER)
library(gtools)
library(glmnet)
library(dplyr)
library(MASS)

# Variables ----
P = 2 # nr of regressors

n <- 5000 # nr of observations

L <- 21; # Nr of IVs
s0 <- 12; # Nr of invalid IVs
pi = c(rep(1,6), rep(0.5,6), rep(0,(L-s0))) # direct effect coefficients
beta <- rep(0, P); # true effects

set.seed(15046)
gamma1 <- runif(L,1,2);
gamma2 <- runif(L,3,4); # First stage coefficient matrix
gamma = cbind(gamma1, gamma2)

sv = which(pi == 0) # valid IV ids
si = which(pi != 0) # invalid IV ids

cgamma <- 1
epsvar <- 1; # error variance
rho1 = 0.25; # rhos govern correlations between first and second stage errors
rho2 = 0.25; 
rho = c(rho1, rho2)

Sigmaz <- matrix(rep(0,L*L), nrow=L) # Covariance matrix
for(i in 1:L){
  for(j in 1:L){
    Sigmaz[i,j] <- 0.5^abs(i-j) 
  }
}

eps <- rnorm(n, 0, epsvar); # Error
Z <- matrix(mvrnorm(n,mu=rep(0,L),Sigma=Sigmaz), ncol=L); # IV matrix
colnames(Z) <- paste("Z", 1:dim(Z)[2],sep=""); # Set colnames of IV matrix
error <- matrix(rnorm(P*n,0,1),nrow=n, ncol=P) + rho*matrix(rep(eps, times = P), nrow = n);
D <- 0.1 + Z %*% gamma * cgamma + error; # First stage
Y <- 0.1 + Z %*% pi + D %*% beta + eps; # Structural equation

# Demean
Y <- Y-mean(Y)
D <- D-mean(D)
Z <- sweep(x=Z, MAR=2, STATS=colMeans(Z), FUN="-")

source("ahciv_function.R")

# Analysis ----
# Multiple
ahc_res <- ahc_iv(Y=Y, D=D, Z=Z, tau=0.1/log(length(Y))) 

# The following results are stored:
ahc_res$wv
ahc_res$wi
summary(ahc_res$ivreg)