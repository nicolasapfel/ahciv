# ahciv

## Description: 
This repository provides the function **ahciv()** that implements the method proposed in "Agglomerative Hierarchical Clustering for Selecting Valid Instrumental Variables" by Nicolas Apfel and Xiaoran Liang. \
The functions is in **ahciv_function.R** and can be sourced from there.\
Synthetic data is created in the file **example.R** This file also illustrates how to run the command.

## Options: 
*Y*: Dependent variable, as vector.\
*D*: Independent variable, as vector.\
*Z*: Matrix of instruments. Columns should have names.\
*tau*: significance level used in downward testing procedure. We propose 0.1/log(n) in the paper. \
All variables should be pretransformed as necessary, in particular if an intercept should be included, the variables should be demeaned and if controls should be included variables should be residualized with respect to observable controls.

## Saved outputs: 
*wv*: Column names of IVs selected as valid.\
*wi*: Column names of IVs selected as invalid.\
*ivreg*: ivreg object of the regression using the IVs selected as invalid as controls and those selected as valid as IVs. Get summary via summary().
