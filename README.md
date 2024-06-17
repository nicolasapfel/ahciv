# ahciv
 
This repository provides the function ahciv that implements the method proposed in "Agglomerative Hierarchical Clustering for Selecting Valid Instrumental Variables" by Nicolas Apfel and Xiaoran Liang. \
The functions is in ahciv_function.R and can be sourced from there.\
Synthetic data is created in the file example.R This file also illustrates how to run the command.\
Options: *Y*: Dependent variable, as vector.\
*D*: Independent variable, as vector.\
*Z*: Matrix of instruments. Columns should have names.\
All variables should be pretransformed, in particular if an intercept should be included, the variables should be demeaned and variables should be residualized with respect to observable controls.\
*tau*: significance level used in downward testing procedure. We propose 0.1/log(n) in the paper. 
