# hc function ---- 
ahc_iv <- function(Y, D, Z, tau){
  # Define:
  n <- dim(D)[1] # Number of obs.
  P <- dim(D)[2] # Number of endogenous regressors
  L <- dim(Z)[2] # Number of IVs
  
  names_Z <- colnames(Z)
  
  # Create all just identifying combinations
  comb <- combinations(n=length(colnames(Z)), v=colnames(Z), r=P, repeats.allowed=F)
  # Number of just-id. combinations:
  Len <- choose(length(colnames(Z)), P)
  
  bvec <- matrix(NA, nrow=Len, ncol=P) # Empty vector, for just-id. estimates
  
  for(j in 1:Len){ # Generate all just-id. estimates
    u <- comb[j,]
    Zt.end <- Z[, colnames(Z)[-which(colnames(Z) %in% u)]]
    bvec[j,] <- ivreg(Y ~ D + Zt.end -1 | Z -1 )$coefficients[1:P]
  }
  plot_bvec <- plot(sort(bvec)) # Sort and plot for visual inspection.
  
  # Distance matrix
  dist <- dist(bvec, method = 'euclidean')
  # Dendrogram
  dend <- hclust(dist, method="complete")
  #Matrix: one col. for step number (,1), one with largest cluster identity (,2), one for largest family size (,3), on for IVs involved (,4), one for p-value (,5)
  steps <- matrix(NA, nrow=Len-1, ncol=5)
  steps[,1] <- 1:(Len-1)
  # < Fill matrix with values
  
  for (b in 1:(Len-1)){
    # Pick a step of the dendrogram, the one with k=b clusters
    dstep <- cutree(dend, k=b)
    # Select the largest cluster(s)
    (dstep2 <- which(table(dstep) == max(table(dstep))))
    # Temporary table which includes cluster number (1), clustersize (2), number of IVs involved (3), p-value (4)
    tempstep <- matrix(NA, nrow=length(dstep2), ncol=4)
    tempstep[,1] <- dstep2
    tempstep[,2] <- max(table(dstep))
    tempstep
    
    tryCatch({ # In case one step fails, skip. Should display error.
      for(k in 1:length(dstep2)){ # For each maximal cluster
        dtemp <- dstep2[k] # Sequentially select cluster ids with maximal cluster number
        dtemp <- which(dstep == dtemp) # Which estimations are part of the maximal cluster
        (dtemp <- unique(as.vector(comb[dtemp,]))) # All IVs involved in these estimations
        tempstep[,3] <- length(dtemp) # Save number of IVs
        wv <- unique(dtemp) # Which IVs are involved in the estimation of the estimates selected as consistent?
        wi <- which(!(names_Z %in% wv)) # Which are invalid
        Ze <- Z[,wi] # Matrix of invalid IVs
        X <- cbind(Ze, D) # Invalid IVs and treatment variable(s)
        res_FirstStep <- residuals(AER::ivreg(Y ~ X - 1 | Z - 1)) # First step 2SLS regression
        
        Weight_SecondStep <- crossprod(res_FirstStep * Z); # Weight matrix
        
        Coef_SecondStep <- solve(
          t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
        ) %*% t(X) %*% Z %*% solve(Weight_SecondStep) %*% t(Z) %*% Y; # Coefficients from second step
        
        res_SecondStep <- as.vector(Y - X %*% Coef_SecondStep); 
        
        sd_SecondStep <- sqrt(diag(solve(
          t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
        ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)%*%crossprod(res_SecondStep * Z)%*%t(
          solve(
            t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
          ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)
        ))); # SE from second step
        
        HansenJ_Stat <- t(res_SecondStep) %*% Z %*% solve(Weight_SecondStep) %*%
          t(Z) %*% res_SecondStep; # Hansen Statistic
        
        # Compute p-value and save it
        pv <- pchisq(HansenJ_Stat, df=ncol(Z) - ncol(X), lower.tail = FALSE)
        tempstep[k,4] <- pv
      }
      # For the case of a tie, select the step with the most IVs involved.
      tempstep <- matrix(tempstep[which(tempstep[,3]==max(tempstep[,3])),], ncol=4)
      dstep2 <- tempstep[which(tempstep[,4]==max(tempstep[,4])), ] 
      # < If tie subsists, select the step with the highest p-value
      if(class(dstep2)[1]=="matrix"){
        print("Warning: Multiple results!")
        dstep2 <- matrix(dstep2, ncol=4)[1,]
      } # If tie subsists, print warning and select the first row.
      steps[b,2:5] <- dstep2 
    }, error=function(e){})
  }
  # In steps, we now have: step number (1), cluster number (2), clustersize (3), number of IVs involved (4), p-value (5)
  #Select the largest family for which HS-test not rejected
  maxover <- steps[which(steps[,5]> tau), , drop =F] # Select steps which do not violate HS criterion
  #Change to matrix if this is a vector already
  
  # If there are no models with a p-value larger the cutoff
  # return error and break
  if(dim(maxover)[1]<1){
    print("There was no single model which fulfilled the HS criterion. Use naive model.")
    summary(ivreg(Y ~ D - 1| Z - 1))
  } else { # If there are one or more models with p-value larger the cutoff:
    maxover <- maxover[which(maxover[,3]==max(maxover[,3])), , drop=F] # Select steps with maximal number of clusters
    maxover <- maxover[which(maxover[,4]==max(maxover[,4])), , drop=F] # Select steps with maximal number of IVs
    maxover <- maxover[which(maxover[,5]==max(maxover[,5])), ,drop=F] # Select step(s) with maximal p-value.
    
    selfam <- maxover[1,] # Select first row, for case that multiple results generated (this is effectively excluded by earlier calculations)
    finalK <- selfam[1] # Final step
    finalId <- selfam[2] # Final cluster ID
    findstep <- cutree(dend, k=finalK)
    findstep <- which(findstep==finalId)
    #Select IVs
    wv_names <- unique(as.vector(comb[findstep,]))
    wi <- which(!(names_Z %in% wv_names))
    wi_names <- colnames(Z)[wi]
    
    Ze <- matrix(Z[,wi], nrow=n)
    # Valid and invalid instruments, corrected shift-share IV
    if(sum(wi)==0){
      bc <- ivreg(Y ~ D -1 | Z -1)
    } else {
      bc <- ivreg(Y ~ D + Ze -1 | Z -1)
    }
    
    #Save nInv
    res_list <- list(wi=wi_names, wv=wv_names, ivreg=bc)
    
    print("IVs selected as invalid:")
    print(wi_names)
    print("IVs selected as valid:")
    print(wv_names)
    return(res_list)
  }
}