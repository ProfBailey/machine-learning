# art_Lin-Luo_2019__Multilevel-CART-binary-outcomes.pdf

# algorithm

M.CART <- function (formula, data, random = "SID", 
                   rslope = "+X1+X2", lgmodel = "int", 
                   initialRandomEffects = rep(0, TotalObs),
                   ErrorTolerance = 0.00001, MaxIterations = 10000,
                   verbose = FALSE, cpmin = 0.01) {
  TotalObs <- dim(data)[1]
  originaldata <- data
  Predictors <- paste(attr(terms(formula), "term.labels"), 
                      collapse = "+")
  TargetName <- formula[[2]]
  if (length(TargetName) > 1) TargetName <- TargetName[3]
  if (verbose) print(paste("Target variable: ", TargetName))
  Target <- data[, toString(TargetName)]
  newdata <- data
  
  #COULD BE ERROR BELOW HERE (with "-" (minus or page cutoff?)
  originaldata$random <- rep(0, TotalObs)
  for (i in 1:length(summary(originaldata$TID))) {
    originaldata$random[originaldata$TID == i] <- 
      mean(originaldata[originaldata$TID == i, toString(TargetName)]) - 
      mean(originaldata[, toString(TargetName)])
  }
  
  AdjustedTarget <- data[, toString(TargetName)]- originaldata$random
  ContinueCondition <- TRUE
  iterations <- 0
  oldlik <- 0
  
  while (ContinueCondition) {
    iterations <- iterations + 1
    newdata[, "AdjustedTarget"] <- AdjustedTarget
    tree <- rpart(formula(paste(c("AdjustedTarget", Predictors), 
                                collapse = "~")), 
                  data = newdata, method = "anova", 
                  control = rpart.control(cp = cpmin, xval = 10))
    if (verbose)
      print(tree)
    
    newdata[, "nodeInd"] <- tree$where
    if (min(tree$where) == max(tree$where)) {
      glmerfit2 <-glmer(
        formula(paste(c(toString(TargetName), 
                        paste("1+(1|", random,")", sep = "")), 
                      collapse = "~")), data = newdata, 
        family = binomial, nAGQ = 0, na.action = na.exclude
      )
    }
    else {
      if (lgmodel == "int") {
        glmerfit2 <- glmer(
          formula(paste(c(toString(TargetName), 
                          paste("as.factor(nodeInd)+(1|", random,")",
                                sep = "")), 
                        collapse = "~")), data = newdata,
          family = binomial, nAGQ = 0, na.action = na.exclude)
      }
      else {
        glmerfit2 <- glmer(
          formula(paste(c(toString(TargetName), 
                          paste("as.factor(nodeInd)", 
                                paste("(1", rslope, "j", random, ")",
                                      sep = ""), sep = "+")),
                        collapse = "~")), data = newdata, 
          family = binomial, nAGQ = 0, na.action = na.exclude)
      }
    }
    
    newlik <- logLik(glmerfit2)
    ContinueCondition <- (abs(newlik - oldlik) >
                            ErrorTolerance & iterations < MaxIterations)
    oldlik <- newlik 
    AdjustedTarget2 <- (as.matrix(getME(glmerfit2, name = "X"))) %*%
      (as.matrix(getME(glmerfit2, name = "beta")))
    AdjustedTarget <- exp(AdjustedTarget2) / (1 + exp(AdjustedTarget2))
    AdjustedTarget[is.na(AdjustedTarget[,1]),1] <- 
      (data[,toString(TargetName)]-
         originaldata$random)[is.na(AdjustedTarget[,1])]
  }
  
  if (lgmodel != "int") {
    Between <- cbind(VarCorr(glmerfit2)[[1]][1:2], 
                     VarCorr(glmerfit2)[[1]][3:4])
  }
  else {
    Between <- VarCorr(glmerfit2)[[1]][1]
  }
  preditFinal<-predict(glmerfit2, type = "response")
  result <- list(data = data, IterationsUsed = iterations, 
                 Random = random, Rslope = rslope, 
                 ErrorTolerance = ErrorTolerance, prune_cv = prune_cv,
                 predMtree = preditFinal, Tree = tree, 
                 Treewhere = tree$where, EffectModel = glmerfit2,
                 RandomEffects = ranef(glmerfit2),
                 BetweenMatrix = as.matrix(Between)*sigma(glmerfit2)^2,
                 logLik = newlik, AIC = AIC(glmerfit2), 
                 BIC = BIC(glmerfit2), deviance = deviance(glmerfit2),
                 df = as.numeric(summary(glmerfit2)$AICtab[5]),
                 Formula = formula, Totalinter = interations)
  class(result) <- "M.CART"
  return(result)
}

# empirical demonstration
# (there may be errors in spacing btwn "+", unsure if it matters?)
M.CART(Y~ReadScore+MathScore+ScienceScore+EnjoySch+
         CloseT + CloseS + GradeImp + Home + FlGood +
       ParentEdu + ParentTlk + ParentHlp + Public + PctHsp +
       PctBlk + SchFund, 
       data = train, random = "schoolID", lgmodel = "slope", 
       rslope = "+ ReadScore + MathScore + ScienceScore + CloseT +
       ParentEdu * SchFund", cpmin  = 0.001)
