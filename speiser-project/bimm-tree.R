# dis_Speiser_2017__Decision-tree-and-random-forest-longitudinal-binary.pdf
# BiMM-tree functions

#BiMM tree with one iteration

BiMMtree1 <- function(traindata,testdata,formula,random){
  
  #initialize parameters
  minsize = round(length(traindata[, 1]) / 10, 0)
  data = traindata
  initialRandomEffects = rep(0, length(data[, 1]))
  ErrorTolerance = 0.001
  MaxIterations = 1000
  tree.control = rpart.control(minbucket=minsize)
  
  #parse formula
  Predictors <- paste(attr(terms(formula), "term.labels"), 
                      collapse = "+")
  TargetName <- formula[[2]]
  Target <- data[, toString(TargetName)]
  
  #set up variables for loop
  ContinueCondition <- TRUE
  iterations <- 0
  
  #initial values
  AdjustedTarget <- as.numeric(Target) - initialRandomEffects
  oldlik <- -Inf
  
  # Make a new data frame to include all the new variables
  newdata <- data
  
  #run 1 iteration of algorithm
  newdata[, "AdjustedTarget"] <- AdjustedTarget
  iterations <- iterations + 1
  
  #build tree
  tree <- rpart(formula(paste(c("AdjustedTarget", Predictors),
                              collapse = "~")), 
                data = data, method = "class", control = tree.control)
  
  ## Estimate New Random Effects and Errors using BLMER
  # Get variables that identify the node for each observation
  data[, "nodeInd"] <- 0
  data["nodeInd"] <- tree$where
  
  # Fit linear model with nodes as predictors (we use the original
  # target so likelihoods are comparable)
  # Check that the fitted tree has at least two nodes.
  if(min(tree$where) == max(tree$where)){
    lmefit <- tryCatch(bglmer(
      formula(c(paste(paste(c(toString(TargetName),1), collapse = "~"),
                                    "+(1|random)",sep = ""))), 
      data = data, family = binomial, 
      control = glmerControl(optCtrl = list(maxfun = 20000))), 
      error = function(cond)"skip"
      )
  }
  else{
    lmefit <- tryCatch(bglmer(
    formula(c(paste(paste(c(toString(TargetName),"as.factor(nodeInd)"), 
                          collapse = "~"), "+(1|random)", sep = ""))),
    data = data, family = binomial, 
    control = glmerControl(optimizer = "bobyqa",
                           optCtrl = list(maxfun = 20000000000))),
    error = function(cond)"skip")
  }
  
  #if GLMM did not converge, return NA's for accuracy statistics
  if(class(lmefit)[1] == "character"){
    #return train and test confusion matrices
    return(list(
      c(NA, NA, NA, NA),
      c(NA, NA, NA, NA),
      NA
    ))
  }
  else if(!(class(lmefit)[1] == "character")){
    
    #train dataset predictions
    train.preds.ave <- AdjustedTarget
    train.preds <- predict(tree, traindata, type = "class")
    
    #test dataset predictions
    test.preds <- predict(tree, testdata, type = "class")
    
    #format table to make sure it always has 4 entries, even if it is
    #only 2 by 1 (0's in other spots)
    t1 <- table(data$comagradelow, train.preds.ave)
    t4 <- table(testdata$comagradelow, test.preds)
    
    #code if table for train or test data if all predictions 
    #are for same group
    if(ncol(t1) == 1 & train.preds.ave[1] == 1){
      t1 <- c(0, 0, t1[1, 1], t1[2, 1])
    }
    else if(ncol(t1) == 1 & train.preds.ave[1] == 0){
      t1 <- c(t1[1, 1], t1[2, 1], 0, 0)
    }
    if(ncol(t4) == 1 & test.preds[1] == 1){
      t4 <- c(0, 0, t4[1, 1], t4[2, 1])
    }
    else if(ncol(t4) == 1 & test.preds[1] == 0){
      t4 <- c(t4[1, 1], t4[2, 1], 0, 0)
    }
    
    #return train and test confusion matrices, # iterations
    return(list(
      c(t1),
      c(t4),
      iterations
    ))
  }
}

#BiMM tree with H1 updates

BiMMtreeH1<-function(traindata, testdata, formula, random, seed){
  
  #set up variables for Bimm method
  data = traindata1
  initialRandomEffects = rep(0, length(data[, 1]))
  ErrorTolerance = 0.006
  MaxIterations = 1000
  
  #parse formula
  Predictors <- paste(attr(terms(formula), "term.labels"), collapse="+")
  TargetName <- formula[[2]]
  Target <- data[, toString(TargetName)]
  
  #set up variables for loop
  ContinueCondition <- TRUE
  iterations <- 0
  
  #initial values
  AdjustedTarget <- as.numeric(Target) - initialRandomEffects
  oldlik <- -Inf
  
  # Make a new data frame to include all the new variables
  newdata <- data
  
  while(ContinueCondition){
    
    # Current values of variables
    newdata[, "AdjustedTarget"] <- AdjustedTarget
    iterations <- iterations + 1
    
    #build tree
    tree <- rpart(formula(paste(c("AdjustedTarget", Predictors),
                                collapse = "~")),
                  data = data, method = "class", control = tree.control)
    
    ## Estimate New Random Effects and Errors using BLMER
    # Get variables that identify the node for each observation
    data[,"nodeInd"] <- 0
    data["nodeInd"] <- tree$where
    
    # Fit linear model with nodes as predictors (we use the original
    #target so likelihoods are comparable)
    # Check that the fitted tree has at least two nodes.
    if(min(tree$where) == max(tree$where)){
      lmefit <- tryCatch(bglmer(
        formula(c(paste(paste(c(toString(TargetName), 1), collapse = "~"),
                        "+(1|random)", sep = ""))),
        data = data, family = binomial, 
        control = glmerControl(optCtrl=list(maxfun=20000))),
        error=function(cond)"skip")
    }
    else {
      lmefit <- tryCatch(bglmer(
      formula(c(paste(paste(c(toString(TargetName),
      "as.factor(nodeInd)"), collapse = "~"), "+(1|random)",sep = ""))),
      data = data, family = binomial, 
      control = glmerControl(optimizer = "bobyqa",
                             optCtrl = list(maxfun = 20000000000))),
      error = function(cond)"skip")
    }
    # Get the likelihood to check on convergence
    if(!(class(lmefit)[1] == "character")){
      newlik <- logLik(lmefit)
      ContinueCondition <- (((newlik - oldlik) > ErrorTolerance) & 
                              (iterations < MaxIterations))
      oldlik <- newlik
      
      # Extract random effects to make the new adjusted target
      logit <- predict(tree, type = "prob")[, 2]
      logit2 <- exp(predict(lmefit, re.form = NA)) / 
        (1 + exp(predict(lmefit, re.form = NA))) #population level effects
      AllEffects <- (logit + logit2) / 2 #average them
      #split function h1
      AdjustedTarget <- ifelse(as.numeric(AdjustedTarget) +
                                 AllEffects > 0.5, 1, 0)
    }
    else{
      ContinueCondition <- FALSE
    }
  }
  
  if(class(lmefit)[1] == "character"){
    #return train and test confusion matrices
    return(list(
      c(NA,NA,NA,NA),
      c(NA,NA,NA,NA),
      NA
    ))
  }
  else if(!(class(lmefit)[1] == "character")){
    #average effects
    train.preds.ave <- AdjustedTarget
    
    #test dataset predictions-same for all 3 updating methods for the
    #1 iteration model
    test.preds <- predict(tree, testdata, type = "class")
    
    #format table to make sure it always has 4 entries, even if it is
    #only 2 by 1 (0's in other spots)
    t1 <- table(data$ys, train.preds.ave)
    t4 <- table(testdata$ys, test.preds)
    if(ncol(t1) == 1 & train.preds.ave[1] == 1){
      t1 <- c(0, 0, t1[1, 1], t1[2, 1])
    }
    else if(ncol(t1) == 1 & train.preds.ave[1] == 0){
      t1 <- c(t1[1, 1], t1[2, 1], 0, 0)
    }
    if(ncol(t4) == 1 & test.preds[1] == 1){
      t4 <- c(0, 0, t4[1, 1], t4[2, 1])
    }
    else if(ncol(t4) == 1 & test.preds[1] == 0){
      t4 <- c(t4[1, 1], t4[2, 1], 0, 0)
    }
    
    #return train and test confusion matrices
    return(list(
      c(t1),
      c(t4),
      iterations
    ))
  }
}

#BiMM tree with H3 updates

BiMMtreeH3<-function(traindata, testdata, formula, random, seed){
  
  #set up variables for Bimm method
  data = traindata1
  initialRandomEffects = rep(0, length(data[, 1]))
  ErrorTolerance = 0.006
  MaxIterations = 1000
  
  #parse formula
  Predictors <- paste(attr(terms(formula), "term.labels"), collapse = "+")
  TargetName <- formula[[2]]
  Target <- data[ ,toString(TargetName)]
  
  #set up variables for loop
  ContinueCondition <- TRUE
  iterations <- 0
   
  #initial values
  AdjustedTarget <- as.numeric(Target) - initialRandomEffects
  oldlik <- -Inf
  
  # Make a new data frame to include all the new variables
  newdata <- data
  
  while(ContinueCondition){
    # Current values of variables
    newdata[,"AdjustedTarget"] <- AdjustedTarget
    iterations <- iterations + 1
    
    #build tree
    tree <- rpart(
      formula(paste(c("AdjustedTarget", Predictors),
                    collapse = "~")),
      data = data, method = "class", control = tree.control)
    
    ## Estimate New Random Effects and Errors using BLMER
    # Get variables that identify the node for each observation
    data[, "nodeInd"] <- 0
    data["nodeInd"] <- tree$where
    
    # Fit linear model with nodes as predictors (we use the original
    #target so likelihoods are comparable)
    #Check that the fitted tree has at least two nodes.
    if(min(tree$where) == max(tree$where)){
      lmefit <- tryCatch(bglmer(
        formula(c(paste(paste(c(toString(TargetName), 1), collapse = "~"),
                        "+(1|random)", sep = ""))),
        data = data, family = binomial, 
        control = glmerControl(optCtrl = list(maxfun = 20000))),
        error = function(cond)"skip")
    } 
    else {
      lmefit <-tryCatch(bglmer(
      formula(c(paste(paste(c(toString(TargetName), "as.factor(nodeInd)"), 
                            collapse = "~"), "+(1|random)", sep = ""))),
      data = data, family = binomial, 
      control = glmerControl(optimizer = "bobyqa",
                             optCtrl = list(maxfun = 20000000000))),
      error = function(cond)"skip")
    }
    
    # Get the likelihood to check on convergence
    if(!(class(lmefit)[1] == "character")){
      newlik <- logLik(lmefit)
      ContinueCondition <- (newlik - oldlik > ErrorTolerance & 
                            iterations < MaxIterations)
      oldlik <- newlik
    # Extract random effects to make the new adjusted target
      logit <- predict(tree, type = "prob")[, 2]
      logit2 <- exp(predict(lmefit, re.form = NA)) / 
        (1 + exp(predict(lmefit, re.form = NA))) #population level effects
      AllEffects <- (logit + logit2) / 2 #average them
      AdjustedTarget <- ifelse(as.numeric(AdjustedTarget) + 
                                 AllEffects > 0.5, 1, 0)
      #new split function h3
      for(k in 1:length(AllEffects)){
        if(as.numeric(AdjustedTarget[k]) + AllEffects[k] < 0.5){
          AdjustedTarget[k] = 0
        }
        else if(as.numeric(AdjustedTarget[k]) + AllEffects[k] > 1.5){
          AdjustedTarget[k] = 1}
        else{
          #generate random probability coin flip basedon AllEffects 
          #(q notation in paper)
          AdjustedTarget[k]<-rbinom(1,1,AllEffects[k])
        }
      }
    }
    else{ 
      ContinueCondition <- FALSE
    }
  }
  if(class(lmefit)[1] == "character"){
    #return train and test confusion matrices
    return(list(
      c(NA,NA,NA,NA),
      c(NA,NA,NA,NA),
      NA
    ))
  }
  else if(!(class(lmefit)[1] == "character")){
    #average effects
    train.preds.ave <- AdjustedTarget
    
    #test dataset predictions-same for all 3 updating methods for the
    #1 iteration model
    test.preds <- predict(tree, testdata, type = "class")
    
    #format table to make sure it always has 4 entries, even if it is
    #only 2 by 1 (0's in other spots)
    t1 <- table(data$ys, train.preds.ave)
    t4 <- table(testdata$ys, test.preds)
    if(ncol(t1) == 1 & train.preds.ave[1] == 1){
      t1 <- c(0, 0, t1[1, 1], t1[2, 1])
    }
    else if(ncol(t1) == 1 & train.preds.ave[1] == 0){
      t1 <- c(t1[1, 1], t1[2, 1], 0, 0)
    }
    if(ncol(t4) == 1 & test.preds[1] == 1){
      t4 <- c(0, 0, t4[1, 1], t4[2, 1])
    }
    else if(ncol(t4) == 1 & test.preds[1] == 0){
      t4 <- c(t4[1, 1], t4[2, 1], 0, 0)
    }
    #return train and test confusion matrices, # iterations
    return(list(
      c(t1),
      c(t4),
      iterations
    ))
  }
}
