

library(glmnet)
library(survival)
library(survivalsvm)
library(tidyverse)
library(survminer)
library(CoxBoost)
library(superpc)
library(survcomp)
library(randomForestSRC)
library(plsRcox)
library(gbm)


cindex <- function(data,RS){
  rs=as.data.frame(cbind(data[,1:2],RS))
  aa=coxph(Surv(OS.time,OS)~RS,rs)
  cindex=as.numeric(summary(aa)$concordance[1])
  return(cindex)
}

####LASSO
set.seed(seed = 123)
my_lasso <- function(training, val, name){
  y <- data.matrix(Surv(time = training$OS.time, event = training$OS))
  x <- training[, -c(1,2)]
  
  cv.fit = cv.glmnet(x = as.matrix(x), y = y, family = "cox", alpha = 1, nfolds = 10)
  fit = glmnet(x = x,y = y, family = "cox", alpha = 1, lambda = cv.fit$lambda.min)
  coef.min = coef(cv.fit, s = "lambda.min")  
  active.min = which(coef.min@i != 0) ## 找出那些回归系数没有被惩罚为0的
  lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] ## 提取基因名称
  fit$subFeature <- lasso_geneids
  RS <- as.numeric(predict(fit, as.matrix(training[, -c(1,2)])))
  ci <- cindex(training,RS)
  fit$cindex <- ci
  
  cs <- c(ci)
  for (i in 1:length(val)) {
    data <- val[[i]]
    y <- data.matrix(Surv(time = data$OS.time, event = data$OS))
    x <- data[, -c(1,2)]
    RS <- as.numeric(predict(fit, as.matrix(data[, -c(1,2)])))
    ci <- cindex(data, RS)
    cs <- c(cs, ci)
  }
  fit$name <- name
  fit$cs <- cs
  return(list(fit))
}


  


#stepcox
my_stepcox <- function(training, val, name){
  y <- data.matrix(Surv(time = training$OS.time, event = training$OS))
  x <- training[, -c(1,2)]
  direction = c("both", "backward", "forward")
  models <- list()
  for (i in direction) {
    fit <- step(coxph(Surv(OS.time,OS)~., training), direction = i, trace = 0)
    fit$subFeature = names(coef(fit))
    
    RS <- as.numeric(predict(fit, training))
    ci <- cindex(training,RS)
    fit$cindex <- ci
    
    
    name2 <- paste0(name, "[", i, "]")
    cs <- c(ci)
    for (j in 1:length(val)) {
      data <- val[[j]]
      y <- data.matrix(Surv(time = data$OS.time, event = data$OS))
      x <- data[, -c(1,2)]
      RS <- as.numeric(predict(fit, data))
      ci <- cindex(data, RS)
      cs <- c(cs, ci)
    }
    fit$name <- name2
    fit$cs <- cs
    
    models[[i]] <- fit
  }
  return(models)
}



#survival-SVM
my_survivalsvm <- function(training, val, name){
  y <- data.matrix(Surv(time = training$OS.time, event = training$OS))
  x <- training[, -c(1,2)]
  fit <- survivalsvm(Surv(OS.time,OS)~., training,type = "regression", gamma.mu = 1,
                     opt.meth = "quadprog", kernel = "rbf_kernel")
  fit$subFeature = colnames(x)
  RS <- as.numeric(predict(fit, training)$predicted)
  ci <- cindex(training,RS)
  fit$cindex <- ci
  
  
  cs <- c(ci)
  for (i in 1:length(val)) {
    data <- val[[i]]
    y <- data.matrix(Surv(time = data$OS.time, event = data$OS))
    x <- data[, -c(1,2)]
    RS <- as.numeric(predict(fit, data)$predicted)
    ci <- cindex(data, RS)
    cs <- c(cs, ci)
  }
  fit$name <- name
  fit$cs <- cs
  return(list(fit))
}




#CoxBoost
my_CoxBoost <- function(training, val, name){
  pen <- optimCoxBoostPenalty(training[,'OS.time'],
                              training[,'OS'],
                              as.matrix(training[,-c(1,2)]),
                              trace=TRUE,
                              #start.penalty=500,
                              parallel = T)
  
  cv.res <- cv.CoxBoost(training[,'OS.time'],
                        training[,'OS'],
                        as.matrix(training[,-c(1,2)]),
                        maxstepno=500,
                        K=5,
                        type="verweij",
                        penalty= pen$penalty,
                        multicore=1)
  fit <- CoxBoost(training[,'OS.time'],
                  training[,'OS'],
                  as.matrix(training[,-c(1,2)]),
                  stepno=cv.res$optimal.step,
                  penalty=pen$penalty)
  fit$subFeature <- names(coef(fit)[abs(coef(fit)) > 0])
  RS <- as.numeric(predict(fit, training[,-c(1,2)]))
  ci <- cindex(training,RS)
  fit$cindex <- ci
  
  cs <- c(ci)
  for (i in 1:length(val)) {
    data <- val[[i]]
    y <- data.matrix(Surv(time = data$OS.time, event = data$OS))
    x <- data[, -c(1,2)]
    RS <- as.numeric(predict(fit, x))
    ci <- cindex(data, RS)
    cs <- c(cs, ci)
  }
  fit$name <- name
  fit$cs <- cs
  return(list(fit))
}



#Rideg
my_Rideg <- function(training, val, name){
  y <- data.matrix(Surv(time = training$OS.time, event = training$OS))
  x <- training[, -c(1,2)]
  
  cv.fit = cv.glmnet(x = as.matrix(x), y = y, family = "cox", alpha = 0, nfolds = 10)
  fit = glmnet(x = x,y = y, family = "cox", alpha = 0, lambda = cv.fit$lambda.min)
  coef.min = coef(cv.fit, s = "lambda.min")  
  active.min = which(coef.min@i != 0) ## 找出那些回归系数没有被惩罚为0的
  lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] ## 提取基因名称
  fit$subFeature <- lasso_geneids
  RS <- as.numeric(predict(fit, as.matrix(training[, -c(1,2)])))
  ci <- cindex(training,RS)
  fit$cindex <- ci
  
  
  cs <- c(ci)
  for (i in 1:length(val)) {
    data <- val[[i]]
    y <- data.matrix(Surv(time = data$OS.time, event = data$OS))
    x <- data[, -c(1,2)]
    RS <- as.numeric(predict(fit, as.matrix(x)))
    ci <- cindex(data, RS)
    cs <- c(cs, ci)
  }
  fit$name <- name
  fit$cs <- cs
  return(list(fit))
}

#superpc
my_superpc <- function(training, val, name){
  y <- data.matrix(Surv(time = training$OS.time, event = training$OS))
  x <- training[, -c(1,2)]
  data <- list(x = as.data.frame(t(scale(x))),
               y = training$OS.time,
               censoring.status = training$OS,
               featurenames = colnames(x))
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) # default
  cv.fit <- suppressWarnings(superpc.cv(fit, data,
                                        n.threshold = 20, # default
                                        n.fold = 10,
                                        n.components = 3,
                                        min.features = 1,
                                        max.features = nrow(data$x),
                                        compute.fullcv = TRUE,
                                        compute.preval = TRUE))
  fit$threshold <- cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])]
  fit$subFeature = names(fit$feature.scores)[abs(fit$feature.scores) > 0.5]
  
  
  ff <- superpc.predict(fit, data, data, threshold = fit$threshold, n.components = 1)
  RS <- as.numeric(ff$v.pred) 
  ci <- cindex(training,RS)
  fit$cindex <- ci
  
  
  cs <- c(ci)
  for (i in 1:length(val)) {
    vals <- val[[i]]
    y <- data.matrix(Surv(time = vals$OS.time, event = vals$OS))
    x <- vals[, -c(1,2)]
    data2 <- list(x = as.data.frame(t(scale(x))),
                 y = vals$OS.time,
                 censoring.status = vals$OS,
                 featurenames = colnames(x))
    
    ff <- superpc.predict(fit, data, data2, threshold = fit$threshold, n.components = 1)
    RS <- as.numeric(ff$v.pred) 
    ci <- cindex(vals, RS)
    cs <- c(cs, ci)
  }
  
  fit$name <- name
  fit$cs <- cs
  return(list(fit))
}


#RSF
my_RSF <- function(training, val, name){
  fit <- rfsrc(Surv(OS.time,OS)~.,data = training,
               ntree = 200,
               nodesize = 5,
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T)
  fit$subFeature = var.select(fit, verbose = F)$topvars
  
  
  
  ff <- predict(fit, as.data.frame(training))
  RS <- as.numeric(ff$predicted) 
  ci <- cindex(training,RS)
  fit$cindex <- ci
  
  cs <- c(ci)
  for (i in 1:length(val)) {
    vals <- val[[i]]
    ff <- predict(fit, as.data.frame(vals))
    RS <- as.numeric(ff$predicted) 
    ci <- cindex(vals,RS)
    cs <- c(cs, ci)
  }
  
  fit$name <- name
  fit$cs <- cs
  return(list(fit))
}


#Enet

my_Enet <- function(training, val, name){
  y <- data.matrix(Surv(time = training$OS.time, event = training$OS))
  x <- training[, -c(1,2)]
  models <- list()
  for (i in (1:9)/10) {
    cv.fit = cv.glmnet(x = as.matrix(x), y = y, family = "cox", alpha = i, nfolds = 10)
    fit = glmnet(x = x,y = y, family = "cox", alpha = i, lambda = cv.fit$lambda.min)
    coef.min = coef(cv.fit, s = "lambda.min")  
    active.min = which(coef.min@i != 0) ## 找出那些回归系数没有被惩罚为0的
    lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] ## 提取基因名称
    fit$subFeature <- lasso_geneids
    RS <- as.numeric(predict(fit, as.matrix(training[, -c(1,2)])))
    ci <- cindex(training,RS)
    fit$cindex <- ci
    
    name2 <- paste0(name, "[", i, "]")
    cs <- c(ci)
    for (j in 1:length(val)) {
      data <- val[[j]]
      RS <- as.numeric(predict(fit, as.matrix(data[, -c(1,2)])))
      ci <- cindex(data, RS)
      cs <- c(cs, ci)
    }
    fit$name <- name2
    fit$cs <- cs
    
    models[[i*10]] <- fit
  }
  return(models)
}



#plsRcox
my_plsRcox <- function(training, val, name){
  y <- data.matrix(Surv(time = training$OS.time, event = training$OS))
  x <- training[, -c(1,2)]
  data <- list(x = x,
               time = training$OS.time,
               status = training$OS)
  cv.plsRcox.res = cv.plsRcox(data = data,
                              nt=10, verbose = FALSE)
  fit <- plsRcox(Xplan = data$x,
                 time = data$time,
                 event = data$status,
                 nt = as.numeric(cv.plsRcox.res[5]),
                 verbose = F, sparse = T)
  fit$subFeature = rownames(fit$Coeffs)[fit$Coeffs != 0]
  
  RS <- as.numeric(predict(fit, as.matrix(x)))
  ci <- cindex(training,RS)
  fit$cindex <- ci
  
  
  cs <- c(ci)
  for (i in 1:length(val)) {
    vals <- val[[i]]
    y <- data.matrix(Surv(time = vals$OS.time, event = vals$OS))
    x <- vals[, -c(1,2)]
    data <- list(x = x,
                 time = vals$OS.time,
                 status = vals$OS)
    RS <- as.numeric(predict(fit, as.matrix(x)))
    ci <- cindex(vals,RS)
    cs <- c(cs, ci)
  }
  
  fit$name <- name
  fit$cs <- cs
  return(list(fit))
}


#GBM
my_GBM <- function(training, val, name){
  y <- data.matrix(Surv(time = training$OS.time, event = training$OS))
  x <- training[, -c(1,2)]
  
  fit <- gbm(formula = Surv(training$OS.time, training$OS) ~ .,
             data = as.data.frame(x),
             distribution = 'coxph',
             n.trees = 250,
             interaction.depth = 5,
             n.minobsinnode = 5,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  best <- which.min(fit$cv.error)
  fit <- gbm(formula = Surv(training$OS.time, training$OS) ~ .,
             data = as.data.frame(x),
             distribution = 'coxph',
             n.trees = best,
             interaction.depth = 5,
             n.minobsinnode = 10,
             shrinkage = 0.001, n.cores = 8)
  fit$subFeature = rownames(summary.gbm(fit, plotit = FALSE))[summary.gbm(fit, plotit = FALSE)$rel.inf > 0]
  
  RS <- as.numeric(predict(fit, as.data.frame(x)))
  ci <- cindex(training,RS)
  fit$cindex <- ci
  
  cs <- c(ci)
  for (i in 1:length(val)) {
    vals <- val[[i]]
    y <- data.matrix(Surv(time = vals$OS.time, event = vals$OS))
    x <- vals[, -c(1,2)]
    RS <- as.numeric(predict(fit, as.data.frame(x)))
    ci <- cindex(vals,RS)
    cs <- c(cs, ci)
  }
  
  fit$name <- name
  fit$cs <- cs
  return(list(fit))
}









SimpleHeatmap <- function(Cindex_mat = NULL, 
                          avg_Cindex = NULL, 
                          CohortCol = NULL, 
                          barCol = NULL,
                          col = c("#4195C1", "#FFFFFF", "#CB5746"), # 红蓝配色
                          cellwidth = 1, cellheight = 0.5, 
                          cluster_columns, cluster_rows){
  col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                            col = list("Cohort" = CohortCol),
                            show_annotation_name = F)
  
  row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                            gp = gpar(fill = barCol, col = NA),
                                            add_numbers = T, numbers_offset = unit(-10, "mm"),
                                            axis_param = list("labels_rot" = 0),
                                            numbers_gp = gpar(fontsize = 9, col = "white"),
                                            width = unit(3, "cm")),
                         show_annotation_name = F)
  
  Heatmap(as.matrix(Cindex_mat), name = "C-index",
          right_annotation = row_ha, 
          top_annotation = col_ha,
          col = col, 
          rect_gp = gpar(col = "black", lwd = 1), # 边框设置为黑色
          cluster_columns = cluster_columns, cluster_rows = cluster_rows, # 不进行聚类，无意义
          show_column_names = FALSE, 
          show_row_names = TRUE,
          row_names_side = "left",
          width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
          height = unit(cellheight * nrow(Cindex_mat), "cm"),
          column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
          column_title = NULL,
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                      x, y, gp = gpar(fontsize = 10))
          }
  )
}













