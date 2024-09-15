


Training_expr <- "Training_expr.txt"
Training_surv <- "Training_surv.txt"

Testing_expr <- "Testing_expr.txt"
Testing_surv <- "Testing_surv.txt"



train <- read.table(Training_expr, row.names = 1, check.names = F)
train_s <- read.table(Training_surv, row.names = 1, check.names = F)




training <- cbind(train_s, t(train))
training <- training[which(training$OS.time >= 1), ]
#training[, c(3:ncol(training))] <- log2(training[, c(3:ncol(training))] + 1)
#training$OS.time <- training$OS.time/30

val <- read.table(Testing_expr, row.names = 1, check.names = F)
val_s <- read.table(Testing_surv, row.names = 1, check.names = F)
colnames(val_s)[1:2] <- c("OS.time", "OS")
val_s <- val_s[which(val_s$OS.time >=1), ]
val_s <- val_s[which(val_s$OS <=1), ]
val <- val[, rownames(val_s)]
#val_s$OS.time <- val_s$OS.time/30

val_n <- list()
Cohort <- unique(val_s$Cohort)
for (i in 1:length(Cohort)) {
  val1 <- rownames(val_s)[which(val_s$Cohort == Cohort[i])]
  val1 <- cbind(val_s[val1, c("OS", "OS.time")], t(val[, val1]))
  val_n[[i]] <- val1
}


val <- val_n


m <- c("stepcox","GBM","RSF", "superpc", "lasso", "survivalsvm", "CoxBoost", "Rideg", "Enet", "plsRcox")
#m <- c("lasso","RSF")
#m <- c("lasso", "superpc")



cindexs <- as.data.frame(matrix(nrow = length(Cohort)+1, ncol = 1))
colnames(cindexs) <- m[1]
rownames(cindexs) <- c("Training", unique(val_s[,3]))
models <- list()


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

library(RColorBrewer)
library(ComplexHeatmap)

source("ML.R")



for(i in m){
  name <- i
  mm <- m[-which(m==i)]
  eval(parse(text = paste0("fit1 <- my_", i, "(training, val, name)")))

  for (j in 1:length(fit1)) {
    fit <- fit1[[j]] 
    models[[fit$name]] <- fit
    
    
    gene <- fit$subFeature
    for (o in 1:length(gene)) {
      gene[o] <- gsub('[`"`]', '', gene[o])
    }
    name1 <- fit$name
    print(name1)
    cindexs[, fit$name] <- fit$cs
    
    training_2 <- cbind(training[,1:2], training[,gene[which(gene %in% colnames(training))]])
    if(ncol(training_2)<=4){
      next
    }
    val12 <- cbind(val[[1]][,1:2], val[[1]][,gene])
    val22 <- cbind(val[[2]][,1:2], val[[2]][,gene])
    val32 <- cbind(val[[3]][,1:2], val[[3]][,gene])
    
    
    
    val_2 <- list(val12, val22, val32)
    
    for(k in mm){
      name2 <- paste0(name1, " + ", k)
      if(name2 %in% names(models)){
        next
      }
      print(name2)
      
      #fit2 <- my_superpc(training_2, val_2, name2)
      eval(parse(text = paste0("fit2 <- my_", k, "(training_2, val_2, name2)")))
      for (u in 1:length(fit2)) {
        fit <- fit2[[u]] 
        models[[fit$name]] <- fit
        
        cindexs[, fit$name] <- fit$cs
      }
      
    }
  }
  
}







save(models, file = "models.rda")
cindexs <- t(cindexs)
write.csv(cindexs, "cindexs.csv")

Cindex_mat <- read.csv("cindexs.csv", check.names = F, row.names = 1, header = T)
colnames(Cindex_mat)[1] <- "TCGA"

ts <- Cindex_mat[,-1]
ts <- Cindex_mat
avg_Cindex <- sort(apply(ts, 1, mean), decreasing = T) # 计算每种算法在所有队列中平均C-index，并降序排列


Cindex_mat <- Cindex_mat[names(avg_Cindex), ] # 对C-index矩阵排序
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数
#fea_sel <- fea_list[[rownames(Cindex_mat)[1]]] # 最优模型（即测试集[或者训练集+测试集]C指数均值最大）所筛选的特征

CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired")# 设置绘图时的队列颜色
names(CohortCol) <- colnames(Cindex_mat)

# 调用简易绘图函数
cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = Cindex_mat, # 主矩阵0..
                    avg_Cindex = avg_Cindex, # 侧边柱状图
                    CohortCol = CohortCol, # 列标签颜色
                    barCol = "steelblue", # 右侧柱状图颜色
                    col = c("#1CB8B2", "#FFFFFF", "#EEB849",
                            "red", "blue","green"), # 热图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是否对行列进行聚类

pdf(file.path("heatmap of cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 7, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right") # 热图注释均放在右侧
invisible(dev.off())









