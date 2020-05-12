### clear envirinment
rm(list=ls())

# 1. Packages- install -----------------------------------------------
list.of.packages <- c("BGLR","grid", "futile.logger", "glmnet", "tidyverse", "dplyr","ggfortify", "ggthemes", "tableone", "formattable", "gplots", "heatmap", "pheatmap") ### add all the packages that are used in your script here
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])] ### check if they are installed
if(length(new.packages)) install.packages(new.packages) ### install any that need installing

lapply(list.of.packages, require, character.only = TRUE)

### Load BiocManager Libraries

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

### Install Venn Diagram package from Bio of required
if(!"VennDiagram" %in% installed.packages()){
  BiocManager::install("VennDiagram")
}

library("VennDiagram")

### Save output file path 
out.path<-"~/Documents/shati/" #insert your file pathway here

# 2. Data import and clean --------------------------------------------------

data(mice)
#please note data exploration was conducted earlier on- i.e. missing data checked for and accounted for/exploratory visualisations conducted
#we have 4 litters- proxy for different ethnic groups
#we could have denoised the litter (but for this case was not necessary)

# check for missing values in mice.X
sum(is.na(mice.X)) == 0 
# no missing vlaues in mice.X

# check for missing values in phenotype
sum(is.na(mice.pheno)) == 0
# no missing vlaues in phenotype


# 3. Loop to create PRS for each litter -----------------------------------------------------

for (lit in 1:8){
  
  print(paste("Analysing litter number", lit, "..."))
  
  Litter1 <- subset.data.frame(mice.pheno, Litter==lit)
  
  # select litter subjects in miceX
  mice.X1 <- mice.X[match(Litter1$SUBJECT.NAME, rownames(mice.X)),]
  
  # keeping the variables for the phenotype (aka obesity.BMI)
  df1 = Litter1[c("SUBJECT.NAME", "Obesity.BMI")]
  
  # check the Y is the same order as the mice.X
  all(df1$SUBJECT.NAME==rownames(mice.X1))
  
  # # turn them as matrices and vectors
  mice.X1 <- as.matrix(mice.X1)
  typeof(mice.X1)
  
  
  
  ### Split training and testing ###
  
  # Split mice.X (Litter1) (snps) into training data
  set.seed(1)
  train1 <- sample(1:nrow(mice.X1), 0.7*nrow(mice.X1), replace=FALSE)
  snps_train1 <- mice.X1[train1,]
  
  # Split pheno.type (obesity) into training data
  phenotypes_train1 <- df1[train1,]
  
  # split mice.X(snps) into test data
  snps_test1 <- mice.X1[-train1,]
  
  # split pheno,type (obesity) into test data
  phenotypes_test1 <- df1[-train1,]
  
  
  #check
  stopifnot(length(train1)==nrow(phenotypes_train1))
  stopifnot(nrow(snps_train1)==nrow(phenotypes_train1))
  stopifnot(nrow(snps_test1)==nrow(phenotypes_test1))
  
  
  ### Cross validating glmnet ###
  
  # Run penalized linear regression (EN, alpha=0.5) for obesity
  
  X_train1 <- snps_train1
  Y_train1 <- phenotypes_train1[,"Obesity.BMI"]
  X_test1 <- snps_test1
  Y_test1 <- phenotypes_test1[,"Obesity.BMI"]
  
  
  ### cannot build or crossvalidate a model with < 20 obs. Next litter in loop if so
  if(nrow(X_train1) < 20){
    print(paste("Insufficient observations in litter", lit, "for modelling"))
    print(paste("Number of obs ==", nrow(X_train1), "Moving to next litter"))
    
    next
  }
  
  set.seed(2)
  model.EN.obesityBMI.Litter1 <- cv.glmnet(x=X_train1, y=Y_train1, 
                                           nfolds = min(ceiling(nrow(X_train1)/5),10), ### to ensure 5 obs per fold
                                           family="gaussian", alpha=0.5, type.measure="deviance")
  

  save(model.EN.obesityBMI.Litter1, file=paste0(out.path,paste0("model.EN.obesityBMI.Litter_",lit,".rdata")))
  
  pdf(paste0(out.path,paste0("litter_",lit,"on",lit,"_LitterENplot.pdf")) )
  plot(model.EN.obesityBMI.Litter1)
  dev.off()
  
  # coef(model.EN.obesityBMI.Litter1, s = "lambda.min")
  
  optimal.lambda1<-model.EN.obesityBMI.Litter1$lambda.min
  
  # Working on test data ###
  EN.obesityBMI.pred.Litter1 <- predict(model.EN.obesityBMI.Litter1, newx = X_test1, s = optimal.lambda1,type="response")
  EN.obesityBMI.pred.Litter1.min <- predict(model.EN.obesityBMI.Litter1, newx = X_test1, s = "lambda.min",type="response")
  EN.obesityBMI.pred.Litter1.1se <- predict(model.EN.obesityBMI.Litter1, newx = X_test1, s = "lambda.1se",type="response")
  
  
  #MSE
  MSELitter1EN<-mean((EN.obesityBMI.pred.Litter1.min - Y_test1)^2)
  save(MSELitter1EN, file=paste0(out.path,paste0("MSElitter",lit,"EN.RData")))
  
  
  ### Calculating PRS for litter ###

  
  # getting Betas
  plot(model.EN.obesityBMI.Litter1)
  beta1 <- coef(model.EN.obesityBMI.Litter1,s="lambda.min") # NB. coef is giving the wrong values without s=
  summary(beta1)

  # converting s4 into dataframe
  beta1_selected<-data.frame(name = beta1@Dimnames[[1]][beta1@i + 1], coefficient = beta1@x)
  
  
  # remove intercept from beta
  beta1_selectedt=beta1_selected[-1,]
  
  ### discover number of nonzero betas. If zero, break look
  if(dim(beta1_selectedt)[1] == 0){
    print(paste("No SNPS selected in litter", lit, ". Moving to next litter"))
    next
  }
  
  # make sure it is a dataframe
  b1<-as.data.frame(beta1_selectedt)
  
  # some plots
  pdf(paste0 (out.path, paste0("histbetaselectedENlitter_",lit,".pdf")))
  hist(b1$coefficient)
  dev.off()
  
  pdf(paste0 (out.path,paste0("boxplotselectedvariableslitter_",lit,"EN.pdf")))
  boxplot(b1$coefficient, main="boxplot of beta coefficients of SNPs selected in litter 1")
  dev.off()
  
  # transpose mice.X1
  mice.X1t<-t(mice.X1)
  dim(mice.X1t)
  
  # turn mice.X1t into a dataframe
  mice.X1t<-data.frame(mice.X1t)
  
  # select the snps and alleles from b1 in mice.X1, aka make sure the same snps have been selected
  selected_snps_mice1 <- mice.X1t[match(b1$name, rownames(mice.X1t)),]
  dim(selected_snps_mice1)#we want matrix, 86rows by 613mice
  
  # check to see if b2 and selected_snps_mice2 are in the same order
  all(b1$name==rownames(selected_snps_mice1)) #they are in same order
  
  # sanity checks
  stopifnot((b1$name)==rownames(selected_snps_mice1))
  nrow(selected_snps_mice1) == nrow(b1)
  
  b1coefficients<-b1$coefficient
  b1coefficients<-as.matrix(b1coefficients)
  
  # sweep to calculate PRS scores
  resultsl1<-sweep(selected_snps_mice1, MARGIN=1,b1coefficients, `*`)
  dim (resultsl1)
  
  # getting PRS score per individual mouse- summing each row
  
  PRSscoresforlitter1EN<-colSums(resultsl1)
  save(PRSscoresforlitter1EN, file= paste0(out.path, paste0("ENPRSscoresforlitter_",lit,".RData")))
  hist(PRSscoresforlitter1EN)
  
  print(paste("PRS on litter", lit, "complete"))
  
}



# 4. PRS on all litters -------------------------------------------------------------


# Sorting all litters 

data(mice)
mice.Xx<-mice.X


# keeping the variables for the phenotype (aka obesity.BMI)
dfalllitter = mice.pheno[c("SUBJECT.NAME", "Obesity.BMI")]

# check the Y is the same order as the mice.X
all(dfalllitter$SUBJECT.NAME==rownames(mice.Xx))

# # turn them as matrices and vectors
mice.Xx <- as.matrix(mice.Xx)
typeof(mice.Xx)


# Splitting training and testing data -------------------------------------

set.seed(1)
trainx <- sample(1:nrow(mice.Xx), 0.7*nrow(mice.Xx), replace=FALSE)
snps_trainx<- mice.Xx[trainx,]

# Split pheno.type (obesity) into training data
phenotypes_trainx <- dfalllitter[trainx,]

# split mice.X(snps) into test data
snps_testx <- mice.Xx[-trainx,]

# split pheno,type (obesity) into test data
phenotypes_testx <- dfalllitter[-trainx,]


stopifnot(length(trainx)==nrow(phenotypes_trainx))
stopifnot(nrow(snps_trainx)==nrow(phenotypes_trainx))
stopifnot(nrow(snps_testx)==nrow(phenotypes_testx))


# Cross validating glmnet -------------------------------------------------

# Run penalized linear regression for obesity

X_trainx <- snps_trainx
Y_trainx <- phenotypes_trainx[,"Obesity.BMI"]
X_testx <- snps_testx
Y_testx <- phenotypes_testx[,"Obesity.BMI"]


set.seed(2)

model.EN.obesityBMI.Litterx <- cv.glmnet(x=X_trainx, y=Y_trainx, family="gaussian", alpha=0.5, type.measure="deviance")


save(model.EN.obesityBMI.Litterx, file= paste0(out.path,"obesityonlyalllitersxEN.RData"))

pdf(paste0(out.path,"MSc_Project/alllittersENlot.pdf")) 
plot(model.EN.obesityBMI.Litterx)
dev.off()


coef(model.EN.obesityBMI.Litterx, s = "lambda.min")

optimal.lambdax<-model.EN.obesityBMI.Litterx$lambda.min


# Working on test data ----------------------------------------------------

EN.obesityBMI.pred.Litterx<- predict(model.EN.obesityBMI.Litterx, newx = X_testx, s = optimal.lambdax,type="response")
EN.obesityBMI.pred.Litterx.min <- predict(model.EN.obesityBMI.Litterx, newx = X_testx, s = "lambda.min",type="response")
EN.obesityBMI.pred.Litterx.1se <- predict(model.EN.obesityBMI.Litterx, newx = X_testx, s = "lambda.1se",type="response")

#MSE
MSELitterx<-mean((EN.obesityBMI.pred.Litterx.min - Y_testx)^2)
save(MSELitterx, file=paste0(out.path, "MSElitterxEN.RData"))



# Calculating PRS ---------------------------------------------------------

# getting Beta's
plot(model.EN.obesityBMI.Litterx)
betax <- coef(model.EN.obesityBMI.Litterx,s="lambda.min") # NB. coef is giving the wrong values without s=
summary(betax)

# converting s4 into datafram
betax_selected<-data.frame(name = betax@Dimnames[[1]][betax@i + 1], coefficient = betax@x)


# remove intercept from beta
betax_selectedt=betax_selected[-1,]
dim(betax_selectedt)

# make sure it is a dataframe
bx<-as.data.frame(betax_selectedt)

# some plots
pdf(paste0(out.path, "histbetaselectedENlitterx.pdf")) 
hist(bx$coefficient)
dev.off()

pdf(paste0(out.path, "boxplotselectedvariableslitterxEN.pdf"))
boxplot(bx$coefficient, main="boxplot of beta coefficients of SNPs selected in litter x")
dev.off()

# now select snps from mice.Xx with b
dim(mice.Xx)

# transpose mice.X1
mice.Xxt<-t(mice.Xx)
dim(mice.Xxt)

# turn mice.Xt into a dataframe
mice.Xxt<-data.frame(mice.Xxt)

# select the snps and alleles from b in mice.X, aka make sure the same snps have been selected
selected_snps_micex <- mice.Xxt[match(bx$name, rownames(mice.Xxt)),]
dim(selected_snps_micex)


# check to see if b and selected_snps_mice are in the same order
all(bx$name==rownames(selected_snps_micex)) #they are in same order


# sanity check
stopifnot((bx$name)==rownames(selected_snps_micex))

dim(selected_snps_micex)
dim(bx)


bxcoefficients<-bx$coefficient
bxcoefficients<-as.matrix(bxcoefficients)

# sweep to calculate PRS scores
resultslx<-sweep(selected_snps_micex, MARGIN=1,bxcoefficients, `*`)
dim (resultslx)


# getting PRS score per individual mouse- summing each row
PRSscoresforalllittersEN<-colSums(resultslx)
save(PRSscoresforalllittersEN, file=paste0(out.path,"PRSscoresforalllittersEN.RData"))
hist(PRSscoresforalllittersEN)


# Plot PRS scores ---------------------------------------------------------

scriptName <- "*PRS calculated by EN, no score calculated for Litters 5-8"

footnote <- paste(scriptName)

# default footnote is today's date, cex=.7 (size) and color
# is a kind of grey


makeFootnote <- function(footnoteText=
                           format(scriptName),
                         size= .7, color= grey(.5))
{
  require(grid)
  pushViewport(viewport())
  grid.text(label= footnoteText ,
            x = unit(1,"npc") - unit(2, "mm"),
            y= unit(2, "mm"),
            just=c("right", "bottom"),
            gp=gpar(cex= size, col=color))
  popViewport()
}

makeFootnote(footnote)



png(paste0(out.path,"DistributionofPRSscoresEN.png")) 
plot(density(PRSscoresforlitter1EN, adjust=1),col="dark blue", main="Distribution of PRS scores- Elastic Net", xlab="PRS score", ylab="Density", xlim=c(-0.09,0.12), ylim=c(0,25),)
lines(density(PRSscoresforlitter2EN, adjust=1),col="dark green")
lines(density(PRSscoresforlitter3EN, adjust=1),col="black")
lines(density(PRSscoresforlitter4EN, adjust=1),col="light blue")
lines(density(PRSscoresforalllittersEN, adjust-1),col="red",lty=6,)
makeFootnote(footnote)



# add a legend
legend(0.070, 25,
       legend = c("PRS for litter 1", "PRS for litter 2", "PRS for litter 3", "PRS for litter 4", "PRS for all litters"), 
       col = c("dark blue", "dark green", "black", "light blue", "red"),
       bty= "n",
       lty=1:1,
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1),
       title ="Legend",
       cex=0.6
)   




dev.off()




pdf(paste0(out.path,"DistributionofPRS-EN.pdf")) 
plot(density(PRSscoresforlitter1EN, adjust=1),col="dark blue", main="Distribution of PRS - Elastic Net", xlab="PRS score", ylab="Density", xlim=c(-0.09,0.12), ylim=c(0,25),)
lines(density(PRSscoresforlitter2EN, adjust=1),col="dark green")
lines(density(PRSscoresforlitter3EN, adjust=1),col="black")
lines(density(PRSscoresforlitter4EN, adjust=1),col="light blue")
lines(density(PRSscoresforalllittersEN, adjust=1),col="red",lty=6,)
makeFootnote(footnote)

# add a legend
legend(0.070, 25,
       legend = c("PRS for litter 1", "PRS for litter 2", "PRS for litter 3", "PRS for litter 4", "PRS for all litters"), 
       col = c("dark blue", "dark green", "black", "light blue", "red"),
       bty= "n",
       lty=1:1,
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1),
       title ="Legend",
       cex=0.6
)  

dev.off()

hist(PRSscoresforlitter1EN)
hist(PRSscoresforalllittersEN)


dev.off()


# Plotting for Elastic Net Lambdas ----------------------------------------


# basic information at the beginning of each script
scriptName <- "*Elastic Net failed on Litter 8"

footnote <- paste(scriptName)

# default footnote is today's date, cex=.7 (size) and color
# is a kind of grey

makeFootnote <- function(footnoteText=
                           format(scriptName),
                         size= .7, color= grey(.5))
{
  require(grid)
  pushViewport(viewport())
  grid.text(label= footnoteText ,
            x = unit(1,"npc") - unit(2, "mm"),
            y= unit(2, "mm"),
            just=c("right", "bottom"),
            gp=gpar(cex= size, col=color))
  popViewport()
}

makeFootnote(footnote)


png(paste0(out.path,"plotsforEN1.png")) 
par(mfrow=c(4,2))
plot(model.EN.obesityBMI.Litter1) 
title("Litter 1", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litter2) 
title("Litter 2", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litter3) 
title("Litter 3", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litter4) 
title("Litter 4", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litter5) 
title("Litter 5", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litter6) 
title("Litter 6", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litter7) 
title("Litter 7", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litterx)
title("All litters combined", adj=0.5, line=2.5)
par(oma=c(0,0,0.6,0))
title("Elastic Net Regularization" , outer=TRUE)
makeFootnote(footnote)
dev.off()


pdf(paste0(out.path,"plotsforElasticNetfinal.pdf")) 
par(mfrow=c(4,2))
plot(model.EN.obesityBMI.Litter1) 
title("Litter 1", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litter2) 
title("Litter 2", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litter3) 
title("Litter 3", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litter4) 
title("Litter 4", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litter5) 
title("Litter 5", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litter6) 
title("Litter 6", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litter7) 
title("Litter 7", adj=0.5, line=2.5)
plot(model.EN.obesityBMI.Litterx)
title("All Litters combined", adj=0.5, line=2.5)
par(oma=c(0,0,0.5,0))
title("MSE for Elastic Net per Litter " , outer=TRUE)
makeFootnote(footnote)
dev.off()





# Venn Diagrams -----------------------------------------------------------



rm(list=ls())

getwd()

snps_selected_EN<- read.csv(paste0(out.path, 'selected_snps_dataset.csv'))

attach(snps_selected)

#setting up Venn diagram for mRNA, for EN
All<-snps_selected_EN$All.litters
Litter1<-snps_selected_EN$Litter1
Litter2<-snps_selected_EN$Litter2
Litter3<-snps_selected_EN$Litter3
Litter4<-snps_selected_EN$Litter4

mysets= list(All=All, Litter1=Litter1, Litter2=Litter2, Litter3=Litter3, Litter4=Litter4) 

colors = c("tomato", "forestgreen", "skyblue", "gold",
           "navy")
names = paste0(c("All litters", "Litter 1", "Litter 2", "Litter 3", "Litter 4"))

venn.diagram(mysets, filename = paste0(out.path, "MSc_Project/Venn_EN.png"),
             fill = colors, category.names = names, imagetype = "png",
             cat.just = list(c(0.5, 1), c(-0.5, -5), c(0, 0),
                             c(0.5, 0), c(1, -4)), main="Selected SNPs-Elastic Net"
)

dev.off()


# Pheatmap for Elastic Net ------------------------------------------------

# creating pheatmap

pdf(paste0(out.path,"pheatmapsforsnpsElasticNet.pdf"))

beta1[1:length(beta1)]
heatmap.df <- cbind(beta1[1:length(beta1)],beta2[1:length(beta2)],beta3[1:length(beta3)],beta4[1:length(beta4)],betax[1:length(betax)])
cormat <- cor(heatmap.df)
rownames(cormat)<- c("Litter1","Litter2","Litter3","Litter4","All_litters")
colnames(cormat)<- c("Litter1","Litter2","Litter3","Litter4","All_litters")

pheatmap(cormat, 
         main= "Heatmaps for Litters- Elastic Net", 
         cluster_rows=FALSE, cluster_cols=FALSE,
         #turn off clustering so rows and cols don't reorder
         show_rownames=TRUE,
         show_colnames=TRUE,
)
dev.off()


# Running T-tests ---------------------------------------------------------


write.table(PRSscoresforalllittersEN, "~/documents/Msc_Project/scoreforallls.txt", sep="\t")
write.table(PRSscoresforlitter1EN, "~/documents/Msc_Project/score1lsEN.txt", sep="\t")
write.table(PRSscoresforlitter2EN, "~/documents/Msc_Project/score2lsEN.txt", sep="\t")
write.table(PRSscoresforlitter3EN, "~/documents/Msc_Project/score3lsEN.txt", sep="\t")
write.table(PRSscoresforlitter4EN, "~/documents/Msc_Project/score4lsEN.txt", sep="\t")


# RMSE loop ---------------------------------------------------------------


# Adding subject name col to mice.X
mice.X <- as.data.frame(mice.X)
mice.X$SUBJECT.NAME <- rownames(mice.X)
mice.X <- merge(mice.X, mice.pheno[,c("Litter","SUBJECT.NAME")], by =  "SUBJECT.NAME")
rownames(mice.X) <- mice.X$SUBJECT.NAME
rownames(mice.pheno) <- mice.pheno$SUBJECT.NAME

models.list <- list()
predictions.matrix <- as.data.frame(matrix(nrow=nrow(mice.pheno), ncol = 8))
predictions.matrix$SUBJECT.NAME <- mice.X$SUBJECT.NAME

set.seed(2)


for (i in 1:max(mice.X$Litter)){
  X_train <- as.matrix(mice.X[mice.X$Litter == i,!names(mice.X) %in% c("Litter","SUBJECT.NAME")])
  Y_train <- as.matrix(mice.pheno[mice.pheno$Litter == i, "Obesity.BMI"])
  X <- as.matrix(mice.X[,!names(mice.X) %in% c("Litter","SUBJECT.NAME")])
  Y <- as.matrix(mice.pheno[, "Obesity.BMI"])
  mod <- cv.glmnet(x=X_train, y=Y_train, family="gaussian", alpha=0.5, type.measure="deviance")
  optimal.lambda <- mod$lambda.min
  models.list[[i]] <- mod
  predictions.matrix[,i] <- predict(mod, newx = X,s = optimal.lambda,type="response")
}

predictions.matrix <- merge(predictions.matrix, mice.X[,c("Litter", "SUBJECT.NAME")], by = "SUBJECT.NAME")
predictions.matrix <- predictions.matrix[,c(2:9,1,10)]


RMSE.matrix <- as.data.frame(matrix(nrow= 8, ncol = 8, dimnames = list(sprintf("test[%s]", seq(1:8)), sprintf("train[%s]", seq(1:8)))))

# i is the litter; j is the model
# Columns are models; rows are test data sets (by litter)

for (i in 1:8){
  for (j in 1:8){
    RMSE <- (mean((predictions.matrix[predictions.matrix$Litter == i, j] - mice.pheno[mice.pheno$Litter == i, "Obesity.BMI"])^2))^0.5
    RMSE.matrix[i,j] <- RMSE
  }
}


# pheatmap RMSE -----------------------------------------------------------


library(gplots)
pdf("~/documents/MSc_project/heatmapEN.pdf")
heatmap.2(as.matrix(RMSE.matrix),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none', main= "Root Mean Square Error-Elastic Net")
dev.off()


table(mice.X$Litter)


# RMSE as excel -----------------------------------------------------------

write.table(RMSE.matrix, "~/documents/Msc_Project/RMSEmatrix1_8.txt", sep="\t")



# PCA ---------------------------------------------------------------------
#conducting a PCA
pca.out <- prcomp(mice.X[,1:1000], scale=TRUE)
var.explained <- cumsum(pca.out$sdev^2.0)/sum(pca.out$sdev^2.0)
pdf(paste0(out.path,"CVEPCA.pdf")) 
plot(var.explained, main="Cumulative Variance Explained by Principal Components")
dev.off()

### Exploring PCA by litter
pca.predictor.df <-data.frame(pca.out$x)[1:100]
pca.predictor.df <- cbind(as.factor(mice.pheno$Litter), pca.predictor.df)
names(pca.predictor.df)[1] <- "Litter"

#Plot how specific litters vary from the rest of the litters
ggplot(data = pca.predictor.df, aes(y = PC1, x = PC2, col = Litter, fill = Litter)) + geom_point(alpha = 0.5) + theme_clean() #they should all be factors

# plot all litters on PCA graph with grouping frames
pdf(paste0(out.path,"plotsforPCA.pdf")) 
autoplot(pca.out, data = pca.predictor.df, col = "Litter", frame =T) + theme_clean() + ggtitle ("Visualizations of SNPs data by Litter, using PCA")
dev.off()


# Table 1 for write up ----------------------------------------------------

# need to turn litter into factors
mice.pheno$Litter<-as.factor(mice.pheno$Litter)

#We need the list of ALL variables we want (in your initial script you had just included the continuous ones in listVars)
listVars <- c("Obesity.BMI", "Litter", "GENDER", "DATE.MONTH", "CageDensity", "CoatColour","Date.season")
#Factors
catVars <- c("Gender", "Litter","CoatColour")


#Note that since you made sure all your factors are actually factors beforehand,
# I don't think the factorVars argument is necessary
table1 <-CreateTableOne( vars = listVars, data = mice.pheno, strata="GENDER")
print(table1)

tab1Mat <- print(table1, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
## Save to a CSV file
write.csv(tab1Mat, file = paste0(out.path, "myTable1.csv"))

