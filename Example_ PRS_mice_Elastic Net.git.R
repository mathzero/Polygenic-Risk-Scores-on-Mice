
# Packages- install -----------------------------------------------
list.of.packages <- c("BGLR", "glmnet", "tidyverse", "dplyr","ggfortify", "ggthemes", "tableone", "formattable", "gplots", "heatmap", "pheatmap") ### add all the packages that are used in your script here
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])] ### check if they are installed
if(length(new.packages)) install.packages(new.packages) ### install any that need installing

lapply(list.of.packages, require, character.only = TRUE)


#please note while here I have saved everything to my personal files you can reproduce the plots etc: by using
#out.path<-"~/Users/shatirahman/PRS"
#save(paste00(out.path,"nameofile.rdata")



# Load Libraries ----------------------------------------------------------
source("http://bioconductor.org/biocLite.R")
install.packages("biocLite")
install.packages("BiocManager")
BiocManager::install("VennDiagram")
library(grid)
library(futile.logger)
library("VennDiagram")



# Select the data to use --------------------------------------------------
data(mice)
#please note data exploration was conducted earlier on- i.e. missing data checked for and accounted for/exploratory visualisations conducted
#we have 7 litters- proxy for different ethnic groups
#we could have denoised the litter (but for this case was not necessary)


# Sorting by litter 1 -----------------------------------------------------

litter1 <- subset.data.frame(mice.pheno, Litter=='1')

# select litter1 subjects in miceX
mice.X1 <- mice.X[match(Litter1$SUBJECT.NAME, rownames(mice.X)),]
#View(mice.X1)
dim(mice.X1)

# keeping the variables for the phenotype (aka obesity.BMI)
df1 = Litter1[c("SUBJECT.NAME", "Obesity.BMI")]

# check the Y is the same order as the mice.X
all(df1$SUBJECT.NAME==rownames(mice.X1))

# # turn them as matrices and vectors
mice.X1 <- as.matrix(mice.X1)
typeof(mice.X1)



# Split training and testing ----------------------------------------------


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



# Cross validating glmnet -------------------------------------------------

# Run penalized linear regression (EN, alpha=0.5) for obesity

X_train1 <- snps_train1
Y_train1 <- phenotypes_train1[,"Obesity.BMI"]
X_test1 <- snps_test1
Y_test1 <- phenotypes_test1[,"Obesity.BMI"]


set.seed(2)
model.EN.obesityBMI.Litter1 <- cv.glmnet(x=X_train1, y=Y_train1, family="gaussian", alpha=0.5, type.measure="deviance")

save(model.EN.obesityBMI.Litter1, file="~/documents/Msc_Project/obesityonlyLitter1EN.RData")

pdf("~/documents/MSc_Project/litter1on1LitterENplot.pdf") 
plot(model.EN.obesityBMI.Litter1)
dev.off()

coef(model.EN.obesityBMI.Litter1, s = "lambda.min")

optimal.lambda1<-model.EN.obesityBMI.Litter1$lambda.min

# Working on test data  ---------------------------------------------------
EN.obesityBMI.pred.Litter1 <- predict(model.EN.obesityBMI.Litter1, newx = X_test1, s = optimal.lambda1,type="response")
EN.obesityBMI.pred.Litter1.min <- predict(model.EN.obesityBMI.Litter1, newx = X_test1, s = "lambda.min",type="response")
EN.obesityBMI.pred.Litter1.1se <- predict(model.EN.obesityBMI.Litter1, newx = X_test1, s = "lambda.1se",type="response")


#MSE
MSELitter1EN<-mean((EN.obesityBMI.pred.Litter1.min - Y_test1)^2)
save(MSELitter1EN, file="~/documents/Msc_Project/MSElitter1EN.RData")


# Calculating PRS for litter ----------------------------------------------

# getting Beta's
plot(model.EN.obesityBMI.Litter1)
beta1 <- coef(model.EN.obesityBMI.Litter1,s="lambda.min") # NB. coef is giving the wrong values without s=
summary(beta1)

# converting s4 into dataframe
beta1_selected<-data.frame(name = beta1@Dimnames[[1]][beta1@i + 1], coefficient = beta1@x)


# remove intercept from beta
beta1_selectedt=beta1_selected[-1,]
dim(beta1_selectedt)

# make sure it is a dataframe
b1<-as.data.frame(beta1_selectedt)

# some plots
pdf("~/documents/MSc_Project/histbetaselectedENlitter1.pdf") 
hist(b1$coefficient)
dev.off()

pdf("~/documents/MSc_Project/boxplotselectedvariableslitter1EN.pdf") 
boxplot(b1$coefficient, main="boxplot of beta coefficients of SNPs selected in litter 1")
dev.off()

# now select snps from mice.X1 with b1
mice.X1
dim(mice.X1)

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

# sanity check
stopifnot((b1$name)==rownames(selected_snps_mice1))

dim(selected_snps_mice1)
dim(b1)

b1coefficients<-b1$coefficient
b1coefficients<-as.matrix(b1coefficients)

# sweep to calculate PRS scores
resultsl1<-sweep(selected_snps_mice1, MARGIN=1,b1coefficients, `*`)
dim (resultsl1)

# getting PRS score per individual mouse- summing each row

PRSscoresforlitter1EN<-colSums(resultsl1)
save(PRSscoresforlitter1EN, file="~/documents/Msc_Project/ENPRSscoresforlitter1.RData")
hist(PRSscoresforlitter1EN)

# Litter 2 ----------------------------------------------------------------

Litter2 <- subset.data.frame(mice.pheno, Litter=='2')


# Sorting the litter  -----------------------------------------------------

# Select the same mice (subject id ) from mice X, for Litter 2
# make header the same
# mice.X <- data.frame("SUBJECT.NAME"=rownames(mice.X), mice.X)

# select litter subjects in miceX
mice.X2 <- mice.X[match(Litter2$SUBJECT.NAME, rownames(mice.X)),]
#View(mice.X1)
dim(mice.X2)

# keeping the variables for the phenotype (aka obesity.BMI)
df2 = Litter2[c("SUBJECT.NAME", "Obesity.BMI")]


# check the Y is the same order as the mice.X
all(df2$SUBJECT.NAME==rownames(mice.X2))

# # turn them as matrices and vectors
mice.X2 <- as.matrix(mice.X2)
typeof(mice.X2)


# Split training and testing data -----------------------------------------


# Split mice.X (Litter) (snps) into training data
set.seed(1)
train2 <- sample(1:nrow(mice.X2), 0.7*nrow(mice.X2), replace=FALSE)
snps_train2 <- mice.X2[train2,]

# Split pheno.type (obesity) into training data
phenotypes_train2 <- df2[train2,]

# split mice.X(snps) into test data
snps_test2 <- mice.X2[-train2,]

# split pheno,type (obesity) into test data
phenotypes_test2 <- df2[-train2,]


stopifnot(length(train2)==nrow(phenotypes_train2))
stopifnot(nrow(snps_train2)==nrow(phenotypes_train2))
stopifnot(nrow(snps_test2)==nrow(phenotypes_test2))


# Cross validating glment -------------------------------------------------


# cross validating glmnet

# Run penalized linear regression for obesity

X_train2 <- snps_train2
Y_train2 <- phenotypes_train2[,"Obesity.BMI"]
X_test2 <- snps_test2
Y_test2 <- phenotypes_test2[,"Obesity.BMI"]

set.seed(2)
model.EN.obesityBMI.Litter2 <- cv.glmnet(x=X_train2, y=Y_train2, family="gaussian", alpha=0.5, type.measure="deviance")


save(model.EN.obesityBMI.Litter2, file="~/documents/Msc_Project/obesityonlyLitter2EN.RData")

pdf("~/documents/MSc_Project/litter2on2ENplot.pdf") 
plot(model.EN.obesityBMI.Litter2)
dev.off()

coef(model.EN.obesityBMI.Litter2, s = "lambda.min")

optimal.lambda2<-model.EN.obesityBMI.Litter2$lambda.min


# Working on test data ----------------------------------------------------

EN.obesity.BMI.pred.Litter2 <- predict(model.EN.obesityBMI.Litter2, newx = X_test2, s = optimal.lambda2,type="response")
EN.obesity.BMI.pred.Litter2.min <- predict(model.EN.obesityBMI.Litter2, newx = X_test2, s = "lambda.min",type="response")
EN.obesityBMI.pred.Litter2.1se <- predict(model.EN.obesityBMI.Litter2, newx = X_test2, s = "lambda.1se",type="response")

#MSE
MSELitter2EN<-mean((EN.obesityBMI.pred.Litter2.min - Y_test2)^2)
save(MSELitter2EN, file="~/documents/Msc_Project/MSElitter2EN.RData")


# Calculating PRS ---------------------------------------------------------


# getting Beta's
plot(model.EN.obesityBMI.Litter2)
beta2 <- coef(model.EN.obesityBMI.Litter2,s="lambda.min") # NB. coef is giving the wrong values without s=
summary(beta2)

# converting s4 into datafram
beta2_selected<-data.frame(name = beta2@Dimnames[[1]][beta2@i + 1], coefficient = beta2@x)

# remove intercept from beta
beta2_selectedt=beta2_selected[-1,]

# make sure it is a dataframe
b2<-as.data.frame(beta2_selectedt)

# some plots
pdf("~/documents/MSc_Project/histbetaselectedlassolitter2on2EN.pdf") 
hist(b2$coefficient)
dev.off()

pdf("~/documents/MSc_Project/boxplotselectedvariableslitter2EN.pdf") 
boxplot(b2$coefficient, main="boxplot of beta coefficients of SNPs selected in litter 2")
dev.off()

# now select snps from mice.X with b
mice.X2
dim(mice.X2)

# transpose mice.X
mice.X2t<-t(mice.X2)
dim(mice.X2t)

# turn mice.Xt into a dataframe
mice.X2t<-data.frame(mice.X2t)

# select the snps and alleles from b in mice.X, aka make sure the same snps have been selected
selected_snps_mice2 <- mice.X2t[match(b2$name, rownames(mice.X2t)),]
dim(selected_snps_mice2)#we want matrix, 128 rows by 515 col

# check to see if b and selected_snps_mice are in the same order
all(b2$name==rownames(selected_snps_mice2)) #they are in same order

# sanity check
stopifnot((b2$name)==rownames(selected_snps_mice2))

dim(selected_snps_mice2)
dim(b2)

b2coefficients<-b2$coefficient
b2coefficients<-as.matrix(b2coefficients)

# sweep to calculate PRS scores
resultsl2EN<-sweep(selected_snps_mice2, MARGIN=1,b2coefficients, `*`)
dim (resultsl2EN)


# getting PRS score per individual mouse- summing each row

PRSscoresforlitter2EN<-colSums(resultsl2EN)
save(PRSscoresforlitter2EN, file="~/documents/Msc_Project/PRSscoresforlitter2EN.RData")
hist(PRSscoresforlitter2EN)


# Litter 3 ----------------------------------------------------------------

# Sorting by litter -------------------------------------------------------


Litter3 <- subset.data.frame(mice.pheno, Litter=='3')

# Select the same mice (subject id ) from mice X, for Litter 3

# make header the same
# mice.X <- data.frame("SUBJECT.NAME"=rownames(mice.X), mice.X)

# select litter subjects in miceX
mice.X3 <- mice.X[match(Litter3$SUBJECT.NAME, rownames(mice.X)),]
dim(mice.X3)

# keeping the variables for the phenotype (aka obesity.BMI)
df3 = Litter3[c("SUBJECT.NAME", "Obesity.BMI")]

# check the Y is the same order as the mice.X
all(df3$SUBJECT.NAME==rownames(mice.X3))

# # turn them as matrices and vectors
mice.X3 <- as.matrix(mice.X3)
typeof(mice.X3)



# Split training and testing data -----------------------------------------

# Split mice.X (Litter) (snps) into training data
set.seed(1)
train3 <- sample(1:nrow(mice.X3), 0.7*nrow(mice.X3), replace=FALSE)
snps_train3 <- mice.X3[train3,]

# Split pheno.type (obesity) into training data
phenotypes_train3 <- df3[train3,]

# split mice.X(snps) into test data
snps_test3 <- mice.X3[-train3,]

# split pheno,type (obesity) into test data
phenotypes_test3 <- df3[-train3,]

stopifnot(length(train3)==nrow(phenotypes_train3))
stopifnot(nrow(snps_train3)==nrow(phenotypes_train3))
stopifnot(nrow(snps_test3)==nrow(phenotypes_test3))


# Cross validating glmnet -------------------------------------------------

# Run penalized linear regression for obesity

X_train3 <- snps_train3
Y_train3 <- phenotypes_train3[,"Obesity.BMI"]
X_test3 <- snps_test3
Y_test3 <- phenotypes_test3[,"Obesity.BMI"]

set.seed(2)
model.EN.obesityBMI.Litter3 <- cv.glmnet(x=X_train3, y=Y_train3, family="gaussian", alpha=0.5, type.measure="deviance")

save(model.EN.obesityBMI.Litter3, file="~/documents/Msc_Project/obesityLitter3EN.RData")

pdf("~/documents/MSc_Project/litter3ENplot.pdf") 
plot(model.EN.obesityBMI.Litter3)
dev.off()

coef(model.EN.obesityBMI.Litter3, s = "lambda.min")

optimal.lambda3<-model.EN.obesityBMI.Litter3$lambda.min


# Working on test data ----------------------------------------------------

EN.obesityBMI.pred.Litter3 <- predict(model.EN.obesityBMI.Litter3, newx = X_test3, s = optimal.lambda3,type="response")
EN.obesityBMI.pred.Litter3.min <- predict(model.EN.obesityBMI.Litter3, newx = X_test3, s = "lambda.min",type="response")
EN.obesityBMI.pred.Litter3.1se <- predict(model.EN.obesityBMI.Litter3, newx = X_test3, s = "lambda.1se",type="response")

#MSE
MSELitter3EN<-mean((EN.obesityBMI.pred.Litter3.min - Y_test3)^2)
save(MSELitter3, file="~/documents/Msc_Project/MSElitter3EN.RData")


# Calculating PRS ---------------------------------------------------------

# getting Beta's
plot(model.EN.obesityBMI.Litter3)
beta3 <- coef(model.EN.obesityBMI.Litter3,s="lambda.min") # NB. coef is giving the wrong values without s=
summary(beta3)

# converting s4 into datafram
beta3_selected<-data.frame(name = beta3@Dimnames[[1]][beta3@i + 1], coefficient = beta3@x)

# remove intercept from beta
beta3_selectedt=beta3_selected[-1,]
dim(beta3_selectedt)

# make sure it is a dataframe
b3<-as.data.frame(beta3_selectedt)

# some plots
pdf("~/documents/MSc_Project/histbetaselectedENlitter3.pdf") 
hist(b3$coefficient)
dev.off()

pdf("~/documents/MSc_Project/boxplotselectedvariableslitter3EN.pdf") 
boxplot(b3$coefficient, main="boxplot of beta coefficients of SNPs selected in litter 3")
dev.off()

# now select snps from mice.X with b
dim(mice.X3)

# transpose mice.X
mice.X3t<-t(mice.X3)
dim(mice.X3t)

# turn mice.Xt into a dataframe
mice.X3t<-data.frame(mice.X3t)

# select the snps and alleles from b in mice.X, aka make sure the same snps have been selected
selected_snps_mice3 <- mice.X3t[match(b3$name, rownames(mice.X3t)),]
dim(selected_snps_mice3)#we want matrix, 86rows by 613mice

# check to see if b and selected_snps_mice are in the same order
all(b3$name==rownames(selected_snps_mice3)) #they are in same order

# sanity check
stopifnot((b3$name)==rownames(selected_snps_mice3))

dim(selected_snps_mice3)
dim(b3)

b3coefficients<-b3$coefficient
b3coefficients<-as.matrix(b3coefficients)

# sweep to calculate PRS scores
resultsl3<-sweep(selected_snps_mice3, MARGIN=1,b3coefficients, `*`)
dim (resultsl3)

# getting PRS score per individual mouse- summing each row
PRSscoresforlitter3EN<-colSums(resultsl3)
save(PRSscoresforlitter3EN, file="~/documents/Msc_Project/PRSscoresforlitter3EN.RData")
hist(PRSscoresforlitter3EN)



# Litter 4 ----------------------------------------------------------------


# Sorting litter 4 --------------------------------------------------------



# train litter 4

Litter4 <- subset.data.frame(mice.pheno, Litter=='4')


# Select the same mice (subject id ) from mice X, for Litter 4
# make header the same
# mice.X <- data.frame("SUBJECT.NAME"=rownames(mice.X), mice.X)

# select litter subjects in miceX
mice.X4 <- mice.X[match(Litter4$SUBJECT.NAME, rownames(mice.X)),]
#View(mice.X4)
dim(mice.X4)

# keeping the variables for the phenotype (aka obesity.BMI)
df4 = Litter4[c("SUBJECT.NAME", "Obesity.BMI")]

# check the Y is the same order as the mice.X
all(df4$SUBJECT.NAME==rownames(mice.X4))

# # turn them as matrices and vectors
mice.X4 <- as.matrix(mice.X4)
typeof(mice.X4)

# Split training and testing data -----------------------------------------

# Split mice.X (Litter) (snps) into training data
set.seed(1)
train4 <- sample(1:nrow(mice.X4), 0.7*nrow(mice.X4), replace=FALSE)
snps_train4 <- mice.X4[train4,]

# Split pheno.type (obesity) into training data
phenotypes_train4 <- df4[train4,]

# split mice.X(snps) into test data
snps_test4 <- mice.X4[-train4,]

# split pheno,type (obesity) into test data
phenotypes_test4 <- df4[-train4,]

stopifnot(length(train4)==nrow(phenotypes_train4))
stopifnot(nrow(snps_train4)==nrow(phenotypes_train4))
stopifnot(nrow(snps_test4)==nrow(phenotypes_test4))


# Cross-validation of glmnet ----------------------------------------------


# Run penalized linear regression for obesity
X_train4 <- snps_train4
Y_train4 <- phenotypes_train4[,"Obesity.BMI"]
X_test4 <- snps_test4
Y_test4 <- phenotypes_test4[,"Obesity.BMI"]

set.seed(2)
model.EN.obesityBMI.Litter4 <- cv.glmnet(x=X_train4, y=Y_train4, family="gaussian", alpha=0.5, type.measure="deviance")

save(model.EN.obesityBMI.Litter4, file="~/documents/Msc_Project/obesityLitter4EN.RData")

pdf("~/documents/MSc_Project/litter4ENplot.pdf") 
plot(model.EN.obesityBMI.Litter4)
dev.off()

coef(model.EN.obesityBMI.Litter4, s = "lambda.min")

optimal.lambda4<-model.EN.obesityBMI.Litter4$lambda.min

# working on test data set (prediction)
EN.obesityBMI.pred.Litter4 <- predict(model.EN.obesityBMI.Litter4, newx = X_test4, s = optimal.lambda4,type="response")
EN.obesityBMI.pred.Litter4.min <- predict(model.EN.obesityBMI.Litter4, newx = X_test4, s = "lambda.min",type="response")
EN.obesityBMI.pred.Litter4.1se <- predict(model.EN.obesityBMI.Litter4, newx = X_test4, s = "lambda.1se",type="response")

#MSE
MSELitter4<-mean((EN.obesityBMI.pred.Litter4.min - Y_test4)^2)
save(MSELitter4, file="~/documents/Msc_Project/MSElitter4EN.RData")



# Calculating PRS ---------------------------------------------------------


# getting Beta's
plot(model.EN.obesityBMI.Litter4)
beta4 <- coef(model.EN.obesityBMI.Litter4,s="lambda.min") # NB. coef is giving the wrong values without s=
summary(beta4)


# converting s4 into datafram
beta4_selected<-data.frame(name = beta4@Dimnames[[1]][beta4@i + 1], coefficient = beta4@x)


# remove intercept from beta
beta4_selectedt=beta4_selected[-1,]
dim(beta4_selectedt)

# make sure it is a dataframe
b4<-as.data.frame(beta4_selectedt)

# some plots
pdf("~/documents/MSc_Project/histbetaselectedENlitter4.pdf") 
hist(b4$coefficient)
dev.off()

pdf("~/documents/MSc_Project/boxplotselectedvariableslitter4EN.pdf") 
boxplot(b4$coefficient, main="boxplot of beta coefficients of SNPs selected in litter 4")
dev.off()

# now select snps from mice.X with b
dim(mice.X4)

# transpose mice.X
mice.X4t<-t(mice.X4)
dim(mice.X4t)

# turn mice.Xt into a dataframe
mice.X4t<-data.frame(mice.X4t)

# select the snps and alleles from b in mice.X, aka make sure the same snps have been selected
selected_snps_mice4 <- mice.X4t[match(b4$name, rownames(mice.X4t)),]
dim(selected_snps_mice4)#we want matrix, 51rows by 1763mice

# check to see if b and selected_snps_mice are in the same order
all(b4$name==rownames(selected_snps_mice4)) #they are in same order

# sanity check
stopifnot((b4$name)==rownames(selected_snps_mice4))

dim(selected_snps_mice4)
dim(b4)

b4coefficients<-b4$coefficient
b4coefficients<-as.matrix(b4coefficients)

# sweep to calculate PRS scores
resultsl4<-sweep(selected_snps_mice4, MARGIN=1,b4coefficients, `*`)
dim (resultsl4)

# getting PRS score per individual mouse- summing each row
PRSscoresforlitter4EN<-colSums(resultsl4)
save(PRSscoresforlitter4EN, file="~/documents/Msc_Project/PRSscoresforlitter4EN.RData")
hist(PRSscoresforlitter4EN)


# Litter 5 ----------------------------------------------------------------


# Sorting by litter 5 -----------------------------------------------------

Litter5 <- subset.data.frame(mice.pheno, Litter=='5')

# Select the same mice (subject id ) from mice X, for Litter 
# make header the same
# mice.X <- data.frame("SUBJECT.NAME"=rownames(mice.X), mice.X)

# select litter subjects in miceX
mice.X5 <- mice.X[match(Litter5$SUBJECT.NAME, rownames(mice.X)),]
#View(mice.X5)
dim(mice.X5)

# keeping the variables for the phenotype (aka obesity.BMI)
df5 = Litter5[c("SUBJECT.NAME", "Obesity.BMI")]

# check the Y is the same order as the mice.X
all(df5$SUBJECT.NAME==rownames(mice.X5))

# # turn them as matrices and vectors
mice.X5 <- as.matrix(mice.X5)
typeof(mice.X5)

# Split training and testing data -----------------------------------------

# Split mice.X (Litter) (snps) into training data
set.seed(1)
train5 <- sample(1:nrow(mice.X5), 0.7*nrow(mice.X5), replace=FALSE)
snps_train5 <- mice.X5[train5,]

# Split pheno.type (obesity) into training data
phenotypes_train5 <- df5[train5,]

# split mice.X(snps) into test data
snps_test5 <- mice.X5[-train5,]

# split pheno,type (obesity) into test data
phenotypes_test5 <- df5[-train5,]


stopifnot(length(train5)==nrow(phenotypes_train5))
stopifnot(nrow(snps_train5)==nrow(phenotypes_train5))
stopifnot(nrow(snps_test5)==nrow(phenotypes_test5))


# Cross validating glmnet -------------------------------------------------


# Run penalized linear regression for obesity

X_train5 <- snps_train5
Y_train5 <- phenotypes_train5[,"Obesity.BMI"]
X_test5 <- snps_test5
Y_test5 <- phenotypes_test5[,"Obesity.BMI"]


set.seed(2)
model.EN.obesityBMI.Litter5 <- cv.glmnet(x=X_train5, y=Y_train5, family="gaussian", alpha=0.5, type.measure="deviance")


save(model.EN.obesityBMI.Litter5, file="~/documents/Msc_Project/obesityLitter5EN.RData")

pdf("~/documents/MSc_Project/litter5ENplot.pdf") 
plot(model.EN.obesityBMI.Litter5)
dev.off()

coef(model.EN.obesityBMI.Litter5, s = "lambda.min")

optimal.lambda5<-model.EN.obesityBMI.Litter5$lambda.min


# Working on test data set ------------------------------------------------

EN.obesityBMI.pred.Litter5 <- predict(model.EN.obesityBMI.Litter5, newx = X_test5, s = optimal.lambda5,type="response")
EN.obesityBMI.pred.Litter5.min <- predict(model.EN.obesityBMI.Litter5, newx = X_test5, s = "lambda.min",type="response")
EN.obesityBMI.pred.Litter5.1se <- predict(model.EN.obesityBMI.Litter5, newx = X_test5, s = "lambda.1se",type="response")

#MSE
MSELitter5<-mean((EN.obesityBMI.pred.Litter5.min - Y_test5)^2)
save(MSELitter5, file="~/documents/Msc_Project/MSElitter5EN.RData")


# Calculating PRS  ------------------------------------------------------------

# getting Beta's
plot(model.EN.obesityBMI.Litter5)
beta5 <- coef(model.EN.obesityBMI.Litter5,s="lambda.min") # NB. coef is giving the wrong values without s=
summary(beta5)


# converting s4 into datafram
beta5_selected<-data.frame(name = beta5@Dimnames[[1]][beta5@i + 1], coefficient = beta5@x)


# remove intercept from beta
beta5_selectedt=beta5_selected[-1,]
dim(beta5_selectedt)

# make sure it is a dataframe
b5<-as.data.frame(beta5_selectedt)

# some plots
pdf("~/documents/MSc_Project/histbetaselectedENlitter5.pdf") 
plot(b5$coefficient)
dev.off()

pdf("~/documents/MSc_Project/boxplotselectedvariableslitter5EN.pdf") 
boxplot(b5$coefficient, main="boxplot of beta coefficients of SNPs selected in litter 5")
dev.off()

# now select snps from mice.X with b
dim(mice.X5)

# transpose mice.X
mice.X5t<-t(mice.X5)
dim(mice.X5t)

# turn mice.Xt into a dataframe
mice.X5t<-data.frame(mice.X5t)

# select the snps and alleles from b in mice.X, aka make sure the same snps have been selected
selected_snps_mice5 <- mice.X5t[match(b5$name, rownames(mice.X5t)),]
dim(selected_snps_mice5)#no results as no snps selected


# check to see if b and selected_snps_mice are in the same order
all(b5$name==rownames(selected_snps_mice5)) #they are in same order


# sanity check
stopifnot((b5$name)==rownames(selected_snps_mice5))

dim(selected_snps_mice5)
dim(b5)

b5coefficients<-b5$coefficient
b5coefficients<-as.matrix(b5coefficients)

# sweep to calculate PRS scores
resultsl5<-sweep(selected_snps_mice5, MARGIN=1,b5coefficients, `*`)
dim (resultsl5)


# getting PRS score per individual mouse- summing each row
PRSscoresforlitter5EN<-colSums(resultsl5)
save(PRSscoresforlitter5EN, file="~/documents/Msc_Project/PRSscoresforlitter5EN.RData")
hist(PRSscoresforlitter5EN)

##elastic net failed on litter 5 as well!!!!!



# Litter 6 ----------------------------------------------------------------


# Sorting litter 6 --------------------------------------------------------


Litter6 <- subset.data.frame(mice.pheno, Litter=='6')

# Select the same mice (subject id ) from mice X, for Litter 6
# make header the same
# mice.X <- data.frame("SUBJECT.NAME"=rownames(mice.X), mice.X)

# select litter subjects in miceX
mice.X6 <- mice.X[match(Litter6$SUBJECT.NAME, rownames(mice.X)),]
#View(mice.X6)
dim(mice.X6)

# keeping the variables for the phenotype (aka obesity.BMI)
df6 = Litter6[c("SUBJECT.NAME", "Obesity.BMI")]

# check the Y is the same order as the mice.X
all(df6$SUBJECT.NAME==rownames(mice.X6))

# # turn them as matrices and vectors
mice.X6 <- as.matrix(mice.X6)
typeof(mice.X6)


# Split into training and testing data ------------------------------------


set.seed(1)
train6 <- sample(1:nrow(mice.X6), 0.7*nrow(mice.X6), replace=FALSE)
snps_train6 <- mice.X6[train6,]

# Split pheno.type (obesity) into training data
phenotypes_train6 <- df6[train6,]

# split mice.X(snps) into test data
snps_test6 <- mice.X6[-train6,]

# split pheno,type (obesity) into test data
phenotypes_test6 <- df6[-train6,]

stopifnot(length(train6)==nrow(phenotypes_train6))
stopifnot(nrow(snps_train6)==nrow(phenotypes_train6))
stopifnot(nrow(snps_test6)==nrow(phenotypes_test6))


# Cross-validating glmnet -------------------------------------------------

# Run penalized linear regression for obesity

X_train6 <- snps_train6
Y_train6 <- phenotypes_train6[,"Obesity.BMI"]
X_test6 <- snps_test6
Y_test6 <- phenotypes_test6[,"Obesity.BMI"]


set.seed(2)
model.EN.obesityBMI.Litter6 <- cv.glmnet(x=X_train6, y=Y_train6, family="gaussian", alpha=0.5, type.measure="deviance")


save(model.EN.obesityBMI.Litter6, file="~/documents/Msc_Project/obesityLitter6lasso.RData")

pdf("~/documents/MSc_Project/litter6ENplot.pdf") 
plot(model.EN.obesityBMI.Litter6)
dev.off()

coef(model.EN.obesityBMI.Litter6, s = "lambda.min")

optimal.lambda6<-model.EN.obesityBMI.Litter6$lambda.min


# Working on test data prediction -----------------------------------------

EN.obesityBMI.pred.Litter6 <- predict(model.EN.obesityBMI.Litter6, newx = X_test6, s = optimal.lambda6,type="response")
EN.obesityBMI.pred.Litter6.min <- predict(model.EN.obesityBMI.Litter6, newx = X_test6, s = "lambda.min",type="response")
EN.obesityBMI.pred.Litter6.1se <- predict(model.EN.obesityBMI.Litter6, newx = X_test6, s = "lambda.1se",type="response")

#MSE
MSELitter6<-mean((EN.obesityBMI.pred.Litter6.min - Y_test6)^2)
save(MSELitter6, file="~/documents/Msc_Project/MSElitter6EN.RData")


# Calculating PRS ---------------------------------------------------------


# getting Beta's
plot(model.EN.obesityBMI.Litter6)
beta6 <- coef(model.EN.obesityBMI.Litter6,s="lambda.min") # NB. coef is giving the wrong values without s=
summary(beta6)

# converting s4 into datafram
beta6_selected<-data.frame(name = beta6@Dimnames[[1]][beta6@i + 1], coefficient = beta6@x)

# remove intercept from beta
beta6_selectedt=beta6_selected[-1,]
dim(beta6_selectedt)

# make sure it is a dataframe
b6<-as.data.frame(beta6_selectedt)

# some plots
pdf("~/documents/MSc_Project/histbetaselectedENlitter6.pdf") 
hist(b6$coefficient)
dev.off()

pdf("~/documents/MSc_Project/boxplotselectedvariableslitter6EN.pdf") 
boxplot(b6$coefficient, main="boxplot of beta coefficients of SNPs selected in litter 6")
dev.off()

# now select snps from mice.X with b
dim(mice.X6)

# transpose mice.X
mice.X6t<-t(mice.X6)
dim(mice.X6t)

# turn mice.Xt into a dataframe
mice.X6t<-data.frame(mice.X6t)

# select the snps and alleles from b in mice.X, aka make sure the same snps have been selected
selected_snps_mice6 <- mice.X6t[match(b6$name, rownames(mice.X6t)),]
dim(selected_snps_mice6)

# check to see if b and selected_snps_mice are in the same order
all(b6$name==rownames(selected_snps_mice6)) #they are in same order

# sanity check
stopifnot((b6$name)==rownames(selected_snps_mice6))

dim(selected_snps_mice6)
dim(b6)

b6coefficients<-b6$coefficient
b6coefficients<-as.matrix(b6coefficients)

# sweep to calculate PRS scores
resultsl6<-sweep(selected_snps_mice6, MARGIN=1,b6coefficients, `*`)
dim (resultsl6)

# getting PRS score per individual mouse- summing each row
PRSscoresforlitter6EN<-colSums(resultsl6)
save(PRSscoresforlitter6EN, file="~/documents/Msc_Project/PRSscoresforlitter6EN.RData")
hist(PRSscoresforlitter6EN)

##failed on litter 6 EN


# Litter 7 ----------------------------------------------------------------

# Sorting Litter 7 --------------------------------------------------------

Litter7 <- subset.data.frame(mice.pheno, Litter=='7')

# Select the same mice (subject id ) from mice X, for Litter 7
# make header the same
# mice.X <- data.frame("SUBJECT.NAME"=rownames(mice.X), mice.X)

# select litter subjects in miceX
mice.X7 <- mice.X[match(Litter7$SUBJECT.NAME, rownames(mice.X)),]
#View(mice.X7)
dim(mice.X7)

# keeping the variables for the phenotype (aka obesity.BMI)
df7 = Litter7[c("SUBJECT.NAME", "Obesity.BMI")]

# check the Y is the same order as the mice.X
all(df7$SUBJECT.NAME==rownames(mice.X7))

# # turn them as matrices and vectors
mice.X7 <- as.matrix(mice.X7)
typeof(mice.X7)


# Split into training and testing data ------------------------------------

# Split mice.X (Litter1) (snps) into training data
set.seed(1)
train7 <- sample(1:nrow(mice.X7), 0.7*nrow(mice.X7), replace=FALSE)
snps_train7 <- mice.X7[train7,]

# Split pheno.type (obesity) into training data
phenotypes_train7 <- df7[train7,]

# split mice.X(snps) into test data
snps_test7 <- mice.X7[-train7,]

# split pheno,type (obesity) into test data
phenotypes_test7 <- df7[-train7,]

stopifnot(length(train7)==nrow(phenotypes_train7))
stopifnot(nrow(snps_train7)==nrow(phenotypes_train7))
stopifnot(nrow(snps_test7)==nrow(phenotypes_test7))



# Cross validating glmnet -------------------------------------------------

X_train7 <- snps_train7
Y_train7 <- phenotypes_train7[,"Obesity.BMI"]
X_test7 <- snps_test7
Y_test7 <- phenotypes_test7[,"Obesity.BMI"]


set.seed(2)
model.EN.obesityBMI.Litter7 <- cv.glmnet(x=X_train7, y=Y_train7, family="gaussian", alpha=0.5, type.measure="deviance")


save(model.EN.obesityBMI.Litter7, file="~/documents/Msc_Project/obesityLitter7EN.RData")

pdf("~/documents/MSc_Project/litter7ENplot.pdf") 
plot(model.EN.obesityBMI.Litter7)
dev.off()


coef(model.EN.obesityBMI.Litter7, s = "lambda.min")

optimal.lambda7<-model.EN.obesityBMI.Litter7$lambda.min


# Working on test data ----------------------------------------------------



EN.obesityBMI.pred.Litter7 <- predict(model.EN.obesityBMI.Litter7, newx = X_test7, s = optimal.lambda7,type="response")
EN.obesityBMI.pred.Litter7.min <- predict(model.EN.obesityBMI.Litter7, newx = X_test7, s = "lambda.min",type="response")
EN.obesityBMI.pred.Litter7.1se <- predict(model.EN.obesityBMI.Litter7, newx = X_test7, s = "lambda.1se",type="response")

#MSE
MSELitter7<-mean((EN.obesityBMI.pred.Litter7.min - Y_test7)^2)
save(MSELitter7, file="~/documents/Msc_Project/MSElitter7EN.RData")


# Calculating PRS  --------------------------------------------------------


# getting Beta's
plot(model.EN.obesityBMI.Litter7)
beta7 <- coef(model.EN.obesityBMI.Litter7,s="lambda.min") # NB. coef is giving the wrong values without s=
summary(beta7)

# converting s4 into dataframe
beta7_selected<-data.frame(name = beta7@Dimnames[[1]][beta7@i + 1], coefficient = beta7@x)

# remove intercept from beta
beta7_selectedt=beta7_selected[-1,]
dim(beta7_selectedt)

# make sure it is a dataframe
b7<-as.data.frame(beta7_selectedt)

# some plots
pdf("~/documents/MSc_Project/histbetaselectedENlitter7.pdf") 
hist(b7$coefficient)
dev.off()

pdf("~/documents/MSc_Project/boxplotselectedvariableslitter7EN.pdf") 
boxplot(b7$coefficient, main="boxplot of beta coefficients of SNPs selected in litter 7")
dev.off()


# transpose mice.X
mice.X7t<-t(mice.X7)
dim(mice.X7t)

# turn mice.Xt into a dataframe
mice.X7t<-data.frame(mice.X7t)

# select the snps and alleles from b in mice.X, aka make sure the same snps have been selected
selected_snps_mice7 <- mice.X7t[match(b7$name, rownames(mice.X7t)),]
dim(selected_snps_mice7)

# check to see if b and selected_snps_mice are in the same order
all(b7$name==rownames(selected_snps_mice7)) #they are in same order

# sanity check
stopifnot((b7$name)==rownames(selected_snps_mice7))

dim(selected_snps_mice7)
dim(b7)

b7coefficients<-b7$coefficient
b7coefficients<-as.matrix(b7coefficients)

# sweep to calculate PRS scores
resultsl7<-sweep(selected_snps_mice7, MARGIN=1,b7coefficients, `*`)
dim (resultsl7)

# getting PRS score per individual mouse- summing each row
PRSscoresforlitter7EN<-colSums(resultsl7)
save(PRSscoresforlitter7EN, file="~/documents/Msc_Project/PRSscoresforlitter7EN.RData")
hist(PRSscoresforlitter7EN)

#####failed on litter 7 - EN (lambda min)


# Litter 8 ----------------------------------------------------------------

# Sorting litter 8 --------------------------------------------------------

Litter8 <- subset.data.frame(mice.pheno, Litter=='8')



# Select the same mice (subject id ) from mice X, for Litter 8
# make header the same
# mice.X <- data.frame("SUBJECT.NAME"=rownames(mice.X), mice.X)

# select litter subjects in miceX
mice.X8 <- mice.X[match(Litter8$SUBJECT.NAME, rownames(mice.X)),]
#View(mice.X7)
dim(mice.X8)

# keeping the variables for the phenotype (aka obesity.BMI)
df8 = Litter8[c("SUBJECT.NAME", "Obesity.BMI")]


# check the Y is the same order as the mice.X
all(df8$SUBJECT.NAME==rownames(mice.X8))

# # turn them as matrices and vectors
mice.X8 <- as.matrix(mice.X8)
typeof(mice.X8)



# Split into training and testing data ------------------------------------


# Split mice.X (Litter1) (snps) into training data
set.seed(1)
train8 <- sample(1:nrow(mice.X8), 0.7*nrow(mice.X8), replace=FALSE)
snps_train8 <- mice.X8[train8,]

# Split pheno.type (obesity) into training data
phenotypes_train8 <- df8[train8,]

# split mice.X(snps) into test data
snps_test8 <- mice.X8[-train8,]

# split pheno,type (obesity) into test data
phenotypes_test8 <- df8[-train8,]


stopifnot(length(train8)==nrow(phenotypes_train8))
stopifnot(nrow(snps_train8)==nrow(phenotypes_train8))
stopifnot(nrow(snps_test8)==nrow(phenotypes_test8))


# Cross-validation of glmnet ----------------------------------------------

# Run penalized linear regression for obesity

X_train8 <- snps_train8
Y_train8 <- phenotypes_train8[,"Obesity.BMI"]
X_test8 <- snps_test8
Y_test8 <- phenotypes_test8[,"Obesity.BMI"]


set.seed(2)
model.EN.obesityBMI.Litter8 <- cv.glmnet(x=X_train8, y=Y_train8, family="gaussian", alpha=0.5, type.measure="deviance")


save(model.EN.obesityBMI.Litter8, file="~/documents/Msc_Project/obesityLitter8EN.RData")

pdf("~/documents/MSc_Project/litter8ENplot.pdf") 
plot(model.EN.obesityBMI.Litter8)
dev.off()


coef(model.EN.obesityBMI.Litter8, s = "lambda.min")

optimal.lambda8<-model.EN.obesityBMI.Litter8$lambda.min


# Working on test data ----------------------------------------------------

EN.obesityBMI.pred.Litter8 <- predict(model.EN.obesityBMI.Litter8, newx = X_test8, s = optimal.lambda8,type="response")
EN.obesityBMI.pred.Litter8.min <- predict(model.EN.obesityBMI.Litter8, newx = X_test8, s = "lambda.min",type="response")
EN.obesityBMI.pred.Litter8.1se <- predict(model.EN.obesityBMI.Litter8, newx = X_test8, s = "lambda.1se",type="response")

#MSE
MSELitter8<-mean((EN.obesityBMI.pred.Litter8.min - Y_test8)^2)
save(MSELitter8, file="~/documents/Msc_Project/MSElitter8EN.RData")


# Calculating PRS  --------------------------------------------------------

# getting Beta's
plot(model.EN.obesityBMI.Litter8)
beta8 <- coef(model.EN.obesityBMI.Litter8,s="lambda.min") # NB. coef is giving the wrong values without s=
summary(beta8)

# converting s4 into datafram
beta8_selected<-data.frame(name = beta8@Dimnames[[1]][beta8@i + 1], coefficient = beta8@x)

# remove intercept from beta
beta8_selectedt=beta8_selected[-1,]
dim(beta8_selectedt)

# make sure it is a dataframe
b8<-as.data.frame(beta8_selectedt)

# some plots
pdf("~/documents/MSc_Project/histbetaselectedENlitter8.pdf") 
hist(b8$coefficient)
dev.off()

pdf("~/documents/MSc_Project/boxplotselectedvariableslitter8EN.pdf") 
boxplot(b8$coefficient, main="boxplot of beta coefficients of SNPs selected in litter 8")
dev.off()

# now select snps from mice.X with b
dim(mice.X8)

# transpose mice.X
mice.X8t<-t(mice.X8)
dim(mice.X8t)

# turn mice.Xt into a dataframe
mice.X8t<-data.frame(mice.X8t)

# select the snps and alleles from b in mice.X, aka make sure the same snps have been selected
selected_snps_mice8 <- mice.X8t[match(b8$name, rownames(mice.X8t)),]
dim(selected_snps_mice8)

# check to see if b and selected_snps_mice are in the same order
all(b8$name==rownames(selected_snps_mice8)) #they are in same order

# sanity check
stopifnot((b8$name)==rownames(selected_snps_mice8))

dim(selected_snps_mice8)
dim(b8)

b8coefficients<-b8$coefficient
b8coefficients<-as.matrix(b8coefficients)

# sweep to calculate PRS scores
resultsl8<-sweep(selected_snps_mice8, MARGIN=1,b8coefficients, `*`)
dim (resultsl8)

# getting PRS score per individual mouse- summing each row
PRSscoresforlitter8EN<-colSums(resultsl8)
save(PRSscoresforlitter8EN, file="~/documents/Msc_Project/PRSscoresforlitter8EN.RData")
hist(PRSscoresforlitter8EN)


# All litters -------------------------------------------------------------


# Sorting all litters -----------------------------------------------------

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


save(model.EN.obesityBMI.Litterx, file="~/documents/Msc_Project/obesityonlyalllitersxEN.RData")

pdf("~/documents/MSc_Project/alllittersENlot.pdf") 
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
save(MSELitterx, file="~/documents/Msc_Project/MSElitterxEN.RData")



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
pdf("~/documents/MSc_Project/histbetaselectedENlitterx.pdf") 
hist(bx$coefficient)
dev.off()

pdf("~/documents/MSc_Project/boxplotselectedvariableslitterxEN.pdf") 
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
save(PRSscoresforalllittersEN, file="~/documents/Msc_Project/PRSscoresforalllittersEN.RData")
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



png("~/documents/MSc_Project/DistributionofPRSscoresEN.png") 
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




pdf("~/documents/MSc_Project/DistributionofPRS-EN.pdf") 
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

## Example ##
# plot(1:10)
# makeFootnote(footnote)


# MSE plots 

# 4 figures arranged in 2 rows and 2 columns


png("~/documents/MSc_project/plotsforEN1.png") 
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


pdf("~/documents/MSc_project/plotsforElasticNetfinal.pdf") 
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

snps_selected_EN<- read.csv('~/Documents/Msc_Project/selected_snps_dataset.csv')

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

venn.diagram(mysets, filename = "~/Documents/MSc_Project/Venn_EN.png",
             fill = colors, category.names = names, imagetype = "png",
             cat.just = list(c(0.5, 1), c(-0.5, -5), c(0, 0),
                             c(0.5, 0), c(1, -4)), main="Selected SNPs-Elastic Net"
)

dev.off()


# Pheatmap for Elastic Net ------------------------------------------------

# creating pheatmap

pdf("~/documents/MSc_project/pheatmapsforsnpsElasticNet.pdf") 

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
pdf("~/documents/MSc_project/CVEPCA.pdf") 
plot(var.explained, main="Cumulative Variance Explained by Principal Components")
dev.off()

### Exploring PCA by litter
pca.predictor.df <-data.frame(pca.out$x)[1:100]
pca.predictor.df <- cbind(as.factor(mice.pheno$Litter), pca.predictor.df)
names(pca.predictor.df)[1] <- "Litter"

#Plot how specific litters vary from the rest of the litters
ggplot(data = pca.predictor.df, aes(y = PC1, x = PC2, col = Litter, fill = Litter)) + geom_point(alpha = 0.5) + theme_clean() #they should all be factors

# plot all litters on PCA graph with grouping frames
pdf("~/documents/MSc_project/plotsforPCA.pdf") 
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
write.csv(tab1Mat, file = "~/documents/MSc_project/myTable1.csv")

