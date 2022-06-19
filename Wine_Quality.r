#---
#title: "Modeling and Predicting Red Wind Quality"
#author: "Shayna Gaulden and Tien Nguyen"
#date: "11/29/2021"
#---

### Loading the Data
WineQuality <- read.table("winequality-red.csv",header=TRUE,sep=",")
summary(WineQuality)

################################################################################
# Splitting data into train and test ###########################################
################################################################################
library(sampling)
N = nrow(WineQuality)
N # definitely enough data points to split the data!
n = round(.20*N) # we will use 20% test data and 80% train data
set.seed(10) # making the results the same
tr = srswor(N-n,N)

# train data
x1.tr = WineQuality$fixed.acidity[tr==1]
x2.tr = WineQuality$volatile.acidity[tr==1]
x3.tr = WineQuality$citric.acid[tr==1]
x4.tr = WineQuality$residual.sugar[tr==1]
x5.tr = WineQuality$chlorides[tr==1]
x6.tr = WineQuality$free.sulfur.dioxide[tr==1]
x7.tr = WineQuality$density[tr==1]
x8.tr = WineQuality$pH[tr==1]
x9.tr = WineQuality$sulphates[tr==1]
x10.tr = WineQuality$alcohol[tr==1]
y.tr = WineQuality$quality[tr==1]

# put it into an easy to use data frame
Wine_df_train = data.frame(x1.tr,x2.tr,x3.tr,x4.tr,x5.tr,x6.tr,x7.tr,x8.tr,
                           x9.tr,x10.tr,y.tr)
colnames(Wine_df_train) = colnames(WineQuality)[-7]

# test data
x1.te = WineQuality$fixed.acidity[tr==0]
x2.te = WineQuality$volatile.acidity[tr==0]
x3.te = WineQuality$citric.acid[tr==0]
x4.te = WineQuality$residual.sugar[tr==0]
x5.te = WineQuality$chlorides[tr==0]
x6.te = WineQuality$free.sulfur.dioxide[tr==0]
x7.te = WineQuality$density[tr==0]
x8.te = WineQuality$pH[tr==0]
x9.te = WineQuality$sulphates[tr==0]
x10.te = WineQuality$alcohol[tr==0]
y.te = WineQuality$quality[tr==0]

# put it into an easy to use data frame
Wine_df_test = data.frame(x1.te,x2.te,x3.te,x4.te,x5.te,x6.te,x7.te,x8.te,
                          x9.te,x10.te,y.te)
colnames(Wine_df_test) = colnames(WineQuality)[-7]

################################################################################
# Fitting Model 1 ##############################################################
################################################################################
# fitting the model
res <- lm(y.tr~x1.tr+x2.tr+x3.tr+x4.tr+x5.tr+x6.tr+x7.tr+x8.tr+x9.tr+x10.tr)
res$coefficients

summary(res)
anova(res)

## Checking the Model Assumptions
#We will standardize the residuals to the R-student residuals.
stud = rstudent(res)
### Residuals vs. Time Order
par(mfrow=c(1,1))
plot(stud,ylab="Student Residuals",main = "Residuals vs. Time Order")
line = lm(formula = stud ~ seq(length(stud))) # fitting regression line through residuals
abline(0,0,col="red")
abline(coef=line$coefficients,col="blue") # showing line on the plot

# Residuals in this plot are fairly spread out confirming that variance of the
# model is constant and the mean of the variance looks to be approximately 0.

### Residuals vs. Fitted Values
plot(stud,res$fitted.values)
#In this plot the residuals look spread out as well and there are no assumptions
# that appear to be violated. We can see a possible outlier in the top left hand
# corner. The lines are seen because the response variable is an ordinal scale.

### QQ-plot
qqnorm(stud)
qqline(stud)
#Data drifts apart from the line at the end points but overall looks close to a
# normal distribution.

### Residuals vs. all Regressor Variables in the Model
# large plots
plot(x1.tr,stud) # looks fine
plot(x2.tr,stud) # possible outlier
plot(x3.tr,stud) # possible outlier
plot(x4.tr,stud) # fan shape / not looking great
plot(x5.tr,stud) # not good
plot(x6.tr,stud) # fan shape
plot(x7.tr,stud) # looks good I think
plot(x8.tr,stud) # looks good
plot(x9.tr,stud) # okay maybe some outliers/fan shape
plot(x10.tr,stud) # looks okay
# try different transformations
# little plots all together
library(car) # needed for residualPlots() function
residualPlots(res, layout = c(3,4), type = "rstudent")

################################################################################
# Trying Transformations on Model 1 ############################################
################################################################################
# trying box cox method for finding optimal transformation
library(MASS)
bc = boxcox(res)
max_lambda = bc$x[order(bc$y,decreasing=TRUE)[1]]
max_lambda # this is basically 1 so I would probably not transform y using y^{lambda}

# trying to do transformation on y with log
y_log = log(y.tr)
res_y_log = lm( y_log~ x1.tr + x2.tr + x3.tr + x4.tr + x5.tr + x6.tr + x7.tr
                + x8.tr + x9.tr + x10.tr)
summary(res_y_log)
plot(res_y_log$residuals)
residualPlots(res_y_log, layout = c(3,4), type = "rstudent")
# transformation with log y does not look good so we will not use it

# checking for variables that may be near linear combinations
## centering
x=sweep(WineQuality[,c(-7,-12)],2, FUN='-',apply(WineQuality[,c(-7,-12)],2,mean))
VIF=diag(solve(cor(x)))
VIF # none are over 15

### Transformations
# trying log transformations on variables that seemed to violate model assumptions
x4_log = log(x4.tr)
x5_log = log(x5.tr)
x6_log = log(x6.tr)
x9_log = log(x9.tr)
res_log <- lm(y.tr~x1.tr+x2.tr+x3.tr+x4_log+x5_log+x6_log+x7.tr+x8.tr+x9_log+x10.tr)
summary(res_log)

plot(res_log$residuals)
residualPlots(res_log, layout = c(3,4), type = "rstudent")
# from the plots it appears the log transformations help to not violate model assumptions

################################################################################
# Method A #####################################################################
################################################################################
### Polynomial Full Model
x1.tr.2 = x1.tr^2
x2.tr.2 = x2.tr^2
x3.tr.2 = x3.tr^2
x4.tr.2 = x4_log^2
x5.tr.2 = x5_log^2
x6.tr.2 = x6_log^2
x7.tr.2 = x7.tr^2
x8.tr.2 = x8.tr^2
x9.tr.2 = x9_log^2
x10.tr.2 = x10.tr^2

# creating x matrix
x = cbind(x1.tr,x1.tr.2,x2.tr,x2.tr.2,x3.tr,x3.tr.2,x4_log,x4.tr.2,x5_log,
          x5.tr.2,x6_log,x6.tr.2,x7.tr,x7.tr.2,x8.tr,x8.tr.2,x9_log,x9.tr.2,
          x10.tr,x10.tr.2)
# checking correlation
# too many to check by eye so this function will sort out which ones have a
# correlation higher than 0.7
correlated = c()
for (i in seq(ncol(x))) {
  check = which(abs(cor(x)[i,])>0.7)
  for (j in check) {
    if (j != i) { # avoiding diagonals
      correlated = append(correlated, paste(i,j))
      # pasting on the index of correlated variables
    }
  }
}
correlated # going to need to center

# Performing the centering
cx = x
for (i in seq(to=ncol(x),from=1,by=2)) { # going by un-squared terms first
  cx[,i] = x[,i] - mean(x[,i])
}
# now adding back in the centered squared terms
for (i in seq(to=ncol(x),from=2,by=2)) { # going by un-squared terms first
  cx[,i] = cx[,i-1]^2
}
### Checking correlation again
correlated = c()
for (i in seq(ncol(x))) {
  check = which(abs(cor(cx)[i,])>0.7)
  for (j in check) {
    if (j != i) {
      correlated = append(correlated, paste(i,j))
    }
  }
}
correlated # now is empty so centering helped!

### Variable Selection 

# Method A
# Forward Selection added by highest F
res0 <- lm(y.tr~1)
all <- lm(y.tr~cx[,1]+cx[,2]+cx[,3]+cx[,4]+cx[,5]+cx[,6]+cx[,7]
          +cx[,8]+cx[,9]+cx[,10]+cx[,11]+cx[,12]+cx[,13]+cx[,14]+cx[,15]
          +cx[,16]+cx[,17]+cx[,18]+cx[,19]+cx[,20])
# step function is based on AIC value but we would like to stop adding or remove
# a variable if the p-value for the F test is above 0.05 so we modify the k
# parameter
k = qchisq(0.05,1,lower.tail=FALSE)
both <- step(res0, direction='forward', scope=formula(all), trace=0,k=k)
both

# Backwards Selection dropped by lowest F
both <- step(all, direction='backward', scope=formula(all), trace=0,k=k)
both

both <- step(res0, direction='both', scope=formula(all), trace=0,k=k)
both
# forwards, backwards, and step selection gave the same model
res_log_poly <- lm(y.tr ~ cx[,3]+cx[,5]+cx[,9]+cx[,15]+cx[,17]+cx[,18]+cx[,19])
# Key for variables that made the cut
# cx[,2] is fixed.acidity ^ 2
# cx[,3] is volatile.acidity
# cx[,5] is citric acidity
# cx[,9] is log(chlorides)
# cx[,15] is ph
# cx[,17] is log(sulphates)
# cx[,18] is log(sulphates) ^2 
# cx[,19] is alcohol

summary(res_log_poly)
anova(res_log_poly)

plot(res_log_poly$residuals)
residualPlots(res_log_poly, layout = c(3,3), type = "rstudent")
# everything looks better there are definitely some outliers
model_a <- res_log_poly
# need to fit with different variable names for prediction later on
x2c = cx[,3]
x3c = cx[,5]
x5c= cx[,9]
x8c = cx[,15]
x9c = cx[,17]
x9csq = cx[,18]
x10c = cx[,19]
model_a <- lm(y.tr~x2c+x3c+x5c+x8c+x9c+x9csq+x10c)

################################################################################
# Method B #####################################################################
################################################################################
# Method B
df_Wine_log_train = data.frame(x1.tr, x2.tr, x3.tr, log(x4.tr),log(x5.tr),
                               log(x6.tr),x7.tr,x8.tr,log(x9.tr),x10.tr,y.tr)
colnames(df_Wine_log_train) = c("x1", "x2", "x3", "x4", "x5", "x6", "x7","x8",
                                "x9", "x10", "y")
#define intercept-only model
intercept_only <- lm(y ~ 1, data= df_Wine_log_train)
#define model with all predictors
all <- lm(y ~ ., data= df_Wine_log_train)
#perform backward stepwise regression
both <- step(intercept_only, direction='both', scope=formula(all), trace=0,k=k)
both # this is the stepwise regression model
both$coefficients
summary(both)

# trying a polynomial model with previously selected variables only
res_stepwise_poly = lm(y ~ poly(x2,2) + poly(x5,2)
                       + poly(x8,2) + poly(x9,2) + poly(x10,2),
                       data = df_Wine_log_train )
summary(res_stepwise_poly)
anova(res_stepwise_poly)

# based on summary only using poly() function on the variables that were significant
res_stepwise_poly_reduce = lm(y ~ x2+ x5 + poly(x8,2) + poly(x9,2) + x10,
                              data = df_Wine_log_train)
summary(res_stepwise_poly_reduce)
anova(res_stepwise_poly_reduce)

plot(res_stepwise_poly_reduce$residuals)
residualPlots(res_stepwise_poly_reduce, layout = c(3,3), type = "rstudent")

model_b <- res_stepwise_poly_reduce 

################################################################################
# Outliers #####################################################################
################################################################################
# looking for leverage points
n = length(y.tr)
k=7
pivot = 2*(k+1)/n
# checking model a
# x is the variables used in the model
x = cbind(cx[,3],cx[,5],cx[,9],cx[,15],cx[,17],cx[,18],cx[,19])
H = x%*% solve(t(x)%*%x)%*%t(x)

influence.pts = which(pivot < diag(H))
# this is potentially many influential/leverage points

length(influence.pts) # this is so many

# this for loop will print out some summary information on the model with the
# individual potential outlier removed
for (i in influence.pts) {
  t.res= lm(y.tr[-i]~cx[,3][-i]+cx[,5][-i]+cx[,9][-i]+cx[,15][-i]
            +cx[,17][-i]+cx[,18][-i]+cx[,19][-i])
  sumry = summary(t.res)
  anv = anova(t.res)
  print(paste("without point ",i,"adjrsq:", sumry$adj.r.squared,"rsqr:",
              sumry$r.squared,"MSRES",anv$'Mean Sq'[9]))
  print(t.res$coefficients)
  print("-------------------------------------------------------------------")
} # 

# a model with all influence points removed
res.none <- lm(y.tr[-influence.pts] ~ +cx[,3][-influence.pts]
               +cx[,5][-influence.pts]+cx[,9][-influence.pts]
               +cx[,15][-influence.pts]+cx[,17][-influence.pts]
               +cx[,18][-influence.pts]+cx[,19][-influence.pts])
summary(res.none) # decreases r square
anova(res.none) # decreases msres
# remove both and compare msres rsqrd and coeff estimations to see if model is
# better without

residualPlots(res.none, layout = c(3,3), type = "rstudent")

# going to remove only if square goes up and msres goes down
msres_compare = anova(model_a)$'Mean Sq'[8]
rsqr_compare = summary(model_a)$r.squared
points_to_remove = c() # empty list to put the points in
for (i in influence.pts) {
  t.res= lm(y.tr[-i]~cx[,3][-i]+cx[,5][-i]+cx[,9][-i]+cx[,15][-i]
            +cx[,17][-i]+cx[,18][-i]+cx[,19][-i])
  msres_dif = msres_compare - anova(t.res)$'Mean Sq'[8]
  # we want msres dif to be negative meaning msres has gone down
  rsqr_dif = rsqr_compare - summary(t.res)$r.squared
  # we want rsqr_dif to be positive meaning rsqr has gone up
  if (msres_dif<0 || rsqr_dif>0){
    points_to_remove = append(points_to_remove,i)
  }
} 

# model with only the influence points removed that when removed by themselves
# made the rsqr and msres decrease
points_to_remove
res.some <- lm(y.tr[-points_to_remove] ~ +cx[,3][-points_to_remove]
               +cx[,5][-points_to_remove]+cx[,9][-points_to_remove]
               +cx[,15][-points_to_remove]+cx[,17][-points_to_remove]
               +cx[,18][-points_to_remove]+cx[,19][-points_to_remove])
summary(res.some) # decreases r square by so much!!!
anova(res.some)
# possible cluster of differences the graphs look better but maybe it is best
# not to remove the outliers
################################################################################
# Model Testing ################################################################
################################################################################
# Testing model b

# calculate h_hat on model b
y_hat_test_b = predict(model_b, newdata= data.frame(x2=x2.te, x5=log(x5.te),
                                                    x8=x8.te, x9=log(x9.te),
                                                    x10=x10.te))
range(y_hat_test_b) # checking for a reasonable range

# New R^2
R2_test_b =1- (sum((y.te-y_hat_test_b)^2)/sum((y.te-mean(y.te))^2))
R2_test_b

# average square prediction error
avgsqpred_b = sum((y.te-y_hat_test_b)^2)/length(y.te)
avgsqpred_b

# confusion matrix
y_hat_test_b_round = round(y_hat_test_b)
library(caret)
expected_value_test_b <- factor(y.te)
predicted_value_test_b <- factor(y_hat_test_b_round)
confusionMatrix(data= predicted_value_test_b, reference = expected_value_test_b)

# Testing model a
# need to transform testing data so that it fits the model
x1.a = x2.te-mean(x2.te)# this is cx[,3] which is volatile acidity centered
x2.a = x3.te-mean(x3.te)# citric acidity
x3.a = log(x5.te)-mean(log(x5.te))
x4.a = x8.te - mean(x8.te)
x5.a = log(x9.te)-mean(log(x9.te))
x6.a =  (log(x9.te)-mean(log(x9.te)))^2
x7.a = x10.te - mean(x10.te)
df_Wine_log_test_a = data.frame(x2c=x1.a, x3c=x2.a,x5c=x3.a,x8c=x4.a,x9c=x5.a,
                                x9csq=x6.a,x10c=x7.a)

# calculate h_hat on model a
y_hat_test_a = predict(model_a, newdata=data.frame(x2c=x1.a, x3c=x2.a,x5c=x3.a,
                                                   x8c=x4.a,x9c=x5.a,x9csq=x6.a,
                                                   x10c=x7.a))
range(y_hat_test_a) # reasonable

# R sqr on new values
R2_test_a =1- (sum((y.te-y_hat_test_a)^2)/sum((y.te-mean(y.te))^2))
R2_test_a

# average square prediction error
avgsqpred_a = sum((y.te-y_hat_test_a)^2)/length(y.te)
avgsqpred_a

# confusion matrix
y_hat_test_a_round = round(y_hat_test_a)
library(caret)
expected_value_test_a <- factor(y.te)
predicted_value_test_a <- factor(y_hat_test_a_round)
confusionMatrix(data= predicted_value_test_a, reference = expected_value_test_a)
