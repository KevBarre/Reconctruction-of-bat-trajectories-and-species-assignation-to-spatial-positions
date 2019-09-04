rm(list=ls())

library(data.table) # version 1.12.2
library(glmmTMB) # version 0.2.3
library(sjstats) # version 0.17.5
library(DHARMa) # version 0.2.4
library(piecewiseSEM) # version 2.0.2
library(MuMIn) # version 1.43.6
library(lme4) # version 1.1-21

setwd("./")
data <- fread("./dataFinal_ProbForest.csv")

##########################################################################################################################
##########################################################################################################################
#               Probability to be inside forest
##########################################################################################################################
##########################################################################################################################

EdgeSideBinom = rep(999, nrow(data))
data = cbind(data, EdgeSideBinom)
data$EdgeSideBinom[which(data$EdgeSide == "Forest")] = "1"
data$EdgeSideBinom[which(data$EdgeSide == "Open")] = "0"
nrow(subset(data, data$EdgeSideBinom == 1))
nrow(subset(data, data$EdgeSideBinom == 0)) 

data$Precision = 1/data$Imp^2

vif.mer <- function (fit) {
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
} # to test for collinearity

# PIPISTRELLES ###########################################################################################################
dataPip = subset(data, data$SpeciesGroup == "Pip")
mPipHab1 <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mPipHab1) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataPip, weights = Precision)
AIC(m)
mPipHab1 <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
AIC(mPipHab1)
summary(mPipHab1)
r2(mPipHab1)
dataPip <- within(dataPip, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mPipHab1)

# Indivual R2
mPipHab1 <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataPip, weights = Precision)
r2(mPipHab1)
mPipHab1 <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
r2(mPipHab1)
res = simulateResiduals(mPipHab1)
plot(res, rank = T) # check for residuals

# PLEC/MYO ###########################################################################################################
dataPle = subset(data, data$SpeciesGroup == "MyoPlec")
mPleHab1 <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mPleHab1) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataPle, weights = Precision)
AIC(m)
mPleHab1 <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPle, weights = Precision)
AIC(mPleHab1)
summary(mPleHab1)
r2(mPleHab1)
dataPle <- within(dataPle, spetrum <- relevel(as.factor(spetrum), ref = "white"))
summary(mPleHab1) 

# Indivual R2
dataPle <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataPle, weights = Precision)
r2(dataPle)
dataPle <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataPle, weights = Precision)
r2(dataPle)
res = simulateResiduals(mPleHab1)
plot(res, rank = T) # check for residuals

# EPT/NYC ###########################################################################################################
dataEpt = subset(data, data$SpeciesGroup == "EptNyc")
mEptHab1 <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mEptHab1) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataEpt, weights = Precision)
AIC(m)
mEptHab1 <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataEpt, weights = Precision)
AIC(mEptHab1)
summary(mEptHab1)
r2(mEptHab1)
dataEpt <- within(dataEpt, spetrum <- relevel(as.factor(spetrum), ref = "white"))
summary(mEptHab1)

# Indivual R2
mEptHab1 <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataEpt, weights = Precision)
r2(mEptHab1)
mEptHab1 <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataEpt, weights = Precision)
r2(mEptHab1)
res = simulateResiduals(mEptHab1)
plot(res, rank = T)

# Predictions in relation with distance to the lamp ###################################################################
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataPip$Precision), Angle = "All")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataPip$Precision), Angle = "All"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataPip$Precision), Angle = "All"))
pred = predict(mPipHab1, newdata = newdata, type = "response", se.fit = TRUE)
pred2 = cbind(newdata, pred, species = "Pipistrellus")
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataPle$Precision), Angle = "All")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataPle$Precision), Angle = "All"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataPle$Precision), Angle = "All"))
pred = predict(mPleHab1, newdata = newdata, type = "response", se.fit = TRUE)
pred3 = cbind(newdata, pred, species = "Myotis/Plecotus")
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataEpt$Precision), Angle = "All")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataEpt$Precision), Angle = "All"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataEpt$Precision), Angle = "All"))
pred = predict(mEptHab1, newdata = newdata, type = "response", se.fit = TRUE)
pred4 = cbind(newdata, pred, species = "Eptesicus/Nyctalus")
pred = rbind(pred2, pred3, pred4)
vline = rep(999, nrow(pred))
pred = cbind(pred, vline)

# Average probabilities predictions #####################################################################################
newdata = data.frame(distLight = mean(data$distLight), spetrum = c("control", "red", "white", "control", "red", "white", "control", "red", "white"), 
                     species = c("Pipistrellus","Pipistrellus","Pipistrellus","Myotis/Plecotus","Myotis/Plecotus","Myotis/Plecotus",
                                 "Eptesicus/Nyctalus","Eptesicus/Nyctalus","Eptesicus/Nyctalus"),
                     ID = NA, Precision = mean(dataPip$Precision), Angle = "All")
pred1 = as.data.frame(predict(mPipHab1, newdata = newdata[c(1:3),], type = "response",se.fit = TRUE))
pred2 = as.data.frame(predict(mPleHab1, newdata = newdata[c(4:6),], type = "response",se.fit = TRUE))
pred3 = as.data.frame(predict(mEptHab1, newdata = newdata[c(7:9),], type = "response", se.fit = TRUE))
pred = rbind(pred1,pred2,pred3)
pred= cbind(newdata, pred)
pred$conf.low <- pred$fit - 1.96 * pred$se.fit
pred$conf.high <- pred$fit + 1.96 * pred$se.fit


# IN RELATION WITH VERTICAL POSITION (UNDER OR ABOVE LIGHT HEIGHT)
# PIPISTRELLES ###########################################################################################################
# under light
dataPipNegAngle = subset(dataPip, dataPip$angle < 0)
mPipHab2Neg <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mPipHab2Neg) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataPipNegAngle, weights = Precision)
summary(m) 
mPipHab2Neg <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPipNegAngle, weights = Precision)
summary(mPipHab2Neg) 
dataPipNegAngle <- within(dataPipNegAngle, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mPipHab2Neg)
r2(mPipHab2Neg)

# Indivual R2
m <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataPipNegAngle, weights = Precision)
r2(m)
m <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataPipNegAngle, weights = Precision)
r2(m)

# above light
dataPipPosAngle = subset(dataPip, dataPip$angle > 0)
mPipHab2Pos <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mPipHab2Pos) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataPipPosAngle, weights = Precision)
summary(m) 
mPipHab2Pos <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPipPosAngle, weights = Precision)
summary(mPipHab2Pos)
dataPipPosAngle <- within(dataPipPosAngle, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mPipHab2Pos)
r2(mPipHab2Pos)

# Indivual R2
m <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataPipPosAngle, weights = Precision)
r2(m)
m <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataPipPosAngle, weights = Precision)
r2(m)

# MYOTIS/PLECOTUS ###########################################################################################################
# under light
dataPleNegAngle = subset(dataPle, dataPle$angle < 0)
mPleHab2Neg <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mPleHab2Neg) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataPleNegAngle, weights = Precision)
summary(m) 
mPleHab2Neg <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPleNegAngle, weights = Precision)
summary(mPleHab2Neg) 
dataPleNegAngle <- within(dataPleNegAngle, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mPleHab2Neg)
r2(mPleHab2Neg) 

# Indivual R2
m <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataPleNegAngle, weights = Precision)
r2(m)
m <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataPleNegAngle, weights = Precision)
r2(m)

# above light
dataPlePosAngle = subset(dataPle, dataPle$angle > 0)
mPleHab2Pos <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mPleHab2Pos) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataPlePosAngle, weights = Precision)
summary(m) 
mPleHab2Pos <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPlePosAngle, weights = Precision)
summary(mPleHab2Pos) 
dataPlePosAngle <- within(dataPlePosAngle, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mPleHab2Pos)
r2(mPleHab2Pos)

# Indivual R2
m <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataPlePosAngle, weights = Precision)
r2(m)
m <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataPlePosAngle, weights = Precision)
r2(m)

# EPT/NYC ###########################################################################################################
# under light
dataEptNegAngle = subset(dataEpt, dataEpt$angle < 0)
mEptHab2Neg <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mEptHab2Neg) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataEptNegAngle, weights = Precision)
summary(m) 
mEptHab2Neg <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataEptNegAngle, weights = Precision)
summary(mEptHab2Neg) 
dataEptNegAngle <- within(dataEptNegAngle, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mEptHab2Neg)
r2(mEptHab2Neg) 

# Indivual R2
m <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataEptNegAngle, weights = Precision)
r2(m)
m <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataEptNegAngle, weights = Precision)
r2(m)

# above light
dataEptPosAngle = subset(dataEpt, dataEpt$angle > 0)
mEptHab2Pos <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mEptHab2Pos) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataEptPosAngle, weights = Precision)
summary(m)
mEptHab2Pos <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataEptPosAngle, weights = Precision)
summary(mEptHab2Pos)
dataEptPosAngle <- within(dataEptPosAngle, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mEptHab2Pos)
r2(mEptHab2Pos) 

# Indivual R2
m <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataEptPosAngle, weights = Precision)
r2(m)
m <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataEptPosAngle, weights = Precision)
r2(m)

# Predictions in relation with distance to the lamp ###################################################################
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataPipNegAngle$Precision), Angle = "Under")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataPipNegAngle$Precision), Angle = "Under"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataPipNegAngle$Precision), Angle = "Under"))
pred = predict(mPipHab2Neg, newdata = newdata, type = "response", se.fit = TRUE)
pred2 = cbind(newdata, pred, species = "Pipistrellus")
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataPipPosAngle$Precision), Angle = "Above")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataPipPosAngle$Precision), Angle = "Above"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataPipPosAngle$Precision), Angle = "Above"))
pred = predict(mPipHab2Pos, newdata = newdata, type = "response", se.fit = TRUE)
pred3 = cbind(newdata, pred, species = "Pipistrellus")
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataPleNegAngle$Precision), Angle = "Under")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataPleNegAngle$Precision), Angle = "Under"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataPleNegAngle$Precision), Angle = "Under"))
pred = predict(mPleHab2Neg, newdata = newdata, type = "response", se.fit = TRUE)
pred4 = cbind(newdata, pred, species = "Myotis/Plecotus")
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataPlePosAngle$Precision), Angle = "Above")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataPlePosAngle$Precision), Angle = "Above"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataPlePosAngle$Precision), Angle = "Above"))
pred = predict(mPleHab2Pos, newdata = newdata, type = "response", se.fit = TRUE)
pred5 = cbind(newdata, pred, species = "Myotis/Plecotus")
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataEptNegAngle$Precision), Angle = "Under")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataEptNegAngle$Precision), Angle = "Under"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataEptNegAngle$Precision), Angle = "Under"))
pred = predict(mEptHab2Neg, newdata = newdata, type = "response", se.fit = TRUE)
pred6 = cbind(newdata, pred, species = "Eptesicus/Nyctalus")
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataEptPosAngle$Precision), Angle = "Above")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataEptPosAngle$Precision), Angle = "Above"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataEptPosAngle$Precision), Angle = "Above"))
pred = predict(mEptHab2Pos, newdata = newdata, type = "response", se.fit = TRUE)
pred7 = cbind(newdata, pred, species = "Eptesicus/Nyctalus")
pred8 = rbind(pred2, pred3, pred4, pred5, pred6, pred7)



# AND IN RELATION WITH HORIZONTAL POSITION (LEFT OR RIGHT LIGHT SIDE)
# PIPISTRELLUS ###########################################################################################################
# left side
dataPipLeft = subset(dataPip, dataPip$x..m.-4 > 0)
mPipLeft <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mPipLeft) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataPipLeft, weights = Precision)
summary(m)
mPipLeft <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPipLeft, weights = Precision)
summary(mPipLeft)
dataPipLeft <- within(dataPipLeft, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mPipLeft)
r2(mPipLeft) 

# Indivual R2
m <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataPipLeft, weights = Precision)
r2(m)
m <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataPipLeft, weights = Precision)
r2(m)

# right side
dataPipRight = subset(dataPip, dataPip$x..m.-4 < 0)
mPipRight <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mPipRight) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataPipRight, weights = Precision)
summary(m) 
mPipRight <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPipRight, weights = Precision)
summary(mPipRight)
dataPipRight <- within(dataPipRight, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mPipRight)
r2(mPipRight) 

# Indivual R2
m <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataPipRight, weights = Precision)
r2(m)
m <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataPipRight, weights = Precision)
r2(m)

# MYOTIS/PLECOTUS ###########################################################################################################
# left side
dataPleLeft = subset(dataPle, dataPle$x..m.-4 > 0)
mPleLeft <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mPleLeft) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataPleLeft, weights = Precision)
summary(m)
mPleLeft <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPleLeft, weights = Precision)
summary(mPleLeft) 
dataPleLeft <- within(dataPleLeft, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mPleLeft)
r2(mPleLeft) 

# Indivual R2
m <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataPleLeft, weights = Precision)
r2(m)
m <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataPleLeft, weights = Precision)
r2(m)

# right side
dataPleRight = subset(dataPle, dataPle$x..m.-4 < 0)
mPleRight <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mPleRight) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataPleRight, weights = Precision)
summary(m) 
mPleRight <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPleRight, weights = Precision)
summary(mPleRight) 
dataPleRight <- within(dataPleRight, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mPleRight)
r2(mPleRight)

# Indivual R2
m <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataPleRight, weights = Precision)
r2(m)
m <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataPleRight, weights = Precision)
r2(m)

# EPT/NYC ###########################################################################################################
# left side
dataEptLeft = subset(dataEpt, dataEpt$x..m.-4 > 0)
mEptLeft <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mEptLeft) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataEptLeft, weights = Precision)
summary(m) 
mEptLeft <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataEptLeft, weights = Precision)
summary(mEptLeft) 
dataEptLeft <- within(dataEptLeft, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mEptLeft)
r2(mEptLeft) 

# Indivual R2
m <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataEptLeft, weights = Precision)
r2(m)
m <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataEptLeft, weights = Precision)
r2(m)

# right side
dataEptRight = subset(dataEpt, dataEpt$x..m.-4 < 0)
mEptRight <- glmer(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataPip, weights = Precision)
vif.mer(mEptRight) # check for collinearity
m <- glmmTMB(factor(EdgeSideBinom) ~ 1 + (1|ID), family = binomial, data=dataEptRight, weights = Precision)
summary(m) 
mEptRight <- glmmTMB(factor(EdgeSideBinom) ~ spetrum * distLight + (1|ID), family = binomial, data=dataEptRight, weights = Precision)
summary(mEptRight) 
dataEptRight <- within(dataEptRight, spetrum <- relevel(as.factor(spetrum), ref = "white")) # intercept permutation
summary(mEptRight)
r2(mEptRight) 

# Indivual R2
m <- glmmTMB(factor(EdgeSideBinom) ~ spetrum + (1|ID), family = binomial, data=dataEptRight, weights = Precision)
r2(m)
m <- glmmTMB(factor(EdgeSideBinom) ~ distLight + (1|ID), family = binomial, data=dataEptRight, weights = Precision)
r2(m)

# Predictions in relation with distance to the lamp ###################################################################
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataPipLeft$Precision), Side = "Left")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataPipLeft$Precision), Side = "Left"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataPipLeft$Precision), Side = "Left"))
pred = predict(mPipLeft, newdata = newdata, type = "response", se.fit = TRUE)
pred2 = cbind(newdata, pred, species = "Pipistrellus")
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataPipRight$Precision), Side = "Right")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataPipRight$Precision), Side = "Right"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataPipRight$Precision), Side = "Right"))
pred = predict(mPipRight, newdata = newdata, type = "response", se.fit = TRUE)
pred3 = cbind(newdata, pred, species = "Pipistrellus")
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataPleLeft$Precision), Side = "Left")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataPleLeft$Precision), Side = "Left"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataPleLeft$Precision), Side = "Left"))
pred = predict(mPleLeft, newdata = newdata, type = "response", se.fit = TRUE)
pred4 = cbind(newdata, pred, species = "Myotis/Plecotus")
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataPleRight$Precision), Side = "Right")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataPleRight$Precision), Side = "Right"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataPleRight$Precision), Side = "Right"))
pred = predict(mPleRight, newdata = newdata, type = "response", se.fit = TRUE)
pred5 = cbind(newdata, pred, species = "Myotis/Plecotus")
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataEptLeft$Precision), Side = "Left")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataEptLeft$Precision), Side = "Left"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataEptLeft$Precision), Side = "Left"))
pred = predict(mEptLeft, newdata = newdata, type = "response", se.fit = TRUE)
pred6 = cbind(newdata, pred, species = "Eptesicus/Nyctalus")
newdata = data.frame(distLight = seq(0,30,0.1), spetrum = "control", ID = NA, Precision = mean(dataEptRight$Precision), Side = "Right")
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "red", ID = NA, Precision = mean(dataEptRight$Precision), Side = "Right"))
newdata = rbind(newdata, data.frame(distLight = seq(0,30,0.1), spetrum = "white", ID = NA, Precision = mean(dataEptRight$Precision), Side = "Right"))
pred = predict(mEptRight, newdata = newdata, type = "response", se.fit = TRUE)
pred7 = cbind(newdata, pred, species = "Eptesicus/Nyctalus")
pred8 = rbind(pred2, pred3, pred4, pred5, pred6, pred7)




