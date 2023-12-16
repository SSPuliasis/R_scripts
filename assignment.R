#QUESTION 1
attach(polychor)
plot(age, conc)
model <- glm(conc~age, poisson)
summary(model)
coef(model)

xpred <- seq(0,12,0.01)
ypred <- predict(model, list(age=xpred), type="response", se.fit=T)

lines(ypred$fit~xpred, col="blue")
lines(ypred$fit+1.96*ypred$se.fit~xpred, col="red", lty=2)
lines(ypred$fit-1.96*ypred$se.fit~xpred, col="red", lty=2)

simplemodel<- glm(conc~1)

anova(model, simplemodel, test="Chisq")
summary(simplemodel)

yval<- predict(model, list(age=5), type="response", se.fit=T)
yval
abline(h=137)
detach(polychor)

#QUESTION 2
attach(radiation)

prop <- (Survived/400)

plot(prop~lab)

y<-cbind(Survived, 400-Survived)

model<- glm(y~lab, binomial)
summary(model, type="response")
coef(model)
confint(model)

model2 <- glm(y~1, binomial)
anova (model, model2, test="Chisq")

predict(model, list(lab="g"), type="response")

detach(radiation)

#QUESTION 3
attach(wall)
plot(Weight~Age, pch=19, col="lightblue")

#ASYMTOTIC MODEL
asymodel<- nls(Weight~SSasymp(Age,a,b,c), data=wall)
summary(asymodel)
coef(asymodel)
xpred<- seq(0,1300,1)
ypred<- predict(asymodel, list(Age=xpred), type="response")
lines(xpred,ypred, col="blue", lty=1)
confint(asymodel)

yasyval <-predict(asymodel, list(Age=250), type="response")
abline(h=17204)
abline(v=250)

#LOGISTIC MODEL
logmodel<- nls(Weight~SSlogis(Age,a,b,c))
summary(logmodel)
coef(logmodel)
bpred<- predict(logmodel, list(Age=xpred), type="response")
lines(xpred,bpred, col="darkgreen", lty=1)
confint(logmodel)
curve(55701.1/(1+exp((334.6-x)/71.1)), from=0, to=1300, add=T, col="green")
curve(58094.3/(1+exp((348.8-x)/80.8)), from=0, to=1300, add=T, col="green")
#Asym/(1+exp((xmid-input)/scal))
#this model is a better fit

ylogval <- predict(logmodel, list(Age=500), type="response")
abline(h=50617)
abline(v=500)

detach(wall)
