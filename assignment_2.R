attach(actg)

#PART 1
summary(glm(cd4_count~karnof, family= poisson))



plot(cd4_count~karnof, col=c("indianred", "deepskyblue", "palegreen", "purple"),
     main= "Karnofsky - CD4 relationship", ylab="cd4count", 
     xlab="Karnofsky Stage")

#PART 2
plot(cd4_count~sex)
#plot shows that mean is slightly higher for women
women<- cd4_count[sex=="Female"]
men <- cd4_count[sex=="Male"]

#CONVENTIONAL APPROACH
mean(women)
sd(women)
se <- sd(women)/sqrt(length(women))
se
confint <- c(qt(0.975, (length(women)-1)), qt(0.025, length(women)-1))
confint
standardconf <- mean(women)+se*confint
standardconf

#BOOTSTRAPPING APPROACH
bootmean <- replicate(10000, mean(sample(cd4_count, replace=T)))
hist(bootmean, breaks=100, main="Histogram of bootstrapped  cd4 count mean",
     xlab="bootstrapped mean", xlim=c(80, 93))
bootconf<-quantile(bootmean,p=c(0.975,0.025))
mean(bootmean)
bootconf


#PART 3
require(survival)

#saturated/maximal model (involving all factors, time & censorship)
saturated <- survreg(Surv(time_d, censor)~ idv+sex+ivdrug+hemo+karnof+cd4_count+zdv_months+age)
summary(saturated)
        
#exponential hazard model
expo <- survreg(Surv(time_d, censor)~ idv+sex+ivdrug+hemo+karnof+cd4_count+zdv_months+age, dist="exponential")
summary(expo)#

#using anova, which model is the better fit?
anova(saturated, expo)

#simplification of chosen model
simp <- step(survreg(Surv(time_d, censor)~ idv+sex+ivdrug+hemo+karnof+cd4_count+zdv_months+age, dist="exponential"))
summary(simp)

# mam = minimal adequate model
mam <- survreg(Surv(time_d, censor) ~ idv+karnof+cd4_count+age, dist ="exponential")
summary(mam)
coef(mam)

#compare the mam to the chosen saturated model
anova(expo, mam)


predict(mam, list(age=40, cd4_count=50, karnof="Stage IV", ivdrug="user", idv="untreated"), type="response")

#km plot of time to death with karnof score as a factor
model<- (survfit(Surv(time_d, censor)~karnof))
plot(model, conf.int="both", ylab="Survivorship", xlab="Time to death", col=c("red", "blue", "green", "black"))
legend("bottomleft", legend=c("Stage I", "Stage II", "Stage III", "Stage IV"), col=c("red", "blue", "green", "black"), lwd=2)
