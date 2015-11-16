grp1 <- rep(NULL, 2539)
grp2 <- rep(NULL, 2539)
grp3 <- rep(NULL, 2539)
grp4 <- rep(NULL, 2539)
grp5 <- rep(NULL, 2539)
grp6 <- rep(NULL, 2539) 
grp7 <- rep(NULL, 2539)
grp8 <- rep(NULL, 2539)
grp9 <- rep(NULL, 2539)
grp10 <- rep(NULL, 2539)
grp11 <- rep(NULL, 2539)
grp12 <- rep(NULL, 2539)
grp1 <- ifelse(ccp2[,3]==0, 0, 1)
grp2 <- ifelse(ccp2[,4]==0, 0, 1)
grp3 <- ifelse(ccp2[,5]==0, 0, 1)
grp4 <- ifelse(ccp2[,6]==0, 0, 1)
grp5 <- ifelse(ccp2[,7]==0, 0, 1)
grp6 <- ifelse(ccp2[,8]==0, 0, 1)
grp7 <- ifelse(ccp2[,9]==0, 0, 1)
grp8 <- ifelse(ccp2[,10]==0, 0, 1)
grp9 <- ifelse(ccp2[,11]==0, 0, 1)
grp10 <- ifelse(ccp2[,12]==0, 0, 1)
grp11 <- ifelse(ccp2[,13]==0, 0, 1)
grp12 <- ifelse(ccp2[,14]==0, 0, 1)
grp <- data.frame(grp1, grp2, grp3, grp4, grp5, grp6, grp7, grp8, grp9, grp10, grp11, grp12)
grp <- as.matrix(grp)

library(nnet)
mod <- multinom(ccp2[,21]~grp1+grp2+grp3+grp4+grp5+grp6+grp7+grp8+grp9+grp10+grp11+grp12, model=T)
summary(mod)

mod <- multinom(ccp2[,22]~KL_biogrp1+KL_biogrp2+KL_biogrp3+KL_biogrp4+KL_biogrp5+KL_biogrp6+KL_biogrp7+KL_biogrp8+KL_biogrp9+KL_biogrp10+KL_biogrp11+KL_biogrp12, data=ccp2, model=T)

library(MASS)
aic <- stepAIC(mod, direction="backward")
aic$AIC

library(leaps)
modall <- regsubsets(ccp2[,21]~grp1+grp2+grp3+grp4+grp5+grp6+grp7+grp8+grp9+grp10+grp11+grp12, data=ccp2, nbest=1, nvmax=12)

modall <- regsubsets(ccp2[,21]~KL_biogrp1+KL_biogrp2+KL_biogrp3+KL_biogrp4+KL_biogrp5+KL_biogrp6+KL_biogrp7+KL_biogrp8+KL_biogrp9+KL_biogrp10+KL_biogrp11+KL_biogrp12, data=ccp2, nbest=1, nvmax=8)

plot(modall, scale="adjr2")


sum.mod <- summary(modall)
str(sum.mod)
pt(abs(sum.mod$coefficients/sum.mod$standard.errors), df=nrow(ccp2)-16, lower=F)


#alternative function: leaps() is compatibility wrapper for regsubsets()
library(leaps)
leapmod <- leaps(x = grp, y = ccp2[,21], names=names(ccp2)[3:14], method=c("Cp", "adjr2", "r2"), nbest=1)

library(car)
subsets(modall, statistic="cp")
abline(1, 1, lty=2, col="red")


#formular in regression models
new<-"rs4713243_G"
for (i in 22:1187) new <- append(new ,paste("+", names(ccp2)[i]))
formul <- reformulate(termlabels=new, response=names(ccp2)[3])

#ANOVA
library(multcomp)
K <- diag(length(coef(aov)))[-1,]
rownames(K) <- names(coef(aov))[-1]
K
ci <- glht(aov, linfct=K)
plot(ci)
