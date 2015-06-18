# MNL

rm(list = ls())

########## Libraries to be loaded ##############

library(nnet)
library(scales)

########## Functions to be used ################

rename <- function(df1,u1,v1){
  oldnames=u1
  newnames=v1
  for(i in 1:length(u1))names(df1)[names(df1)==oldnames[i]]=newnames[i]
  return(df1)
}

############## Applying MNL for initial data ##################

data <- read.csv("USeBay.csv")
data$y.obs <- ifelse(data$a1_1 <= 6,-1,ifelse((data$a1_1 == 7 | data$a1_1 == 8),0,
                                              ifelse(data$a1_1 > 8,1,NA)))

data$y.obs.fact <- as.factor(ifelse(data$y.obs == 1,"Promoter",ifelse(data$y.obs == 0,"Passive",
                                                   ifelse(data$y.obs == -1,"Detractor",NA))))

# attach(data)
# detach(data)

data$y.obs.rl <- relevel(data$y.obs.fact, ref = "Passive")

init.mnl <- multinom(y.obs.rl ~ a4_a19_1 + a4_a19_2 + a4_a19_2 + a4_a19_3 + a4_a19_4 + a4_a19_5
                     + a4_a19_6, data)
summary(init.mnl)

z <- summary(init.mnl)$coefficients/summary(init.mnl)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2

init.pred <- data.frame(init.mnl$fitted.values)
init.data <- cbind(data,init.mnl$fitted.values)

init.data$y.pred <- (ifelse(apply(init.pred,1,max)==init.data$Promoter,1,
                            ifelse(apply(init.pred,1,max)==init.data$Passive,0,
                                   ifelse(apply(init.pred,1,max)==init.data$Detractor,-1,NA))))

a <- mean(data$y.obs)
b <- mean(init.data$y.pred)

init.a4_a19_1.per <- percent(sum(table(init.data$a4_a19_1)[9:10])/dim(init.data)[1])
init.a4_a19_2.per <- percent(sum(table(init.data$a4_a19_2)[9:10])/dim(init.data)[1])
init.a4_a19_3.per <- percent(sum(table(init.data$a4_a19_3)[9:10])/dim(init.data)[1])
init.a4_a19_4.per <- percent(sum(table(init.data$a4_a19_4)[9:10])/dim(init.data)[1])
init.a4_a19_5.per <- percent(sum(table(init.data$a4_a19_5)[9:10])/dim(init.data)[1])
init.a4_a19_6.per <- percent(sum(table(init.data$a4_a19_6)[9:10])/dim(init.data)[1])

write.csv(init.data,"init.data.csv",row.names = F)

init.nps.per <- percent(mean(init.data$y.pred))


per <- dim(data)[1]*0.01
random <- sample(1:500,per,replace = F)
data$a <- ifelse(data$respid %in% random,data$respid,0)
data$x1.n <- ifelse((data$a == 0 | data$x1 == max(data$x1)),data$x1,data$x1 + 1)
data$x2.n <- ifelse((data$a == 0 | data$x2 == max(data$x2)),data$x2,data$x2 + 1)

mnl.n <- multinom(data$y.obs ~ data$x1.n + data$x2.n, data)

pred.n <- data.frame(mnl.n$fitted.values)
data.n <- cbind(data,mnl.n$fitted.values)

data.n <- rename(data.n,c("1","0","-1"), 
                 c("pred.prob.prom.n","pred.prob.passi.n","pred.prob.detra.n"))

data.n$y.pred.n <- (ifelse(apply(pred.n,1,max)==data.n$pred.prob.prom.n,1,
                           ifelse(apply(pred.n,1,max)==data.n$pred.prob.passi.n,0,
                                  ifelse(apply(pred.n,1,max)==data.n$pred.prob.detra.n,-1,NA))))

x1.n.per <- percent(sum(table(data.n$x1.n)[9:10])/dim(data.n)[1])
x2.n.per <- percent(sum(table(data.n$x2.n)[9:10])/dim(data.n)[1])
nps.n.per <- percent(mean(data.n$y.pred.n))




