library(tidyr)

df <- data.frame(
  TreatmentA = c(175, 168, 168, 190, 156, 181, 182, 175, 174, 179),
  TreatmentB = c(185, 169, 173, 173, 188, 186, 175, 174, 179, 180) 
)

summary(df)

sample.size <- length(df$TreatmentA)
reshape_df <- data.frame(
  Treatment = rep(c("A","B"),each=sample.size),
  Mass = c(df$TreatmentA,df$TreatmentB)
)
plot(Mass~Treatment, data=reshape_df)

# boxplot(df$GroupA,df$GroupB,names=c("GroupA","GroupB"))  # (alternative method!)

observed_dif <- mean(df$GroupA) - mean(df$GroupB)
observed_dif

t.test(df$TreatmentA,df$TreatmentB, var.equal=TRUE, paired=FALSE)

lots <- 1000000  # large number approximating infinity 

popMean_null <- mean(reshape_df$Mass)        # assume groups A and B come from a population with common mean 
popSD_null <- sd(reshape_df$Mass)                      # and common standard deviation... 
popData_null <- rnorm(n=lots,mean=popMean_null,sd=popSD_null)    # the statistical "population" of interest (under null model w no treatment effect)

sampleA <- sample(popData_null,size=sample.size)    # use R's native "sample()" function
sampleB <- sample(popData_null,size=sample.size)

round(sampleA)
difference <- mean(sampleA)-mean(sampleB)   # sample statistic = difference between sample means
difference

reps <- 1000
null_difs <- numeric(reps)

for(i in 1:reps){
  sampleA <- sample(popData_null,size=sample.size)   
  sampleB <- sample(popData_null,size=sample.size)
  null_difs[i] <- mean(sampleA)-mean(sampleB)
}

hist(null_difs)
abline(v=observed_dif,col="green",lwd=3)

  ordered_difs <- sort(abs(null_difs))   
  higher_anomaly <- length(which(ordered_difs>=abs(observed_dif)))
  p_value <- higher_anomaly/reps  
  p_value

reps <- 5000
null_difs <- numeric(reps)
for (i in 1:reps){					
	newindices <- sample(c(1:nrow(reshape_df)))   			
  newGroup <- reshape_df$Group[newindices]			
	dif <- mean(reshape_df$Mass[newGroup=="A"])	- mean(reshape_df$Mass[newGroup=="B"])	
	null_difs[i] <- dif	
}
hist(null_difs)
abline(v=observed_dif,col="green",lwd=3)


higher_anomaly <- length(which(abs(null_difs)>=abs(observed_dif)))
p_value <- higher_anomaly/reps  
p_value

head(trees)   # use help(trees) for more information

plot(trees$Volume~trees$Height, main = 'Black Cherry Tree Height/Volume Relationship', xlab = 'Height', ylab = 'Volume', pch = 16, col ='blue')
plot(trees$Volume~trees$Girth, main = 'Black Cherry Tree Girth/Volume Relationship', xlab = 'Girth', ylab = 'Volume', pch = 16, col ='red')

Rsquared <- function(df,responsevar="Volume"){    # interactions not yet implemented
  response <- df[,responsevar]
  names <- names(df)
  rsq <- numeric(length(names))
  names(rsq) <- names(df)
  rsq <- rsq[names(rsq)!=responsevar]
  for(i in names(rsq)){         # loop through predictors
      predictor <- df[,i]
      model <- lm(response~predictor)  
      rsq[i] <- summary(model)$r.square 
  }
  return(rsq)
}

stat <- Rsquared(trees,"Volume")
stat


boot_sample <- function(df,statfunc,n_samples,n_stats){
  indices <- c(1:nrow(df))
  output <- matrix(NA,nrow=n_samples,ncol=n_stats)
  for(i in 1:n_samples){
    boot_rows <- sample(indices,size=nrow(df),replace=T)
    newdf <- df[boot_rows,]
    output[i,] <- statfunc(newdf)
  }
  return(output)
}

boot <- boot_sample(trees,Rsquared,10,length(stat))
colnames(boot) <- names(stat)

boot
stat

boot <- boot_sample(trees,Rsquared,1000,length(stat))   # 1000 bootstrap samples
confint <- apply(boot,2,function(t)  quantile(t,c(0.025,0.5,0.975)))
colnames(confint) <- names(stat)
t(confint)
quantile(boot,c(0.025,0.975))

