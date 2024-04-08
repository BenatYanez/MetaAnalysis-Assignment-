library(metafor)

slope <- -0.25
intercept <- 0
predictor <- rnorm(n=100, mean=10, sd=10)
response <- intercept+ slope*predictor+ rnorm(n=100, mean=0, sd=40)
plot(predictor, response)
model <- lm(response ~predictor)
model
summary(model)
#Create a loop generating many datasets with diffrent size and store outputs in a matrix
store <-matrix(nrow=200,ncol=4)
for(x in 1:200) {
  samplesize<- ceiling(exp(rnorm(1,4.5,1.5)))+3
  #we're using this code to select sample sizes at random from a log normal distribution, so that small sample sizes are common and large sample sizes are rare. And n is always > 3.                  
  
  predictor <-rnorm(n=samplesize,mean=10,sd=10)
  response <- intercept+slope*predictor+rnorm(n=samplesize, 0,40)
  
  model<-lm(response ~predictor)
  store[x,] <-c(samplesize, summary(model)$coefficients[2,1:2],summary(model)$coefficients[2,4])
}
store <- as.data.frame(store)
names(store)<- c("n","slope","standard.error","p.value")

#Make a funnel plot from this simulated dataset
par(mfrow=c(1,2))
plot(store$slope,store$n,xlab="Slope",ylab="Sample size")
plot(store$slope, (1/store$standard.error), xlab="Slope", ylab="Precision, (1/SE)")

#Indicate the slope we used in the simulation and ahow which points have a p-value bellow 0.05
signslope <- which(store$p.value < 0.05)
par(mfrow=c(1,2))
plot(store$slope,store$n,xlab="Slope",ylab="Sample size")
points(store$slope[signslope],store$n[signslope],pch=16,col="red")
abline(v=slope,lty=2)
plot(store$slope, (1/store$standard.error), xlab="Slope", ylab="Precision, (1/SE)")
points(store$slope[signslope],(1/store$standard.error[signslope]),pch=16,col="red")
abline(v=slope,lty=2)

#Q1 What happens to effect size as sample size and precision increse?
  #They tend to be closer to the true value which is -0.25
#Q2 When precision is low what can you say about the estimates of effect size that tend to be significant?
  #They tend to have a negative values, they have quite a large magnitude (extreme outcomes) not close to the true value 
#Q3 If there is a tendency to only publish results when the effect is significant and or sample size is large, which points would be missing from the plot?
  #The points with low precision /sample size that are near a slope value of 0.

#Estimating the mean effect size
model2 <- lm(slope~1,data=store)
summary(model2)
#The problem with this analysis is that it ignores sampling variance , some slopes are estimated more acuratelly than others

?rma
meta <-rma(yi=store$slope,sei=store$standard.error, data=store)
summary(meta)
funnel(meta)
forest(meta, cex.lab=0.8, cex.axis=0.8, addfit=TRUE, shade="zebra")

#Simulate a new dataset with changing slope similar to change in temperature with latitude
latitude <- runif(100,0,90)#Randomly sample latitude from 0,90

slope <-0+latitude*-0.1+rnorm(100,0,3)
plot(latitude, slope)

store2 <-matrix(nrow=200,ncol=7)
#Generate species random effects
species<-rep(1:20,each=10)
specieseffect <-rep(rnorm(20,0,2),each=10)
for(x in 1:200) {
  #we're using this code to select sample sizes at random from a log normal distribution, so that small sample sizes are common and large sample sizes are rare. And n is always > 3.                  
  latitude <- runif(1,0,90)#Randomly sample latitude from 0,90
  slope <- specieseffect[x]+0+latitude*-0.1+rnorm(1,0,3)
  
  samplesize<- ceiling(exp(rnorm(1,4.5,1.5)))
  if(samplesize >3) {
    #Include this so that we dont run analyses on sample sizes that are too small
    predictor <-rnorm(n=samplesize,mean=10,sd=10)
    response <- intercept+ slope*predictor+ rnorm(n=samplesize, 0,40)
    model<-lm(response ~predictor)
    store2[x,] <-c(samplesize, summary(model)$coefficients[2,1:2],summary(model)$coefficients[2,4],latitude,species[x],x)
    
  }
}
store2 <- as.data.frame(store2)
names(store2)<- c("n","slope","standard.error","p.value","latitude","species","ID")

par(mfrow=c(1,2))
plot(store2$slope,store2$n,xlab="Slope",ylab="Sample size")
plot(store2$slope, (1/store2$standard.error), xlab="Slope", ylab="Precision, (1/SE)")

meta2<-rma(yi=store2$slope,sei=store2$standard.error, data=store2)
summary(meta2)
funnel(meta2)
forest(meta2, cex.lab=0.8, cex.axis=0.8, addfit=TRUE, shade="zebra")
#Q1 Why does the slope not estimate a funnel so much
  #Because each species has a different slope, there is no one true value for a slope for all the samples.

#Run a meta analysis controlling for latitude
meta3<-rma(yi=store2$slope,sei=store2$standard.error,mods=~latitude ,data=store2)
meta3
funnel(meta3)

#Add random term to the analysis
store2$se2 <-store2$standard.error^2
store3 <- store2[-which(is.na(store2$slope)==TRUE),]
meta4 <- rma.mv(yi=slope, V=se2,mods=~latitude, random=~1|species/ID,data=store3)
meta4

#Use real data
birdbroods <- read.csv("birdbroods.csv",sep=",",header=TRUE)
par(mfrow=c(1,2))
plot(birdbroods$slope, (1/birdbroods$slope.SE), xlab="Slope", ylab="Precision, (1/SE)")
birdbroods$se2 <-birdbroods$slope.SE^2

meta5 <- rma.mv(yi=slope, V=se2,random=~1|Species/id.pop,data=birdbroods)
meta5
funnel(meta5)
forest(meta5, cex.lab=0.8, cex.axis=0.8, addfit=TRUE, shade="zebra",order="obs")
#Q1 HAs brood size of the average bird species declined?
  #The estimate is very slightly negative (although not significant), so there is no significant change in average clutch size over time
#Q2 IS more of the variation in slope estimates distributed among or between species?
  #There is more variation within species 0.0010 than among species 0.0003 although this difference is small
#Q3 IS the trend in brood size more positive for populations in protected areas?
meta6 <- rma.mv(yi=slope, V=se2,mods=~protected.area,random=~1|Species/id.pop,data=birdbroods)
meta6
#There is no significant evidence that brood size is more positive in protected areas, in fact it would be slightly negative -0.008
#Q4 If the information was available what other terms would you include?
  #Habitat?, Study as a random effects, length of study, average year of study
