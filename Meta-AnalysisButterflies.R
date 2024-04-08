library(metafor)
eval(metafor:::.MuMIn)
library(MuMIn)
library(tidyverse)
library(kableExtra)
#Open the raw data and adjust it to only include studies with Standard errors
data <-read.csv("data/Meta_AnalysisData.csv",sep=";")
colnames(data)<- c("Study","Taxon","Slope","SE","SampleYears","SampleSites","Samples","Country","Measurement","AnalysisType","Years","pvalue","MeanLatitude","NºYears","Voltinism")
data$Slope<-as.numeric(data$Slope)
data$Study<-as.factor(data$Study)
MissingSE<-is.na(data$SE)
dataWithSE<-subset(data,subset=!MissingSE)
dataWithSE$ID<-as.factor(c(1:nrow(dataWithSE)))


#Generate the model to obtain the overall effect (Just include the random effect)
EstimateModel<-rma.mv(yi=dataWithSE$Slope,V=dataWithSE$SE2,random=list(~1|Taxon/ID,~1|Study), data=dataWithSE)
summary(EstimateModel)
#Obtain heterogeneity in the model
W <- diag(1/EstimateModel$vi)
X <- model.matrix(EstimateModel)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(EstimateModel$sigma2) / (sum(EstimateModel$sigma2) + (EstimateModel$k-EstimateModel$p)/sum(diag(P)))

sum(EstimateModel$fit.stats)

#Generate the maximal model and generate a table to save the results
metaMuMin<-rma.mv(yi=dataWithSE$Slope,V=dataWithSE$SE2,mods=~Voltinism+MeanLatitude+AnalysisType+Measurement,random=list(~1|Taxon/ID,~1|Study),slab=Taxon, data=dataWithSE)
summary(metaMuMin)
candidate_models<-dredge(metaMuMin)
candidates_aic2 <- subset(candidate_models, delta<=2) # Display all models within 2 values of AICc
importance <- sw(model.avg(candidate_models, subset=delta<=2))# relative importance (sum of weights) of the moderators
mod.avg <- summary(model.avg(candidate_models, subset=delta<=2)) # Generate model-averaged estimates 
confidence <- confint(mod.avg, full=TRUE) # Generate confidence intervals for the estimates averaged using full-averages procedures
#Create a table to store the results
tibble(`Fixed effect` = c("Intercept", "Temporal", "Temporal+Spatial", "Mean Latitude",
                           "Mean GBIF Record","First 20% Appearance", "Univoltine"), Estimate = mod.avg$coefmat.full[,
                                                                                                                              1], `Lower CI [0.025]` = confidence[, 1], `Upper CI  [0.975]` = confidence[,
                                                                                                                                                                                                         2], `p-value` = c("0.0008",0.7104,"<0.0001","0.0006","<0.0001",0.1341,0.8623),`Significance`=c("***","","***","***","***","","")) %>%
  kable("html", digits = 3) %>%
  kable_styling("striped", position = "left")

#Generate forest plot of effect sizes grouped by species
EscalcedData <-escalc(dataWithSE,measure="SMD",yi=Slope,vi=SE,data=dataWithSE)
agg <-aggregate(EscalcedData,cluster=Taxon,V=vcov(metaMuMin,type="obs"),addK=T)
EscalcedMeta <-rma(yi,vi,method="EE",data=agg,slab=Taxon)
pdf(file="ForestPlot.pdf",height=16,width=12)

  forest(EscalcedMeta, addpred=T,xlim=c(-40,40),at=seq(-20,30,by=1),order="obs",  shade="zebra",
         ilab=Voltinism, ilab.xpos=c(-25),
         cex=0.75, header="Species",tabfig=2, col="#CA625A",)
  addpoly.default(EstimateModel$beta,sei=EstimateModel$se,rows=c(-3),col="#5CC8A7")
  text(c(-25), EscalcedMeta$k+3, c("Voltinism"), cex=0.75, font=2 )
  text(c(-36.5),-3, c("Random effect Model"),cex=0.75)
  dev.off()
#Generate a plot of latitude against effect size
regplot(metaMuMin,mod="MeanLatitude",bg="#5CC8A7",ylab="Slope")
#Generate an orchard plot of the impact of different measurement types on effect size
install.packages("remotes")
remotes::install_github("daniel1noble/orchaRd")
library(orchaRd)
mr2 <- mod_results(metaMuMin, mod = "Measurement",group="Measurement")
mr2_mod_table <- mr2$mod_table

mr2_data <- mr2$data
mr2_data$moderator <- factor(mr2_data$moderator, levels = mr2_mod_table$name, labels = mr2_mod_table$name)
mr2_data$scale <- (1/sqrt(mr2_data[,"vi"]))
mr2_mod_table$K <- as.vector(by(mr2_data, mr2_data[,"moderator"], function(x) length(x[,"yi"])))
mr2_group_no <- nrow(mr2_mod_table)

fig_method <- ggplot(data = mr2_mod_table, aes(x = estimate, y = name)) + 
  scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, by = 1), labels = c("-10", "-9", "-8", "-7", "-6", "-5", "-4", "-3","-2","-1","0","1", "2", "3", "4", "5", "6", "7", "8","9","10")) + 
  ggbeeswarm::geom_quasirandom(data = mr2_data, aes(x = yi, y = moderator, size = scale, colour=moderator), groupOnX=FALSE, alpha = 0.4) +
  geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR), height = 0, show.legend = F, size = 0.5, alpha = 0.6) + # CI
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL), height = 0, show.legend = F, size = 1.2) + 
  geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.5) +
  geom_point(aes(fill=name), size = 3, shape = 21) + 
  annotate('text', x = 0.93, y = (seq(1, mr2_group_no, 1)+0.35), label= paste("italic(k)==", mr2_mod_table$K), parse=TRUE, hjust = "left", size=3.5) +
  scale_color_manual(values = c("Appearance Mean First" = "#CC79A7",  "Observed 20% first appearances " = "#D55E00",  "Mean Day of Record in GBIF"= "#0072B2")) +
  scale_fill_manual(values = c("Appearance Mean First" = "#CC79A7",  "Observed 20% first appearances " = "#D55E00",  "Mean Day of Record in GBIF"= "#0072B2")) +
  ggplot2::theme_bw() +
  ggplot2::guides(fill="none", colour="none") +
  ggplot2::theme(axis.text.y = element_text(size=9, colour="black", hjust=0.5, angle=90),
                 legend.position = "none", axis.title.y = element_blank(), axis.text.x=element_text(size=9)) +
  labs(x="Effect size (Slope)",y=c("Mean First Appearance","Mean GBIF Record","First 20% Observed"))

#Check for bias in publication
regtest(x=metaMuMin, sei=SE, data=dataWithSE,
        model="rma", predictor="sei", ret.fit=FALSE)
regtest(x=EstimateFixedModel, sei=SE, data=dataWithSE,
        model="rma", predictor="sei", ret.fit=FALSE)
#Generate a funnel plot
funnel(EstimateModel,level=c(90, 95, 99),lty=1 ,shade=c("#595959", "#7F7F7F", "#A6A6A6"),back="white" , legend=TRUE)
