data2<-read.csv("rnd_eff_temp.csv",sep=",")

library("rgbif")
library(tidyverse)
library(dplyr)
spp_check_ok  <- name_backbone("Felis catus",verbose=T,strict=T) #We initialize a data frame for a species we know is in backbone taxonomy
spp_check_ok  <- spp_check_ok[-1,] #Remove this now
spp_check_bad  <- name_backbone("xxx",verbose = T,strict = T) #We start a dataframe with a species that in not in the backbone
spp_check_bad  <- spp_check_bad[-1,]

for(i in 1:nrow(data2)) {
  toto  <- name_backbone(data2[i,1],verbose=T,strict = T) #Check species I againts the abckbone
  if(length(which(names(toto)=="acceptedUsageKey"))==1) { #IF there is a colum acceptedusagekey we remove it and will not be included for all species
    toto  <- toto[,-which(names(toto)=="acceptedUsageKey")]
  }
  if(ncol(toto)==ncol(spp_check_ok)) { #If there are 23 colums the species name was recognised
    if(length(which(toto$status=="ACCEPTED")) >0) {
      #If there is a species with status ACCEPTED in the returned datframe
      spp_check_ok  <- rbind(spp_check_ok,toto[which(toto$status=="ACCEPTED"),]) #If so we only keep this name
    } else if(length(which(toto$status=="SYNONYM"))>0 ){
      #If there is no species with the status ACCEPTED in the returend datframe is there a species with the name SYNONYM instead?
      warning(paste("Species",data2[i,1],"is a synonym")) #We print a warning
      spp_check_ok   <- rbind(spp_check_ok,toto[which(toto$status=="SYNONYM")[1],])
    } else if(length(which(toto$status=="DOUBTFUL"))>0) {
      warning(paste("Species",data2[i,1],"is doubtful"))
      spp_check_ok  <- rbind(spp_check_ok,toto[which(toto$status=="DOUBTFUL")[1],])
    } else {
      stop("Status unknown") #The status is neither of the others 
    }
    
  }
  else if(ncol(toto)==ncol(spp_check_bad)){
    spp_check_bad  <- rbind(spp_check_bad,toto)
  }
  #else{
  # stop("Unknown length") #If we have  a data frame witha  different size we wnat to check why
  #}
}
Lepidoptera <-spp_check_ok[which(spp_check_ok$order == "Lepidoptera"),c(3,11)]
LepidopteraResults<-cbind(1,1,1,1)
LepidopteraResults<- data.frame(LepidopteraResults)
LepidopteraResults<-LepidopteraResults[-1,]
colnames(LepidopteraResults)<-c("Species","Slope","Standard Error", "Family")
for (i in 1:nrow(Lepidoptera) ) {
 LepidopteraResults[i,]<-  data2[which(data2[,1] == paste(Lepidoptera[i,1])),c(1,4,6)]
 LepidopteraResults[i,4]<-  Lepidoptera[i,2]
}
  ButterflyResults<-  LepidopteraResults[LepidopteraResults$Family != "Noctuidae" & LepidopteraResults$Family != "Geometridae" & LepidopteraResults$Family != "Zygaenidae" & LepidopteraResults$Family != "Saturniidae" & LepidopteraResults$Family != "Limacodidae" & LepidopteraResults$Family != "Erebidae"  & LepidopteraResults$Family != "Sphingidae" & LepidopteraResults$Family != "Drepanidae" & LepidopteraResults$Family != "Lasiocampidae",]
write.csv(ButterflyResults, file="ButterflyRelations_Paper32.csv")
  