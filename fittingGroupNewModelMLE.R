#clean working space for a good start.
rm(list=ls())

#read the  model in
source(paste(getwd(),"ACT-M_model.R",sep="/"))
source(paste(getwd(),"createStudySequence.R",sep="/"))

#error calculation.
MODlnl <-function(parms, d) {
  predictions<-getregpred(parms, d)
  if ((min(parms) < 0)){
    return(.Machine$double.xmax)
  } else {
    #lnL <- -sum(log(predictions)
    lnL <- -sum(d[,4]*log(predictions) + (1-data[,4])*log(1-predictions))
      #-sum(d[,4]*log(predictions)) #+ (240-d)*log(1-predictions))
return (lnL)}
}

#the model goes here.
getregpred <- function(parms,data) {
  weight_shared_same <- parms[1]
  weight_shared_different<- parms[2]
  weight_distinctive_same <-parms[3]
  weight_distinctive_different <-parms[4]
  c <- parms[5]
  gamma <- parms[6]
  predictions = c()
  for(i in 1:nrow(data)){
    predictions[i] <- model(data[i,1],data[i,2],data[i,3],c,2,weight_shared_same,weight_shared_different,weight_distinctive_same,weight_distinctive_different,gamma,2,0,studysequenceIntHighSim = studyItemsIHS,studysequenceBlkHighSim = studyItemsBHS ,studysequenceIntLowSim = studyItemsILS,studysequenceBlkLowSim = studyItemsBLS,randChoice = 0.25)}
  
  getregpred = predictions
  #print(getregpred)
  return(predictions)
}

#read subject data in
d<-read.csv(paste(getwd(),"data_all.csv",sep="/"), header=TRUE, row.names=NULL)
library(plyr)
data = as.matrix(d)
rm(d)

iterations = 50
fittingResults = matrix(NA, nrow = iterations, ncol = 20)

for(iteration in 1:iterations){
  IHSname = paste("studyItemsIHS",iteration,sep="")
  BHSname = paste("studyItemsBHS",iteration,sep="")
  ILSname = paste("studyItemsILS",iteration,sep="")
  BLSname = paste("studyItemsBLS",iteration,sep="")
  resultsname = paste("fittingIteration",iteration,sep="")
  
  studyItemsIHS = createStudySequence(probability_of_cat_repetition = .25,stimuliNumber = 0)
  studyItemsBHS = createStudySequence(probability_of_cat_repetition = .75,stimuliNumber = 0)
  studyItemsILS = createStudySequence(probability_of_cat_repetition = .25,stimuliNumber = 1)
  studyItemsBLS = createStudySequence(probability_of_cat_repetition = .75,stimuliNumber = 1)
  
  assign(IHSname,studyItemsIHS)
  assign(BHSname,studyItemsBHS)
  assign(ILSname,studyItemsILS)
  assign(BLSname,studyItemsBLS)
  
  startParms <- c(0.90, 0.10, 0.40, 0.90, 1,1)
  lowerParm <- c(0.0000000001,0.00000000001,0.0000000001,0.00000000001,0,0.000000001)
  upperParm <- c(1,1,1,1,100,100)
  xout <- nlminb(startParms ,MODlnl,gradient=NULL,d=data,control=list(trace=1),lower=lowerParm, upper=upperParm) 
  
  assign(resultsname,xout)
  
  fittingResults[iteration,1] = xout$par[1]
  fittingResults[iteration,2] = xout$par[2]
  fittingResults[iteration,3] = xout$par[3]
  fittingResults[iteration,4] = xout$par[4]
  fittingResults[iteration,5] = xout$par[5]
  fittingResults[iteration,6] = xout$par[6]
  fittingResults[iteration,7] = xout$objective
  fittingResults[iteration,8] = xout$iterations
  fittingResults[iteration,9] = xout$evaluations[1]
  fittingResults[iteration,10] = xout$convergence
  
  highsim_old_int = model(stimuliNumber = 0,typeItem = 0,probability_of_cat_repetition = 0.25,c = xout[[1]][5],r = 2,weight_shared_same = xout[[1]][1],weight_shared_different = xout[[1]][2],weight_distinctive_same = xout[[1]][3],weight_distinctive_different = xout[[1]][4],gamma = xout[[1]][6],p = 2,L=0,studysequenceIntHighSim = studyItemsIHS,studysequenceBlkHighSim = studyItemsBHS ,studysequenceIntLowSim = studyItemsILS,studysequenceBlkLowSim = studyItemsBLS,randChoice = 0.25)
  
  highsim_old_blk = model(stimuliNumber = 0,typeItem = 0,probability_of_cat_repetition = 0.75,c = xout[[1]][5],r = 2,weight_shared_same = xout[[1]][1],weight_shared_different = xout[[1]][2],weight_distinctive_same = xout[[1]][3],weight_distinctive_different = xout[[1]][4],gamma = xout[[1]][6],p = 2,L=0,studysequenceIntHighSim = studyItemsIHS,studysequenceBlkHighSim = studyItemsBHS ,studysequenceIntLowSim = studyItemsILS,studysequenceBlkLowSim = studyItemsBLS,randChoice = 0.25)
  
  highsim_new_int = model(stimuliNumber = 0,typeItem = 1,probability_of_cat_repetition = 0.25,c = xout[[1]][5],r = 2,weight_shared_same = xout[[1]][1],weight_shared_different = xout[[1]][2],weight_distinctive_same = xout[[1]][3],weight_distinctive_different = xout[[1]][4],gamma = xout[[1]][6],p = 2,L=0,studysequenceIntHighSim = studyItemsIHS,studysequenceBlkHighSim = studyItemsBHS ,studysequenceIntLowSim = studyItemsILS,studysequenceBlkLowSim = studyItemsBLS,randChoice = 0.25)
  
  highsim_new_blk = model(stimuliNumber = 0,typeItem = 1,probability_of_cat_repetition = 0.75,c = xout[[1]][5],r = 2,weight_shared_same = xout[[1]][1],weight_shared_different = xout[[1]][2],weight_distinctive_same = xout[[1]][3],weight_distinctive_different = xout[[1]][4],gamma = xout[[1]][6],p = 2,L=0,studysequenceIntHighSim = studyItemsIHS,studysequenceBlkHighSim = studyItemsBHS ,studysequenceIntLowSim = studyItemsILS,studysequenceBlkLowSim = studyItemsBLS,randChoice = 0.25)
  
  lowsim_old_int = model(stimuliNumber = 1,typeItem = 0,probability_of_cat_repetition = 0.25,c = xout[[1]][5],r = 2,weight_shared_same = xout[[1]][1],weight_shared_different = xout[[1]][2],weight_distinctive_same = xout[[1]][3],weight_distinctive_different = xout[[1]][4],gamma = xout[[1]][6],p = 2,L=0,studysequenceIntHighSim = studyItemsIHS,studysequenceBlkHighSim = studyItemsBHS ,studysequenceIntLowSim = studyItemsILS,studysequenceBlkLowSim = studyItemsBLS,randChoice = 0.25)
  
  lowsim_old_blk = model(stimuliNumber = 1,typeItem = 0,probability_of_cat_repetition = 0.75,c = xout[[1]][5],r = 2,weight_shared_same = xout[[1]][1],weight_shared_different = xout[[1]][2],weight_distinctive_same = xout[[1]][3],weight_distinctive_different = xout[[1]][4],gamma = xout[[1]][6],p = 2,L=0,studysequenceIntHighSim = studyItemsIHS,studysequenceBlkHighSim = studyItemsBHS ,studysequenceIntLowSim = studyItemsILS,studysequenceBlkLowSim = studyItemsBLS,randChoice = 0.25)
  
  lowsim_new_int = model(stimuliNumber = 1,typeItem = 1,probability_of_cat_repetition = 0.25,c = xout[[1]][5],r = 2,weight_shared_same = xout[[1]][1],weight_shared_different = xout[[1]][2],weight_distinctive_same = xout[[1]][3],weight_distinctive_different = xout[[1]][4],gamma = xout[[1]][6],p = 2,L=0,studysequenceIntHighSim = studyItemsIHS,studysequenceBlkHighSim = studyItemsBHS ,studysequenceIntLowSim = studyItemsILS,studysequenceBlkLowSim = studyItemsBLS,randChoice = 0.25)
  
  lowsim_new_blk = model(stimuliNumber = 1,typeItem = 1,probability_of_cat_repetition = 0.75,c = xout[[1]][5],r = 2,weight_shared_same = xout[[1]][1],weight_shared_different = xout[[1]][2],weight_distinctive_same = xout[[1]][3],weight_distinctive_different = xout[[1]][4],gamma = xout[[1]][6],p = 2,L=0,studysequenceIntHighSim = studyItemsIHS,studysequenceBlkHighSim = studyItemsBHS ,studysequenceIntLowSim = studyItemsILS,studysequenceBlkLowSim = studyItemsBLS,randChoice = 0.25)
  
  fittingResults[iteration,11] = highsim_old_int
  fittingResults[iteration,12] = highsim_old_blk
  fittingResults[iteration,13] = highsim_new_int
  fittingResults[iteration,14] = highsim_new_blk
  fittingResults[iteration,15] = lowsim_old_int
  fittingResults[iteration,16] = lowsim_old_blk
  fittingResults[iteration,17] = lowsim_new_int
  fittingResults[iteration,18] = lowsim_new_blk
  
  BIC = 2*log(fittingResults[iteration,7])+6*log(8)
  
  fittingResults[iteration,19] = BIC
  
  AIC =  2*log(fittingResults[iteration,7])+2*8
  
  fittingResults[iteration,20] = AIC
  
}

print("finished iterating")
library(data.table)
fwrite(as.data.table(fittingResults),"FittingResults.txt",sep = "\t")

mean(fittingResults[,19])
mean(fittingResults[,20])

mean(fittingResults[,1]);sd(fittingResults[,1])
mean(fittingResults[,2]);sd(fittingResults[,2])
mean(fittingResults[,3]);sd(fittingResults[,3])
mean(fittingResults[,4]);sd(fittingResults[,4])
mean(fittingResults[,5]);sd(fittingResults[,5])
mean(fittingResults[,6]);(sd(fittingResults[,6]/sqrt(30)))
mean(fittingResults[,7]);(sd(fittingResults[,7]/sqrt(30)))

modelPred = matrix(data = c(mean(fittingResults[,14]),mean(fittingResults[,13]),mean(fittingResults[,12]),mean(fittingResults[,11]),mean(fittingResults[,18]),mean(fittingResults[,17]),mean(fittingResults[,16]) ,mean(fittingResults[,15])),nrow = 4,ncol = 2,byrow = TRUE)

modelPred = data.frame(modelPred)
names(modelPred) = c("Blocked","Interleaved")
modelPred$typeItem = c("New Items","Old Items","New Items","Old Items")
modelPred$catStructure = c("High Similarity","High Similarity","Low Similarity","Low Similarity")

library(reshape2)
modelPred = melt(data = modelPred,id.vars=c("typeItem","catStructure"),variable.name = "sequence",value.name="Performance")

### confidence intervals

n = nrow(fittingResults)
modelSDs = matrix(data = c(sd(fittingResults[,14]),sd(fittingResults[,13]),sd(fittingResults[,12]),sd(fittingResults[,11]),sd(fittingResults[,18]),sd(fittingResults[,17]),sd(fittingResults[,16]) ,sd(fittingResults[,15])),nrow = 4,ncol = 2,byrow = TRUE)

modelError = matrix(data = c((qnorm(0.975)*modelSDs[1,1]/sqrt(n)),(qnorm(0.975)*modelSDs[1,2]/sqrt(n)),(qnorm(0.975)*modelSDs[2,1]/sqrt(n)),(qnorm(0.975)*modelSDs[2,2]/sqrt(n)),(qnorm(0.975)*modelSDs[3,1]/sqrt(n)),(qnorm(0.975)*modelSDs[3,2]/sqrt(n)),(qnorm(0.975)*modelSDs[4,1]/sqrt(n)),(qnorm(0.975)*modelSDs[4,2]/sqrt(n))),nrow = 4,ncol = 2,byrow = TRUE)

modelError = data.frame(modelError)
names(modelError) = c("Blocked","Interleaved")
modelError$typeItem = c("New Items","Old Items","New Items","Old Items")
modelError$catStructure = c("High Similarity","High Similarity","Low Similarity","Low Similarity")

library(reshape2)
modelError = melt(data = modelError,id.vars=c("typeItem","catStructure"),variable.name = "sequence",value.name="ConfInt")

modelPred = merge(modelPred,modelError)

library(plyr)
#sem <- function(x) sd(x)/sqrt(length(x))
conf.int <-function(x) qnorm(0.975)*sd(x)/sqrt(length(x))
d<-read.csv(paste(getwd(),"data_all.csv",sep="/"), header=TRUE, row.names=NULL)
library(plyr)
d = ddply(.data = d,.variables = .(Similarity,TypeItem,Sequence),.fun = summarize,performance=mean(Performance),sem=conf.int(Performance))

empdata = data.frame(d)
names(empdata) = c("catStructure","typeItem","sequence","Performance","ConfInt")
empdata$catStructure = ifelse(empdata$catStructure==0,"High Similarity","Low Similarity")
empdata$typeItem = ifelse(empdata$typeItem==0,"Old Items","New Items")
empdata$sequence = ifelse(empdata$sequence==0.25, "Interleaved","Blocked")
rm(d)

save.image("FitMLE.RData")

library(ggplot2)
library(grid)
library(gridExtra)


plot <- ggplot(empdata, aes(x=catStructure, y=Performance,fill = sequence)) + geom_bar(position="dodge",stat="identity",color="black",alpha=0.4) + coord_cartesian(ylim = c(0,1), expand = TRUE) + geom_errorbar(data = empdata,aes(ymax = Performance + ConfInt, ymin= Performance - ConfInt, width=0.2),position=position_dodge(width = 0.90)) + geom_point(data = modelPred, aes(x=catStructure, y=Performance,group = sequence,colour=sequence),position=position_dodge(width = 0.90),size=4) + geom_errorbar(data = modelPred, aes(ymax = Performance + ConfInt, ymin= Performance - ConfInt, width=0.2,colour=sequence),position=position_dodge(width = 0.90)) + scale_colour_manual(values=c("#e66101","#b2abd2")) + labs(x="Category Structure",y="Probability of Correct Classification") + theme_bw() +  theme(legend.position="right",legend.title = element_blank(),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank()) + scale_fill_manual(values=c("#e66101","#b2abd2")) + theme(legend.position = c(0.1,0.8))
plot + facet_grid(. ~ typeItem,as.table = T)

ggsave(filename = "Fitting results.png",plot = plot + facet_grid(. ~ typeItem,as.table = T),width = 8,height = 5)
