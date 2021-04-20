# model of sequential learning category learning
# this model is a function with several variables passed (most of them parameters of the model)
# this version was updated and optimized for fitting (quicker, more efficient). Carvalho 3/7/2016

model = function(stimuliNumber,typeItem,probability_of_cat_repetition,c,r,weight_shared_same,weight_shared_different,weight_distinctive_same,weight_distinctive_different,gamma,p,L,studysequenceIntHighSim,studysequenceIntLowSim,studysequenceBlkHighSim,studysequenceBlkLowSim,randChoice){
  
  #get the stimuli. if stimuliNumber == 0, highsim, if == 1, lowsim
  if(stimuliNumber==0){
    stimuli = "highsim"
  }else{stimuli="lowsim"}
  file_name<-paste(paste(getwd(),"exemplarsStudy_",sep="/"),stimuli,".csv",sep="") #stimuli used to train the model
  file_name2<-paste(paste(getwd(),"exemplarsNovel_",sep="/"),stimuli,".csv",sep="") #stimuli used to test the model
  features<-as.matrix(read.table(file_name, header = TRUE,sep=',')) #read in data from file with headers for the study items
  
  #get the test items. if typeItem == 0, old items, if == 1 new items.
  if(typeItem == 0){
    featuresTest <- features
  }else{featuresTest<-as.matrix(read.table(file_name2, header = TRUE,sep=','))} #read in data from file with headers for the test items)
  
  
  studyitems = rbind(features,features,features) #repeat the data as many times as needed to simulate one block of study.
  
  randchoice = randChoice
  L = L ## this is the constant for how much frequency matters!
  num_features<-8 #how many features do the stimuli have?
  num_trials <- 72 #how many trials were there?
  
  #these are parameters of the model.
  weight_shared_same<-weight_shared_same
  weight_shared_different<-weight_shared_different
  weight_distinctive_same <-weight_distinctive_same
  weight_distinctive_different<-weight_distinctive_different
  
  #are we running the interleaved or the blocked condition?
  probability_of_cat_repetition <- probability_of_cat_repetition #interleaved = 0.25; blocked =0.75
  if(probability_of_cat_repetition == .25){cond="interleaved"}else{cond="blocked"}
  number_of_examples=nrow(features) #number of stimuli in file
  
  #parameters of GCM
  c<-c # sensitivity parameter of GCM
  r<-r # Minkowski distance formula parameter. r=1 for city-block; r=2 for Euclidean
  
  #functions needed later.
  softmax<-function(activation,power,gamma) #not used, but could later on.  Really, varying C acts very similar to changing gamma in soft-max
    
    #power determines how much more accurate people are for objects similar to both categories than dissimilar from both
  {
    return (exp((activation^power)*gamma))
    #another way of interpreting gamma is that it is the threshold of a non-linear sigmoidal function, with gamma controlling how non-linear
  }
  
  similarity <-function(probe,stored)
  {
    distance<-0
    for (i in 1:num_features)
    {
      if (probe[i]==stored[i])
      {difference<-0}
      else
      {difference<-1}
      distance<-distance+((difference^r)*stored[num_features+1+i]) #multiply by the weight to the feature in stored memory
      #Later on, it might be interesting to also multiply by probed weight of feature, which reflects its diagnosticity
      
    }
    distance<-distance^(1/r)
    #print(exp((-c*distance)^p))
    return(exp(-c*(distance^p))) #returns similarity ranging from 0 to 1
  }
  
  summed_sim_to_cat <- function (probe,categ)
  {
    sim<-0
    right_category_exemplars<-subset(exemplars,exemplars[,9]==categ)
    for (i in 1:nrow(right_category_exemplars))
    {
      sim<-sim+similarity(probe,right_category_exemplars[i,])
    }
    sum(sim)
    return(sim)
  }
  
  probability_categories <-function(probe)
  {
    correct_prob<-0
    categories<-factor(exemplars[,9])
    category_names<-levels(categories)
    sim<-numeric(length(category_names)) #start with a bunch of zeros to initialize similarities to each category
    total<-0
    for (i in 1:length(category_names))
    {
      sim[i]<-summed_sim_to_cat(probe,i)
      total<-total+sim[i]
    }
    for (i in 1:length(category_names))
    {
      #cat("Probability of item belonging to Category ",category_names[i]," is ", (sim[i]^gamma)/total^gamma,"\n")
      if (category_names[i]==probe[9]) {correct_prob<-randchoice*0.5+(1-randchoice)*(sim[i]^gamma)/total^gamma}
    }
    return(correct_prob)
  }
  
  #variable to store the results
  results = c()
  
  #run the model runs times
    ##create a defined order following the repetition probability (and interval around it, defined)
    
  #import the study sequence
  if(stimuliNumber==0  & probability_of_cat_repetition == .25){
    studysequence_items = studysequenceIntHighSim
  }
  
  if(stimuliNumber==0  & probability_of_cat_repetition == .75){
    studysequence_items = studysequenceBlkHighSim
  }
  
  if(stimuliNumber==1  & probability_of_cat_repetition == .75){
    studysequence_items = studysequenceBlkLowSim
  }
  
  if(stimuliNumber==1  & probability_of_cat_repetition == .25){
    studysequence_items = studysequenceIntLowSim
  }

  exemplars <- matrix(0,nrow=num_trials+1, ncol=num_features*2+1)
  exemplars[1,] <- c(rep(0,9),c(rep(1/num_features,num_features)))
  #names(exemplars) <- c(paste("F",1:num_features,sep=""), paste("W",1:num_features,sep=""),"category")
  
  FreqFeatures = as.character(levels(factor(c(unique(studysequence_items[,1:8])))))
  Freqcats = c(rep(unique(studysequence_items[,9])[1],length(FreqFeatures)),rep(unique(studysequence_items[,9])[2],length(FreqFeatures)),rep(unique(studysequence_items[,9])[3],length(FreqFeatures)))
  
  FreqTable = cbind(matrix(0,nrow = length(FreqFeatures)*3,ncol=num_features+2))
  FreqTable[,1] = Freqcats
  FreqTable[,2] = FreqFeatures
  
  previous_exemplar<- 1
  
  for (j in 1:nrow(studysequence_items)){
    sampled_exemplar <- j
    same_category<-(studysequence_items[sampled_exemplar,9]==exemplars[previous_exemplar,9]) #do two examples share category label?
    for(c in 1:ncol(studysequence_items)-1){
      tally = as.numeric(FreqTable[,c+2][FreqTable[,1] == studysequence_items[sampled_exemplar,9] & FreqTable[,2] == studysequence_items[sampled_exemplar,c]])
      FreqTable[,c+2][FreqTable[,1] == studysequence_items[sampled_exemplar,9] & FreqTable[,2] == studysequence_items[sampled_exemplar,c]] = tally + 1
    }
    for (i in 1:num_features){
      exemplars[previous_exemplar+1,i]=studysequence_items[sampled_exemplar,i] #copy over features from sampled example to memory
      freqVal = as.numeric(FreqTable[,i+2][FreqTable[,1]==studysequence_items[sampled_exemplar,9]&FreqTable[,2]==studysequence_items[sampled_exemplar,i]])/sum(as.numeric(FreqTable[,i+2][FreqTable[,1]==studysequence_items[sampled_exemplar,9]]))
      if (studysequence_items[sampled_exemplar,9]==exemplars[previous_exemplar,9]) #if the previous and current examples are in the same category
      {if (exemplars[previous_exemplar+1,i] == exemplars[previous_exemplar,i])
      {exemplars[previous_exemplar+1,num_features+1+i]<-weight_shared_same+L*freqVal} #+i to add weights to second set of columns
        else
        {exemplars[previous_exemplar+1,num_features+1+i]<-weight_distinctive_same+L*freqVal}}
      else
      {if (exemplars[previous_exemplar+1,i] == exemplars[previous_exemplar,i]) 
      {exemplars[previous_exemplar+1,num_features+1+i]<-weight_shared_different+L*freqVal} #+i to add weights to second set of columns
        else
        {exemplars[previous_exemplar+1,num_features+1+i]<-weight_distinctive_different+L*freqVal}}
    }
    
    
    sum_of_weights<-0 #prepare to normalize weights so that they add up to 1
    #making the sum of weights add to 1 is critical because it reward shared/distinctive features that STAND OUT
    for (i in (num_features+2):(2*num_features+1)) #to target only the weights - the second set of items in exemplars
    {
      sum_of_weights<-sum_of_weights+exemplars[previous_exemplar+1,i]
    }
    for (i in (num_features+2):(2*num_features+1)) #to target only the weights - the second set of items in exemplars
    {
      exemplars[previous_exemplar+1,i]<-exemplars[previous_exemplar+1,i]/sum_of_weights #give each weight it scaled value
    }
    
    exemplars[previous_exemplar+1,9]<-studysequence_items[sampled_exemplar,9] #copy appropriate category for exemplar
    previous_exemplar<-previous_exemplar+1} #increment previous_exemplar
  exemplars<-exemplars[-1,] #this will remove the first row of exemplars.  This was just a starting off point, not a real exemplar.
  
    #run the test items classifications and calculate probability of correct categorization for each item.
    average_pct<-0
    for (i in 1:nrow(featuresTest))
    {
      average_pct<-average_pct+probability_categories(featuresTest[i,])
    }
    average_pct<-average_pct/nrow(featuresTest)
    results = average_pct
  
  return(results)
  #return(list(d,results)) #return the results of the simulations, the condition ran and the stimuli type used.
}