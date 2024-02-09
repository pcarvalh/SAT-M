
createStudySequence = function(probability_of_cat_repetition,stimuliNumber,num_trials){

  # five families chosen randomly per block. needs to be changed below.
  
  if(stimuliNumber==0){
    stimuli = "highsim"
  }else{stimuli="lowsim"}
  
  probability_of_cat_repetition <- probability_of_cat_repetition #interleaved = 0.25; blocked =0.75 

  #this needs three values instead of only one .60 probabily of change. there are two "blocked/interleaved" versions: by letter and by font. then there is fully blocked.

  if(probability_of_cat_repetition == .25){cond="interleaved"}else{cond="blocked"}
  
  file_name<-paste(paste(getwd(),"exemplarsStudy_",sep="/"),stimuli,".csv",sep="") #stimuli used to train the model
  features<-as.matrix(read.table(file_name, header = TRUE,sep=',')) #read in data from file with headers for the study items
  studyitems = rbind(features,features,features) #repeat the data as many times as needed to simulate one block of study.
  
  num_trials <- num_trials
  
  #studysequence_items = studyitems
  studysequence_items = matrix(NA, nrow = num_trials, ncol = 9)
  #names(studysequence_items) = c(paste("F",1:8,sep=""),"category")
  success = 0
  
  while(success != 1){
    rep_array = c() #array to put the rep info to check the sequence created at the end.
    itemstostudy = as.matrix(studyitems)
    
    #create an initial sequence by randomly picking from the table with ALL the stimuli one by one
    for(i in 1:num_trials){
      if(i == 1){
        item = round(runif(1,1,nrow(itemstostudy)),0)
        studysequence_items[1,1:9] = itemstostudy[item,1:9]
        itemstostudy = itemstostudy[-item,]
      }
      if(i > 1 & i < max(num_trials)){
        continue = 0
        while(continue != 1){
          item = round(runif(1,1,nrow(itemstostudy)),0)
          #what do we want?
          repeatcat = runif(1) < probability_of_cat_repetition
          same_category = (itemstostudy[item,9]==studysequence_items[i-1,9])
          if(repeatcat == same_category){
            studysequence_items[i,1:9] = itemstostudy[item,1:9]
            itemstostudy = itemstostudy[-item,]
            continue = 1 #only continue if the criterion is satisfied.
            rm(repeatcat,same_category)
          }
        }
      }
      if(i == max(num_trials)){
        studysequence_items[i,1:9] = itemstostudy
        rm(itemstostudy)
      }
    }
    
    #let's check if the final sequence actually has a probability of repetition that falls within a defined interval around the theoritical one defined.
    for(i in 1:nrow(studysequence_items)){
      if(i == 1){
        next
      }
      if(i > 1){
        if(studysequence_items[i,9]==studysequence_items[i-1,9]){
          rep_array[i-1] = 1
        }else{
          rep_array[i-1] = 0
        }
      }
    }
    
    success_min = probability_of_cat_repetition - 0.05
    success_max = probability_of_cat_repetition + 0.05
    
    if(mean(rep_array) > success_min & mean(rep_array) < success_max){
      success = 1 #the sequence created follows the probabilities! success!
    }
  }
  return(studysequence_items)
}

