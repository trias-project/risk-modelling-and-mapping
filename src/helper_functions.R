#This function calculates the number of decimal places in any given numeric value 
# eg., 15.21 has 2 decimal places, 15.2569 has 4 decimal places, 15.25690 also has 4, as 0 in the end doesn't count
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

#Divide a numerical value by 10
divide10<-function(x){
  value<-x/10
  return(value)
}

#Function to return threshold where sens=spec from caret results 
findThresh<-function(df){
  df<-df[c("rowIndex","obs","present")]
  df<-df %>%
    mutate(observed= ifelse(obs == "present",1,0)) %>%
    select(rowIndex,observed,predicted=present)
  result<-PresenceAbsence::optimal.thresholds(df,opt.methods = 2)
  return(result)
}

#Recalculate accuracy for a given model with the threshold that has been optimized
#using fingThresh
accuracyStats<-function(df,y){
  df<-df[c("rowIndex","obs","present")]
  df<-df %>%
    mutate(observed= ifelse(obs == "present",1,0)) %>%
    select(rowIndex,observed,predicted=present)
  result<-PresenceAbsence::presence.absence.accuracy(df,threshold = y,st.dev=FALSE)
  return(result)
}

