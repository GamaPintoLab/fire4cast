#load classifier rules

rulemat=read.csv("rulemat4.csv",header=T,stringsAsFactors = F)


#load climate data for prediction
#csv file with column headers:
# date tmed tmax tmin hr wind rain rad

clim=read.csv("clim_example.csv",header=T,stringsAsFactors = F)


#function to make the predictions based on the rule set
rulepred=function(rules,clim){
  varnames=names(clim)
  rules=rules[rules$xt*rulest$yt!=0,]
  predmat=matrix(data=0,nrow=nrow(clim),ncol=nrow(rules)+1)
  for (i in 1:nrow(rules)){
    v1=which(varnames==rules$x[i])
    v2=which(varnames==rules$y[i])
    predmat[,i]=1*(clim[,v1]>=rules$xt[i])*(clim[,v2]>=rules$yt[i])
  }
  predmat[,nrow(rules)+1]=rowSums(predmat[,1:nrow(rules)])
  predmat
}

#make the predictions
predmat=rulepred(rulemat,clim)

#the output predmat is a matrix with the predictions
#each line corresponds to a date (according to the dates in clim)
# it has one column for each rule and one final column with the sum 
# of the predictions from individual rules, that is, the last column
# shows how many rules are positive/true in the corresponding line/date

# we considered that at least one positive rule is sufficient for a 
# warning of higher risk of transmission and infection, but
# higher numbers of positive rules appear to correlate with more 
# precise predicions