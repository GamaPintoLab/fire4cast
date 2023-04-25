setwd("/Users/franciscopinto/Library/CloudStorage/OneDrive-UniversidadedeLisboa/FRANCISCO/work/working_papers/fire4cast")

#import climate data

library(readxl)
vn_qn_clim <- read_excel("vn_qn_clim.xlsx",
col_types = c("date", "numeric", "numeric",
"numeric", "numeric", "numeric",
"numeric", "numeric", "numeric",
"numeric", "numeric", "text", "numeric",
"numeric", "numeric", "numeric",
"numeric", "numeric", "numeric"))

#import dates of phenotype change 

feno2r <- read_excel("feno2r.xlsx", col_types = c("text",
"text", "date", "numeric", "numeric",
"numeric"))

# import bacterial detection data

sampling2r <- read_excel("sampling2r.xlsx",
col_types = c("text", "text", "date",
"text", "numeric", "numeric", "numeric",
"numeric", "numeric", "numeric",
"numeric", "numeric", "numeric",
"numeric"))

#function to compute day of the year (1 to 365) from day, month, year

yearday=function(dd,mm,yy){
  diasmes=c(31,28,31,30,31,30,31,31,30,31,30,30)
  diasmes2020=c(31,29,31,30,31,30,31,31,30,31,30,30)
  cumdias=c(0,cumsum(diasmes))
  cumdias2020=c(0,cumsum(diasmes2020))
  if (yy==2020) {
    cd=cumdias2020
  } else {
    cd=cumdias
  }
  yd=cd[mm]+dd
  yd
}

#add day of the year to bacterial detection data frame

sampling=as.data.frame(sampling2r)
sampling$yday=0

for (i in 1:nrow(sampling)){
  sampling$yday[i]=yearday(sampling[i,12],sampling[i,13],sampling[i,14])
}

# add day of the year to phenotype dates data frame

feno=as.data.frame(feno2r)
feno$yday=0
for (i in 1:nrow(feno)){
  feno$yday[i]=yearday(feno$day[i],feno$month[i],feno$year[i])
}

#add day of the year to climate data frame

clim=as.data.frame(vn_qn_clim)
clim$yday=0
for (i in 1:nrow(clim)){
  clim$yday[i]=yearday(clim$day[i],clim$month[i],clim$year[i])
}

#missing data imputation

lines2input=(103:120) #lines were all climate variables are missing in vn location

otherlines=setdiff((1:nrow(clim)),lines2input)

# for each date with missing vn data, select the two dates where the qn climate data is more similar to the qn data inthe missing date.
# imput the missing vn data with the average of the vn data in the two most similar days according to qn data 

newclim=clim
for (i in 1:length(lines2input)){
  baitline=clim[lines2input[i],13:19]
  d2bait=vector()
  for (j in 1:nrow(clim)){
    d2bait[j]=sum((clim[j,13:19]-baitline)^2)
  }
  d2bait[lines2input]=max(d2bait)
  dord=order(d2bait)
  modellines=clim[dord[1:2],2:8]
  newclim[lines2input[i],2:8]=colMeans(modellines)
  
}

# data imputation for individual windmed missing values using the average of the two most similar dates arrording to the other climate variables 

lines2input=which(is.na(newclim$windmedvn))
subclim=newclim[,c(2,3,4,5,7,8)]
for (i in 1:length(lines2input)){
  baitline=subclim[lines2input[i],]
  d2bait=vector()
  for (j in 1:nrow(newclim)){
    d2bait[j]=sum((subclim[j,]-baitline)^2)
  }
  d2bait[lines2input]=max(d2bait)
  dord=order(d2bait)
  modellines=newclim[dord[1:2],6]
  newclim[lines2input[i],6]=mean(modellines)
  
}

save(newclim,file="newclim.RData")

load("newclim.RData")

#convert qPCRq to dCq, so that it is positively correlated with bacterial load

sampling$dCq=40-sampling$qPCRCq

#logarithm transformation of flow cytometry estimates, with (_b) and without taking into acount limit of detection

a=sampling$IFCMTotal

a[is.na(a)]=1
sampling$logTotal=log10(a)
a[a<1.24*(10^5)]=1
sampling$logTotal_b=log10(a)

a=sampling$IFCMLive
a[is.na(a)]=1
sampling$logLive=log10(a)
a[a<1.24*(10^5)]=1
sampling$logLive_b=log10(a)

#separate data frames according to location

loc=unique(sampling$Location)
 
loc1=sampling[sampling$Location==loc[1],]
loc2=sampling[sampling$Location==loc[2],]

#sets of unique sampling dates

sdates1=unique(loc1$Date)
sdates2=unique(loc2$Date)

### generate data frames with one row per sampling date
# bacterial data is the sum of all individual trees in each sampling date

i=1
d=loc1[loc1$Date==sdates1[i],]
l=d[1,]
l$dCq=sum(d$dCq)
l$logTotal=sum(d$logTotal)
l$logTotal_b=sum(d$logTotal_b)
l$logLive=sum(d$logLive)
l$logLive_b=sum(d$logLive_b)
loc1sum=l

for (i in 2:length(sdates1)){
  d=loc1[loc1$Date==sdates1[i],]
  l=d[1,]
  l$dCq=sum(d$dCq)
  l$logTotal=sum(d$logTotal)
  l$logTotal_b=sum(d$logTotal_b)
  l$logLive=sum(d$logLive)
  l$logLive_b=sum(d$logLive_b)
  loc1sum=rbind(loc1sum,l)
}


i=1
d=loc2[loc2$Date==sdates2[i],]
l=d[1,]
l$dCq=sum(d$dCq)
l$logTotal=sum(d$logTotal)
l$logTotal_b=sum(d$logTotal_b)
l$logLive=sum(d$logLive)
l$logLive_b=sum(d$logLive_b)
loc2sum=l

for (i in 2:length(sdates2)){
  d=loc2[loc2$Date==sdates2[i],]
  l=d[1,]
  l$dCq=sum(d$dCq)
  l$logTotal=sum(d$logTotal)
  l$logTotal_b=sum(d$logTotal_b)
  l$logLive=sum(d$logLive)
  l$logLive_b=sum(d$logLive_b)
  loc2sum=rbind(loc2sum,l)
}

#### same as above but values per sampling date are counts across all trees, instead of sums

i=1
d=loc1[loc1$Date==sdates1[i],]
l=d[1,]
l$dCq=sum(d$dCq>0)
l$logTotal=sum(d$logTotal>0)
l$logTotal_b=sum(d$logTotal_b>0)
l$logLive=sum(d$logLive>0)
l$logLive_b=sum(d$logLive_b>0)
loc1count=l

for (i in 2:length(sdates1)){
  d=loc1[loc1$Date==sdates1[i],]
  l=d[1,]
  l$dCq=sum(d$dCq>0)
  l$logTotal=sum(d$logTotal>0)
  l$logTotal_b=sum(d$logTotal_b>0)
  l$logLive=sum(d$logLive>0)
  l$logLive_b=sum(d$logLive_b>0)
  loc1count=rbind(loc1count,l)
}


i=1
d=loc2[loc2$Date==sdates2[i],]
l=d[1,]
l$dCq=sum(d$dCq>0)
l$logTotal=sum(d$logTotal>0)
l$logTotal_b=sum(d$logTotal_b>0)
l$logLive=sum(d$logLive>0)
l$logLive_b=sum(d$logLive_b>0)
loc2count=l

for (i in 2:length(sdates2)){
  d=loc2[loc2$Date==sdates2[i],]
  l=d[1,]
  l$dCq=sum(d$dCq>0)
  l$logTotal=sum(d$logTotal>0)
  l$logTotal_b=sum(d$logTotal_b>0)
  l$logLive=sum(d$logLive>0)
  l$logLive_b=sum(d$logLive_b>0)
  loc2count=rbind(loc2count,l)
}


### create vector of phenotype states per sampling dates in both locations

df_1=vector()
df2_1=vector()
dh_1=vector()
fenostate_1=vector()

for (i in 1:length(sdates1)){
  yy=loc1sum$year[i]
  fday=feno$yday[feno$year==yy & feno$local=="QN" & feno$floracao=="F"]
  f2day=feno$yday[feno$year==yy & feno$local=="QN" & feno$floracao=="F2"]
  hday=feno$yday[feno$year==yy & feno$local=="QN" & feno$floracao=="H"]
  df_1[i]=loc1sum$yday[i]-fday
  df2_1[i]=loc1sum$yday[i]-f2day
  dh_1[i]=loc1sum$yday[i]-hday
  if (df_1[i]<0){
    fenostate_1[i]=0
  } else {
    if (df2_1[i]<0){
      fenostate_1[i]=1
    } else {
      if (dh_1[i]<0){
        fenostate_1[i]=2
      } else {
        fenostate_1[i]=3
      }
    }
  }
  
}



df_2=vector()
df2_2=vector()
dh_2=vector()
fenostate_2=vector()

for (i in 1:length(sdates2)){
  yy=loc2sum$year[i]
  fday=feno$yday[feno$year==yy & feno$local=="VN" & feno$floracao=="F"]
  f2day=feno$yday[feno$year==yy & feno$local=="VN" & feno$floracao=="F2"]
  hday=feno$yday[feno$year==yy & feno$local=="VN" & feno$floracao=="H"]
  df_2[i]=loc2sum$yday[i]-fday
  df2_2[i]=loc2sum$yday[i]-f2day
  dh_2[i]=loc2sum$yday[i]-hday
  if (df_2[i]<0){
    fenostate_2[i]=0
  } else {
    if (df2_2[i]<0){
      fenostate_2[i]=1
    } else {
      if (dh_2[i]<0){
        fenostate_2[i]=2
      } else {
        fenostate_2[i]=3
      }
    }
  }
  
}

# create climatic data frames with variables accumulated after n days (n = 1, 4 or 7)

clim1=newclim[,c(12,9,10,11,20,13,14,15,16,17,18,19)]
clim2=newclim[,c(12,9,10,11,20,2,3,4,5,6,7,8)]
names(clim1)= c("date", "day","month","year","yday","tmed","tmax","tmin","hr","wind", "rain","rad")
names(clim2)= c("date", "day","month","year","yday","tmed","tmax","tmin","hr","wind", "rain","rad")


climacum=function(locvar,climvar,acum){
  i=1
  yy=locvar$year[i]
  yday=locvar$yday[i]
  cc=climvar[climvar$year==yy,]
  cc=cc[(cc$yday<=yday),]
  cc=cc[(cc$yday>(yday-acum)),]
  if (nrow(cc)>1){
    l=cc[cc$yday==yday,]
    l[1,6:12]=colSums(cc[,6:12])  
  } else {
    if (nrow(cc)==1){
      l=cc
    }
  }
  
  cvar=l
  for (i in 2:nrow(locvar)){
    yy=locvar$year[i]
    yday=locvar$yday[i]
    cc=climvar[climvar$year==yy,]
    cc=cc[(cc$yday<=yday),]
    cc=cc[(cc$yday>(yday-acum)),]
    if (nrow(cc)>1){
      l=cc[cc$yday==yday,]
      l[1,6:12]=colSums(cc[,6:12]) 
      cvar=rbind(cvar,l)
    } else {
      if (nrow(cc)==1){
        l=cc
        cvar=rbind(cvar,l)
      }
    }
    
  }
  cvar
}


clim1_1=climacum(loc1sum,clim1,1)
clim1_4=climacum(loc1sum,clim1,4)
clim1_7=climacum(loc1sum,clim1,7)

clim2_1=climacum(loc2sum,clim2,1)
clim2_4=climacum(loc2sum,clim2,4)
clim2_7=climacum(loc2sum,clim2,7)



#data frames to model

a=names(loc1count)
a[16:20]=c("n_dCq","nTot","nTot_b","nLive","nLive_b")
names(loc1count)=a
df1_noclim=cbind(loc1sum[,c(2,(12:20))],loc1count[,16:20])
df1_lag1=cbind(df1_noclim,clim1_1[,6:12])
df1_lag1$feno=fenostate_1
df1_lag1$f=df_1
df1_lag1$f2=df2_1
df1_lag1$h=dh_1

df1_lag4=cbind(df1_noclim,clim1_4[,6:12])
df1_lag4$feno=fenostate_1
df1_lag4$f=df_1
df1_lag4$f2=df2_1
df1_lag4$h=dh_1

df1_lag7=cbind(df1_noclim,clim1_7[,6:12])
df1_lag7$feno=fenostate_1
df1_lag7$f=df_1
df1_lag7$f2=df2_1
df1_lag7$h=dh_1

names(loc2count)=a
df2_noclim=cbind(loc2sum[,c(2,(12:20))],loc2count[,16:20])
df2_lag1=cbind(df2_noclim,clim2_1[,6:12])
df2_lag1$feno=fenostate_2
df2_lag1$f=df_2
df2_lag1$f2=df2_2
df2_lag1$h=dh_2

df2_lag4=cbind(df2_noclim,clim2_4[,6:12])
df2_lag4$feno=fenostate_2
df2_lag4$f=df_2
df2_lag4$f2=df2_2
df2_lag4$h=dh_2

df2_lag7=cbind(df2_noclim,clim2_7[,6:12])
df2_lag7$feno=fenostate_2
df2_lag7$f=df_2
df2_lag7$f2=df2_2
df2_lag7$h=dh_2


df_lag1=rbind(df1_lag1,df2_lag1)
df_lag4=rbind(df1_lag4,df2_lag4)
df_lag7=rbind(df1_lag7,df2_lag7)

alarms <- read.delim("alarms.txt", quote="")

alldf=cbind(df_lag4,alarms)

write.csv(alldf,file="alldf.csv",row.names = F)

#functions to choose best thresholds for a pair of variables
threshand=function(x,y,al,b=1,pt=0.8){
  xval=sort(unique(x))
  yval=sort(unique(y))
  outtab=data.frame(xt=rep(0,length(xval)*length(yval)),yt=0,p=0,r=0,f=0,np=0)
  lin=1
  for (i in 1:length(xval)){
    for (j in 1:length(yval)){
      thpos=thresandeval(x,y,xval[i],yval[j])
      tp=sum(thpos*al)
      fp=sum(thpos*(al==0))
      tn=sum((thpos==0)*(al==0))
      fn=sum((thpos==0)*al)
      prec=tp/(tp+fp)
      rec=tp/(tp+fn)
      fm=(1+b^2)*prec*rec/((b^2)*rec+prec)
      outtab[lin,]=c(xval[i],yval[j],prec,rec,fm,tp+fp)
      lin=lin+1
    }
  }
  havena=rowSums(is.na(outtab))
  outtab=outtab[havena==0,]
  if (max(outtab[,3],na.rm=T)<pt){
    bestthresh=c(xt=0,yt=0,p=0,r=0,f=0,np=0)
    outclass=rep(0,length(al))
  } else {
    highp=outtab[outtab$p>=pt,]
    maxnp=max(highp$np,na.rm=T)
    highnp=highp[highp$np==maxnp,]
    if (nrow(highnp)==1){
      bestthresh=highnp[1,]
      outclass=(x>=highnp[1,1])*(y>=highnp[1,2])
    } else {
      unix=unique(highnp$xt)
      uniy=unique(highnp$yt)
      if (length(unix)<length(uniy)){
        highnp=highnp[order(highnp$yt),]
      } else {
        highnp=highnp[order(highnp$xt),]
      }
      bestthresh=highnp[floor(nrow(highnp)/2),]
      outclass=(x>=bestthresh[1,1])*(y>=bestthresh[1,2])
    }
  }
  
  list(tab=outtab,bestt=bestthresh,outclass=outclass)
}
  
thresandeval=function(x,y,xt,yt){
  outvec=(x>=xt)*(y>=yt)
  outvec
}


# determine thresholds for the 11 pairs of predictor variables

predmat1=matrix(0, nrow=40, ncol=11)
predmat4=matrix(0, nrow=40, ncol=11)
predmat7=matrix(0, nrow=40, ncol=11)

rulemat1=data.frame(x="-",y="-",xt=rep(0,11),yt=0,p=0,r=0,f=0,np=0)
rulemat4=data.frame(x="-",y="-",xt=rep(0,11),yt=0,p=0,r=0,f=0,np=0)
rulemat7=data.frame(x="-",y="-",xt=rep(0,11),yt=0,p=0,r=0,f=0,np=0)

  
tmed1andhr=threshand(df_lag1$tmed,df_lag1$hr,alarms$alarm_n,b=1,pt=0.75)
tmed4andhr=threshand(df_lag4$tmed,df_lag4$hr,alarms$alarm_n,b=1,pt=0.75)
tmed7andhr=threshand(df_lag7$tmed,df_lag7$hr,alarms$alarm_n,b=1,pt=0.75)

predmat1[,1]=tmed1andhr$outclass
predmat4[,1]=tmed4andhr$outclass
predmat7[,1]=tmed7andhr$outclass

rulemat4[1,1]="tmed"
rulemat4[1,2]="hr"
rulemat4[1,3:8]=tmed4andhr$bestt

tmed1andrain=threshand(df_lag1$tmed,df_lag1$rain,alarms$alarm_n,b=1,pt=0.75)
tmed4andrain=threshand(df_lag4$tmed,df_lag4$rain,alarms$alarm_n,b=1,pt=0.75)
tmed7andrain=threshand(df_lag7$tmed,df_lag7$rain,alarms$alarm_n,b=1,pt=0.75)


predmat1[,2]=tmed1andrain$outclass
predmat4[,2]=tmed4andrain$outclass
predmat7[,2]=tmed7andrain$outclass

rulemat4[2,1]="tmed"
rulemat4[2,2]="rain"
rulemat4[2,3:8]=tmed4andrain$bestt

tmed1andwind=threshand(df_lag1$tmed,df_lag1$wind,alarms$alarm_n,b=1,pt=0.75)
tmed4andwind=threshand(df_lag4$tmed,df_lag4$wind,alarms$alarm_n,b=1,pt=0.75)
tmed7andwind=threshand(df_lag7$tmed,df_lag7$wind,alarms$alarm_n,b=1,pt=0.75)

predmat1[,3]=tmed1andwind$outclass
predmat4[,3]=tmed4andwind$outclass
predmat7[,3]=tmed7andwind$outclass

rulemat4[3,1]="tmed"
rulemat4[3,2]="wind"
rulemat4[3,3:8]=tmed4andwind$bestt


tmed1andrad=threshand(df_lag1$tmed,df_lag1$rad,alarms$alarm_n,b=1,pt=0.75)
tmed4andrad=threshand(df_lag4$tmed,df_lag4$rad,alarms$alarm_n,b=1,pt=0.75)
tmed7andrad=threshand(df_lag7$tmed,df_lag7$rad,alarms$alarm_n,b=1,pt=0.75)

predmat1[,4]=tmed1andrad$outclass
predmat4[,4]=tmed4andrad$outclass
predmat7[,4]=tmed7andrad$outclass

rulemat4[4,1]="tmed"
rulemat4[4,2]="rad"
rulemat4[4,3:8]=tmed4andrad$bestt

thr1andrain=threshand(df_lag1$hr,df_lag1$rain,alarms$alarm_n,b=1,pt=0.75)
thr4andrain=threshand(df_lag4$hr,df_lag4$rain,alarms$alarm_n,b=1,pt=0.75)
thr7andrain=threshand(df_lag7$hr,df_lag7$rain,alarms$alarm_n,b=1,pt=0.75)

predmat1[,5]=thr1andrain$outclass
predmat4[,5]=thr4andrain$outclass
predmat7[,5]=thr7andrain$outclass

rulemat4[5,1]="hr"
rulemat4[5,2]="rain"
rulemat4[5,3:8]=thr4andrain$bestt

thr1andwind=threshand(df_lag1$hr,df_lag1$wind,alarms$alarm_n,b=1,pt=0.75)
thr4andwind=threshand(df_lag4$hr,df_lag4$wind,alarms$alarm_n,b=1,pt=0.75)
thr7andwind=threshand(df_lag7$hr,df_lag7$wind,alarms$alarm_n,b=1,pt=0.75)

predmat1[,6]=thr1andwind$outclass
predmat4[,6]=thr4andwind$outclass
predmat7[,6]=thr7andwind$outclass

rulemat4[6,1]="hr"
rulemat4[6,2]="wind"
rulemat4[6,3:8]=thr4andwind$bestt

thr1andrad=threshand(df_lag1$hr,df_lag1$rad,alarms$alarm_n,b=1,pt=0.75)
thr4andrad=threshand(df_lag4$hr,df_lag4$rad,alarms$alarm_n,b=1,pt=0.75)
thr7andrad=threshand(df_lag7$hr,df_lag7$rad,alarms$alarm_n,b=1,pt=0.75)

predmat1[,7]=thr1andrad$outclass
predmat4[,7]=thr4andrad$outclass
predmat7[,7]=thr7andrad$outclass

rulemat4[7,1]="hr"
rulemat4[7,2]="rad"
rulemat4[7,3:8]=thr4andrad$bestt

twind1andrain=threshand(df_lag1$wind,df_lag1$rain,alarms$alarm_n,b=1,pt=0.75)
twind4andrain=threshand(df_lag4$wind,df_lag4$rain,alarms$alarm_n,b=1,pt=0.75)
twind7andrain=threshand(df_lag7$wind,df_lag7$rain,alarms$alarm_n,b=1,pt=0.75)

predmat1[,8]=twind1andrain$outclass
predmat4[,8]=twind4andrain$outclass
predmat7[,8]=twind7andrain$outclass

rulemat4[8,1]="wind"
rulemat4[8,2]="rain"
rulemat4[8,3:8]=twind4andrain$bestt


twind1andrad=threshand(df_lag1$wind,df_lag1$rad,alarms$alarm_n,b=1,pt=0.75)
twind4andrad=threshand(df_lag4$wind,df_lag4$rad,alarms$alarm_n,b=1,pt=0.75)
twind7andrad=threshand(df_lag7$wind,df_lag7$rad,alarms$alarm_n,b=1,pt=0.75)

predmat1[,9]=twind1andrad$outclass
predmat4[,9]=twind4andrad$outclass
predmat7[,9]=twind7andrad$outclass

rulemat4[9,1]="wind"
rulemat4[9,2]="rad"
rulemat4[9,3:8]=twind4andrad$bestt


train1andrad=threshand(df_lag1$rain,df_lag1$rad,alarms$alarm_n,b=1,pt=0.75)
train4andrad=threshand(df_lag4$rain,df_lag4$rad,alarms$alarm_n,b=1,pt=0.75)
train7andrad=threshand(df_lag7$rain,df_lag7$rad,alarms$alarm_n,b=1,pt=0.75)

predmat1[,10]=train1andrad$outclass
predmat4[,10]=train4andrad$outclass
predmat7[,10]=train7andrad$outclass

rulemat4[10,1]="rain"
rulemat4[10,2]="rad"
rulemat4[10,3:8]=train4andrad$bestt


tmin1andmax=threshand(df_lag1$tmin,df_lag1$tmax,alarms$alarm_n,b=1,pt=0.75)
tmin4andmax=threshand(df_lag4$tmin,df_lag4$tmax,alarms$alarm_n,b=1,pt=0.75)
tmin7andmax=threshand(df_lag7$tmin,df_lag7$tmax,alarms$alarm_n,b=1,pt=0.75)

predmat1[,11]=tmin1andmax$outclass
predmat4[,11]=tmin4andmax$outclass
predmat7[,11]=tmin7andmax$outclass

rulemat4[11,1]="tmin"
rulemat4[11,2]="tmax"
rulemat4[11,3:8]=tmin4andmax$bestt

#compute global precision and recall when all rules are applied

rp=sum(alarms$alarm_n)

tp1=sum((rowSums(predmat1)>0)*alarms$alarm_n)
pp1=sum((rowSums(predmat1)>0))
prec1=tp1/pp1
rec1=tp1/rp

tp4=sum((rowSums(predmat4)>0)*alarms$alarm_n)
pp4=sum((rowSums(predmat4)>0))
prec4=tp4/pp4
rec4=tp4/rp

tp7=sum((rowSums(predmat7)>0)*alarms$alarm_n)
pp7=sum((rowSums(predmat7)>0))
prec7=tp7/pp7
rec7=tp7/rp

colSums(predmat1)
colSums(predmat4)
colSums(predmat7)

#function to compute classification performance measures
classperf=function(thpos,al,b=1){
  tp=sum(thpos*al)
  fp=sum(thpos*(al==0))
  tn=sum((thpos==0)*(al==0))
  fn=sum((thpos==0)*al)
  prec=tp/(tp+fp)
  rec=tp/(tp+fn)
  fm=(1+b^2)*prec*rec/((b^2)*rec+prec)
  out=c(tp=tp,fp=fp,tn=tn,fn=fn,prec=prec,rec=rec,fm=fm)
  out
}


# compute expected distribution of performance measures with random classifiers

randperf=data.frame(tp=rep(0,50000),fp=0,tn=0,fn=0,prec=0,rec=0,fm=0)
for (i in 1:50000){
  randpred=sample(40,20)
  rpvec=rep(0,40)
  rpvec[randpred]=1
  randperf[i,]=classperf(rpvec,alarms$alarm_n)
}

quantile(randperf$prec,probs=c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(randperf$rec,probs=c(0,0.025,0.25,0.5,0.75,0.975,1))

#compare best classifier with random performance

p_p=sum(randperf$prec>=prec4)/50000
p_r=sum(randperf$rec>=prec4)/50000

#write summary information for best classifier

alldf=cbind(df_lag4,alarms)
alldf=cbind(alldf,predmat4)
alldf$predtot=rowSums(predmat4)

write.csv(alldf,file="alldf.csv",row.names = F)


write.csv(rulemat4,file="rulemat4.csv",row.names = F)

