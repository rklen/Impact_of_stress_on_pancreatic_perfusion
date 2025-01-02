library(readxl)
linearInter<-function(x1,x,y){
  if(x1==x[1]){return(y[1])}
  for(j in 1:(length(x)-1)){
    if(x1>x[j] & x1<=x[j+1]){
      return(y[j]+(x1-x[j])/(x[j+1]-x[j])*(y[j+1]-y[j]))
    }
  }
  if(x1>x[length(x)]){return(y[length(x)])}
}
fit_model<-function(param,input,j){
  k1<-abs(param[1])
  k2<-abs(param[2])
  Vb<-1/(1+exp(-param[3]))
  cT<-0
  for(i in 2:length(input)){
    if(i-j-1<1){c0<-0}else{
      if(i-j-1>length(input)){c0<-input[length(input)]}else{c0<-input[i-j-1]}
    }
    cT<-c(cT,k1*c0+(1-k2)*cT[i-1])
  }
  v<-cT+Vb*input
  return(v)
}
errorfunction<-function(param,input,organCurve,j){
  return(sum((organCurve-fit_model(param,input,j))^2))
}
meanRelError<-function(modelCurve,organCurve){
  r<-0
  for(i in 2:length(organCurve)){
    if(modelCurve[i]!=organCurve[i]){
      r<-r+abs(modelCurve[i]-organCurve[i])/organCurve[i]
    }
  }
  r<-r/(length(organCurve)-1)
  return(r)
}
optimizationAlg<-function(initialValues,input,organCurve,j){
  tryCatch(
    {
      res<-nlm(errorfunction,initialValues,input=input,
               organCurve=organCurve,j=j,stepmax=1000,print.level=1)
      return(res$estimate)
    },
    error=function(e){
      return(initialValues)
    }
  )
}
colMedians<-function(d){
  medians<-c()
  for(i in 1:ncol(d)){
    medians[i]<-median(d[,i])
  }
  return(medians)
}
colSds<-function(d){
  sds<-c()
  for(i in 1:ncol(d)){
    sds[i]<-sd(d[,i])
  }
  return(sds)
}
wilcoxSymbol<-function(x,y){
  s<-''
  p<-wilcox.test(x,y,alternative='less',paired=TRUE)$p.value
  if(p<0.05){s<-'*'}
  if(p<0.01){s<-'**'}
  if(p<0.001){s<-'***'}
  return(s)
}
findInput<-function(studyNumber,stress,times){
  df1<-read_excel('E:/koveri/uusi_metadata_aivot.xlsx')
  codesInExcel<-df1$study_code
  acCodesInExcel<-df1$ac
  for(l in 1:length(codesInExcel)){
    if(!is.na(codesInExcel[l])){
      if(studyNumber==codesInExcel[l]){
        acCode<-acCodesInExcel[l]
      }
    }
  }
  files<-list.files('E:/koveri/inputs')
  for(l in 1:length(files)){
    if(stress){
      if(substr(files[l],1,nchar(acCode)+7)==tolower(paste(acCode,'_stress',sep=''))){
        file1<-files[l]
      }
    }else{
      if(substr(files[l],1,nchar(acCode)+5)==tolower(paste(acCode,'_rest',sep=''))){
        file1<-files[l]
      }
    }
  }
  df<-read.table(paste('E:/koveri/inputs/',file1,sep=''))
  input<-c()
  for(l in 1:length(times)){
    input[l]<-linearInter(times[l],df$V1,df$V2)
  }
  return(input)
}

files<-list.files('C:/Users/Oona/Documents/Tpc/isp')
studyNumbers<-c()
for(i in 1:length(files)){
  if(substr(files[i],1,14)=='array1_koveri0'){
    studyNumbers<-c(studyNumbers,substr(files[i],8,17))
  }
}
t<-c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,80,90,100,120,140,
     160,190,220,250,280)
times<-c(0:280)

#1TCM with time delay and V_b
df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=7)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/Oona/Documents/Tpc/isp/array1_',studyNumbers[i],'.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,5)
  #input<-c()
  organCurve<-c()
  input<-findInput(studyNumbers[i],FALSE,times)
  for(l in 1:length(times)){
    #input[l]<-linearInter(times[l],t,df[2,])
    organCurve[l]<-linearInter(times[l],t,df[3,])
    if(is.na(input[l])){input[l]<-0}
  }
  #plot(times,input,type='l')
  plot(times,organCurve,type='l',ylab=i)
  initialValues<-c(1,1,-log(9))
  errorsForj<-c()
  for(j in 0:15){
    param<-optimizationAlg(initialValues,input,organCurve,j)
    err<-errorfunction(param,input,organCurve,j)
    errorsForj<-c(errorsForj,err)
  }
  j<-c(0:15)[which.min(errorsForj)]
  param<-optimizationAlg(initialValues,input,organCurve,j)
  points(times,fit_model(param,input,j),type='l')
  err<-errorfunction(param,input,organCurve,j)
  df1[i,1:length(param)]<-param
  df1[i,4]<-j
  df1[i,5]<-err/(length(times)-1)
  df1[i,6]<-meanRelError(fit_model(param,input,j),organCurve)
  df1[i,7]<-df[1,3]
}
write.csv(df1,file=paste('isp.csv',sep=''),row.names=F)
df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=6)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/Oona/Documents/Tpc/isp/array1_',studyNumbers[i],'.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,5)
  #input<-c()
  organCurve<-c()
  input<-findInput(studyNumbers[i],TRUE,times)
  for(l in 1:length(times)){
    #input[l]<-linearInter(times[l],t,df[4,])
    organCurve[l]<-linearInter(times[l],t,df[5,])
    if(is.na(input[l])){input[l]<-0}
  }
  #plot(times,input,type='l')
  plot(times,organCurve,type='l',ylab=i)
  initialValues<-c(1,1,-log(9))
  errorsForj<-c()
  for(j in 0:15){
    param<-optimizationAlg(initialValues,input,organCurve,j)
    err<-errorfunction(param,input,organCurve,j)
    errorsForj<-c(errorsForj,err)
  }
  j<-c(0:15)[which.min(errorsForj)]
  param<-optimizationAlg(initialValues,input,organCurve,j)
  points(times,fit_model(param,input,j),type='l')
  err<-errorfunction(param,input,organCurve,j)
  df1[i,1:length(param)]<-param
  df1[i,4]<-j
  df1[i,5]<-err/(length(times)-1)
  df1[i,6]<-meanRelError(fit_model(param,input,j),organCurve)
}
write.csv(df1,file=paste('isp_1.csv',sep=''),row.names=F)

df<-matrix(data=NA,nrow=length(studyNumbers),ncol=20)
df[,1]<-studyNumbers
filepath<-'C:/Users/Oona/Documents/Tpc/lok/Info for Oona.txt'
df1<-read.table(filepath)
for(i in 1:length(studyNumbers)){
  for(j in 1:length(df1$V1)){
    if(df[i,1]==df1$V1[j]){
      df[i,2]<-df1$V2[j]
      df[i,3]<-df1$V3[j]
    }
  }
}
indexes<-c()
for(i in 1:length(studyNumbers)){
  indexes<-c(indexes,as.numeric(substr(studyNumbers[i],7,10)))
}
library(readxl)
df1<-read_excel('E:/koveri/uusi_metadata_aivot.xlsx')
ids<-as.vector(df1$patient_id)[indexes]
for(i in 1:length(ids)){
  if(as.numeric(substr(ids[i],10,10))%%2==0){
    df[i,3]='Female'
  }else{df[i,3]='Male'}
}
for(i in 1:length(ids)){
  if(is.na(df[i,2])){
    df[i,2]=124-as.numeric(substr(ids[i],5,6))
  }
}
df[,4]<-as.numeric(as.vector(df1$weight)[indexes])
df[,5]<-as.numeric(as.vector(df1$height)[indexes])
df[,6]<-as.numeric(df[,4])/(as.numeric(df[,5])/100)^2
df[,7]<-as.numeric(as.vector(df1$'dose/lepo')[indexes])
df[,8]<-as.numeric(as.vector(df1$'dose/rasitus')[indexes])
df1<-read.csv('C:/Users/Oona/Documents/Tpc/isp/isp.csv')
df[,9]<-abs(df1[,1])*60
df[,10]<-abs(df1[,2])*60
df[,11]<-abs(df1[,1]/df1[,2])
df[,12]<-1/(1+exp(-df1[,3]))
df[,13]<-df1[,4]
df[,14]<-df1[,7]
df1<-read.csv('C:/Users/Oona/Documents/Tpc/isp/isp_1.csv')
df[,15]<-abs(df1[,1])*60
df[,16]<-abs(df1[,2])*60
df[,17]<-abs(df1[,1]/df1[,2])
df[,18]<-1/(1+exp(-df1[,3]))
df[,19]<-df1[,4]
df[,20]<-(as.numeric(df[,15])-as.numeric(df[,9]))/as.numeric(df[,9])

plot(1:72,as.numeric(df[,9]))
plot(1:72,as.numeric(df[,15]))
plot(1:72,as.numeric(df[,20]))

#df[,14]<-as.numeric(df[,14])*0.165*0.165*0.3
womendf<-df[df[,3]=="Female",]
mendf<-df[df[,3]=="Male",]

dim(df)
dim(mendf)
dim(womendf)

wilcox.test(as.numeric(df[,9]),as.numeric(df[,15]),paired=TRUE)
wilcox.test(as.numeric(df[,10]),as.numeric(df[,16]),paired=TRUE)
wilcox.test(as.numeric(df[,13]),as.numeric(df[,19]),paired=TRUE)
  
i=20
mean(as.numeric(df[,i]),na.rm=TRUE)
sd(as.numeric(df[,i]),na.rm=TRUE)
min(as.numeric(df[,i]),na.rm=TRUE)
max(as.numeric(df[,i]),na.rm=TRUE)

mean(as.numeric(womendf[,i]),na.rm=TRUE)
sd(as.numeric(womendf[,i]),na.rm=TRUE)
min(as.numeric(womendf[,i]),na.rm=TRUE)
max(as.numeric(womendf[,i]),na.rm=TRUE)

mean(as.numeric(mendf[,i]),na.rm=TRUE)
sd(as.numeric(mendf[,i]),na.rm=TRUE)
min(as.numeric(mendf[,i]),na.rm=TRUE)
max(as.numeric(mendf[,i]),na.rm=TRUE)

wilcox.test(as.numeric(womendf[,i]),as.numeric(mendf[,i]),paired=FALSE)

i=13
mean(as.numeric(womendf[,i]),na.rm=TRUE)
sd(as.numeric(womendf[,i]),na.rm=TRUE)
mean(as.numeric(mendf[,i]),na.rm=TRUE)
sd(as.numeric(mendf[,i]),na.rm=TRUE)
mean(as.numeric(df[,i]),na.rm=TRUE)
sd(as.numeric(df[,i]),na.rm=TRUE)
wilcox.test(as.numeric(womendf[,i]),as.numeric(mendf[,i]),paired=FALSE)

wilcox.test(as.numeric(womendf[,i]),as.numeric(womendf[,i+6]),paired=TRUE)
wilcox.test(as.numeric(mendf[,i]),as.numeric(mendf[,i+6]),paired=TRUE)
wilcox.test(as.numeric(df[,i]),as.numeric(df[,i+6]),paired=TRUE)

wilcox.test(as.numeric(womendf[,9]),as.numeric(mendf[,9]),paired=FALSE)
wilcox.test(as.numeric(womendf[,15]),as.numeric(mendf[,15]),paired=FALSE)
wilcox.test(as.numeric(womendf[,20]),as.numeric(mendf[,20]),paired=FALSE)

cor.test(as.numeric(womendf[,9]),as.numeric(womendf[,15]))
cor.test(as.numeric(mendf[,9]),as.numeric(mendf[,15]))
cor.test(as.numeric(df[,9]),as.numeric(df[,15]))

cor.test(as.numeric(womendf[,9]),as.numeric(womendf[,7]))
cor.test(as.numeric(mendf[,9]),as.numeric(mendf[,7]))
cor.test(as.numeric(df[,9]),as.numeric(df[,7]))

cor.test(as.numeric(womendf[,15]),as.numeric(womendf[,8]))
cor.test(as.numeric(mendf[,15]),as.numeric(mendf[,8]))
cor.test(as.numeric(df[,15]),as.numeric(df[,8]))

i=9
cor.test(as.numeric(df[,i]),as.numeric(df[,2]))
cor.test(as.numeric(df[,i]),as.numeric(df[,5]))
cor.test(as.numeric(df[,i]),as.numeric(df[,4]))
cor.test(as.numeric(df[,i]),as.numeric(df[,6]))
cor.test(as.numeric(df[,i]),as.numeric(df[,14]))

cor.test(as.numeric(womendf[,i]),as.numeric(womendf[,2]))
cor.test(as.numeric(womendf[,i]),as.numeric(womendf[,5]))
cor.test(as.numeric(womendf[,i]),as.numeric(womendf[,4]))
cor.test(as.numeric(womendf[,i]),as.numeric(womendf[,6]))
cor.test(as.numeric(womendf[,i]),as.numeric(womendf[,14]))

cor.test(as.numeric(mendf[,i]),as.numeric(mendf[,2]))
cor.test(as.numeric(mendf[,i]),as.numeric(mendf[,5]))
cor.test(as.numeric(mendf[,i]),as.numeric(mendf[,4]))
cor.test(as.numeric(mendf[,i]),as.numeric(mendf[,6]))
cor.test(as.numeric(mendf[,i]),as.numeric(mendf[,14]))

wilcox.test(as.numeric(df[,9]),as.numeric(df[,15]),paired=TRUE)

df0<-matrix(data=NA,nrow=length(studyNumbers),ncol=281)
df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=281)
df2<-matrix(data=NA,nrow=length(studyNumbers),ncol=281)
df3<-matrix(data=NA,nrow=length(studyNumbers),ncol=281)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/Oona/Documents/Tpc/isp/array1_',studyNumbers[i],'.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,5)
  for(l in 1:length(times)){
    #df0[i,l]<-linearInter(times[l],t,df[2,])/1000
    df1[i,l]<-linearInter(times[l],t,df[3,])/1000
    #df2[i,l]<-linearInter(times[l],t,df[4,])/1000
    df3[i,l]<-linearInter(times[l],t,df[5,])/1000
  }
  df0[i,]<-findInput(studyNumbers[i],FALSE,times)/1000
  df2[i,]<-findInput(studyNumbers[i],TRUE,times)/1000
}
l0<-c(colMeans(df0)-colSds(df0),rev(colMeans(df0)+colSds(df0)))
l1<-c(colMeans(df1)-colSds(df1),rev(colMeans(df1)+colSds(df1)))
l2<-c(colMeans(df2)-colSds(df2),rev(colMeans(df2)+colSds(df2)))
l3<-c(colMeans(df3)-colSds(df3),rev(colMeans(df3)+colSds(df3)))
plot(n,type='l',main='Mean +/- SD of mean TACs in the rest PET images',
     xlim=c(0,280),ylim=c(0,150),
     lwd=2,ylab='Activity concentration (kBq/mL)',xlab='Time (s)',
     cex.axis=1.3,cex.lab=1.3)
for(i in 1:length(l0)){
  if(!is.na(l0[i])){
    if(l0[i]<0){l0[i]<-0}
  }else{
    l0[i]<-0
  }
}
polygon(c(times,rev(times)),l0,col='white',border='black',lwd=2,lty=2)
points(times,colMeans(df0),type='l',lwd=2,col='black')

mycol0 <- rgb(0, 128, 255, max = 255, alpha = 100, names = "blue50")
polygon(c(times,rev(times)),l1,col=mycol0,border='blue',lwd=2,lty=2)
points(times,colMeans(df1),type='l',lwd=2,col='blue')

legend('topright',legend=c('Aorta',
                           'Pancreas'),
       lty=c(1,1),
       col=c('black','blue'),
       lwd=c(2,2),cex=1.4)

plot(n,type='l',main='Mean +/- SD of mean TACs in the stress PET images',
     xlim=c(0,280),ylim=c(0,150),
     lwd=2,ylab='Activity concentration (kBq/mL)',xlab='Time (s)',
     cex.axis=1.3,cex.lab=1.3)
for(i in 1:length(l2)){
  if(!is.na(l2[i])){
    if(l2[i]<0){l2[i]<-0}
  }else{
    l2[i]<-0
  }
}
polygon(c(times,rev(times)),l2,col='white',border='black',lwd=2,lty=2)
points(times,colMeans(df2),type='l',lwd=2,col='black')

mycol0 <- rgb(0, 128, 255, max = 255, alpha = 100, names = "blue50")
polygon(c(times,rev(times)),l3,col=mycol0,border='blue',lwd=2,lty=2)
points(times,colMeans(df3),type='l',lwd=2,col='blue')

legend('topright',legend=c('Aorta',
                           'Pancreas'),
       lty=c(1,1),
       col=c('black','blue'),
       lwd=c(2,2),cex=1.4)


boxplot(as.numeric(df[,9]),as.numeric(df[,15]),col=c('white','white'),
        ylab=c('Blood flow (mL/min/mL)'),
        names=c('Rest','Stress'),cex.axis=1.3,cex.lab=1.3,notch=TRUE)

plot(as.numeric(mendf[,9]),as.numeric(mendf[,15]),
     xlim=c(0,2),ylim=c(0,3),
     lwd=2,ylab='Stress blood flow (mL/min/mL)',
     xlab='Rest blood flow (mL/min/mL)',
     cex.axis=1.3,cex.lab=1.3,col='blue')
points(as.numeric(womendf[,9]),as.numeric(womendf[,15]),
       xlim=c(0,3),ylim=c(0,4.5),
       lwd=2,pch=3,axes=FALSE)
fit<-lm(as.numeric(mendf[,15])~as.numeric(mendf[,9]))
abline(fit,col='blue',lw=2)
fit<-lm(as.numeric(womendf[,15])~as.numeric(womendf[,9]))
abline(fit,lw=2,lty=2)
legend('topleft',legend=c('Female','Male'),lty=c(2,1),pch=c(3,1),
       col=c('black','blue'),lwd=c(2,2),cex=1.4)

plot(as.numeric(mendf[,5]),as.numeric(mendf[,9]),
     xlim=c(150,190),ylim=c(0,2),
     lwd=2,ylab='Rest blood flow (mL/min/mL)',xlab='Height (cm)',
     cex.axis=1.3,cex.lab=1.3,col='blue')
points(as.numeric(womendf[,5]),as.numeric(womendf[,9]),
       xlim=c(150,190),ylim=c(0,3),
       lwd=2,pch=3,axes=FALSE)
fit<-lm(as.numeric(mendf[,9])~as.numeric(mendf[,5]))
abline(fit,col='blue',lw=2)
fit<-lm(as.numeric(womendf[,9])~as.numeric(womendf[,5]))
abline(fit,lw=2,lty=2)
legend('bottomleft',legend=c('Female','Male'),lty=c(2,1),pch=c(3,1),
       col=c('black','blue'),lwd=c(2,2),cex=1.4)

plot(as.numeric(mendf[,14]),as.numeric(mendf[,9]),
     xlim=c(25,145),ylim=c(0,3),
     lwd=2,ylab='Rest blood flow (mL/min/mL)',
     xlab='Pancreatic volume (cm^3)',
     cex.axis=1.3,cex.lab=1.3,col='blue')
points(as.numeric(womendf[,14]),as.numeric(womendf[,9]),
       xlim=c(25,85),ylim=c(0,3),
       lwd=2,pch=3,axes=FALSE)
fit<-lm(as.numeric(mendf[,9])~as.numeric(mendf[,14]))
abline(fit,col='blue',lw=2)
fit<-lm(as.numeric(womendf[,9])~as.numeric(womendf[,14]))
abline(fit,lw=2,lty=2)
legend('topleft',legend=c('Female','Male'),lty=c(2,1),pch=c(3,1),
       col=c('black','blue'),lwd=c(2,2),cex=1.4)
