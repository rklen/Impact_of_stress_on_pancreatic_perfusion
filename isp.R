linearInter<-function(x1,x,y){
  if(x1==x[1]){return(y[1])}
  for(j in 1:(length(x)-1)){
    if(x1>x[j] & x1<=x[j+1]){
      return(y[j]+(x1-x[j])/(x[j+1]-x[j])*(y[j+1]-y[j]))
    }
  }
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
  v<-(1-Vb)*cT+Vb*input
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
v1<-c()
df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=7)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/Oona/Documents/Tpc/isp/array1_',studyNumbers[i],'.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,5)
  input<-c()
  organCurve<-c()
  for(l in 1:length(times)){
    input[l]<-linearInter(times[l],t,df[2,])
    organCurve[l]<-linearInter(times[l],t,df[3,])
  }
  v1<-c(v1,max(input))
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
  input<-c()
  organCurve<-c()
  for(l in 1:length(times)){
    input[l]<-linearInter(times[l],t,df[4,])
    organCurve[l]<-linearInter(times[l],t,df[5,])
  }
  v1<-c(v1,max(input))
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
df[,4]<-c(90,94,79,85,73,82,81,108,63,92,64,117,99,59,89,100,103,
         101,78,86,87,77,114,81,90,78,65,94,121,61,95,64,60,92,
         101,69,62,64,132,98,130,67,123,80,NA,102,93,76,91,72,87,
         88,78,113,110,81,95,69,107,NA,92,100,74,94,85,92,NA,100,
         102,78,83,69,65,86,70,78,98,74,105,90,92,80,91,73,111,120,
         73,81,75,63,86,98,74,98,NA,88,NA,88,53,84,61,97,75,88,103,
         67,69,80,120,78,81,76,84,88,88,72,89,130,98,97,100,63,71,
         127,76,98,98,77,85,68)[indexes]
df[,5]<-c(168,180,172,169,180,172,170,173,166,174,158,185,164,166,159,
         177,172,170,179,176,174,178,184,160,164,167,163,164,176,152,
         183,161,158,192,180,164,158,161,167,171,182,167,190,171,NA,
         180,182,162,167,165,172,163,181,170,192,165,182,168,170,NA,
         186,174,168,172,177,173,NA,180,177,162,183,165,166,170,175,
         178,186,165,176,168,169,162,169,161,178,173,162,164,173,165,
         172,190,172,183,NA,178,NA,184,158,174,162,160,176,189,168,
         157,166,160,185,161,164,160,177,163,175,165,168,185,170,160,
         174,164,158,172,168,188,174,171,165,
         176)[indexes]
df[,6]<-as.numeric(df[,4])/(as.numeric(df[,5])/100)^2
df[,7]<-c(355,356,351,308,316,336,364,354,350,406,357,319,334,350,356,
         341,339,332,324,322,408,402,364,385,363,334,380,364,336,351,
         362,366,396,390,359,364,334,345,379,390,348,366,321,332,295,
         370,356,371,323,307,374,383,333,381,349,338,353,360,362,366,
         369,342,302,322,307,297,353,406,333,362,327,348,373,339,360,
         328,321,329,347,327,330,318,389,340,336,344,330,342,328,347,
         339,346,369,NA,359,384,NA,363,355,360,366,337,325,328,353,
         356,331,350,364,399,361,376,402,357,378,312,322,323,305,339,
         340,315,323,315,306,304,402,386,419,
         337)[indexes]
df[,8]<-c(339,350,353,315,321,333,337,379,345,342,350,347,351,343,372,
        340,319,NA,332,323,407,409,352,362,382,328,369,359,341,367,
        338,343,394,397,358,374,360,370,420,359,351,360,316,336,386,
        355,341,380,332,310,372,388,347,370,319,353,348,350,395,349,
        382,354,315,336,326,306,367,348,341,357,350,351,358,331,374,
        324,350,343,349,337,334,307,374,337,340,339,321,323,328,367,
        360,361,379,NA,357,360,NA,373,360,349,381,345,337,327,356,
        351,327,362,354,382,375,383,374,356,385,320,334,336,319,412,
        328,306,320,325,320,298,392,409,418,
        352)[indexes]
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

df<-df[c(1:70,72:109),]
plot(1:108,as.numeric(df[,9]))
plot(1:108,as.numeric(df[,15]))
plot(1:108,as.numeric(df[,20]))

#df[,14]<-as.numeric(df[,14])*0.165*0.165*0.3
womendf<-df[df[,3]=="Female",]
mendf<-df[df[,3]=="Male",]

dim(df)
dim(mendf)
dim(womendf)

wilcox.test(as.numeric(df[,9]),as.numeric(df[,15]),paired=TRUE)
wilcox.test(as.numeric(df[,10]),as.numeric(df[,16]),paired=TRUE)
wilcox.test(as.numeric(df[,13]),as.numeric(df[,19]),paired=TRUE)
  
i=9
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

i=20
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
    df0[i,l]<-linearInter(times[l],t,df[2,])/1000
    df1[i,l]<-linearInter(times[l],t,df[3,])/1000
    df2[i,l]<-linearInter(times[l],t,df[4,])/1000
    df3[i,l]<-linearInter(times[l],t,df[5,])/1000
  }
}
df0<-df0[c(1:70,72:108),]
df1<-df1[c(1:70,72:108),]
df2<-df2[c(1:70,72:108),]
df3<-df3[c(1:70,72:108),]
l0<-c(colMeans(df0)-colSds(df0),rev(colMeans(df0)+colSds(df0)))
l1<-c(colMeans(df1)-colSds(df1),rev(colMeans(df1)+colSds(df1)))
l2<-c(colMeans(df2)-colSds(df2),rev(colMeans(df2)+colSds(df2)))
l3<-c(colMeans(df3)-colSds(df3),rev(colMeans(df3)+colSds(df3)))
plot(n,type='l',main='Mean +/- SD of mean TACs in the rest PET images',
     xlim=c(0,280),ylim=c(0,93),
     lwd=2,ylab='Activity concentration (kBq/mL)',xlab='Time (s)',
     cex.axis=1.3,cex.lab=1.3)

for(i in 1:length(l2)){
  if(l0[i]<0){l0[i]<-0}
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


boxplot(as.numeric(df[,9]),as.numeric(df[,15]),col=c('white','white'),
        ylab=c('Perfusion (mL/min/mL)'),
        names=c('Rest','Stress'),cex.axis=1.3,cex.lab=1.3,notch=TRUE)

plot(as.numeric(mendf[,9]),as.numeric(mendf[,15]),
     xlim=c(0,3),ylim=c(0,4.5),
     lwd=2,ylab='Stress perfusion (mL/min/mL)',
     xlab='Rest perfusion (mL/min/mL)',
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

plot(as.numeric(mendf[,2]),as.numeric(mendf[,9]),
     xlim=c(25,85),ylim=c(0,3),
     lwd=2,ylab='Rest perfusion (mL/min/mL)',xlab='Age',
     cex.axis=1.3,cex.lab=1.3,col='blue')
points(as.numeric(womendf[,2]),as.numeric(womendf[,9]),
     xlim=c(25,85),ylim=c(0,3),
     lwd=2,pch=3,axes=FALSE)
fit<-lm(as.numeric(mendf[,9])~as.numeric(mendf[,2]))
abline(fit,col='blue',lw=2)
fit<-lm(as.numeric(womendf[,9])~as.numeric(womendf[,2]))
abline(fit,lw=2,lty=2)
legend('topleft',legend=c('Female','Male'),lty=c(2,1),pch=c(3,1),
       col=c('black','blue'),lwd=c(2,2),cex=1.4)

plot(as.numeric(mendf[,2]),as.numeric(mendf[,15]),
     xlim=c(25,85),ylim=c(0,4.5),
     lwd=2,ylab='Stress perfusion (mL/min/mL)',xlab='Age',
     cex.axis=1.3,cex.lab=1.3,col='blue')
points(as.numeric(womendf[,2]),as.numeric(womendf[,15]),
       xlim=c(25,85),ylim=c(0,4.5),
       lwd=2,pch=3,axes=FALSE)
fit<-lm(as.numeric(mendf[,15])~as.numeric(mendf[,2]))
abline(fit,col='blue',lw=2)
fit<-lm(as.numeric(womendf[,15])~as.numeric(womendf[,2]))
abline(fit,lw=2,lty=2)
legend('topleft',legend=c('Female','Male'),lty=c(2,1),pch=c(3,1),
       col=c('black','blue'),lwd=c(2,2),cex=1.4)
