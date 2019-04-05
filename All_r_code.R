


require(R.utils)
#sourceDirectory("C:/Users/rfisher/OneDrive - Australian Institute of Marine Science/Documents/AIMS/Statistics and computing/Custom_functions")
devtools::install_github("beckyfisher/custom_functions")
library(custom.functions.pkg)
require(doBy)

#### go get the wq data --------------------------------------------------------
setwd("C:/Users/rfisher/Australian Institute of Marine Science/Ross Jones - Mortality pathways/Mortality_MS/R_clean")

wq.dat=read.csv(file="wq_data.csv")
dat=read.csv(file="coral_health_data.csv")
dat$color=dat$Dredge.Colour.Score-dat$start.color
                 # negative values = bleaching
                 # positive values = darkening


# scoring scale
score.list=list(c.1=0,
     c.2=c(0.01,0.05),
    c.3=c(0.05,0.33),
    c.4=c(0.33,0.66),
    c.5=c(0.66,0.95),
    c.6=c(0.95,0.99),
    c.7=1)

score.vals=unlist(lapply(score.list,FUN=mean))

# bleaching data
dat$prop.bleach=score.vals[dat$Thermal.Bleaching]
# sediment data
dat$prop.sed=score.vals[dat$Sediment.Cover]
# mucus data
dat$prop.mucus=score.vals[dat$Mucous]
# mortality data
dat$prop.mort=score.vals[dat$Mortality]

# accumulated sediment
dat$prop.smoth=dat$prop.sed-dat$prop.mucus
hist(dat$prop.smoth)

require(car)
plot(jitter(dat$color),jitter(logit(dat$Dredge.Bleach)))
plot(jitter(dat$color),jitter(logit(dat$prop.bleach)))
plot(jitter(logit(dat$Dredge.Bleach)),jitter(logit(dat$prop.bleach)))


### remove corals that were more than 33% dead at the start of dredging -------
use.colonies.index=which(is.na(match(dat$ColID,
                            as.character(unique(dat[which(dat$j.date<0 &
                                          dat$Mortality>3),"ColID"]))))!=F)
dat=dat[use.colonies.index,]
# summary of colonies by site
dat$index=1
x=summaryBy(index~Family+use.group+Site.Code+ColID,FUN=sum,data=dat)
x$index=1
xx=summaryBy(index~Family+use.group+Site.Code,FUN=sum,data=x)
xx[which(is.na(xx$use.group)==F),]


### --- create some lagged colony variables (ie state of last image) -----------
# for each colony, get the lagged values
colonies=unique(as.character(dat$ColID))
dat$lag1.jdate=NA
dat$lag1.mort=NA
dat$lag1.bleach=NA
dat$lag1.sed=NA
dat$lag1.mucus=NA
dat$lag1.Live=NA
dat$lag1.Bleach=NA
dat$lag1.color=NA
dat$lag1.smoth=NA

for(i in 1:length(colonies)){
 dat.c=dat[which(dat$ColID==colonies[i]),]
 index.c=which(dat$ColID==colonies[i])
 if(nrow(dat.c)>1){
  dat.c$lag1.jdate[2:nrow(dat.c)]=dat.c$j.date[1:(nrow(dat.c)-1)]
  dat.c$lag1.mort[2:nrow(dat.c)]=dat.c$prop.mort[1:(nrow(dat.c)-1)]
  dat.c$lag1.bleach[2:nrow(dat.c)]=dat.c$prop.bleach[1:(nrow(dat.c)-1)]
  dat.c$lag1.sed[2:nrow(dat.c)]=dat.c$prop.sed[1:(nrow(dat.c)-1)]
  dat.c$lag1.mucus[2:nrow(dat.c)]=dat.c$prop.mucus[1:(nrow(dat.c)-1)]
  dat.c$lag1.Live[2:nrow(dat.c)]=dat.c$Live[1:(nrow(dat.c)-1)]
  dat.c$lag1.Bleach[2:nrow(dat.c)]=dat.c$Bleach[1:(nrow(dat.c)-1)]
  dat.c$lag1.color[2:nrow(dat.c)]=dat.c$color[1:(nrow(dat.c)-1)]
  dat.c$lag1.smoth[2:nrow(dat.c)]=dat.c$prop.smoth[1:(nrow(dat.c)-1)]

  dat$lag1.jdate[index.c]=dat.c$lag1.jdate
  dat$lag1.mort[index.c]=dat.c$lag1.mort
  dat$lag1.bleach[index.c]=dat.c$lag1.bleach
  dat$lag1.sed[index.c]=dat.c$lag1.sed
  dat$lag1.mucus[index.c]=dat.c$lag1.mucus
  dat$lag1.Live[index.c]=dat.c$lag1.Live
  dat$lag1.Bleach[index.c]=dat.c$lag1.Bleach
  dat$lag1.color[index.c]=dat.c$lag1.color
  dat$lag1.smoth[index.c]=dat.c$lag1.smoth
}}

dat=dat[with(dat,order(ColID,j.date)),]

dat$dif.jdate=dat$j.date-dat$lag1.jdate
# calculate the amount of change between the time steps
dat$dif.mort=dat$prop.mort-dat$lag1.mort
dat$dif.bleach=dat$prop.bleach-dat$lag1.bleach
dat$dif.sed=dat$prop.sed-dat$lag1.sed
dat$dif.Live=dat$Live-dat$lag1.Live
dat$dif.Bleach=dat$Bleach-dat$lag1.Bleach
dat$dif.color=dat$color-dat$lag1.color
dat$dif.smoth=dat$prop.smoth-dat$lag1.smoth

hist(dat$dif.jdate)
hist(dat$dif.mort)
hist(dat$dif.bleach)
hist(dat$dif.Live)
hist(dat$dif.Bleach)
hist(dat$dif.color)
hist(dat$dif.smoth)

# calculate the loss as a proportion of the starting value
dat$prop.mort.ratio=round(dat$dif.mort*100)/round(100-dat$lag1.mort*100)
dat$prop.sed.ratio=round(dat$dif.sed*100)/round(100-dat$lag1.sed*100)
dat$prop.bleach.ratio=round(dat$dif.bleach*100)/round(100-dat$lag1.bleach*100)

# remove Inf as these are NA (they occur becuase the previous score at the highest
# mortality, bleaching or sediment score level (ie they can't die, bleach or get smoethered
# because they were already completely dead/bleached/smothered.
dat$prop.mort.ratio[which(dat$prop.mort.ratio=="Inf")]=NA
dat$prop.bleach.ratio[which(dat$prop.bleach.ratio=="Inf")]=NA
dat$prop.sed.ratio[which(dat$prop.bleach.ratio=="Inf")]=NA

dat$prop.mort.ratio[which(dat$prop.mort.ratio=="-Inf")]=NA
dat$prop.bleach.ratio[which(dat$prop.bleach.ratio=="-Inf")]=NA
dat$prop.sed.ratio[which(dat$prop.bleach.ratio=="-Inf")]=NA

# model only positive events(ie there is more bleaching/sediment than there was in the time before)
dat$prop.sed.ratio[which(dat$prop.sed.ratio<0)]=0
dat$prop.bleach.ratio[which(dat$prop.bleach.ratio<0)]=0

# event frequency metrics
# any event
dat$index=1
dat$mort.event=0
dat$mort.event[dat$dif.mort>0]=1
dat$mort.event[which(dat$lag1.mort==1)]=NA # you can't die if you are already dead

dat$bleach.event=0
dat$bleach.event[dat$dif.bleach>0]=1
dat$bleach.event[which(dat$lag1.bleach==1)]=NA # you can't bleach if you are already bleached

# L.e events
dat$Lmort.event=0
dat$Lmort.event[dat$dif.mort>0.4]=1
dat$Lmort.event[which(dat$lag1.mort>0.6)]=NA # you can't die if you are already dead

dat$Lbleach.event=0
dat$Lbleach.event[dat$dif.bleach>0.4]=1
dat$Lbleach.event[which(dat$lag1.bleach>0.6)]=NA # you can't bleach if you are already bleached

col.dat=dat

#### Extract data to use for analysis ------------------------------------------
# only use observations in the analysis were the gap between observations lies between 1 week  and 1 month
dat.use=dat[which(dat$dif.jdate<30 & dat$dif.jdate>7),]
table(dat.use$dif.jdate)
table(paste(dat.use$Site.Code[which(dat.use$use.group=="Acr.poci.br")],
            dat.use$FSN[which(dat.use$use.group=="Acr.poci.br")]))
table(paste(dat.use$Site.Code[which(dat.use$use.group=="Porit.ms")],
            dat.use$FSN[which(dat.use$use.group=="Porit.ms")]))


### create site level data set -------------------------------------------------
# need to model:
# 1. Morality events (dif.mort) ~ previous live cover
# 2. Bleaching events (dif.bleach)~previous unbleached cover
# 3. Sediment.events (dif.sed)~previous "clean" cover
require(doBy)
site.dat=summaryBy(j.date+Live+Bleach+..Points+color+ # chevron data
                   prop.bleach+prop.sed+prop.mucus+prop.mort+prop.smoth+
                   dif.bleach+dif.Live+dif.Bleach+dif.color+dif.smoth+dif.mort+dif.sed+
                   lag1.mort+lag1.bleach+lag1.sed+lag1.mucus+lag1.Live+lag1.Bleach+lag1.color+lag1.smoth+
                   prop.mort.ratio+prop.bleach.ratio+prop.sed.ratio+Size..cm2.~
                   use.group+Site.Code+FSN+dist.d+WaterHtmean,
                   data=dat.use,FUN=mean,na.rm=T,keep.names=T)
plot(site.dat$j.date,site.dat$prop.mort)

site.dat.index=summaryBy(index~
                   use.group+Site.Code+FSN+dist.d+WaterHtmean,
                   data=dat.use,FUN=sum,na.rm=T,keep.names=F)
summaryBy(index.sum~use.group,
                   data=site.dat.index,FUN=median,na.rm=T,keep.names=F)
# what is the maximum number of colonies per site? Ie - largest site level rep
site.acro.rep=summaryBy(index.sum~Site.Code,FUN=max,data=site.dat.index[
  which(site.dat.index$use.group=="Acr.poci.br"),])
range(site.acro.rep$index.sum.max)
median(site.acro.rep$index.sum.max)
site.porit.rep=summaryBy(index.sum~Site.Code,FUN=max,data=site.dat.index[
  which(site.dat.index$use.group=="Porit.ms"),])
range(site.porit.rep$index.sum.max)
median(site.porit.rep$index.sum.max)

length(unique(dat.use$ColID[which(dat.use$use.group=="Acr.poci.br")]))
length(unique(dat.use$ColID[which(dat.use$use.group=="Porit.ms")]))

dat=site.dat

#mortality                                 bleach.event+Lmort.event+Lbleach.event
mort.freq=summaryBy(index+mort.event ~
                           use.group+Site.Code+FSN,
                           data=dat.use[,c("index","mort.event","use.group","Site.Code","FSN")],
                           FUN=sum,keep.names=T,na.rm=T)
colnames(mort.freq)=c("use.group","Site.Code","FSN","mort.index","mort.event")
dat=merge(dat,mort.freq)
#Lmortality
Lmort.freq=summaryBy(index+Lmort.event ~
                           use.group+Site.Code+FSN,
                           data=dat.use[,c("index","Lmort.event","use.group","Site.Code","FSN")],
                           FUN=sum,keep.names=T,na.rm=T)
colnames(Lmort.freq)=c("use.group","Site.Code","FSN","Lmort.index","Lmort.event")
dat=merge(dat,Lmort.freq)
#bleach
bleach.freq=summaryBy(index+bleach.event ~
                           use.group+Site.Code+FSN,
                           data=dat.use[,c("index","bleach.event","use.group","Site.Code","FSN")],
                           FUN=sum,keep.names=T,na.rm=T)
colnames(bleach.freq)=c("use.group","Site.Code","FSN","bleach.index","bleach.event")
dat=merge(dat,bleach.freq)
#Lbleach
Lbleach.freq=summaryBy(index+Lbleach.event ~
                           use.group+Site.Code+FSN,
                           data=dat.use[,c("index","Lbleach.event","use.group","Site.Code","FSN")],
                           FUN=sum,keep.names=T,na.rm=T)
colnames(Lbleach.freq)=c("use.group","Site.Code","FSN","Lbleach.index","Lbleach.event")
dat=merge(dat,Lbleach.freq)


# find the worst case WQ values during the time period since the last observation
wq.pred.vars=c("log.ntu.wc.14d","light.stress.wc.14d","SAS.index.60d","temp")

# create longer term wq cumulative variables -----------------------------------
time.sequence=c(14,28,56,112,224)

all.prev.dat=list()
for(t in 1:length(time.sequence)){
  wq.worst.dat=list()
  wq.mean.dat=list()
  r=1
  for(r in 1:nrow(dat)){
      dat.r=dat[r,]
      tt=wq.dat[which(wq.dat$Site.Code==as.character(dat.r$Site.Code) &
                       wq.dat$j.date<dat.r$j.date &
                       wq.dat$j.date>=(dat.r$j.date-time.sequence[t])),]
      wq.worst.dat=c(wq.worst.dat,list(apply(tt[,wq.pred.vars],MARGIN=2,FUN=max, na.rm=T)))
      wq.mean.dat=c(wq.mean.dat,list(apply(tt[,wq.pred.vars],MARGIN=2,FUN=mean, na.rm=T)))
  }

  wq.worst.dat=lapply(wq.worst.dat,FUN=function(x){
      out=x
      out[which(out==-Inf)]=NA
      return(out)})
  wq.mean.dat=lapply(wq.mean.dat,FUN=function(x){
      out=x
      out[which(out==-Inf)]=NA
      return(out)})

  ww=do.call("rbind", wq.worst.dat)
  mm=do.call("rbind", wq.mean.dat)
  colnames(ww)=paste(colnames(ww),".max.",time.sequence[t],"days.prev",sep="")
  colnames(mm)=paste(colnames(mm),".mean.",time.sequence[t],"days.prev",sep="")

  all.prev.dat=c(all.prev.dat,list(cbind(mm,ww)))
}

aa=do.call("cbind",all.prev.dat)
colnames(aa)=gsub("log.ntu.wc.14d","ntu",colnames(aa))
colnames(aa)=gsub("light.stress.wc.14d","light",colnames(aa))
colnames(aa)=gsub("SAS.index.60d","sas",colnames(aa))

dat=cbind(dat,aa)
nrow(dat[which(dat$use.group=="Acr.poci.br"),])
nrow(dat[which(dat$use.group=="Porit.ms"),])

######### SEM ##################################################################
library(sem)
library(MASS)
#source("C:/Users/rfisher/Australian Institute of Marine Science/Ross Jones - Mortality pathways/figure_scripts/SEM/Path_functions.R")

dat.sem=dat[,c("j.date","Site.Code","FSN","use.group")]

dat.sem$SSC=dat$ntu.max.28days.prev
dat.sem$light=dat$light.max.28days.prev
dat.sem$sedimentation=dat$sas.max.28days.prev
dat.sem$sediment.cover=dat$lag1.sed
dat.sem$bleaching=dat$lag1.bleach
dat.sem$temperature=dat$temp.max.28days.prev

require(car)
dat.sem$response=dat$prop.mort.ratio

dat.sem=na.omit(dat.sem)#[which(dat.sem$response>0),])
dat.mort=dat.sem[which(dat.sem$response>0),]

mod.vars=c( "SSC","light","sedimentation","response",
            "sediment.cover","bleaching","temperature")

use.groups=na.omit(unique(as.character(dat.sem$use.group)))

tot.vals=list()
out.stats=list()
fit.out=list()
#g=1   # model magnitude of mortality events
for(g in 1:length(use.groups)){
  dat.g=dat.sem[which(dat.sem$use.group==use.groups[g]),]# &
                      #dat.sem$response>0),]
  require(car)
  dat.g$response=logit(dat.g$response*100)
  model.fit=specifyModel(paste("mFull_",use.groups[g],".txt",sep=""))
  options(fit.indices = c("GFI", "AGFI", "RMSEA", "NFI", "NNFI", "CFI", "RNI", "IFI", "SRMR", "AIC", "AICc", "BIC", "CAIC"))
  fit.g=sem(model.fit,data=dat.g[,mod.vars])

  fit.g.summary=summary(fit.g)
  out.stats.g=round(unlist(fit.g.summary[c("chisq","df","chisqNull","dfNull","GFI","AGFI",
   "RMSEA","NFI","NNFI","CFI","RNI","IFI","Rsq")]),3)

  coeff.g=stdCoef(fit.g)
  vals=coeff.g[,"Std. Estimate"]
  names(vals)=rownames(coeff.g)

  #bleaching - 1 path only
  tot.temp=vals["bleaching2response"]*vals["temperature2bleaching"]+
           vals["temperature2response"]
  #light
  l.direct=vals["light2response"]
  #l.bleaching=vals["light2bleaching"]*vals["bleaching2response"]
  tot.light=l.direct#+l.bleaching
  # sedimentation
  s.direct=vals["sedimentation2response"]
  s.sed.cover=vals["sedimentation2sediment.cover"]*vals["sediment.cover2response"]
  tot.sed=sum(s.direct,s.sed.cover,na.rm=T)
  # ssc
  tot.ssc=vals["SSC2response"]
  # residual total
  tot.res=vals["response2response"]

  all.totals=c(tot.temp,tot.light,tot.sed,tot.ssc,tot.res)

  names(all.totals)=c("Thermal stress","Light","Sedimentation","SSC","residual")
  tot.vals=c(tot.vals,list(all.totals))
  out.stats=c(out.stats,list(out.stats.g))
  fit.out=c(fit.out,list(fit.g))

}

names(tot.vals)=use.groups
names(out.stats)=use.groups
names(fit.out)=use.groups

plotvals=do.call("rbind",tot.vals)[,1:4]

pdf("SEM_barplot_site_lvl.pdf",width=6.3,height=3.5,pointsize=10)
par(mar=c(4,6,1,1))
barplot(plotvals,beside=T,horiz=T,
        xlab="Total standardized effect size",
        legend.text=use.groups,args.legend=list(x="topright"),las=1)
dev.off()
out.tab=do.call("cbind",lapply(out.stats,FUN=function(x){x[c("chisq","df","chisqNull","dfNull","GFI","AGFI",
                                      "NFI","NNFI","CFI","RNI","IFI","Rsq.response")]}))
write.csv(out.tab,file="sem_results_site_lvl.csv")

#Porites
pathDiagram(fit.out$Porit.ms,edge.labels ="values",file = paste("Porit.ms","pathDiagram",sep=""),standardize=T)
#Acr.poci.br
pathDiagram(fit.out$Acr.poci.br,edge.labels ="values",file = paste("Acr.poci.br","pathDiagram",sep=""),standardize=T)

######### Basic site level summary #############################################
use.groups=na.omit(unique(as.character(dat$use.group)))
#phases=c("p1","p2","p3")
#phase.dates.list=list(p1=c(1,189),p2=c(189,377),p3=c(377,565))
dat.use=dat.use[which(is.na(dat.use$use.group)!=TRUE),]
dd=summaryBy(Mortality+dif.mort+prop.mort~use.group+Site.Code+ColID,FUN=c(min,max),data=dat.use)
dd$dif.val=dd$prop.mort.max-dd$prop.mort.min
dd$check.sudden.death=dd$dif.val-dd$dif.mort.max

dd[which(dd$Mortality.max==7 & dd$check.sudden.death==0),]
write.csv(dd[which(dd$Mortality.max==7),],file="whole_colony_deaths.csv")

dd$index=1
temp.range=range(round(dat$temp.max.28days.prev), na.rm=T)
temp.vec=temp.range[1]:temp.range[2]
temp.cols=rev(colorRampPalette(c("red","yellow","blue"))(length(temp.vec)))#rainbow(length(dist.vec),start=0,end=.7)#rgb( ramp(seq(0, 1, length = length(dist.vec))), max = 255)
dat$temp.cols=temp.cols[match(round(dat$temp.max.28days.prev),temp.vec)]  # heat.colors(length(dist.vec))#rainbow(length(dist.vec))#
dat$Site.Code=reorder(dat$Site.Code, dat$dist.d, FUN = mean,
        order = is.ordered(dat$Site.Code))

dat.use$Site.Code=reorder(dat.use$Site.Code, dat.use$dist.d, FUN = mean,
        order = is.ordered(dat.use$Site.Code))


pretty.titles=c("Branching","Massive")
require(mgcv)

### by site. code
alpha.val=0.5

pdf(file="raw_mortality_plots.pdf",width=6.4, height=7, pointsize=12, onefile=T)
par(mfcol=c(3,2),oma=c(4,4,6,4),mar=c(0.5,0.5,0.5,0.5),bty="n",xpd=NA)
site.dat.list=list()
for(p in 1:length(use.groups)){
dd.p=dd[which(dd$use.group==use.groups[p]),]
f.dd.p=as.data.frame.matrix(xtabs(index~Site.Code+Mortality.max, data=dd.p))
colnames(f.dd.p)=paste("Mcat.",colnames(f.dd.p),sep="")
f.dd.p$Site.Code=rownames(f.dd.p)
f.dd.p$n.colonies=rowSums(f.dd.p[,paste("Mcat.",1:7,sep="")] )

dat.use.p=dat.use[which(dat.use$use.group==use.groups[p]),]
dat.use.p$index=1
dat.p=dat[which(dat$use.group==use.groups[p]),]
dat.p=dat.p[order(dat.p$temp.max.28days.prev),]#,decreasing=T),]
site.dat.p=summaryBy(dif.mort~Site.Code+dist.d+site.category,  data=dat.p, FUN=c(mean,sd),na.rm=T)
site.dat.p=merge(site.dat.p,f.dd.p)
site.dat.p$Site.Code=reorder(site.dat.p$Site.Code, site.dat.p$dist.d, FUN = mean,
        order = is.ordered(site.dat.p$Site.Code))
site.dat.p=site.dat.p[order(site.dat.p$dist.d),]
plot.mat=t(as.matrix(site.dat.p[,paste("Mcat.",1:7,sep="")]))

# raw scores at end of dredging
barplot.cols=colorRampPalette(c("white","yellow","orange","red","purple","black"))(nrow(plot.mat))
mids=barplot(plot.mat,col=barplot.cols,names.arg=rep(NA,ncol(plot.mat)),ylim=c(0,37),
             xaxt="n",yaxt="n")
legend("topleft", letters[p],bty="n", inset=c(-0.07,-0.03),cex=1.2)
axis(side=3,at=mids,labels=levels(dat.p$Site.Code),las=2)
 mtext(side=3,"Site",outer=F,line=3.5)
mtext(side=3,pretty.titles[p],outer=F,line=4.5)
if(p==1){
 legend("left",legend=
              c("0","1-5","6-33","33-66","66-95","95-99","100"),pch=15,col=barplot.cols,
                title="Mortality (%)",title.adj=0,inset=0.02,bty="n",pt.cex=1.2,ncol=1)
 axis(side=2,at=c(0,10,20,30),labels=c(0,10,20,30))
 mtext(side=2,"Number of colonies",outer=F,line=2)
}

# raw partial mortality between surveys
plot(jitter(as.numeric(dat.use.p$Site.Code),factor=0.5),#jitter(dat.use.p$dist.d,factor=500),
     jitter(dat.use.p$dif.mort,factor=20),
             ylim=c(0,1.1), xaxt="n",yaxt="n", cex=1.2,
             col=adjustcolor("black",alpha=0.4),pch=1,xlab="",ylab="")
if(p==1){
 axis(side=2)
 mtext(side=2,text="Partial mortality events",outer=T,line=2)}
legend("topleft", letters[3:4][p],bty="n", inset=c(-0.07,-0.03),cex=1.2)

# mean partial mortality between surveys
par(new=T)
plot(as.numeric(site.dat.p$Site.Code),
 site.dat.p$dif.mort.mean,
 xlab="",ylab="", xaxt="n",yaxt="n",ylim=c(0,0.04),col="red",pch=17,cex=1.2)
lines(as.numeric(site.dat.p$Site.Code),site.dat.p$dif.mort.mean,lty=3,col=2)
if(p==2){
 axis(side=4)
 mtext(side=4,text="Mean partial mortality rate",outer=T,line=2)}

# proportion of colonies showing some bleaching
plot(jitter(as.numeric(dat.p$Site.Code)),#jitter(dat.p$dist.d,factor=1000),
             xaxt="n",
     jitter(dat.p$bleach.event/dat.p$bleach.index,factor=50),
             ylim=c(0,1.1), cex=(dat.p$bleach.index/max(dat.p$bleach.index))+1,
             col=adjustcolor(dat.p$temp.cols,alpha=alpha.val),pch=16,xlab="",ylab="",yaxt="n")
legend("topleft", letters[5:6][p],bty="n", inset=c(-0.07,-0.03),cex=1.2)

if(p==1){
 axis(side=2)
 legend("topleft",legend=temp.vec[seq(1,12,by=2)],pch=16,
                title="Temperature 0C",bty="n",pt.cex=1.2,inset=c(0.04,0),
                col=adjustcolor(temp.cols[seq(1,12,by=2)],alpha=alpha.val))
 mtext(side=2,"Proportion bleached",outer=F,line=2) }

axis(side=1,at=1:length(levels(site.dat.p$Site.Code)),
              labels=signif(site.dat.p$dist.d,2),las=2)

if(p==2){
legend("topleft",
       legend=round(seq(min(dat.p$bleach.index),max(dat.p$bleach.index),length=5)),
       pt.cex=(round(seq(min(dat.p$bleach.index),max(dat.p$bleach.index),length=5))/
          max(dat.p$bleach.index))+1,pch=16,col="darkgrey",bty="n",inset=c(0.04,0),
          title="n. colonies",xpd=NA)
}

site.dat.list=c(site.dat.list,list(site.dat=site.dat.p,
                total=colSums(site.dat.p[,
                c("Mcat.1","Mcat.2","Mcat.3","Mcat.4","Mcat.5","Mcat.6","Mcat.7","n.colonies")]),
                means=colMeans(site.dat.p[,c(
                "dif.mort.mean","dif.mort.sd")])))
}

dev.off()

#site.dat.list

############ gam analysis ######################################################
dat.gam=dat[,c("j.date","Site.Code","FSN","use.group")]
dat.gam$SSC=dat$ntu.max.28days.prev
dat.gam$light=dat$light.max.28days.prev
dat.gam$sedimentation=dat$sas.max.28days.prev
dat.gam$temperature=dat$temp.max.28days.prev
dat.gam$sediment.cover=dat$lag1.sed

require(car)

dat.gam$bleaching=0 #dat$lag1.bleach
dat.gam$bleaching[round(dat$lag1.bleach,2)>=0.03]=1
dat.gam$Distance=dat$dist.d
dat.gam$Depth=dat$WaterHtmean
dat.gam$cumulative.stress.SSC=dat$ntu.mean.224days.prev
dat.gam$cumulative.stress.light=dat$light.mean.224days.prev
dat.gam$cumulative.stress.sediment=dat$sas.mean.224days.prev
dat.gam$m.trials=dat$mort.index #number of colony observations
dat.gam$mort.event=dat$mort.event
dat.gam$bleach.event=dat$bleach.event
dat.gam$b.trials=dat$bleach.index

dat.gam$mortality=dat$prop.mort.ratio
dat.gam$prop.bleach.ratio=dat$prop.bleach.ratio

# load the FSS function --------------------------------------------------------
library(RCurl)

devtools::install_github("beckyfisher/FSSgam_package")
library(FSSgam)

require(mgcv); require(MuMIn); require(doParallel); require(gamm4)

##### now model bleaching  #####################################################
#dat.c=na.omit(dat.gam[which(dat.gam$use.group==use.groups[1] |
#                              dat.gam$use.group==use.groups[2]),
#                        c("prop.bleach.ratio","bleach.event","b.trials","Site.Code","FSN",
#                        "light","temperature","use.group")])
#dat.c$ID=1:nrow(dat.c)
#dat.c$use.group=factor(dat.c$use.group)
#response.mat=as.matrix(data.frame(successes=dat.c$bleach.event,
#                          failures=dat.c$b.trials-dat.c$bleach.event))
#model.1=uGamm(response.mat~t2(light,temperature, by=use.group,k=4,bs="cr"),
#             random=~(1|Site.Code),
#             family="binomial",data=dat.c,lme4=T)
#summary(model.1$gam)
#
#out.list=full.subsets.gam(use.dat=dat.c,        #       check.model.set
#                             test.fit=model.1,
#                             pred.vars.cont=c("light","temperature"),
#                             pred.vars.fact="use.group",
#                             smooth.smooth.interactions=T,
#                             k=5)
##out.list$n.mods
#save(out.list,file="FSS_out_site_lvlBleach_v1.RData")

##### now model PA mortality  #####################################################
# find the worst case WQ values during the time period since the last observation
wq.pred.vars=c("log.ntu.wc.14d","light.stress.wc.14d","SAS.index.60d","SAS.index.14d","temp")

wq.worst.dat=list()
r=1
for(r in 1:nrow(dat.use)){
    dat.r=dat.use[r,]
    tt=wq.dat[which(wq.dat$Site.Code==as.character(dat.r$Site.Code) &
                     wq.dat$j.date<dat.r$j.date &
                     wq.dat$j.date>=dat.r$j.date-dat.r$dif.jdate),]
    wq.worst.dat=c(wq.worst.dat,list(apply(tt[,wq.pred.vars],MARGIN=2,FUN=max, na.rm=T)))
}

wq.worst.dat=lapply(wq.worst.dat,FUN=function(x){
    out=x
    out[which(out==-Inf)]=NA
    return(out)})
dd=do.call("rbind", wq.worst.dat)

dat.use=cbind(dat.use,dd)

# create lagged wq variables ---------------------------------------------------
dat.use$lag1.log.ntu=NA
dat.use$lag1.light.stress=NA
dat.use$lag1.SAS=NA
dat.use$lag1.temp=NA

for(i in 1:length(colonies)){
 dat.use.c=dat.use[which(dat.use$ColID==colonies[i]),]
 index.c=which(dat.use$ColID==colonies[i])
 if(nrow(dat.use.c)>1){
  dat.use.c$lag1.log.ntu[2:nrow(dat.use.c)]=dat.use.c$log.ntu.wc.14d[1:(nrow(dat.use.c)-1)]
  dat.use.c$lag1.light.stress[2:nrow(dat.use.c)]=dat.use.c$light.stress.wc.14d[1:(nrow(dat.use.c)-1)]
  dat.use.c$lag1.SAS[2:nrow(dat.use.c)]=dat.use.c$SAS.index.60d[1:(nrow(dat.use.c)-1)]
  dat.use.c$lag1.temp[2:nrow(dat.use.c)]=dat.use.c$temp[1:(nrow(dat.use.c)-1)]

  dat.use$lag1.log.ntu[index.c]=dat.use.c$lag1.log.ntu
  dat.use$lag1.light.stress[index.c]=dat.use.c$lag1.light.stress
  dat.use$lag1.SAS[index.c]=dat.use.c$lag1.SAS
  dat.use$lag1.sed[index.c]=dat.use.c$lag1.sed
  dat.use$lag1.temp[index.c]=dat.use.c$lag1.temp
}}

dat.use=dat.use[with(dat.use,order(ColID,j.date)),]


head(unique(dat.use$ColID[which(dat.use$Mortality==7)]))
length(unique(dat.use$ColID[which(dat.use$Mortality==7)]))
colonies.died=as.character(unique(dat.use$ColID[which(dat.use$Mortality==7)]))
dat.use$dead=0
if(length(colonies.died)>0){
for(i in 1:length(colonies.died)){
  colony.i=colonies.died[i]
  dat.use.i=dat.use[which(dat.use$ColID==colony.i),c("j.date","Mortality")]
  dat.use.i=dat.use.i[order(dat.use.i$j.date),]
  date.died=min(dat.use.i$j.date[which(dat.use.i$Mortality==7)])
  index.dead=which(dat.use$ColID==colony.i & dat.use$j.date>date.died)
  if(length(index.dead)>0){
  dat.use$dead[index.dead]=1
}}}

dat.use=dat.use[which(dat.use$dead==0),]

dat.use$SSC=dat.use$log.ntu.wc.14d
dat.use$light=dat.use$light.stress.wc.14d
dat.use$sedimentation=dat.use$SAS.index.60d
dat.use$sediment.cover=0
dat.use$sediment.cover[dat.use$lag1.sed>=0.03]=1
dat.use$bleaching=0
dat.use$bleaching[dat.use$lag1.bleach>=0.03]=1
dat.use$temperature=dat.use$temp
dat.use$colony.size=dat.use$Size..cm2.
dat.use$colony.surroundings=jitter(dat.use$Background.substrate)
dat.use$Distance=dat.use$dist.d
dat.use$site.sediment.cover=dat.use$prop.sed.mean
dat.use$Depth=dat.use$WaterHtmean

dat.use$response=dat.use$dif.mort#/(1-dat.use$lag1.mort)#dat.all$dif.Live
dat.use$trials=round(100-dat.use$lag1.mort*100)  # previous estimated live cover
dat.use$mortality=round(dat.use$response*100)
dat.use$mort.event=0
dat.use$mort.event[dat.use$response>0]=1
cont.preds=c("SSC","light","sedimentation","Depth","Distance","temperature")

n.site.fsn=summaryBy(index~Site.Code+use.group+FSN,FUN=sum,data=dat.use)
median(n.site.fsn$index.sum[which(n.site.fsn$use.group=="Acr.poci.br")])
median(n.site.fsn$index.sum[which(n.site.fsn$use.group=="Porit.ms")])

#g=1
#for(g in 1:length(use.groups)){
#  if(use.groups[g]=="Acr.poci.br"){
#         cat.preds="bleaching"}else{cat.preds=c("bleaching","sediment.cover")}
#  dat.c=na.omit(dat.use[which(dat.use$use.group==use.groups[g]),
#                        c("trials","mortality","mort.event","Site.Code","ColID","FSN",cont.preds,cat.preds)])
#  dat.c$ID=1:nrow(dat.c)
#  model.1=uGamm(mort.event~s(SSC,k=5,bs="cr"),
#             random=~(1|ColID),
#             family="binomial",data=dat.c,lme4=T)
#  summary(model.1$gam)
#  out.list=full.subsets.gam(use.dat=dat.c,
#                             test.fit=model.1,
#                             pred.vars.cont=cont.preds,
#                             pred.vars.fact=cat.preds,
#                             smooth.smooth.interactions=T)
#  #out.list$n.mods
#  save(out.list,file=paste("mortality_predictors_fss_gam_PA_v2",use.groups[g],".RData"))
#  print(g)
#}

##### now model loss given mortality  #####################################################
for(g in 1:length(use.groups)){
  if(use.groups[g]=="Acr.poci.br"){
         cat.preds="bleaching"}else{cat.preds=c("bleaching","sediment.cover")}
  dat.c=na.omit(dat.use[which(dat.use$use.group==use.groups[g]),
                        c("trials","mortality","Site.Code","ColID","FSN",cont.preds,cat.preds)])
  dat.c$ID=1:nrow(dat.c)
  dat.c=dat.c[which(dat.c$mortality>0),]
  dat.c$index=1
  n.site.fsn.c=summaryBy(index~Site.Code+FSN,FUN=sum,data=dat.c)
  print(median(n.site.fsn.c$index.sum))
  print(range(n.site.fsn.c$index.sum))
  print(nrow(dat.c))
}

#g=1
#for(g in 1:length(use.groups)){
#  if(use.groups[g]=="Acr.poci.br"){
#         cat.preds="bleaching"}else{cat.preds=c("bleaching","sediment.cover")}
#  dat.c=na.omit(dat.use[which(dat.use$use.group==use.groups[g]),
#                        c("trials","mortality","Site.Code","ColID","FSN",cont.preds,cat.preds)])
#  dat.c$ID=1:nrow(dat.c)
#  dat.c=dat.c[which(dat.c$mortality>0),]
#  response.mat=as.matrix(data.frame(successes=dat.c$mortality,
#                          failures=dat.c$trials-dat.c$mortality))
#  model.1=uGamm(response.mat~s(SSC,k=5,bs="cr"),
#             random=~(1|ColID)+(1|ID),
#             family="binomial",data=dat.c,lme4=T)
#  summary(model.1$gam)
#
#  out.list=full.subsets.gam(use.dat=dat.c,        #       check.model.set
#                             test.fit=model.1,
#                             pred.vars.cont=cont.preds,
#                             pred.vars.fact=cat.preds,
#                             smooth.smooth.interactions=T,
#                             cov.cutoff=0.4)
#  #out.list$n.mods
#  save(out.list,file=paste("FSS_out_no_zero_v1_wT",use.groups[g],".RData"))
#  print(g)
#}

#### now collate all the output ------------------------------------------------
var.imp=list()
top.all=list()
top.fits.list=list()
all.stats.out=list()
cor.matrix.out=list()

### probablity of bleaching
load("FSS_out_site_lvlBleach_v1.RData")
mod.table=out.list$mod.data.out
mod.table$R.sq=unlist(lapply(out.list$success.models,FUN=function(x){summary(x$gam)$r.sq}))
mod.table=mod.table[order(mod.table$AICc),]
out.i=mod.table
all.stats.out=c(all.stats.out,list(out.i))
var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))
all.less.wAICc=mod.table[which(mod.table$delta.AICc<2),]
top.all=c(top.all,list(all.less.wAICc))
top.fits.g=list()
# plot the all best models
pdf(file=paste("mod_fits_",use.groups[g],".pdf",sep=""),onefile=T)
par(oma=c(1,1,4,1))
for(r in 1:nrow(all.less.wAICc)){
best.model.name=as.character(all.less.wAICc$modname[r])
best.model=out.list$success.models[[best.model.name]]
if(best.model.name!="null"){
 plot(best.model$gam,all.terms=T,pages=1,residuals=T,pch=16)
 mtext(side=3,text=use.groups[g],outer=T)}
top.fits.g=c(top.fits.g,list(out.list$success.models[as.character(all.less.wAICc$modname[r])]))
}
names(top.fits.g)=as.character(all.less.wAICc$modname)
dev.off()
top.fits.list=c(top.fits.list,list(top.fits.g))
cor.matrix.out=c(cor.matrix.out,list(out.list$predictor.correlations))
write.csv(out.i[,c("modname","AICc","BIC","delta.AICc","delta.BIC","wi.AICc","wi.BIC","edf","R.sq")],
           file="all_model_stats_bleaching.csv")
write.csv(out.list$predictor.correlations,
           file="Predictor_correlations_bleachin.csv")

### make one uber plot for sup
pdf(file="mod_fits_all.pdf",
   width=10, height=18, pointsize=8, onefile=T)
par(mfrow=c(9,5),oma=c(1,1,4,1),mar=c(4,2,0,0),oma=c(0,2,1,1))

### probability of mortality
for(g in 1:length(use.groups)){
 load(file=paste("mortality_predictors_fss_gam_PA_v2",use.groups[g],".RData"))

 mod.table=out.list$mod.data.out
 mod.table$R.sq=unlist(lapply(out.list$success.models,FUN=function(x){summary(x$gam)$r.sq}))
 mod.table=mod.table[order(mod.table$AICc),]
 out.i=mod.table
 all.stats.out=c(all.stats.out,list(out.i))
 var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))
 all.less.wAICc=mod.table[which(mod.table$delta.AICc<2),]
 top.all=c(top.all,list(all.less.wAICc))
 top.fits.g=list()
 # plot the all best models
 #pdf(file=paste("mod_fits_",use.groups[g],".pdf",sep=""),onefile=T)
 #par(oma=c(1,1,4,1))
 for(r in 1:nrow(all.less.wAICc)){
 best.model.name=as.character(all.less.wAICc$modname[r])
 best.model=out.list$success.models[[best.model.name]]
 if(best.model.name!="null"){
   plot(best.model$gam,all.terms=T,rug=F)#,residuals=T,pch=16)
   mtext(side=3,text=use.groups[g],outer=T)}
 top.fits.g=c(top.fits.g,list(out.list$success.models[as.character(all.less.wAICc$modname[r])]))
 }
 names(top.fits.g)=as.character(all.less.wAICc$modname)
 #dev.off()
 top.fits.list=c(top.fits.list,list(top.fits.g))
 cor.matrix.out=c(cor.matrix.out,list(out.list$predictor.correlations))
 write.csv(out.i[,c("modname","AICc","BIC","delta.AICc","delta.BIC","wi.AICc","wi.BIC","edf","R.sq")],
           file=paste("all_model_stats",use.groups[g],".csv",sep=""))
 write.csv(out.list$predictor.correlations,
           file=paste("Predictor_correlations",use.groups[g],".csv",sep=""))
}

### proportional loss (given mortality)
for(g in 1:length(use.groups)){
 load(file=paste("FSS_out_no_zero_v1_wT",use.groups[g],".RData"))

  #print(head(out.list$mod.data.out))
 # examine the list of failed models
 mod.table=out.list$mod.data.out
 mod.table$R.sq=unlist(lapply(out.list$success.models,FUN=function(x){summary(x$gam)$r.sq}))
 mod.table=mod.table[order(mod.table$AICc),]
 out.i=mod.table
 all.stats.out=c(all.stats.out,list(out.i))
 var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))
 all.less.wAICc=mod.table[which(mod.table$delta.AICc<2),]
 top.all=c(top.all,list(all.less.wAICc))
 top.fits.g=list()
 # plot the all best models
 #pdf(file=paste("mod_fits_no_zero",use.groups[g],".pdf",sep=""),onefile=T)
 #par(oma=c(1,1,4,1))
 for(r in 1:nrow(all.less.wAICc)){
 best.model.name=as.character(all.less.wAICc$modname[r])
 best.model=out.list$success.models[[best.model.name]]
 if(best.model.name!="null"){
   plot(best.model$gam,all.terms=T,rug=F)#,residuals=T,pch=16)
   mtext(side=3,text=use.groups[g],outer=T)}
 top.fits.g=c(top.fits.g,list(out.list$success.models[as.character(all.less.wAICc$modname[r])]))
 }
 names(top.fits.g)=as.character(all.less.wAICc$modname)
 #dev.off()
 top.fits.list=c(top.fits.list,list(top.fits.g))
 cor.matrix.out=c(cor.matrix.out,list(out.list$predictor.correlations))
 write.csv(out.i[,c("modname","AICc","BIC","delta.AICc","delta.BIC","wi.AICc","wi.BIC","edf","R.sq")],
           file=paste("all_model_stats",use.groups[g],"_PresenceOnly.csv",sep=""))
 write.csv(out.list$predictor.correlations,
           file=paste("Predictor_correlations",use.groups[g],"_PresenceOnly.csv",sep=""))
}


dev.off()


name.vars=c("Probability of bleaching",
            paste(use.groups,"Probability of any mortality"),
            paste(use.groups,"Proportional loss"))
names(var.imp)=name.vars
names(top.all)=name.vars
names(top.fits.list)=name.vars
names(all.stats.out)=name.vars
names(cor.matrix.out)=name.vars

### plot the heatmap -----------------------------------------------------------
use.preds=c("SSC",
            "light",
            "sedimentation",
            "temperature",
            "Depth",
            "Distance",
            "sediment.cover",
            "bleaching")
var.imp.plot=do.call("rbind",lapply(var.imp[2:5],FUN=function(x){x[use.preds]}))
colnames(var.imp.plot)=c("SSC",
                         "Light",
                         "Sedimentation index",
                         "Temperature",
                         "Depth",
                         "Distance",
                         "Sediment cover",
                         "Bleached")

var.imp.plot=var.imp.plot[4:1,]
pdf("Mortality_predictors_heatmap.pdf",height=5,width=6.3,pointsize=10)
par(mfrow=c(1,1),oma=c(0,12,5,3))
importance.heatmap(var.imp.plot,inset=-0.2)
dev.off()

require(boot)
load("FSS_out_site_lvlBleach_v1.RData")
dat.c.bleach=out.list$used.data
out.list.bleach=out.list


load(file="mortality_predictors_fss_gam_PA_v2 Acr.poci.br .RData")
out.list.prob.any.acro=out.list
load(file="mortality_predictors_fss_gam_PA_v2 Porit.ms .RData")
out.list.prob.any.porit=out.list
load(file="FSS_out_no_zero_v1_wT Acr.poci.br .RData")
out.list.prop.loss.acro=out.list
dat.c.prop.loss.acro=out.list$used.data
load(file="FSS_out_no_zero_v1_wT Porit.ms .RData")
dat.c.prop.loss.porit=na.omit(out.list$used.data)
out.list.prop.loss.porit=out.list


### make plots of best models --------------------------------------------------
pdf("Bleaching_and_mortality_main_effects.pdf",height=7.5,width=6.3,pointsize=10,onefile=T)
par(mfrow=c(4,2),mar=c(3,1,2,0),oma=c(1,3,1,7),bty="l")
y.lim=1
axis.seq=seq(0.3,0.9,by=0.1)

### plot bleaching
contour.vec=c(0.01,0.03,0.1,0.25,0.5,0.75,0.9)
col.break.vec=c(10^-7,contour.vec,0.999999999)
col.pal=colorRampPalette(colours()[c(491,554,503,499)])
plot.model=out.list.bleach$success.models$'use.group+light.te.temperature'

new.dat=expand.grid(list(use.group=c("Acr.poci.br","Porit.ms"),
                         light=quantile(dat.c.bleach$light,probs=0.05),
                         temperature=quantile(dat.c.bleach$temperature,probs=0.95)))
eta=pred.back.logit(plot.model$gam,new.dat)
barx <- barplot(eta$P,ylab="",
                names.arg=c("Acrop./Pocillop. Br","Porites Ms"),#NA,
                xpd=NA,
                ylim=c(0,1),col=adjustcolor("black",alpha=0.5))
error.bar(x=barx,y=eta$P,upper=eta$P.up,lower=eta$P.lw)
mtext(side=2,"Probability of bleaching",line=2.5,cex=0.75,xpd=NA)
legend("topleft",legend="a",bty="n",xpd=NA,inset=c(-0.1,-0.15),text.font=2)

myvis.gam(plot.model$gam,
 plot.type="contour",cond=list(use.group="Porit.ms"),
 color="custom",custom.colors=col.pal,alpha=0.5,
 #type="response",
 contour.levels=logit(contour.vec),
 contour.labels=signif(contour.vec,2),
 n.grid=200,
 nCol=length(contour.vec)+1,
 col.breaks=logit(col.break.vec),
 #ylim=c(28.5,31.5),
 view=c("light","temperature"),
 xaxt="n",xlab="",ylab="",yaxt="n",)
axis(side=1,at=axis.seq,labels=signif(dli.xform(axis.seq),2))
legend("right",legend=as.character(round(col.break.vec,2)),
         fill=adjustcolor(col.pal(length(col.break.vec)),alpha=0.5),bg="white",
         border=NA,xpd=NA,pt.cex=3,y.intersp=0.8,bty="n",
         x.intersp=0.5,inset=-0.35,
         title="")
mtext(side=3,"",line=0.5,cex=0.75)
axis(side=4)
mtext(side=1,"Light (DLI)",line=2.5,cex=0.75)
mtext(side=4,"Temperature (oC)",line=2.5,cex=0.75)
legend("topleft",legend="b",bty="n",xpd=NA,inset=c(-0.1,-0.15),text.font=2)

## probability of any mortality
col.pal=colorRampPalette(c(1,"floralwhite",2))   # "magenta"  ,"orange
z.lim=c(-6,0.8)
legend.seq=logit(seq(inv.logit(z.lim[1]),inv.logit(z.lim[2]),by=0.1))
#plot.model=out.list.prob.any.acro$success.models$'bleaching+light.by.bleaching+temperature.by.bleaching'
plot.model=out.list.prob.any.acro$success.models$'light+bleaching+temperature.by.bleaching'
# - Acr.poci.br
myvis.gam(plot.model$gam,
 plot.type="contour",
 #type="response",
 zlim=z.lim,
 nCol=50,
 view=c("light","temperature"),
 cond=list(bleaching=0),
 color="custom",custom.colors=col.pal,alpha=1,
 xaxt="n",n.grid=100,plot.contours=F,
 xlab="",ylab="")
axis(side=1,at=axis.seq,labels=signif(dli.xform(axis.seq),2))
legend("topleft",legend="c",bty="n",xpd=NA,inset=c(-0.1,-0.15),text.font=2)
mtext(side=3,"not bleached",line=0.25,cex=0.75)
mtext(side=1,"Light (DLI)",line=2.5,cex=0.75)
mtext(side=2,"Temperature (oC)",line=2.5,cex=0.75)

myvis.gam(plot.model$gam,
 plot.type="contour",
 #type="response",
 zlim=z.lim,
 nCol=50,
 view=c("light","temperature"),
 cond=list(bleaching=1),
 color="custom",custom.colors=col.pal,alpha=1,
 xaxt="n",yaxt="n", n.grid=100,plot.contours=F,
 xlab="",ylab="",xaxt="n")
axis(side=1,at=axis.seq,labels=signif(dli.xform(axis.seq),2))
legend("right",legend=as.character(round(inv.logit(legend.seq),2)),
         fill=adjustcolor(col.pal(length(legend.seq)),alpha=1),
         title="",border=NA,xpd=NA,bty="n",pt.cex=3,y.intersp=0.8,
         x.intersp=0.5,inset=-0.125)
mtext(side=3,"bleached",line=0.25,cex=0.75)
mtext(side=1,"Light (DLI)",line=2.5,cex=0.75)
mtext(side=2,"",line=2.5,cex=0.75)
legend("topleft",legend="d",bty="n",xpd=NA,inset=c(-0.1,-0.15),text.font=2)

# - Porit.ms
range.sedimentation=range(dat.c.prop.loss.porit$sedimentation)
axis.seq=seq(0,1,by=0.1)
#plot.model=out.list.prob.any.porit$success.models$'sedimentation+bleaching+temperature.by.bleaching'
plot.model=out.list.prob.any.porit$success.models$'sedimentation+temperature+bleaching'

# - Acr.poci.br
myvis.gam(plot.model$gam,
 plot.type="contour",
 #type="response",
 zlim=z.lim,
 nCol=50,
 view=c("sedimentation","temperature"),
 cond=list(bleaching=0),
 color="custom",custom.colors=col.pal,alpha=1,
 xaxt="n",n.grid=100,plot.contours=F,
 xlab="",ylab="")
axis(side=1,at=axis.seq,labels=axis.seq)
legend("topleft",legend="c",bty="n",xpd=NA,inset=c(-0.1,-0.15),text.font=2)
mtext(side=3,"not bleached",line=0.25,cex=0.75)
mtext(side=1,"Sedimentation index",line=2.5,cex=0.75)
mtext(side=2,"Temperature (oC)",line=2.5,cex=0.75)

myvis.gam(plot.model$gam,
 plot.type="contour",
 #type="response",
 zlim=z.lim,
 nCol=50,
 view=c("sedimentation","temperature"),
 cond=list(bleaching=1),
 color="custom",custom.colors=col.pal,alpha=1,
 xaxt="n",yaxt="n", n.grid=100,plot.contours=F,
 xlab="",ylab="",xaxt="n")
axis(side=1,at=axis.seq,labels=axis.seq)
legend("right",legend=as.character(round(inv.logit(legend.seq),2)),
         fill=adjustcolor(col.pal(length(legend.seq)),alpha=1),
         title="",border=NA,xpd=NA,bty="n",pt.cex=3,y.intersp=0.8,
         x.intersp=0.5,inset=-0.125)
mtext(side=3,"bleached",line=0.25,cex=0.75)
mtext(side=1,"Sedimentation index",line=2.5,cex=0.75)
mtext(side=2,"",line=2.5,cex=0.75)
legend("topleft",legend="d",bty="n",xpd=NA,inset=c(-0.1,-0.15),text.font=2)

### proportional loss - Acr.poci.br
plot.model=out.list.prop.loss.acro$success.models$'bleaching'
dat.c.prop.loss.acro$prop.loss=dat.c.prop.loss.acro$mortality/dat.c.prop.loss.acro$trials
dat.c.prop.loss.acro$plot.x=NA
dat.c.prop.loss.acro$plot.x[which(dat.c.prop.loss.acro$bleaching==0)]=0.7
dat.c.prop.loss.acro$plot.x[which(dat.c.prop.loss.acro$bleaching==1)]=1.9
dat.c.prop.loss.acro$plot.col=NA
dat.c.prop.loss.acro$plot.col[which(dat.c.prop.loss.acro$bleaching==0)]="navy"
dat.c.prop.loss.acro$plot.col[which(dat.c.prop.loss.acro$bleaching==1)]="orange"

new.dat=expand.grid(list(bleaching=c(0,1)))
eta=pred.back.logit(plot.model$gam,new.dat)
barx <- barplot(eta$P,ylab="Proportional colony loss",border=c("navy","orange"),
                names.arg=c("not bleached","bleached"),xpd=NA,
                ylim=c(0,1.1),col=adjustcolor(c("navy","orange"),alpha=0.5))
#points(jitter(dat.c.prop.loss.acro$plot.x,factor=0.2),jitter(dat.c.prop.loss.acro$prop.loss,factor=20),
#       col=adjustcolor(dat.c.prop.loss.acro$plot.col,alpha=0.2),pch=16)

error.bar(x=barx,y=eta$P,upper=eta$P.up,lower=eta$P.lw,col=c("navy","orange"))

mtext(side=1,"Bleaching status",line=2.5,cex=0.75)
legend("topleft",legend="f",bty="n",xpd=NA,inset=c(-0.1,-0.15),text.font=2)

### proportional loss - Porit.ms
range.light=range(dat.c.prop.loss.porit$light)
axis.seq=seq(0.3,0.9,by=0.1)
plot.model=out.list.prop.loss.porit$success.models$'light+bleaching'
dat.c.prop.loss.porit$plot.col=NA
dat.c.prop.loss.porit$plot.col[which(dat.c.prop.loss.porit$bleaching==0)]="navy"
dat.c.prop.loss.porit$plot.col[which(dat.c.prop.loss.porit$bleaching==1)]="orange"

plot(dat.c.prop.loss.porit$light,jitter(dat.c.prop.loss.porit$mortality/
                                         dat.c.prop.loss.porit$trials,factor=20),pch=16,
     col=adjustcolor(dat.c.prop.loss.porit$plot.col,alpha=0.3),
     xaxt="n",yaxt="n",
     ylim=c(0,y.lim),xlab="",ylab="")
axis(side=1,at=axis.seq,labels=signif(dli.xform(axis.seq),2))
legend("topleft",legend=c("bleached","not bleached"),cex=1,pt.cex=2,
       pch=15,#fill=adjustcolor(c("red","navy"),alpha=0.5),
       lty=1,col=adjustcolor(c("orange","navy"),alpha=0.5),bty="n")
axis(side=4,xpd=NA)
mtext(side=4,"Proportional colony loss",line=2.5,cex=0.75,xpd=NA)
mtext(side=1,"Light (DLI)",line=2.5,cex=0.75)
# bleached corals
new.dat=expand.grid(list(bleaching=1,
                         light=seq(range.light[1],range.light[2],length=100)))
eta=pred.back.logit(plot.model$gam,new.dat)
lines(eta$light,eta$P,col="orange")
polygon(c(eta$light,rev(eta$light)),
          c(eta$P.lw,rev(eta$P.up)),
          col=adjustcolor("orange",alpha=0.5),border=NA)
# unbleached corals
new.dat=expand.grid(list(bleaching=0,
                         light=seq(range.light[1],range.light[2],length=100)))
eta=pred.back.logit(plot.model$gam,new.dat)
lines(eta$light,eta$P,col="navy")
polygon(c(eta$light,rev(eta$light)),
          c(eta$P.lw,rev(eta$P.up)),
          col=adjustcolor("navy",alpha=0.5),border=NA)
legend("topleft",legend="g",bty="n",xpd=NA,inset=c(-0.1,-0.15),text.font=2)

dev.off()

### get some colony Id's for potential images -----------------------------------
find.scenarios=expand.grid(list(light=c("low","high"),temp=c("low","high"),bleach=c(0,1)))
colony.list=list()
for(g in 1:length(use.groups)){
 load(file=paste("FSS_out_no_zero_v1_wT",use.groups[g],".RData"))
dat.c=na.omit(out.list$used.data)
head(dat.c)
colony.list.g=list()
for(r in 1:nrow(find.scenarios)){
     bleaching.state=find.scenarios[r,"bleach"]
     light.state=find.scenarios[r,"light"]
     temp.state=find.scenarios[r,"temp"]
     if(light.state=="high" & temp.state=="low"){
          index=which((dat.c$mortality/dat.c$trials)>0.1 &
            dat.c$light<0.5 &
            dat.c$temperature<27 &
            dat.c$bleaching==bleaching.state)}
     if(light.state=="low" & temp.state=="low"){
          index=which((dat.c$mortality/dat.c$trials)>0.1 &
            dat.c$light>0.7 &
            dat.c$temperature<27 &
            dat.c$bleaching==bleaching.state)}
     if(light.state=="high" & temp.state=="high"){
          index=which((dat.c$mortality/dat.c$trials)>0.1 &
            dat.c$light<0.5 &
            dat.c$temperature>30 &
            dat.c$bleaching==bleaching.state)}
     if(light.state=="low" & temp.state=="high"){
          index=which((dat.c$mortality/dat.c$trials)>0.1 &
            dat.c$light>0.7 &
            dat.c$temperature>30 &
            dat.c$bleaching==bleaching.state)}
     dat.g=dat.c[index,]
     dat.g=dat.g[order(dat.g$mortality,decreasing=T),]
     colony.list.g=c(colony.list.g,list(paste(dat.g$ColID,dat.g$FSN)))
}

names(colony.list.g)=paste(use.groups[g],
                         "_light.",find.scenarios$light,
                         "_temp.",find.scenarios$temp,
                         "_bleach.",find.scenarios$bleach,sep="")
colony.list=c(colony.list,colony.list.g)
}
colony.list
save(colony.list,file="colony_lists_for_pia.RData")

load("colony_lists_for_pia.RData")

### table of stats output ------------------------------------------------------
numeric.vars=c("delta.AICc","delta.BIC",
            "wi.AICc","wi.BIC","edf")
out.stats.tab=do.call("rbind",lapply(top.all,FUN=function(x){
            out=x[,c("modname",numeric.vars)]}))
out.stats.tab[,numeric.vars]=signif(out.stats.tab[,numeric.vars],2)
out.stats.tab$modname=gsub(".te."," x ",out.stats.tab$modname,fixed=T)
out.stats.tab$modname=gsub("+"," + ",out.stats.tab$modname,fixed=T)
out.stats.tab$modname=gsub("use.group","Taxa/group",out.stats.tab$modname)
out.stats.tab$modname=gsub("light","Light",out.stats.tab$modname)
out.stats.tab$modname=gsub("temperature","Temperature",out.stats.tab$modname)
out.stats.tab$modname=gsub("sedimentation","Sedimentation index",out.stats.tab$modname)
out.stats.tab$modname=gsub("bleaching","Bleached",out.stats.tab$modname)
out.stats.tab
write.csv(out.stats.tab,file="model_table.csv")


########### time series plot ###################################################
require(XLConnect)
cyc.hist <- na.omit(readWorksheet(
 loadWorkbook("X:/Wheatstone_analysis/Cyclone_Data/Cyclone_History.xlsx", create = FALSE),
                  sheet="Sheet1"))
cyc.hist$Start.j.date=julian(as.Date(cyc.hist$Start),origin=dredge.date)
cyc.hist$End.j.date=julian(as.Date(cyc.hist$End),origin=dredge.date)
plot.cyclones=function(cyc.hist,y.lim)
 for(r in 1:nrow(cyc.hist)){
  polygon(c(cyc.hist$Start.j.date[r],
            cyc.hist$Start.j.date[r],
            cyc.hist$End.j.date[r],
            cyc.hist$End.j.date[r]),
        c(min(y.lim),max(y.lim),max(y.lim),min(y.lim)),
               col="grey",border="grey",xpd=F)}

head(dat)
head(wq.dat)

# coral head data
dat$site.category=NA
dat$site.category[dat$dist.d<2]="near"
dat$site.category[dat$dist.d>2 & dat$dist.d<20]="intermediate"
dat$site.category[dat$dist.d>20]="far"
dat$site.category=factor(dat$site.category)


detach("package:boot", unload=TRUE)
require(car)
dat$logit.prop.mort=logit(dat$prop.mort)
dat$cols="red"
dat$cols[which(dat$site.category=="intermediate")]="black"
dat$cols[which(dat$site.category=="far")]="blue"

site.dat=unique(dat[,c("Site.Code","site.category","cols","dist.d")])
wq.new=merge(site.dat,wq.dat)
# make a summary across the three types of sites for plotting
# mean ntu, dli, assd, temp
wq.tt=summaryBy(NTU.wc.14d+DLI.wc.14d+temp+SAS.index.60d ~ j.date+site.category+cols,
                data=wq.new,keep.names=T,na.rm=T,FUN=mean)

wq.all=summaryBy(NTU.wc.14d + DLI.wc.14d + temp+ SAS.index.60d ~ j.date,
                data=wq.new,
                keep.names=T,na.rm=T,FUN=mean)
wq.all$date=as.Date(wq.all$j.date,origin=dredge.date)
#write.csv(wq.all[which(wq.all$j.date>-100 & wq.all$j.date<550),], file="wq_summaries.csv")

#AHC BAT DUG LNG0 LNG1 LNG2 LNG3 LNGA LNGC MOF1 MOF3 MOFA MOFB REFN REFS SBS TR
plot(1:18,1:18,pch=1:18,cex=1.5)
pch.vec=rep(c(1,6,15,17,4,3,2,8,7,5,16),2)
dat$pchs=pch.vec[(dat$Site.Code)]
dat$index=1

col.dat$pchs=pch.vec[(col.dat$Site.Code)]

legend.dat=unique(dat[c("site.category","Site.Code","cols","pchs")])
legend.dat=legend.dat[order(legend.dat$site.category),]

dist.range=round(log10(c(0.1,50)+1)*10)#round(range(log10(dat$dist.d+1),na.rm=T)*10)
dist.vec=dist.range[1]:dist.range[2]

cols=adjustcolor(colorRampPalette(c("black","red","purple","blue","cyan","green"))(length(dist.vec)),
               alpha.f=0.75)#rainbow(length(dist.vec),start=0,end=.7)#rgb( ramp(seq(0, 1, length = length(dist.vec))), max = 255)

pia.dat=dat[which(dat$use.group=="Porit.ms"),]
pia.dat$dist.cols=cols[match(round(log10(pia.dat$dist.d+1)*10),dist.vec)]

pia.dat=merge(pia.dat,legend.dat)

pia.summary=summaryBy(prop.smoth+j.date~Site.Code+cols+FSN,
    data=pia.dat,keep.names=T)
pia.near=summaryBy(prop.smoth+j.date~FSN,
    data=pia.dat[which(pia.dat$site.category=="near"),],keep.names=T)

x.lim=c(0,570)
range(na.omit(dat.gam)$j.date)
date.seq=seq(0,570,by=14)

ntu.ylim=range(wq.tt$NTU.wc.14d,na.rm=T)
dli.ylim=range(wq.tt$DLI.wc.14d,na.rm=T)
SAS.index.60d.ylim=range(wq.tt$SAS.index.60d,na.rm=T)
temp.ylim=range(wq.all$temp,na.rm=T)


table(dat$mort.index)
# only use sites during a survey periods with at least 10 colonies
dat=dat[which(dat$mort.index>10),]


require(gamm4)
time.gams=list()
for(g in 1:length(use.groups)){
 dat.g=dat[which(dat$use.group==use.groups[g]),]
 fit=gamm4(logit.prop.mort~s(j.date,bs="cr",by=site.category)+site.category,
              random=~(1|Site.Code),
              family="gaussian",data=na.omit(dat.g))
  time.gams=c(time.gams,list(fit))}

types=c("near","intermediate","far")#as.character(unique(wq.tt$site.category))
col.types=c("red","black","blue")
names(col.types)=types

dredge.date=as.Date("2010-05-19")

pdf("time_series.pdf",height=7.5,width=6.7,pointsize=10)

par(mfrow=c(5,1),mar=c(0,4,0,0.2),oma=c(4,0,6,6),bty="n",lwd=0.75)
# temperature
plot(wq.all$j.date,wq.all$temp,pch=NA,ylab="",xlab="",axes=F,xlim=x.lim,ylim=temp.ylim)
axis(side=3,at=date.seq,labels=format(as.Date(date.seq,origin=dredge.date),"%d %b %y"),
             las=2,xpd=NA,lwd=0.75)
plot.cyclones(cyc.hist,temp.ylim)
lines(wq.all$j.date,wq.all$temp,col="orange")
axis(side=4,col="orange",lwd=0.75,line=-2)
mtext(side=4,outer=F,text="Temperature (daily mean)",line=0,cex=0.8)
# ntu
par(new=T)
plot(NA,xlim=x.lim,ylim=ntu.ylim,pch=NA,axes=F,ylab="",xlab=NA)
axis(side=2,lwd=0.75)
mtext(side=2,outer=F,text="NTU (14 d mean)",line=2,cex=0.8)
lapply(types,FUN=function(x){
   d=wq.tt[which(wq.tt$site.category==x),]
   d=d[order(d$j.date),]
   lines(jitter(d$j.date),d$NTU.wc.14d,col=col.types[names(col.types)==x])})
# dli
plot(NA,xlim=x.lim,ylim=dli.ylim,pch=NA,axes=F,ylab="",xlab=NA)
axis(side=2,lwd=0.75)
mtext(side=2,outer=F,text="DLI (14 d mean)",line=2,cex=0.8)
lapply(types,FUN=function(x){
   d=wq.tt[which(wq.tt$site.category==x),]
   d=d[order(d$j.date),]
   lines(jitter(d$j.date),d$DLI.wc.14d,col=col.types[names(col.types)==x])})
legend("bottomright",legend=types,col=col.types,lty=1,ncol=1,bty="n",xpd=NA)
# SAS.index.1d
plot(NA,xlim=x.lim,ylim=SAS.index.60d.ylim,pch=NA,axes=F,ylab="",xlab=NA)
axis(side=2,lwd=0.75)
mtext(side=2,outer=F,text="SAS (60 d mean)",line=2,cex=0.8)
lapply(types,FUN=function(x){
   d=wq.tt[which(wq.tt$site.category==x),]
   d=d[order(d$j.date),]
   lines(jitter(d$j.date),d$SAS.index.60d,col=col.types[names(col.types)==x])})
# Mortality
for(g in 1:length(use.groups)){
  plot(NA,xlim=x.lim,ylim=c(0,0.9),pch=NA,axes=F,ylab="",xlab=NA)
  axis(side=2,lwd=0.75)
  mean.dat.g=na.omit(dat[which(dat$use.group==use.groups[g]),])
  fit.g=time.gams[[g]]
  for(t in 1:length(types)){
   new.dat=data.frame(list(j.date=seq(x.lim[1],x.lim[2],length=100),
             site.category=types[t]))
   eta=pred.back.logit(fit.g$gam,new.dat,logit.type="car")
   polygon(c(eta$j.date,rev(eta$j.date)),
          c(eta$P.lw,rev(eta$P.up)),
          col=adjustcolor(col.types[names(col.types)==types[t]],alpha=0.2),border=NA)
   lines(eta$j.date,eta$P,col=col.types[names(col.types)==types[t]])
  }
  points(mean.dat.g$j.date,mean.dat.g$prop.mort,col=mean.dat.g$cols,pch=mean.dat.g$pchs)
  legend("topleft",legend=use.groups[g],bty="n")

  # bleaching
  par(new=T)
  bleach.ts=summaryBy(prop.bleach+j.date~FSN,
                      FUN=mean, keep.names=T,
                      data=dat[which(dat$use.group==use.groups[g]),])
  plot(bleach.ts$j.date,bleach.ts$prop.bleach,pch=NA,ylab="",xlab="",
       axes=F,xlim=x.lim,ylim=c(0,30))
  lines(bleach.ts$j.date,bleach.ts$prop.bleach*100,col="orange",lty=1)
  points(bleach.ts$j.date,bleach.ts$prop.bleach*100,col="orange",pch=16)
  axis(side=4,col="orange",lwd=0.75,line=-2)
  if(g==2){mtext(side=4,outer=F,text="% Bleached",line=0,cex=0.8)}
}
axis(side=1)
mtext(side=1,text="Days since start of dredging",outer=T,line=2.5,cex=0.75)
mtext(side=2,outer=F,text="% Mortality",cex=0.8,line=2)
legend("right",legend=legend.dat$Site.Code,col=legend.dat$cols,pch=legend.dat$pchs,
     ncol=1,bty="n",xpd=NA,inset=-0.1)


dev.off()

################################################################################
###### exploring_additive_v_synergistic_effects_bayesian_v3.R
require(rstanarm)
### Bleaching model ######  -----------------------------------------------------
model.bleaching=out.list.bleach$success.models$'use.group+light.te.temperature'
dat.c=out.list.bleach$used.data
dat.c$scale.light=scale(dat.c$light)
dat.c$scale.temperature=scale(dat.c$temperature)
scale.vals.light=attributes(scale(dat.c$light))
scale.vals.temperature=attributes(scale(dat.c$temperature))
scale.vals=list(
                   scale.vals.light,
                   scale.vals.temperature)
scale.vec=unlist(lapply(scale.vals,FUN=function(x){x$'scaled:scale'}))
centre.vec=unlist(unlist(lapply(scale.vals,FUN=function(x){x$'scaled:center'})))
#
#stan.model.bleaching <- stan_gamm4(cbind(bleach.event,b.trials-bleach.event)~
#             s(scale.light,scale.temperature,k=5)+use.group,
#                 data = dat.c,
#                 family="binomial",
#                 random = ~ (1 | Site.Code),
#                 chains = 1, iter = 10000)
##plot_nonlinear(stan.mod)
#x1=predict(model.bleaching$gam,type="response",re.form=NA)
#x3=apply(posterior_predict(stan.model.bleaching,re.form=NA),MARGIN=2,FUN=mean) / dat.c$b.trials
#plot(x1,x3,col=dat.c$use.group)
#abline(c(0,1))
#
### mortality model ######  -----------------------------------------------------
#model.acro.mort=out.list.prob.any.acro$success.models$'bleaching+light.by.bleaching+temperature.by.bleaching'
##'light+bleaching+temperature.by.bleaching'
##'bleaching+light.by.bleaching+temperature.by.bleaching'
#dat.c=out.list.prob.any.acro$used.data
#stan.model.acro.mort <- stan_gamm4(mort.event ~
#    s(light, by = bleaching, k = 5, bs = "cr") +
#    s(temperature, by = bleaching, k = 5, bs = "cr") +
#    bleaching,
#                 data = dat.c,
#                 family="binomial",
#                 random = ~ (1 | ColID),
#                 chains = 1, iter = 10000)
##plot_nonlinear(stan.mod)
#x1=predict(model.acro.mort$gam,type="response",re.form=NA)
#x3=apply(posterior_predict(stan.model.acro.mort,re.form=NA),MARGIN=2,FUN=mean)
#plot(x1,x3)
#abline(c(0,1))
#
### proportional loss model ###### -----------------------------------------------
#model.acro.prop=out.list.prop.loss.acro$success.models$'bleaching'
#dat.c=out.list.prop.loss.acro$used.data
#stan.model.acro.prop <- stan_glmer(cbind(mortality, trials-mortality) ~
#     bleaching + (1|ColID) + (1|ID),
#                 data = dat.c,
#                 family="binomial",
#                 chains = 1, iter = 10000)
##plot_nonlinear(stan.mod)
#x1=predict(model.acro.prop$gam,type="response",re.form=NA)
#x3=apply(posterior_predict(stan.model.acro.prop,re.form=NA) / dat.c$trials,MARGIN=2,FUN=mean)
#plot(x1,x3)
#abline(c(0,1))
#
###### Massive Porites  #########################################################
## for massive corals things are a bit more complicated as have to link sedimentation and
## light. For simplicity assume:
## low light stress = low sedimentation
## moderate light stress = moderate sedmentation
## high light stress = high sedimentation
#
### mortality model ######  -----------------------------------------------------
#model.porit.mort=out.list.prob.any.porit$success.models$'sedimentation+bleaching+temperature.by.bleaching'
##'sedimentation+temperature+bleaching'
##'sedimentation+bleaching+temperature.by.bleaching'
#
#dat.c=out.list.prob.any.porit$used.data
#stan.model.porit.mort <- stan_gamm4(mort.event ~
#    s(sedimentation, k = 5, bs = "cr") +
#    s(temperature, by = bleaching, k = 5, bs = "cr") +
#    bleaching,
#                 data = dat.c,
#                 family="binomial",
#                 random = ~ (1 | ColID),
#                 chains = 1, iter = 10000)
##plot_nonlinear(stan.mod)
#x1=predict(model.porit.mort$gam,type="response",re.form=NA)
#x3=apply(posterior_predict(stan.model.porit.mort,re.form=NA),MARGIN=2,FUN=mean)
#plot(x1,x3)
#abline(c(0,1))
#
### proportional loss model ###### -----------------------------------------------
#model.porit.prop=out.list.prop.loss.porit$success.models$'light+bleaching'
#dat.c=out.list.prop.loss.porit$used.data
#stan.model.porit.prop <- stan_gamm4(cbind(mortality, trials-mortality) ~
#        s(light, k = 5, bs = "cr") + bleaching,
#                 data = dat.c,
#                 family="binomial",
#                 random = ~ (1 | ColID) + (1|ID),
#                 chains = 1, iter = 10000)
##plot_nonlinear(stan.mod)
#x1=predict(model.porit.prop$gam,type="response",re.form=NA)
#x3=apply(posterior_predict(stan.model.porit.prop,re.form=NA) / dat.c$trials,MARGIN=2,FUN=mean)
#plot(x1,x3)
#abline(c(0,1))
#
#### model sedimentation as a function of light stress --------------------------
#head(dat)
#sm.dat=summaryBy(light.max.14days.prev+sas.max.14days.prev~Site.Code+FSN,
#                 FUN=mean,data=dat,keep.names=T)
#colnames(sm.dat)=c("Site.Code","FSN","light","sedimentation")
#sm.dat=na.omit(sm.dat)
#sm.dat=sm.dat[order(sm.dat$light),]
#
#sm.dat$logit.light=logit(sm.dat$light)
#sm.dat$logit.sedimentation=logit(sm.dat$sedimentation)
#
##require(INLA)
##sm.fit=inla(logit.sedimentation~logit.light+f(FSN,model="iid"),data=sm.dat)
#sm.fit=lmer(logit.sedimentation~logit.light+(1|FSN),data=sm.dat)
#
#new.data.sm=expand.grid(list(logit.light=logit(seq(0.3,0.9,length=100))))
#eta=predict(sm.fit,newdata=new.data.sm,re.form=NA)
#pdf("sed_v_light_gam.pdf")
#plot(sm.dat$light,sm.dat$sedimentation,pch=16,
#     ylab="Sedimentation index", xlab="Light stresss")
#lines(inverse.logit(new.data.sm$logit.light,car.rescale=F),
#      inverse.logit(eta,car.rescale=F),col=2)
#axis(side=3,at=seq(0.3,0.9,0.1),labels=signif(dli.xform(seq(0.3,0.9,0.1)),2))
#mtext(side=3,text="DLI",line=2.5,xpd=NA,cex=0.75)
#dev.off()
#
#sm.fit.stan=stan_lmer(logit.sedimentation~logit.light+(1|FSN),
#                      data=sm.dat,
#                      chains = 1, iter = 10000)
#x1=predict(sm.fit,type="response",re.form=NA)
#x3=apply(posterior_predict(sm.fit.stan,re.form=NA),MARGIN=2,FUN=mean)
#plot(x1,x3)
#abline(c(0,1))
#

#head(wq.dat)
#plot(wq.dat$temp,wq.dat$light.stress.wc.14d)
mean(wq.dat$light.stress.wc.14d[which(wq.dat$temp>=31 & wq.dat$distance>19.6)],na.rm=T)
max(wq.dat$light.stress.wc.14d[which(wq.dat$temp>=31 & wq.dat$distance>19.6)],na.rm=T)
#light.mean.dat=wq.dat[which(wq.dat$temp>=31 & wq.dat$distance>19.6),]
#                    # find the mean light under a high thermal stress environment, just for control sites
#light.mean.dat$logit.light=logit(light.mean.dat$light.stress.wc.14d)
#light.mean.dat$Site.Code
#light.mean.fit.stan=stan_glm(logit.light~1,
#                      data=light.mean.dat,
#                      chains = 1, iter = 10000)
#
#save(stan.model.bleaching,
#       stan.model.acro.mort,stan.model.acro.prop,
#       stan.model.porit.mort,stan.model.porit.prop,
#       sm.fit.stan,
#       light.mean.fit.stan,
#       file="stan_model_fits.Rdata")
require(rstanarm)
load("stan_model_fits.Rdata")

##### Now model for a continuous scale of light (dredging) pressure  ###########
models.mort=list(stan.model.acro.mort,stan.model.porit.mort)
names(models.mort)=use.groups
models.prop=list(stan.model.acro.prop,stan.model.porit.prop)
names(models.prop)=use.groups
mod.list=list(bleaching=stan.model.bleaching,mort=models.mort,prop=models.prop,sm=sm.fit.stan)

# set up the prediction data frame
light<<-seq(0.2,0.9,by=0.01)
temperature.conditions=c(25,31)
new.dat=expand.grid(list(light=light,temperature=temperature.conditions))
tt=scale(new.dat[,c("light","temperature")],centre.vec,scale.vec)
colnames(tt)=c("scale.light","scale.temperature")
new.dat=cbind(new.dat,tt)

tt=unique(dat.use[,c("ColID","use.group","index")])
n.colonies=summaryBy(index~use.group, data=tt, FUN=sum)
#
##g=1
##niter=5000
## run 10 iterations to check
#ptm <- proc.time()
#new.dat.all.iter.all=list(list(),list())
#for(g in 1:length(use.groups)){
# new.dat.all.iter=run.iters(new.dat.all.iter.all[[g]],
#                            niter=10,
#                            n.colony.reps=5000,
#                            new.dat.base=new.dat,
#                            mod.list=mod.list,
#                            group.g=use.groups[g],
#                            ncores=4)
# new.dat.all.iter.all[[g]]=new.dat.all.iter
#}
#(ptm - proc.time())["elapsed"]/60
##new.dat.all.iter.all
#
## update with another 100 iterations
#load("new_dat_prediction_iterations.RData")
#ptm <- proc.time()
#for(g in 1:length(use.groups)){
# new.dat.all.iter=run.iters(new.dat.all.iter.all[[g]],
#                            niter=1770,
#                            n.colony.reps=5000,
#                            new.dat.base=new.dat,
#                            mod.list=mod.list,
#                            group.g=use.groups[g],
#                            ncores=4)
# new.dat.all.iter.all[[g]]=new.dat.all.iter
#}
#(ptm - proc.time())["elapsed"]/60
##new.dat.all.iter.all
#save(new.dat.all.iter.all,file="new_dat_prediction_iterations.RData")

# now simulation for bleaching only impacts
natural.light.val=mean(wq.dat$light.stress.wc.14d[which(wq.dat$temp>=31 & wq.dat$distance>19.6)],na.rm=T)
#posterior_predict(light.mean.fit.stan,WTonly.new.dat,draws=1000)

WTonly.new.dat=expand.grid(list(light=natural.light.val,temperature=31))
tt=scale(WTonly.new.dat[,c("light","temperature")],centre.vec,scale.vec)
colnames(tt)=c("scale.light","scale.temperature")
WTonly.new.dat=cbind(WTonly.new.dat,tt)

#g=1
#niter=5000
# run 10 iterations to check
#ptm <- proc.time()
#WTonly.new.dat.all.iter.all=list(list(),list())
#for(g in 1:length(use.groups)){
# WTonly.new.dat.all.iter=run.iters(WTonly.new.dat.all.iter.all[[g]],
#                            niter=10,
#                            n.colony.reps=5000,
#                            new.dat.base=WTonly.new.dat,
#                            mod.list=mod.list,
#                            group.g=use.groups[g],
#                            ncores=4)
# WTonly.new.dat.all.iter.all[[g]]=WTonly.new.dat.all.iter
#}
#(ptm - proc.time())["elapsed"]/60
#save(WTonly.new.dat.all.iter.all,file="WTonly.new_dat_prediction_iterations.RData")
# update with another 100 iterations
#load("WTonly.new_dat_prediction_iterations.RData")
#ptm <- proc.time()
#for(g in 1:length(use.groups)){
# WTonly.new.dat.all.iter=run.iters(WTonly.new.dat.all.iter.all[[g]],
#                            niter=2000,
#                            n.colony.reps=5000,
#                            new.dat.base=WTonly.new.dat,
#                            mod.list=mod.list,
#                            group.g=use.groups[g],
#                            ncores=4)
# WTonly.new.dat.all.iter.all[[g]]=WTonly.new.dat.all.iter
#}
#(ptm - proc.time())["elapsed"]/60
#save(WTonly.new.dat.all.iter.all,file="WTonly.new_dat_prediction_iterations.RData")

### now process the simulation outputs------------------------------------------
load("new_dat_prediction_iterations.RData")
load("WTonly.new_dat_prediction_iterations.RData")
 # natural light is the lowest value o
plot.dat.light=expand.grid(list(light=light))
plot.dat.out=list()
 for(g in 1:length(use.groups)){
 new.dat.all.iter=new.dat.all.iter.all[[g]]
 out.mat=do.call("cbind",lapply(new.dat.all.iter,FUN=function(x){x$total}))

 WTonly.new.dat.all.iter=WTonly.new.dat.all.iter.all[[g]]
 ttWTonlly.out.mat=do.call("cbind",lapply(WTonly.new.dat.all.iter,FUN=function(x){x$total}))
 vec.WT.only=sample(ttWTonlly.out.mat,ncol(out.mat),replace=F)
 # now separate the posterior samples for the two thermal stress regimes
 out.mat.WT=out.mat[which(new.dat$temperature==max(temperature.conditions)),]
 out.mat.NT=out.mat[which(new.dat$temperature==min(temperature.conditions)),]
 # thermal stress in the absence of dredging stress

 out.mat.WT.only=matrix(rep(NA,length(out.mat.NT)),ncol=ncol(out.mat.WT))
 for(f in 1:nrow(out.mat.WT.only)){out.mat.WT.only[f,]=vec.WT.only}
 # additive effect of dredging and thermal stress
 out.mat.additive=out.mat.NT+out.mat.WT.only
 # difference between additive and observed
 out.mat.synergistic=out.mat.WT-out.mat.additive
 prob.antagonistic=rowSums(apply(out.mat.synergistic,MARGIN=2,FUN=function(x){
                                 out=rep(0,length(x))
                                 out[which(x<0)]=1
                                 return(out)}))/ncol(out.mat.synergistic)
 prob.synergistic=rowSums(apply(out.mat.synergistic,MARGIN=2,FUN=function(x){
                                 out=rep(0,length(x))
                                 out[which(x>0)]=1
                                 return(out)}))/ncol(out.mat.synergistic)
 # collate the output matrices
 out.mats.all=list(WT=out.mat.WT,NT=out.mat.NT,WT.only=out.mat.WT.only,
                   additive=out.mat.additive,synergistic=out.mat.synergistic)

 # Now calculate summary stats from the posterior
 mean.vals=do.call("cbind",lapply(out.mats.all,FUN=rowMeans))
 up.vals=do.call("cbind",lapply(out.mats.all,FUN=function(x){
    apply(x,MARGIN=1,FUN=quantile,probs=0.975)}))
 lw.vals=do.call("cbind",lapply(out.mats.all,FUN=function(x){
    apply(x,MARGIN=1,FUN=quantile,probs=0.025)}))
 colnames(up.vals)=paste("up",colnames(up.vals),sep=".")
 colnames(lw.vals)=paste("lw",colnames(lw.vals),sep=".")

 plot.dat=cbind(plot.dat.light,mean.vals,up.vals,lw.vals,
                prob.antagonistic,prob.synergistic)
 plot.dat.out=c(plot.dat.out,list(plot.dat))

}

trans=(lapply(plot.dat.out,FUN=function(x){
  light=x$light
  synergistic=x$synergistic

  trans.val=numeric()
  if(max(synergistic)>0 & min(synergistic)<0){ # if the curves goes either side of 0
   trans.val=max(light[which(diff(abs(synergistic)==synergistic)!=0)+1])}
  if(min(synergistic)>0){# if the curve is always above 0
   trans.val=max(light)}
  if(max(synergistic)<0){# if the curve is always below 0
   trans.val=min(light)}

  inflection=mean(c(light[which.min(x$prob.synergistic)],
                    light[which.max(x$prob.antagonistic)]))

  min.synergistic=light[which.min(x$prob.synergistic)]
  return(list(trans.val=trans.val,
              inflection=inflection))}))

pdf("additive_v_synergistic_impacts_continuous_bayesian.pdf",height=5,width=6.3,pointsize=10,onefile=T)
par(mfcol=c(3,2),mar=c(0,0,1,1),oma=c(4,4,1,1),bty="l")

for(g in 1:length(use.groups)){
 plot.dat.g=plot.dat.out[[g]]
plot(plot.dat.g$light,plot.dat.g$NT,pch=NA,
     ylab="",
     xlab="",lty=1,ylim=c(0,0.25),xaxt="n",col="black")
if(g==1){mtext(side=2,outer=F,text="Total combined live coral loss",cex=0.75,line=1.5)}
axis(side=1,at=seq(0.2,0.9,0.2),labels=signif(dli.xform(seq(0.2,0.9,0.2)),2))
# additive effect lines
lines(plot.dat.g$light,plot.dat.g$additive,col="blue")
polygon(c(plot.dat.g$light,rev(plot.dat.g$light)),
          c(plot.dat.g$lw.additive,rev(plot.dat.g$up.additive)),
          col=adjustcolor("blue",alpha=0.2),border=NA)
# Dredging only mortality (No thermal stress)
lines(plot.dat.g$light,plot.dat.g$NT,col="black")
polygon(c(plot.dat.g$light,rev(plot.dat.g$light)),
          c(plot.dat.g$lw.NT,rev(plot.dat.g$up.NT)),
          col=adjustcolor("black",alpha=0.2),border=NA)
# Thermal and dredging mortality
lines(plot.dat.g$light,plot.dat.g$WT,col="red")
polygon(c(plot.dat.g$light,rev(plot.dat.g$light)),
          c(plot.dat.g$lw.WT,rev(plot.dat.g$up.WT)),
          col=adjustcolor("red",alpha=0.2),border=NA)
if(g==1){legend("topleft",legend=c("Dredging only","Cumulative","Additive"),bty="n",
           fill=adjustcolor(c("black","red","blue"),alpha=0.2),
           border=c("black","red","blue"))}

# synergistic effects
plot(plot.dat.g$light,plot.dat.g$synergistic,pch=NA,
     ylim=c(-0.15,0.15),yaxt="n",
     xaxt="n",ylab="",xlab="")
axis(side=1,at=seq(0.2,0.9,0.2),labels=signif(dli.xform(seq(0.2,0.9,0.2)),2))
axis(side=2,at=seq(-0.5,0.15,by=0.05))

lines(plot.dat.g$light,plot.dat.g$synergistic,col="black")
polygon(c(plot.dat.g$light,rev(plot.dat.g$light)),
          c(plot.dat.g$lw.synergistic,rev(plot.dat.g$up.synergistic)),
          col=adjustcolor("black",alpha=0.2),border=NA)
abline(h=0,col="blue")
#abline(v=trans[[g]],lty=3)
if(g==1){mtext(side=2,outer=F,text="Observed cumulative loss - Additive loss",
         line=1.5,xpd=NA,cex=0.75)}

x.pos=0.25
text(x=x.pos,y=0.15,labels="Synergistic",xpd=NA,srt=0,cex=0.75,col="orange")
text(x=0.975,y=0.0,labels="Additive",xpd=NA,srt=0,cex=0.75,col="blue")
text(x=x.pos,y=-0.15,labels="Antagonistic",xpd=NA,srt=0,cex=0.75,col="cyan")
arrows(x0=rep(x.pos,2),y0=c(-0.01,0.01),
       x1=rep(x.pos,2),y1=c(-0.14,0.14),
       xpd=NA,length=0.1,col=c("cyan","orange"))
#par(new=T)
plot(plot.dat.g$light,plot.dat.g$prob.antagonistic,pch=NA,
     ylim=c(0,1),yaxt="n",
     xaxt="n",ylab="",xlab="",bty="l")
lines(plot.dat.g$light,plot.dat.g$prob.synergistic,col="orange")
lines(plot.dat.g$light,plot.dat.g$prob.antagonistic,col="cyan")
axis(side=2)
axis(side=1,at=seq(0.2,0.9,0.2),labels=signif(dli.xform(seq(0.2,0.9,0.2)),2))
abline(v=trans[[g]],lty=3)
abline(h=0.5,col="blue")

if(g==1){mtext(side=2,outer=F,text="Probability of effect type",line=1.5,cex=0.75)}

}
mtext(side=1,outer=T,text="Increasing dredging stress (DLI)",line=1.5,cex=0.75)
dev.off()


###############################################################################
############ distance_decay_partitioned_mortality_v3.R

### Bleaching model ######  -----------------------------------------------------
model.bleaching=out.list.bleach$success.models$'use.group+light.te.temperature'
dat.c=out.list.bleach$used.data
dat.c$scale.light=scale(dat.c$light)
dat.c$scale.temperature=scale(dat.c$temperature)
scale.vals.light=attributes(scale(dat.c$light))
scale.vals.temperature=attributes(scale(dat.c$temperature))
scale.vals=list(
                   scale.vals.light,
                   scale.vals.temperature)
scale.vec=unlist(lapply(scale.vals,FUN=function(x){x$'scaled:scale'}))
centre.vec=unlist(unlist(lapply(scale.vals,FUN=function(x){x$'scaled:center'})))


require(rstanarm)
load("stan_model_fits.Rdata")

##### Now predict mortality as a function of the observed conditions ###########
##### controlling for dredging and thermal stress respectively       ###########
models.mort=list(stan.model.acro.mort,stan.model.porit.mort)
names(models.mort)=use.groups
models.prop=list(stan.model.acro.prop,stan.model.porit.prop)
names(models.prop)=use.groups
mod.list=list(bleaching=stan.model.bleaching,mort=models.mort,prop=models.prop,sm=sm.fit.stan)

# use the colony data
# take only dredging data
head(dat.use)
head(wq.dat)

dist.dat=wq.dat[which(wq.dat$j.date>0),]
dist.dat.use=dist.dat[,c("Site.Code","j.date","distance")]
dist.dat.use$j.date=dist.dat$j.date
dist.dat.use$fn=as.factor(floor(dist.dat.use$j.date/14))

dist.dat.use$light=dist.dat$light.stress.wc.14d
dist.dat.use$sedimentation=dist.dat$SAS.index.60d
dist.dat.use$SSC=dist.dat$NTU.wc.14d

# generate fortnightly maximim values (to be consistent to with original analysis)
tt=summaryBy(SSC+light+sedimentation~Site.Code+fn+distance,
             FUN=max,keep.names=T,na.rm=T,data=dist.dat.use)
tt <- tt[!is.infinite(rowSums(tt[,3:5])),]
tt$SSC
tt$light
tt$sedimentation[which(tt$sedimentation=="-Inf")]=NA
tt=na.omit(tt)

plot(tt$distance,tt$light)
plot(tt$distance,tt$sedimentation)
plot(tt$distance,tt$SSC)

temperature.conditions=c(25,31)
### predictions without thermal stress
# run 10 iterations to check
tt.NT=tt
tt.NT$temperature=25
tt.WT=tt
tt.WT$temperature=31
tt.all=rbind(tt.WT,tt.NT)
tt.1=scale(tt.all[,c("light","temperature")],centre.vec,scale.vec)
colnames(tt.1)=c("scale.light","scale.temperature")
tt.all=cbind(tt.all,tt.1)

model.bleaching
model.acro.mort=out.list.prob.any.acro$success.models$'bleaching+light.by.bleaching+temperature.by.bleaching'
model.acro.prop=out.list.prop.loss.acro$success.models$'bleaching'
model.porit.mort=out.list.prob.any.porit$success.models$'sedimentation+bleaching+temperature.by.bleaching'
model.porit.prop=out.list.prop.loss.porit$success.models$'light+bleaching'

mort.models.freq=list(model.acro.mort,model.porit.mort)
prop.models.freq=list(model.acro.prop,model.porit.prop)


g=1
tt.out=list()
for(g in 1:length(use.groups)){
  tt.g=tt.all
  tt.g$use.group=use.groups[g]
  model.mort.g=mort.models.freq[[g]]
  model.prop.g=prop.models.freq[[g]]

  tt.g$prob.bleaching=predict(model.bleaching,
                              newdata=tt.g,
                              type="response",re.form=NA)
  # probability of mortality given no bleaching
  tt.g$bleaching=0
  tt.g$prob.mortNB=predict(model.mort.g,
                           newdata=tt.g,
                           type="response",re.form=NA)
  # probability of mortality given bleaching
  tt.g$bleaching=1
  tt.g$prob.mortWB=predict(model.mort.g,
                           newdata=tt.g,
                           type="response",re.form=NA)

  # proportional loss given no bleaching
  tt.g$bleaching=0
  tt.g$prob.propNB=predict(model.prop.g,
                           newdata=tt.g,
                           type="response",re.form=NA)
  # probability of mortality given bleaching
  tt.g$bleaching=1
  tt.g$prob.propWB=predict(model.prop.g,
                           newdata=tt.g,
                           type="response",re.form=NA)
  # calculate the total coral loss for the bleached v unbleached pathways
  tt.g$Total.WBbleaching=tt.g$prob.bleaching*tt.g$prob.mortWB*tt.g$prob.propWB
  tt.g$Total.NBbleaching=(1-tt.g$prob.bleaching)*tt.g$prob.mortNB*tt.g$prob.propNB
  # total mortality is the sum of both pathways
  tt.g$Total.loss=tt.g$Total.WBbleaching+tt.g$Total.NBbleaching

  # now add some stuff to the predictor data frame
  tt.g$j.date=as.numeric(as.character(tt.g$fn))*14
  tt.g$date=as.Date(tt.g$j.date,origin=dredge.date)
  require(lubridate)
  tt.g$month=month(tt.g$date)
  tt.g$season=NA
  tt.g$season[which(tt.g$month==12|tt.g$month==1|tt.g$month==2)]="summer"
  tt.g$log.dist=log(tt.g$distance)
  tt.g$trials=100
  tt.g$ID=as.factor(1:nrow(tt.g))
  tt.g$mortality=round(tt.g$Total.loss*100)

  # separete into with and without thermal stress
  NTtt.g=tt.g[which(tt.g$temperature==min(temperature.conditions)),]
  WTtt.g=tt.g[which(tt.g$temperature==max(temperature.conditions)),]
  WTtt.g=WTtt.g[which(WTtt.g$season=="summer"),]

#  # fit distance decay to each
#  stan.gamm.WT=stan_gamm4(cbind(mortality, trials-mortality)~s(log.dist,bs="cr"),
#                   data = WTtt.g,random=~(1|fn)+(1|Site.Code)+(1|ID),
#                   family="binomial",adapt_delta=0.999,
#                   chains = 3, iter = 10000)
#  stan.gamm.NT=stan_gamm4(cbind(mortality, trials-mortality)~s(log.dist,bs="cr"),
#                   data = NTtt.g,random=~(1|fn)+(1|Site.Code)+(1|ID),
#                   family="binomial",adapt_delta=0.999,
#                   chains = 3, iter = 10000)
#
#  save(stan.gamm.WT,stan.gamm.NT,file=paste("stan_gamm4_distance",use.groups[g],".RData",sep="."))
  tt.out=c(tt.out,list(tt.g))
}

dist.range=range(site.dat$dist.d)
x.vec=exp(seq(log(dist.range[1]),log(dist.range[2]),length=2000))

half.dist.vals=list()
pdf("distance_decay.pdf",height=3.5,width=6.3,pointsize=10,onefile=T)
par(mfrow=c(1,2),mar=c(0,0,1,1),oma=c(4,4,1,1),bty="l")
y.lim=c(0,0.14)
for(g in 1:length(use.groups)){
  tt.g=tt.out[[g]]
  # separete into with and without thermal stress
  NTtt.g=tt.g[which(tt.g$temperature==min(temperature.conditions)),]
  WTtt.g=tt.g[which(tt.g$temperature==max(temperature.conditions)),]
  WTtt.g=WTtt.g[which(WTtt.g$season=="summer"),]
  load(paste("stan_gamm4_distance",use.groups[g],".RData",sep="."))
  pred.tt=data.frame(distance=x.vec)
  pred.tt$mortality=0
  pred.tt$trials=5000
  pred.tt$log.dist=log(pred.tt$distance)
  # with thermal stress
  WTout=posterior_predict(stan.gamm.WT,re.form=NA,newdata=pred.tt)/5000
  WTmean=apply(WTout,MARGIN=2,FUN=median)
  WTup=apply(WTout,MARGIN=2,FUN=quantile,probs=0.975)
  WTlw=apply(WTout,MARGIN=2,FUN=quantile,probs=0.025)
  mean.min=min(WTmean)
  mean.max=max(WTmean)
  mean.dif=mean.max-mean.min
  mean.50=(mean.dif*0.5)+mean.min
  mean.10=(mean.dif*0.1)+mean.min
  # ED50
  mean.ED50=(x.vec[which.min(abs(WTmean-mean.50))])
  up.ED50=(x.vec[which.min(abs(WTup-mean.50))])
  lw.ED50=(x.vec[which.min(abs(WTlw-mean.50))])
  # ED10
  mean.ED10=x.vec[which.min(abs(WTmean-mean.10))]
  up.ED10=x.vec[which.min(abs(WTup-mean.10))]
  lw.ED10=x.vec[which.min(abs(WTlw-mean.10))]
  WT.EDvals=c(mean.50,mean.ED50,up.ED50,lw.ED50,
              mean.10,mean.ED10,up.ED10,lw.ED10)
  # without thermal stress
  NTout=posterior_predict(stan.gamm.NT,re.form=NA,newdata=pred.tt)/5000
  NTmean=apply(NTout,MARGIN=2,FUN=median)
  NTup=apply(NTout,MARGIN=2,FUN=quantile,probs=0.975)
  NTlw=apply(NTout,MARGIN=2,FUN=quantile,probs=0.025)
  mean.min=min(NTmean)
  mean.max=max(NTmean)
  mean.dif=mean.max-mean.min
  mean.50=(mean.dif*0.5)+mean.min
  mean.10=(mean.dif*0.1)+mean.min
  # ED50
  mean.ED50=x.vec[which.min(abs(NTmean-mean.50))]
  up.ED50=x.vec[which.min(abs(NTup-mean.50))]
  lw.ED50=x.vec[which.min(abs(NTlw-mean.50))]
  # ED10
  mean.ED10=x.vec[which.min(abs(NTmean-mean.10))]
  up.ED10=x.vec[which.min(abs(NTup-mean.10))]
  lw.ED10=x.vec[which.min(abs(NTlw-mean.10))]
  NT.EDvals=c(mean.50,mean.ED50,up.ED50,lw.ED50,
              mean.10,mean.ED10,up.ED10,lw.ED10)

  # collate results
  ED.vals=rbind(WT.EDvals,NT.EDvals)
  colnames(ED.vals)=c("total.loss50","ED50","ED50.up","ED50.lw",
                      "total.loss10","ED10", "ED10.up","ED10.lw")
  rownames(ED.vals)=c("WT","NT")
  half.dist.vals=list(half.dist.vals,list(ED.vals))

  # make the plot
  plot.seq=round(seq(1,length(x.vec),length=100))
  plot(jitter(WTtt.g$distance,factor=1),WTtt.g$Total.loss,
                            #log="x",
                            pch=16, ylim=y.lim, yaxt="n",
                            col=adjustcolor("red",alpha=0.3))
  legend("topleft",legend=paste(letters[g],use.groups[g]),bty="n")
  points(jitter(NTtt.g$distance,factor=1),NTtt.g$Total.loss,
             pch=16,
             col=adjustcolor("black",alpha=0.3))
  lines(x.vec[plot.seq],WTmean[plot.seq],col="red")
  polygon(c(pred.tt$distance,rev(pred.tt$distance)),
            c(WTup,rev(WTlw)),
            col=adjustcolor("red",alpha=0.2),border=NA)
  lines(x.vec[plot.seq],NTmean[plot.seq],col="black")
  polygon(c(pred.tt$distance,rev(pred.tt$distance)),
            c(NTup,rev(NTlw)),
            col=adjustcolor("black",alpha=0.2),border=NA)

  abline(h=ED.vals[,"total.loss50"],col=c("red","black"),lty=3)
  abline(v=ED.vals[,"ED50"],col=c("red","black"),lty=3)
 if(g==1){
   axis(side=2)
   mtext(side=2,outer=T,line=2,text="Proportional live coral loss")
 }
}
mtext(side=1,outer=T,line=2,text="Distance (km)")
dev.off()

names(half.dist.vals)=use.groups
half.dist.vals


