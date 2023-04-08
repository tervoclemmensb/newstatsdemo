###Simulation and animation for Sampling Variability of a Correlation####
##library loads###
require(ggplot2)
require(gganimate)
require(dplyr)
require(rockchalk)
require(recipes)
require(emdbook)
require(moments)
plotclean<-function (p, ajust = 0.05) 
{
  p$labels$y <- paste0(p$labels$y, "\n")
  p$labels$x <- paste0("\n", p$labels$x)
  p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                         text = element_text(size = 20), axis.text.y = element_text(hjust = ajust, 
                                                                                    color = "black"), axis.text.x = element_text(vjust = ajust, 
                                                                                                                                 color = "black"))
}
########
######
set.seed(12346) #set seed for reproducibility
rho=.15

myCov=lazyCov(Rho=rho, Sd=1, d=2) #define covariance matrix 
popData=data.frame(mvrnorm(n=1000, mu=c(0,0), Sigma=myCov)) #create two variables with specified covariance matrix 
ggpop<-ggplot(popData,aes(x=X1,y=X2))+geom_point()+xlab("")+ylab("")
ggpop<-plotclean(ggpop)
hist<-ggplot(popData,aes(x=X1))+geom_histogram()+scale_y_continuous(expand=c(0,0))

pop_cor=cor(popData[,1], popData[,2]) 
bins<-round(lseq(25,nrow(popData),20))
iterations=5000
siglevel=.05
outputsubsample=TRUE

subsample<-lapply(bins,function(sb){
  print(sb)
  sbcorrs<-lapply(1:iterations,function(i){
    binsample<-base::sample(nrow(popData),sb,replace=TRUE)
    bincor<-cor(popData[binsample,"X1"],popData[binsample,"X2"])
    pval<-cor.test(popData[binsample,"X1"],popData[binsample,"X2"])$p.value
    datai<-cbind(popData[binsample,"X1"],popData[binsample,"X2"])
    return(cbind(bincor,pval,i,datai))
  })
  sbcofssmat<-do.call(rbind,sbcorrs)
  sbcofssmat_full<-data.frame(sbcofssmat)
  names(sbcofssmat_full)<-c("bincor","pval","iteri","x","y")
  sbcofssmat_full$samplesize<-sb
  sbcofssmat<-sbcofssmat_full[,c("bincor","pval","iteri")] %>% group_by(iteri) %>% dplyr::summarize(bincor=unique(bincor),pval=unique(pval))
  
  sbmean<-mean(sbcofssmat$bincor)
  sbCI=quantile(sbcofssmat$bincor,probs=c(.025,.975))
  
  sigmean<-mean(abs(sbcofssmat$bincor[sbcofssmat$pval<siglevel]))
  sigminimum<-min(abs(sbcofssmat$bincor[sbcofssmat$pval<siglevel]))
  
  sigCI=quantile(sbcofssmat$bincor[sbcofssmat$pval<siglevel],probs=c(.025,.975))
  kurtosis=moments::kurtosis(sbcofssmat$bincor)
  
  returninfo<-c(sb,sbmean,sbCI[1],sbCI[2],sigmean,sigCI[1],sigCI[2],sigminimum,kurtosis)
  names(returninfo)<-c("samplesize","mean","low","high","sigmean","siglow","sighigh","sigminimum","kurtosis")
  
  if (outputsubsample){
    returndat<-list(sbcofssmat_full=sbcofssmat_full,returninfo=returninfo)
  }else{
    returndat<-returninfo
  }
  return(returndat)
})
fullsubsample<-subsample
#########PLOTS#############
subsamplefulldat<-lapply(subsample,"[[",1)
subsamplefulldatdf<-data.frame(do.call(rbind,subsamplefulldat))

ggpop<-ggplot(popData,aes(x=X1,y=X2))+geom_point(alpha=.1)+geom_smooth(se=FALSE,method="lm",colour="black")+xlab("rand var 1")+ylab("rand var 2")+
  geom_text(data=popData,aes(x=-2.25, y=4, label=sprintf("r=%s",rho)),size=8,check_overlap = T)
ggpop<-plotclean(ggpop)+xlab("rand var 1")+ylab("rand var 2")
ggsave(ggpop,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Presentations/SpringJobTalks/PittPsychology/Figures/popeffect.pdf",height = 5.56, width =7.4)

subsamplefulldatdfsim<-subsamplefulldatdf[subsamplefulldatdf$samplesize==30 & subsamplefulldatdf$iteri %in% sample(seq(1,1000),size=20),]
subsamplefulldatdfsim<-subsamplefulldatdfsim %>% group_by(iteri) %>% mutate(cor=cor(x,y))
subsamplefulldatdfsim$cor<-sprintf("r=%s",as.character(round(subsamplefulldatdfsim$cor,3)))

gganimatesubsample<-ggplot()+geom_point(data=popData,aes(x=X1,y=X2),colour="grey33",alpha=.1)+
  geom_point(data=subsamplefulldatdfsim,aes(x=x,y=y),colour="#CE5148")+
  geom_smooth(data=subsamplefulldatdfsim,aes(x=x,y=y),method="lm",colour="#cb4154",se=FALSE)+
  xlab("")+ylab("")+
  geom_text(data = subsamplefulldatdfsim, aes(x=-2.5, y=4, label=cor, group=iteri),size=8,check_overlap = T) + 
  transition_states(iteri, transition_length = 1, state_length = 1, wrap = TRUE)


gganimatesubsample<-plotclean(gganimatesubsample)+xlab("rand var 1")+ylab("rand var 2")

animobj<-animate(gganimatesubsample, height = 400, width =533.333)

anim_save("~/Library/Mobile\ Documents/com~apple~CloudDocs/Presentations/SpringJobTalks/PittPsychology/Figures/samplingvariabilitywithpopeffect.gif")

#################
subsampleplot<-lapply(subsample,"[[",2)
subsampledf<-data.frame(do.call(rbind,subsampleplot))
names(subsampledf)<-c("samplesize","mean","low","high","sigmean","siglow","sighigh","sigminimum","kurtosis")
yscale<-c((range(subsampledf$high)[2]*1.12),(range(subsampledf$low)[1]*1.12))

ggcorbysamplesize<-ggplot(subsampledf,aes(x=samplesize,y=mean))+geom_ribbon(aes(ymin=low,ymax=high),fill="black",alpha=.2)+
  geom_line(data=subsampledf,aes(x=samplesize,y=mean),colour="red",size=2)+xlab("Sample Size")+ylab("Correlation (r)")+
  geom_line(data=subsampledf,aes(x=samplesize,y=low),colour="black",size=1)+
  geom_line(data=subsampledf,aes(x=samplesize,y=high),colour="black",size=1)+
  scale_y_continuous(limits = c(yscale[2],yscale[1]))+
  scale_x_continuous(breaks=c(30,250,500,750,1000))
ggcorbysamplesize<-plotclean(ggcorbysamplesize)+theme(text = element_text(size = 32))
ggsave(ggcorbysamplesize,file=sprintf("~/Library/Mobile\ Documents/com~apple~CloudDocs/Presentations/SpringJobTalks/PittPsychology/Figures/visual1.r%s.n30.pdf",rho),height=8,width=10)


ggcorbysamplesize<-ggplot(subsampledf,aes(x=samplesize,y=mean))+geom_ribbon(aes(ymin=low,ymax=high),fill="black",alpha=.2)+
  geom_line(data=subsampledf,aes(x=samplesize,y=mean),colour="red",size=2)+xlab("Sample Size")+ylab("Correlation (r)")+
  geom_line(data=subsampledf,aes(x=samplesize,y=sigminimum),linetype="dashed",size=1.5)+
  geom_line(data=subsampledf,aes(x=samplesize,y=low),colour="black",size=1)+
  geom_line(data=subsampledf,aes(x=samplesize,y=high),colour="black",size=1)+
  scale_y_continuous(limits = c(yscale[2],yscale[1]))+
  scale_x_continuous(breaks=c(30,250,500,750,1000))
ggcorbysamplesize<-plotclean(ggcorbysamplesize)+theme(text = element_text(size = 32))
ggsave(ggcorbysamplesize,file=sprintf("~/Library/Mobile\ Documents/com~apple~CloudDocs/Presentations/SpringJobTalks/PittPsychology/Figures/visual2.r%s.pdf",rho),height=8,width=10)

ggcorbysamplesize<-ggplot(subsampledf,aes(x=samplesize,y=mean))+geom_ribbon(aes(ymin=low,ymax=high),fill="black",alpha=.2)+
  geom_line(data=subsampledf,aes(x=samplesize,y=mean),colour="red",size=2)+xlab("Sample Size")+ylab("Correlation (r)")+
  #geom_line(data=subsampledf,aes(x=samplesize,y=sigminimum),linetype="dashed",size=1.5)+
  geom_line(data=subsampledf,aes(x=samplesize,y=sigmean),color="blue",size=1.5)+
  geom_line(data=subsampledf,aes(x=samplesize,y=low),colour="black",size=1)+
  geom_line(data=subsampledf,aes(x=samplesize,y=high),colour="black",size=1)+
  scale_y_continuous(limits = c(yscale[2],yscale[1]))
ggcorbysamplesize<-plotclean(ggcorbysamplesize)+theme(text = element_text(size = 32))
ggsave(ggcorbysamplesize,file=sprintf("~/Library/Mobile\ Documents/com~apple~CloudDocs/Presentations/SpringJobTalks/PittPsychology/Figures/visual3.r%s.pdf",rho),height=8,width=10)

ggcorbysamplesize<-ggplot(subsampledf,aes(x=samplesize,y=mean))+geom_ribbon(aes(ymin=low,ymax=high),fill="black",alpha=.2)+
  geom_line(data=subsampledf,aes(x=samplesize,y=mean),colour="red",size=2)+xlab("Sample Size")+ylab("Correlation (r)")+
  geom_line(data=subsampledf,aes(x=samplesize,y=sigminimum),linetype="dashed",size=1.5)+
  geom_line(data=subsampledf,aes(x=samplesize,y=sigmean),color="blue",size=1.5)+
  geom_line(data=subsampledf,aes(x=samplesize,y=low),colour="black",size=1)+
  geom_line(data=subsampledf,aes(x=samplesize,y=high),colour="black",size=1)+
  scale_y_continuous(limits = c(yscale[2],yscale[1]))
ggcorbysamplesize<-plotclean(ggcorbysamplesize)+theme(text = element_text(size = 32))
ggsave(ggcorbysamplesize,file=sprintf("~/Library/Mobile\ Documents/com~apple~CloudDocs/Presentations/SpringJobTalks/PittPsychology/Figures/visual4.r%s.pdf",rho),height=8,width=10)

###################




