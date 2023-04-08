###Simulation and animation for Sampling Variability of a Correlation####
##library loads###
require(ggplot2)
require(gganimate)
require(dplyr)
require(rockchalk)
require(recipes)
require(emdbook)
require(moments)
require(metafor)
require(psych)
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
ggpop<-LNCDR::lunaize(ggpop)
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

subsamplefulldatdfmetasim<-subsamplefulldatdf %>% dplyr::group_by(samplesize,iteri) %>% dplyr::summarize(bincor=unique(bincor),pval=unique(pval)) ###reduce down to just cors, n, p
nstudiesmeta<-seq(5,100,5)

####unbiased meta#####

theseiters<-sample(unique(subsamplefulldatdfmetasim$iteri),size=100)

metabysamplesize<-lapply(bins,function(samplesize){
  #print(samplesize)
metasout<-lapply(nstudiesmeta,function(nmeta){
  #print(nmeta)
  ###note this is a demo for visualization of meta-analytic results. 
  #Consult documentation (e.g.,metafor-package {metafor}) for specific use cases
  ##rma used for general instruction
  metaiters<-theseiters[seq(1,nmeta)]
  metastudies<-subsamplefulldatdfmetasim[subsamplefulldatdfmetasim$samplesize %in% samplesize &
                                           subsamplefulldatdfmetasim$iteri %in% metaiters,]
  
  escalcdat <- metafor::escalc(measure="ZCOR", ri=bincor, ni=samplesize, data=metastudies)
  metaout<-rma(escalcdat$yi,escalcdat$vi)
  returndata<-escalcdat
  returndata$metacoef<-coef(metaout)
  returndata$ci.lb<-metaout$ci.lb
  returndata$ci.ub<-metaout$ci.ub
  returndata$nmeta<-nmeta
  returndata$samplesize<-samplesize

  return(returndata)
})
metasoutdf<-do.call(rbind,metasout)
})

metabysamplesizedf<-do.call(rbind,metabysamplesize)
metabysamplesizedf[,c("metacoef","ci.lb","ci.ub")]<-lapply(metabysamplesizedf[,c("metacoef","ci.lb","ci.ub")],function(x){fisherz2r(x)}) ##back transform to r
metabysamplesizedf$pop_cor<-pop_cor

gganimatemeta<-ggplot()+geom_point(data=metabysamplesizedf[!(metabysamplesizedf$samplesize %in% c(30,47,45,66,90)),],aes(x=samplesize,y=bincor),alpha=.2,size=4)+
  geom_pointrange(data=metabysamplesizedf[!(metabysamplesizedf$samplesize %in% c(30,47,45,66,90)),],aes(x=samplesize,y=metacoef,ymin=ci.lb, ymax=ci.ub),fatten=2.25,colour="#CE5148",size=1.5)+
  geom_hline(data=metabysamplesizedf,aes(yintercept = pop_cor),colour="#CE5148",alpha=.3,linetype="dotted")+
  transition_states(nmeta, transition_length = 1, state_length = 10)+
  labs(title= 'Studies at each Sample Size: {closest_state}')
  
gganimatemeta<-plotclean(gganimatemeta)+xlab("\n Sample Size")+ylab("Correlation (r)\n")+theme(text = element_text(size = 32))

animobj<-animate(gganimatemeta, height = 600, width =750)

anim_save("~/Library/Mobile\ Documents/com~apple~CloudDocs/Presentations/SpringJobTalks/PittPsychology/Figures/metasim.gif")

####publication biased meta#####

####unbiased meta#####
significantitersbysamplesize<-subsamplefulldatdfmetasim[subsamplefulldatdfmetasim$pval <.05,]
nonsignificantitersbysamplesize<-subsamplefulldatdfmetasim[subsamplefulldatdfmetasim$pval >.05,]

theseiters<-sample(unique(subsamplefulldatdfmetasim$iteri),size=100)
nstudiesmeta<-seq(5,100,5)

biasedmetabysamplesize<-lapply(bins,function(samplesize){
  #print(samplesize)
  metasout<-lapply(nstudiesmeta,function(nmeta){
    #print(nmeta)
    ###note this is a demo for visualization of meta-analytic results. 
    #Consult documentation (e.g.,metafor-package {metafor}) for specific use cases
    ##rma used for general instruction
    metaiters<-theseiters[seq(1,nmeta)]
    metastudies<-subsamplefulldatdfmetasim[subsamplefulldatdfmetasim$samplesize %in% samplesize &
                                             subsamplefulldatdfmetasim$iteri %in% metaiters,]
    escalcdat <- metafor::escalc(measure="ZCOR", ri=bincor, ni=samplesize, data=metastudies)
    returndata<-escalcdat
    returndata$samplesize<-samplesize
    returndata$nmeta<-nmeta
    if(length(which(metastudies$pval<.05))>1){
      escalcdatsig<-escalcdat[escalcdat$pval<.05,]
      metaout<-rma(escalcdatsig$yi,escalcdatsig$vi)
      returndata$metacoef<-coef(metaout)
      returndata$ci.lb<-metaout$ci.lb
      returndata$ci.ub<-metaout$ci.ub
    }else{
      returndata$metacoef<-NA
      returndata$ci.lb<-NA
      returndata$ci.ub<-NA
    }
    return(returndata)
  })
  metasoutdf<-do.call(rbind,metasout)
})

biasedmetabysamplesizedf<-do.call(rbind,biasedmetabysamplesize)
biasedmetabysamplesizedf[,c("metacoef","ci.lb","ci.ub")]<-lapply(biasedmetabysamplesizedf[,c("metacoef","ci.lb","ci.ub")],function(x){fisherz2r(x)}) ##back transform to r
biasedmetabysamplesizedf$pop_cor<-pop_cor
biasedmetabysamplesizedf$sigstudycor<-dplyr::if_else(biasedmetabysamplesizedf$pval<.05 ,"significant","non-significant")

gganimatebiasedmeta<-ggplot()+geom_point(data=biasedmetabysamplesizedf[!(biasedmetabysamplesizedf$samplesize %in% c(30,47,45,66,90)),],aes(x=samplesize,y=bincor,shape=sigstudycor),alpha=.2,size=4)+
  geom_pointrange(data=biasedmetabysamplesizedf[!(biasedmetabysamplesizedf$samplesize %in% c(30,47,45,66,90)),],aes(x=samplesize,y=metacoef,ymin=ci.lb, ymax=ci.ub),fatten=2.25,colour="blue",shape=17,size=1.5)+
  geom_hline(data=biasedmetabysamplesizedf,aes(yintercept = pop_cor),colour="#CE5148",alpha=.3,linetype="dotted")+
  transition_states(nmeta, transition_length = 1, state_length = 10)+
  labs(title= 'Studies at each Sample Size: {closest_state}')

gganimatebiasedmeta<-plotclean(gganimatebiasedmeta)+xlab("\n Sample Size")+ylab("Correlation (r)\n")+theme(text = element_text(size = 32))+theme(legend.position = "none")

animobj<-animate(gganimatebiasedmeta, height = 600, width =750)

anim_save("~/Library/Mobile\ Documents/com~apple~CloudDocs/Presentations/SpringJobTalks/PittPsychology/Figures/biasedmetasim.gif")


