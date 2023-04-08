###Simulation and animation for Figure 1 Cumming: "The New Statistics" 2014
##library loads###
require(ggplot2)
require(gganimate)
require(dplyr)
##########
###Simulation Set up#####
set.seed(12346) #set seed for reproducibility
###populations##
Pop1<-rnorm(1000, mean = 10, sd = 20)
Pop2<-rnorm(1000, mean = 20, sd = 20)
mean(Pop2-Pop1)### ~10 unit difference

####create 25 replications at N=32 (per set up from Cumming)
repfunction<-function(j=1){
repsout<-do.call(rbind,lapply(1:25,function(i){
  samplepop1<-sample(Pop1,size=32)
  samplepop2<-sample(Pop2,size=32)
  
  #thisrepdiff<-mean(samplepop2-samplepop1)
  thisrepp<-t.test(c(samplepop1,samplepop2)~c(rep(2,length(samplepop2)),rep(1,length(samplepop1))))
  returndf<-data.frame(i=i,diff=as.numeric(thisrepp$estimate[1]-thisrepp$estimate[2]),
                       pval=thisrepp$p.value,
                       CIlow=thisrepp$conf.int[1],
                       CIhigh=thisrepp$conf.int[2])
  return(returndf)
}))
repsout$groupnumb<-j
return(repsout)
}
#####repeat 25 replications for N groups for animation
allrepsout<-do.call(rbind,lapply(1:10,function(j){repfunction(j)}))
allrepsout$CI_include_pop_est<-as.factor(dplyr::if_else(allrepsout$CIlow<=10 & allrepsout$CIhigh>=10,1,0))
####plot#####
##plot cleanup func##
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
###plots##
ggCIfig<-ggplot()+geom_pointrange(data=allrepsout,aes(x=diff,xmin=CIlow, xmax=CIhigh,y=i,
                                                      colour=CI_include_pop_est,fill=CI_include_pop_est),shape=21)+
  scale_fill_manual(values=c("#cb4154","grey77"))+scale_colour_manual(values=c("#cb4154","black"))+
  geom_vline(xintercept = 0)+geom_vline(xintercept = 10)+
  scale_x_continuous(breaks=c(-20,-10,0,10,20,30))+
  transition_states(groupnumb, transition_length = 0, state_length = 100, wrap = TRUE)
ggCIfig<-plotclean(ggCIfig)+theme(legend.position = "none",
                                  axis.title.y=element_blank(),
                                  axis.text.y=element_blank(),
                                  axis.ticks.y=element_blank(),
                                  axis.line.y = element_blank())+
  xlab("Difference Between Means")

animate(ggCIfig, nframes = 200)

ggCIfigbg<-ggplot()+geom_pointrange(data=allrepsout,aes(x=diff,xmin=CIlow, xmax=CIhigh,y=i,
                                                      colour=CI_include_pop_est,fill=CI_include_pop_est),shape=21)+
  scale_fill_manual(values=c("#cb4154","grey77"))+scale_colour_manual(values=c("#cb4154","black"))+
  geom_vline(xintercept = 0,linetype="dotted",colour="grey66")+geom_vline(xintercept = 10,linetype="dotted",colour="grey66")+
  scale_x_continuous(breaks=c(-20,-10,0,10,20,30))+
  transition_states(groupnumb, transition_length = 0, state_length = 600, wrap = TRUE)
ggCIfigbg<-plotclean(ggCIfigbg)+theme(legend.position = "none",
                                  axis.title.y=element_blank(),
                                  axis.text.y=element_blank(),
                                  axis.ticks.y=element_blank(),
                                  axis.line.y = element_blank())+
  xlab("Difference Between Means")

animate(ggCIfigbg, nframes = 400)

anim_save("~/Library/Mobile\ Documents/com~apple~CloudDocs/Presentations/SpringJobTalks/PittPsychology/Newstatsfig1.sim.gif")










