###Simulation and animation for Figure 1 Cumming: "The New Statistics" 2014
##library loads###
require(ggplot2)
require(gganimate)
require(dplyr)
require(rockchalk)
require(recipes)
require(emdbook)
require(magick)
require(tidyr)
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
##########
###Simulation Set up#####
set.seed(12346) #set seed for reproducibility
rho=.15

myCov=lazyCov(Rho=rho, Sd=1, d=2) #define covariance matrix 
popData=data.frame(mvrnorm(n=1000, mu=c(0,0), Sigma=myCov)) #create two variables with specified covariance matrix 
pop_cor=cor(popData[,1], popData[,2]) 

####create 25 replications at N=30 (~per set up from Cumming)
repfunction<-function(j=1){
repsout<-do.call(rbind,lapply(1:25,function(i){
  #thisrepdiff<-mean(samplepop2-samplepop1)
  theseindices<-sample(1:nrow(popData),30)
  thisrepp<-cor.test(popData$X1[theseindices],popData$X2[theseindices],conf.level = .95)
  returndf<-data.frame(i=i,cor=thisrepp$estimate,
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
allrepsout$CI_include_pop_est<-as.factor(dplyr::if_else(allrepsout$CIlow<=pop_cor & allrepsout$CIhigh>=pop_cor,1,0))
####plot#####
##plot cleanup func##

###plots##
ggCIfig<-ggplot()+geom_pointrange(data=allrepsout,aes(x=cor,xmin=CIlow, xmax=CIhigh,y=i,
                                                      colour=CI_include_pop_est,fill=CI_include_pop_est),shape=21)+
  scale_fill_manual(values=c("#4198ca","grey77"))+scale_colour_manual(values=c("#4198ca","black"))+
  geom_vline(xintercept = 0,alpha=.3)+geom_vline(xintercept = pop_cor,colour="#CE5148",linetype="dotted",size=1)+
  transition_states(groupnumb, transition_length = 0, state_length = 100, wrap = TRUE)
ggCIfig<-plotclean(ggCIfig)+theme(legend.position = "none",
                                  axis.title.y=element_blank(),
                                  axis.text.y=element_blank(),
                                  axis.ticks.y=element_blank(),
                                  axis.line.y = element_blank())+
  xlab("Correlation (r)")

aggCIfig<-animate(ggCIfig, nframes = 200)

anim_save("~/Library/Mobile\ Documents/com~apple~CloudDocs/Presentations/SpringJobTalks/PittPsychology/Figures/Newstatscorrelation.sim.gif")

allrepsoutsummary<-allrepsout %>% dplyr::group_by(groupnumb) %>% 
  dplyr::summarize(percentincludepop=sum(as.numeric(as.character(CI_include_pop_est)))/n(),percentsig=length(which(pval<.05))/n())

allrepsoutsummary_long<-tidyr::pivot_longer(allrepsoutsummary,cols=c("percentincludepop","percentsig"),names_to="type",values_to="percent")
allrepsoutsummary_long$typef<-dplyr::if_else(allrepsoutsummary_long$type=="percentincludepop","CI includes Pop. r","p's < .05")
allrepsoutsummary_long$typef<-factor(allrepsoutsummary_long$typef,levels=c("p's < .05","CI includes Pop. r"))

allrepsoutsummary_long$percentpercent<-allrepsoutsummary_long$percent*100

ggbars<-ggplot()+geom_bar(data=allrepsoutsummary_long,aes(x=typef,y=percent*100),stat='identity')+
  facet_wrap(~typef,scales = "free_x")+scale_y_continuous(limits = c(0,100),expand = c(0,0))+
  transition_states(groupnumb, transition_length = 0, state_length = 100, wrap = TRUE)

ggbars<-plotclean(ggbars)+theme(legend.position = "none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                                strip.background = element_blank())+
  ylab("% of Studies")+xlab("")

aggbars<-animate(ggbars, nframes = 200)

a_mgif <- image_read(aggCIfig)
b_mgif <- image_read(aggbars)

new_gif <- image_append(c(a_mgif[1], b_mgif[1]))
for(i in 2:10){
  combined <- image_append(c(a_mgif[i], b_mgif[i]))
  new_gif <- c(new_gif, combined)
}

animate(new_gif, nframes = 200)

image_write(new_gif,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Presentations/SpringJobTalks/PittPsychology/Figures/Newstatscorrelation.sim.withbars.gif")






