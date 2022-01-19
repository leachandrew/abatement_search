####### LEACH-AD Numeric simulations Nov 2021 #######
#Mac
if(R.version$platform ==  "x86_64-apple-darwin15.6.0")
  setwd("/Users/aleach/Dropbox/Papers/A-D")
#PC
if(R.version$platform ==  "x86_64-w64-mingw32")
  setwd("C:/Users/aleach/Dropbox/Papers/A-D")
print(getwd())

print(getwd())
source("../../andrew_base.R")

#specific to this project

library(cowplot)


#model parameters

c_s<-100
c_a<-100
e_m<-50

#mean of log normal distribution for search success
mean_log<-log(180)


#pi is a constant in r

#lognormal densities dlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)
#plnorm(q, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE)
#qlnorm(p, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE)
#rlnorm(n, meanlog = 0, sdlog = 1)

#MAC

mac<-function(e,theta=0.5){
  theta-theta*e/e_m
}

#integral from e_m to e of mac
tac<-function(e,theta=50){
  theta*e_m-theta/2*e_m-(theta*e-theta/2*e^2/e_m)
}





fig_tac<-function(file_sent="tac_plot.png",width_sent=8,height_sent=5,file_exp=1,theta_high=200,theta_low=50){
  tac_plot<-
    ggplot()+
    geom_function(aes(color="Higher \u03B8",lty="Higher \u03B8"),fun=tac,args=list(theta=theta_high))+
    geom_function(aes(color="Lower \u03B8",lty="Lower \u03B8"),fun=tac,args=list(theta=theta_low))+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0),limits =c(0,e_m))+
    expand_limits(y=theta_high*1.1)+
    scale_color_manual("",values=c("black","grey50","grey70"))+
    scale_linetype_manual("",values=c("solid","22"))+
    labs(x="Emissions",y="Total Abatement Costs ($/emissions unit)")+
    theme_classic()+
    guides(color=guide_legend(nrow=1,byrow=TRUE,keywidth = 6))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      
      legend.position = "bottom",
    )
  
  if(file_exp==1)
    ggsave(tac_plot,file=file_sent, width = width_sent, height = height_sent,dpi=300)
}

#fig_tac()


fig_mac<-function(file_sent="mac_plot.png",width_sent=8,height_sent=5,file_exp=1,theta_high=200,theta_low=50){
  mac_plot<-
    ggplot()+
    geom_function(aes(color="Higher \u03B8",lty="Higher \u03B8"),fun=mac,args=list(theta=theta_high),size=1.25)+
    geom_function(aes(color="Lower \u03B8",lty="Lower \u03B8"),fun=mac,args=list(theta=theta_low),size=1.25)+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0),limits =c(0,e_m))+
    expand_limits(y=theta_high*1.1)+
    scale_color_manual("",values=c("black","grey30","grey50"))+
    scale_linetype_manual("",values=c("solid","22"))+
    labs(x="Emissions",y="Marginal Abatement Costs ($/emissions unit)")+
    theme_classic()+
    guides(color=guide_legend(nrow=1,byrow=TRUE,keywidth = 6))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      
      legend.position = "bottom",
    )
  
  if(file_exp==1)
    ggsave(file=file_sent, width = width_sent, height = height_sent,dpi=300)
}

fig_mac()

library(nleqslv)

equil_emissions<-function(tau_sent=30,theta_sent=50){
  equi_cond<-function(tau_sent,e,theta_sent){mac(e,theta_sent)-tau_sent}
  x_test<-0.5
  test_soln<-nleqslv(x_test,equi_cond,jacobian=NULL,control=list(maxit=500),tau_sent=tau_sent,theta_sent=theta_sent)
  max(test_soln$x,0)
}


e_solve<-equil_emissions(tau_sent=100,theta_sent=150)


tcp<-function(e_sent,theta_sent,tau_sent,alloc=0){ #total cost of pollution net of allocated permits
  tac(e=e_sent,theta = theta_sent)+tau_sent*e_sent-tau_sent*alloc
}

tcp(e_sent = e_solve,theta_sent = 150,tau=100)

fig_tcp<-function(file_sent="tcp_plot.png",width_sent=8,height_sent=5,file_exp=1,theta_high=200,theta_low=50,tau_sent=30){
  #solve for optimal e at each theta_sent
  equil_high<-equil_emissions(tau_sent,theta_high)
  equil_low<-equil_emissions(tau_sent,theta_low)
  #print(paste(equil_high,equil_low))
  tcp_plot<-
    ggplot()+
    geom_function(aes(color="Higher \u03B8",lty="Higher \u03B8"),fun=tcp,args=list(theta_sent=theta_high,tau=tau_sent),size=rel(1.25))+
    geom_function(aes(color="Lower \u03B8",lty="Lower \u03B8"),fun=tcp,args=list(theta_sent=theta_low,tau=tau_sent),size=rel(1.25))+
    geom_point(data = tibble(1),aes(x=equil_high,y=tcp(equil_high,theta_sent = theta_high,tau=tau_sent)),size=rel(2.25))+
    geom_text_repel(data = tibble(1),aes(x=equil_high,y=tcp(equil_high,theta_sent = theta_high,tau=tau_sent)),label=str_wrap("Minimum Total Costs of Pollution Control Higher \u03B8",45),nudge_y = 1000,min.segment.length =0,size=rel(2.5) )+
    geom_point(data = tibble(1),aes(x=equil_low,y=tcp(equil_low,theta_sent = theta_low,tau=tau_sent)),size=rel(2.25),color="black")+
    geom_text_repel(data = tibble(1),aes(x=equil_low,y=tcp(equil_low,theta_sent = theta_low,tau=tau_sent)),label=str_wrap("Minimum Total Costs of Pollution Control Lower \u03B8",45),nudge_y = -500,min.segment.length = 0,size=rel(2.5),color="black")+
    
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0),limits =c(0,e_m))+
    expand_limits(y=0)+
    scale_color_manual("",values=c("black","grey30","grey50"))+
    scale_linetype_manual("",values=c("solid","22"))+
    labs(x="Emissions",y="Total Costs of Pollution Control ($)")+
    theme_classic()+
    guides(color=guide_legend(nrow=1,byrow=TRUE,keywidth = 6))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      
      legend.position = "bottom",
    )
  
  if(file_exp==1)
    ggsave(tcp_plot,file=file_sent, width = width_sent, height = height_sent,dpi=300)
}

fig_tcp()

#reservation technology


res_tech<-function(theta,adopt_cost=c_a,tau_sent){
  tau_sent^2*theta*e_m/(tau_sent^2*e_m+2*theta*adopt_cost)
}


res_tech_plot<-function(file_sent="res_tech_plot.png",width_sent=8,height_sent=5,file_exp=1,theta_high=200,theta_low=100,tau_sent=50){
  #solve for optimal e at each theta_sent
  res_tech_plot<-
    ggplot()+
    geom_function(aes(color="Higher \u03B8, low tax",lty="Higher \u03B8, low tax"),fun=res_tech,args=list(theta=theta_high,tau_sent=tau_sent),size=rel(1.25))+
    geom_function(aes(color="Higher \u03B8, high tax",lty="Higher \u03B8, high tax"),fun=res_tech,args=list(theta=theta_high,tau_sent=tau_sent*1.5),size=rel(1.25))+
    geom_function(aes(color="Lower \u03B8, high tax",lty="Lower \u03B8, high tax"),fun=res_tech,args=list(theta=theta_low,tau_sent=tau_sent*1.5),size=rel(1.25))+
    geom_function(aes(color="Lower \u03B8, low tax",lty="Lower \u03B8, low tax"),fun=res_tech,args=list(theta=theta_low,tau_sent=tau_sent),size=rel(1.25))+
    geom_point(aes(x=0,y = theta_high),size=rel(2.25))+
    geom_text_repel(data = tibble(1),aes(x=0,y=theta_high),label=str_wrap("Original value, higher \u03B8",45),nudge_x = 2000,min.segment.length =0,size=rel(2.5) )+
    geom_text_repel(data = tibble(1),aes(x=0,y=theta_low),label=str_wrap("Original value, lower \u03B8",45),nudge_x = 2000,min.segment.length =0,size=rel(2.5),color="black" )+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    expand_limits(x=c(0,10000))+
    expand_limits(y=c(0,theta_high*1.1))+
    scale_color_manual("",values=c("black","black","grey50","grey50"))+
    scale_linetype_manual("",values=c("solid","11","solid","11"))+
    labs(x="Technology Adoption Cost",y="Reservation Technology Quality (\u03B8)")+
    theme_classic()+
    guides(color=guide_legend(nrow=2,byrow=TRUE,keywidth = 6))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      
      legend.position = "bottom",
    )
  
  if(file_exp==1)
    ggsave(res_tech_plot,file=file_sent, width = width_sent, height = height_sent,dpi=300)
}

res_tech_plot()


#search decision 

#need a distribution for the search function first

#probablity of drawing theta less than res_tech(theta,c_a,tau_sent)

#log-normal distribution


d_search_prob<-function(theta_sent){
  dlnorm(theta_sent,meanlog = mean_log)
}

search_prob<-function(theta_sent){
  plnorm(theta_sent,meanlog = mean_log)
}



search_gains<-function(theta_sent,tau_sent,adopt_cost=c_a,search_cost=c_s,inc=10^(-1)){ #uses global search and adoption costs by default, but you can send your own
  #inc is the step in the numerical expected value of search
  #theta_sent<-50 # for testing
  #tau_sent<-80  # for testing
#reservation technology level
#will adopt a tech with theta less than
res_theta<-res_tech(theta=theta_sent,adopt_cost = adopt_cost,tau_sent=tau_sent)
#print(paste("will adopt any tech less than",res_theta))
#likelihood of that draw is
#print(paste("likelihood of that draw is",search_prob(res_tech(theta=theta_sent,c_a=adopt_cost,tau_sent=tau_sent))))
#define existing cost given optimal choices 
old_cost<-tcp(equil_emissions(tau_sent = tau_sent, theta_sent = theta_sent),theta_sent = theta_sent,tau_sent = tau_sent)
#print(paste("Old cost is", old_cost))
 #numerical expected value of savings of draws less than res tech
 #cost_at_res<-tcp(equil_emissions(tau_sent = tau_sent, theta_sent = res_theta),theta_sent = res_theta,tau_sent = tau_sent)
 #cost_at_res
 num_test<-tibble(theta_draw=seq(10^(-6),res_theta,inc))%>%mutate(prob=d_search_prob(theta_draw))%>%
   group_by(theta_draw) %>% mutate(new_cost=tcp(equil_emissions(tau_sent = tau_sent, theta_sent = theta_draw),theta_sent = theta_draw,tau_sent = tau_sent)+adopt_cost,
                                   savings=old_cost-new_cost)%>%ungroup()
 num_test_out<<-num_test
 num_test<-num_test%>%
   summarize(e_savings=sum(prob*savings)/sum(prob))
       
  #print(paste("expected savings conditional on draw below res tech (",res_theta,") is ",as.numeric(num_test)))
  
  search_gains<-search_prob(res_theta)*as.numeric(num_test)-search_cost
  #print(paste("expected savings conditional on searching is",search_gains))
  search_gains
  
}

search_gains(theta_sent = 200,tau_sent = 50,adopt_cost=100,search_cost = 100)



e_gains_plot<-function(file_sent="e_gains_plot.png",width_sent=8,height_sent=5,file_exp=1,tau_high=50,tau_low=25){
  #solve for optimal e at each theta_sent

  test<-tibble(theta_test=seq(10,500,10))%>%group_by(theta_test)%>%mutate(gains_low=search_gains(theta_test,tau_sent=tau_low,search_cost = 0,adopt_cost = 0),
                                                                          gains_high=search_gains(theta_test,tau_sent=tau_high,search_cost = 0,adopt_cost = 0),
                                                                          gains_low_high_a=search_gains(theta_test,tau_sent=tau_low,search_cost = 0,adopt_cost = 500),
                                                                          gains_high_high_a=search_gains(theta_test,tau_sent=tau_high,search_cost = 0,adopt_cost =500)
                                                                  
                                                                          
                                                                          )
  
  e_gains_out<<-test
  e_gains_plot<-
    ggplot(test)+
      geom_line(aes(color="Low Tax, Low Adoption Cost",lty="Low Tax, Low Adoption Cost",x=theta_test,y=gains_low),size=rel(1.25))+
      geom_line(aes(color="Low Tax, High Adoption Cost",lty="Low Tax, High Adoption Cost",x=theta_test,y=gains_low_high_a),size=rel(1.25))+
      geom_line(aes(color="High Tax, Low Adoption Cost",lty="High Tax, Low Adoption Cost",x=theta_test,y=gains_high),size=rel(1.25))+
      geom_line(aes(color="High Tax, High Adoption Cost",lty="High Tax, High Adoption Cost",x=theta_test,y=gains_high_high_a),size=rel(1.25))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    expand_limits(y=c(0,290))+
    scale_color_manual("",values=c("black","black","grey50","grey50"))+
    scale_linetype_manual("",values=c("solid","32","solid","3211"))+
    labs(x="Initial value of \u03B8",y="Critical value of search costs ($)")+
    theme_classic()+
    guides(color=guide_legend(nrow=2,byrow=TRUE,keywidth = 4))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
    )
  #e_gains_plot
  
  if(file_exp==1)
    ggsave(e_gains_plot,file=file_sent, width = width_sent, height = height_sent,dpi=300)
}

e_gains_plot()



#no need to repeat for cap-and-trade since allocations are not affected by changes in emissions or tech

#equil_emissions are the same for permit price equals tau

#set allocations in tcp to account for changes in allocation of credits

#as long as allocations aren't a function of emissions, then the reservation tech doesn't change


#repeat for standards

#equil_emissions are equal to the standard

#set allocations in tcp to account for changes in allocation of credits
e_std<-7
theta_sent<-240
tcp(e_sent = e_std,theta_sent = theta_sent,tau=0,alloc = 0)
#reservation technology under standards changes because you only get TAX savings

res_tech_std<-function(theta,adopt_cost=c_a,e_s){
#calculate shadow value
    #theta<-240
    #e_s<-7
    sv<-mac(e_s,theta=theta)
    max(10^-6,theta-2*adopt_cost*theta^2/sv^2/e_m)
}

theta_res<-res_tech_std(theta_sent,adopt_cost=c_a,e_s=e_std)

tcp(e_sent = 5,theta_sent = 90,tau=0,alloc = 0)


search_gains_std<-function(theta_sent,std_sent,adopt_cost=c_a,search_cost=c_s,inc=10^(-1)){ #uses global search and adoption costs by default, but you can send your own
  #inc is the step in the numerical expected value of search
  #theta_sent<-240 # for testing
  #std_sent<-5  # for testing
  #reservation technology level
  #will adopt a tech with theta less than
  res_theta<-max(10^(-6),res_tech_std(theta_sent,adopt_cost=adopt_cost,e_s=std_sent))
  #fixed it to 10^(-5) for cases where there is a negative
  #print(paste("will adopt any tech less than",res_theta))
  #likelihood of that draw is
  #print(paste("likelihood of that draw is",search_prob(res_theta)))
  #define existing cost given optimal choices 
  old_cost<-tcp(std_sent,theta_sent = theta_sent,tau_sent = 0)
  #print(paste("Old cost is", old_cost))
  #numerical expected value of savings of draws less than res tech
  #cost_at_res<-tcp(equil_emissions(tau_sent = tau_sent, theta_sent = res_theta),theta_sent = res_theta,tau_sent = tau_sent)
  #cost_at_res
  num_test<<-tibble(theta_draw=seq(10^(-6),res_theta,inc))%>%mutate(prob=d_search_prob(theta_draw))%>%
    group_by(theta_draw) %>% mutate(new_cost=tcp(std_sent,theta_sent = theta_draw,tau_sent = 0)+adopt_cost,
                                    savings=old_cost-new_cost)%>%ungroup()%>%
    summarize(e_savings=sum(prob*savings)/sum(prob))
  
  #print(paste("expected savings conditional on draw below res tech (",res_theta,") is ",as.numeric(num_test)))
  
  search_gains<-search_prob(res_theta)*as.numeric(num_test)-search_cost
  #print(paste("expected savings conditional on searching is",search_gains))
  search_gains
  
}

#search_gains_std(10,6,200)


e_gains_plot_std<-function(file_sent="e_gains_plot_std.png",width_sent=8,height_sent=5,file_exp=1,std_high=6,std_low=2){
  #solve for optimal e at each theta_sent
  
  test<-tibble(theta_test=seq(10,500,10))%>%group_by(theta_test)%>%mutate(gains_low=search_gains_std(theta_test,std_sent = std_low,search_cost = 0,adopt_cost = 0),
                                                                          gains_high=search_gains_std(theta_test,std_sent = std_high,search_cost = 0,adopt_cost = 0),
                                                                          gains_low_high_a=search_gains_std(theta_test,std_sent = std_low,search_cost = 0,adopt_cost = 500),
                                                                          gains_high_high_a=search_gains_std(theta_test,std_sent = std_high,search_cost = 0,adopt_cost =500)
                                                                          
                                                                          
  )
  e_gains_plot<-
    ggplot(test)+
    geom_line(aes(color="Stringent Standard, Low Adoption Cost",lty="Stringent Standard, Low Adoption Cost",x=theta_test,y=gains_low),size=rel(1.25))+
    geom_line(aes(color="Stringent Standard, High Adoption Cost",lty="Stringent Standard, High Adoption Cost",x=theta_test,y=gains_low_high_a),size=rel(1.25))+
    geom_line(aes(color="Weak Standard, Low Adoption Cost",lty="Weak Standard, Low Adoption Cost",x=theta_test,y=gains_high),size=rel(1.25))+
    geom_line(aes(color="Weak Standard, High Adoption Cost",lty="Weak Standard, High Adoption Cost",x=theta_test,y=gains_high_high_a),size=rel(1.25))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    expand_limits(y=c(0,290))+
    scale_color_manual("",values=c("black","black","grey50","grey50"))+
    scale_linetype_manual("",values=c("solid","32","solid","3211"))+
    labs(x="Initial value of \u03B8",y="Critical value of search costs ($)")+
    theme_classic()+
    guides(color=guide_legend(nrow=2,byrow=TRUE,keywidth = 4))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
    )
  e_gains_plot
  
  if(file_exp==1)
    ggsave(e_gains_plot,file=file_sent, width = width_sent, height = height_sent,dpi=300)
}

e_gains_plot_std()


e_gains_plot_mix<-function(file_sent="e_gains_plot_mix.png",width_sent=8,height_sent=5,file_exp=1,e_fix=25){
  #gains from search, fixed emissions by tax or standard
    test<-tibble(theta_test=seq(10^(-2),500,10))%>%group_by(theta_test)%>%mutate(e_std=e_fix,
                                                                            tau_equiv=mac(e_fix,theta_test),
                                                                          tcpc_tax=tcp(e_sent = e_fix,theta_sent = theta_test,tau=tau_equiv),
                                                                          tcpc_std=tcp(e_sent = e_fix,theta_sent = theta_test,tau=0),
                                                                          gains_std=search_gains_std(theta_test,std_sent = e_fix,search_cost = 0,adopt_cost = 0),
                                                                          gains_tax=search_gains(theta_test,tau_sent = tau_equiv,search_cost = 0,adopt_cost = 0),
                                                                          #gains_low_high_a=search_gains_std(theta_test,std_sent = std_low,search_cost = 0,adopt_cost = 500),
                                                                          #gains_high_high_a=search_gains_std(theta_test,std_sent = std_high,search_cost = 0,adopt_cost =500
                                                                          )

  e_gains_plot<-
    ggplot(test)+
    geom_line(aes(color="Standard",lty="Standard",x=theta_test,y=gains_std),size=rel(1.25))+
    geom_line(aes(color="Tax",lty="Tax",x=theta_test,y=gains_tax),size=rel(1.25))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    expand_limits(y=c(0,290))+
    scale_color_manual("",values=c("black","black","grey50","grey50"))+
    scale_linetype_manual("",values=c("solid","32","solid","3211"))+
    labs(x="Initial value of \u03B8",y="Critical value of search costs ($)")+
    theme_classic()+
    guides(color=guide_legend(nrow=1,byrow=TRUE,keywidth = 4))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
    )
  e_gains_plot
  
  if(file_exp==1)
    ggsave(e_gains_plot,file=file_sent, width = width_sent, height = height_sent,dpi=300)
}

e_gains_plot_mix()



#conditions for porter hypothesis


fig_porter_tax<-function(file_sent="porter_tax.png",width_sent=10,height_sent=6,file_exp=1,theta_sent,tau_high=150,tau_low=50,adopt_cost=c_a,search_cost=c_s,inc=0.005,porter_max=200){
  #print(paste(equil_high,equil_low))
  #adopt_cost
  #set theta_grid
 #inc<-1
  
  theta_grid<-tibble(theta_draw=seq(1,porter_max,inc))%>% mutate(res_tech_high=res_tech(theta_draw,adopt_cost,tau_high),
                                                        res_tech_low=res_tech(theta_draw,adopt_cost,tau_low))%>%
    group_by(theta_draw)%>% mutate(gains_high=search_gains(theta_sent = theta_draw,tau_sent = tau_high,adopt_cost=adopt_cost,search_cost = search_cost),
                                   gains_low=search_gains(theta_sent = theta_draw,tau_sent = tau_low,adopt_cost=adopt_cost,search_cost = search_cost),
                                   tcp_high=tcp(equil_emissions(tau_sent = tau_high, theta_sent = theta_draw),theta_sent = theta_draw,tau_sent = tau_high),
                                   tcp_low=tcp(equil_emissions(tau_sent = tau_low, theta_sent = theta_draw),theta_sent = theta_draw,tau_sent = tau_low),
                                   )
  
  theta_grid<-theta_grid %>% mutate(search_cond=((gains_low<0)*(gains_high>0)==1),
                                    low_search=((gains_low>=0)),
                                    e_low=tcp_low-(gains_low>=0)*(gains_low+search_cost),
                                    e_high=tcp_high-(gains_high>=0)*(gains_high+search_cost),
                                    one_search=((gains_low>=0)+(gains_high>=0)>=1))%>% group_by(theta_draw)%>%
    mutate(test2=first(which(abs(theta_grid$tcp_high-tcp_low)==min(abs(theta_grid$tcp_high-tcp_low)))))%>%
    ungroup()%>% mutate(theta_min=theta_draw[test2])
  
  theta_grid<-theta_grid %>% mutate(porter_cond=(theta_min<res_tech_high),
                                    porter_likely=search_prob(theta_min))
  #porter_likely is the probability of getting a draw of theta less than theta at which porter would happen
  theta_grid_out<<-theta_grid
  
  min_high_search<-theta_grid %>% filter(gains_high>0) %>% summarize(value=min(theta_draw)) %>% as.numeric()
  min_low_search<-theta_grid %>% filter(gains_low>0) %>% summarize(value=min(theta_draw)) %>% as.numeric()
  
  
  tcp_plot<-
    ggplot(theta_grid)+
    #geom_ribbon(data=theta_grid%>%filter(search_cond,porter_cond),aes(ymax=theta_min,ymin=0,x=theta_draw),alpha=0.5, colour="blue",fill="blue")+
    geom_line(aes(x=theta_draw,y=theta_draw))+
    geom_text_repel(data=tibble(1),aes(y=porter_max-8,x=porter_max-8),nudge_x = -10,label="45\u00B0 line",size=rel(5))+
    geom_ribbon(data=theta_grid%>%filter(porter_cond),aes(ymax=theta_min,ymin=0,x=theta_draw,fill="C"),alpha=1)+
    
    geom_ribbon(data=theta_grid%>%filter(low_search),aes(x=theta_draw,ymax=res_tech_high,ymin=0,fill="B"),alpha=1)+
    geom_ribbon(data=theta_grid%>%filter(one_search),aes(x=theta_draw,ymax=res_tech_high,ymin=0,fill="A2"),alpha=1)+
    geom_ribbon(aes(x=theta_draw,ymax=res_tech_high,ymin=0,fill="A1"),alpha=1,color="black")+
    geom_ribbon(data=theta_grid%>%filter(one_search),aes(x=theta_draw,ymax=res_tech_high,ymin=0,fill="A2"),alpha=1)+
    geom_ribbon(data=theta_grid%>%filter(one_search),aes(ymax=theta_min,ymin=0,x=theta_draw,fill="C"),alpha=1)+
    
    geom_ribbon(data=theta_grid%>%filter(low_search),aes(x=theta_draw,ymax=res_tech_high,ymin=0,fill="B"),alpha=1)+
    geom_line(aes(x=theta_draw,y=res_tech_high))+
    geom_segment(aes(x=min_high_search,xend=min_high_search,y=0,yend=min_high_search),lty="11")+
    geom_segment(aes(x=min_low_search,xend=min_low_search,y=0,yend=min_low_search),lty="11")+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    expand_limits(y=c(0,porter_max),x=c(0,porter_max))+
    scale_fill_manual("",values=c("grey95","grey80","grey50","black"),labels=c(
      "Adoptible ex-post values of \u03B8, but search not optimal with either tax",                                  
      "Search and adoption optimal under high tax, but ex-post values of adopted \u03B8 not consistent with Porter hypothesis",
                                                                 
                                                                   "Search condition not satisfied (search is optimal under both high and low taxes)",
                                                                   "Search optimal only under high tax, ex-post values of adopted \u03B8 consistent with Porter hypothesis"))+
    scale_color_manual("",values=c("black","grey30","grey50","grey70"))+
    #scale_linetype_manual("",values=c("solid"))+
    labs(x="Initial value of \u03B8",y="Ex-post value of \u03B8")+
    theme_classic()+
    guides(fill=guide_legend(nrow=4,byrow=FALSE,keywidth = 2,keyheight = 2))+
    theme(
      #legend.position=c(.25,.9),
      text=element_text(size=20),legend.text = element_text(size=18),
      #text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      
      legend.position = "bottom",
    )
  save(file="porter_plot.rdata",tcp_plot)
  if(file_exp==1)
    ggsave(tcp_plot,file=file_sent, width = width_sent, height = height_sent,dpi=300)
  tcp_plot
}

fig_porter_tax(file_exp=1,inc=0.1,porter_max=300,width_sent = 16,height_sent = 17,adopt_cost = c_a,search_cost = 600,tau_high=100,tau_low=80)



e_porter_plot<-function(theta_data_sent){
  theta_data_sent<-theta_grid_out
  min_high_search<-theta_data_sent %>% filter(gains_high>0) %>% summarize(value=min(theta_draw)) %>% as.numeric()
  min_low_search<-theta_data_sent %>% filter(gains_low>0) %>% summarize(value=min(theta_draw)) %>% as.numeric()
  mid<-(min_high_search+min_low_search)/2
  e_plot<-ggplot(theta_data_sent)+
    geom_line(aes(x=theta_draw,y=e_low,linetype="Low tax"),size=1.5)+
    geom_line(aes(x=theta_draw,y=e_high,linetype="High tax"),size=1.5)+
    geom_vline(aes(xintercept=min_high_search),lty="11")+
    geom_vline(aes(xintercept=min_low_search),lty="11")+
    annotate("text",x=mid,y=1500,
             label=str_wrap("In this range, a firm would only search when facing the high tax and, for these parameter values, this implies lower expected values of total compliance costs with the high tax despite the increased policy stringency",20))+
    annotate("text",x=120,y=1500,
    label=str_wrap("In this range, the firm would not search under either tax so total compliance costs are a deterministic function of the initial value of \u03B8",20))+
    annotate("text",x=300,y=1500,
             label=str_wrap("In this range, the firm would search under either tax, so expected total compliance costs are lower with less stringent policies.",20))+
   scale_color_manual("",values=c("black","grey30","grey50","grey70"))+
    expand_limits(x=325)+
  scale_linetype_manual("",values=c("solid","32"))+
  labs(x="Initial value of \u03B8",y="Expected total cost of pollution control")+
    theme_classic()+
    guides(color=guide_legend(nrow=1,byrow=TRUE,keywidth = 5),lty=guide_legend(nrow=1,byrow=TRUE,keywidth = 5))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
    )
  ggsave(e_plot,file="expected_cost.png", width = 12, height = 6,dpi=300)
  e_plot
  
}
e_porter_plot(theta_grid_out)







fig_porter_std<-function(file_sent="porter_std.png",width_sent=10,height_sent=6,file_exp=1,theta_sent,sv_weak=90,sv_strong=100,adopt_cost=c_a,
                         search_cost=c_s,inc=0.005,porter_max=200){
  #print(paste(equil_weak,equil_strong))
  #adopt_cost
  #inc<-.5
  #set theta_grid
  #solve for emissions limits
  e_weak<-equil_emissions(tau_sent=sv_weak,theta_sent=theta_sent)
  e_strong<-equil_emissions(tau_sent=sv_strong,theta_sent=theta_sent)
  
  
  
  theta_grid<-tibble(theta_draw=seq(10^-3,porter_max,inc))%>% group_by(theta_draw)%>%
    mutate(emit_imp_weak=e_weak,emit_imp_strong=e_strong)%>%
    mutate(res_tech_weak=res_tech_std(theta = theta_draw,adopt_cost = adopt_cost,e_s=emit_imp_weak),
           res_tech_strong=res_tech_std(theta_draw,adopt_cost,e_s=emit_imp_strong))%>%
     mutate(gains_weak=search_gains_std(theta_draw,emit_imp_weak,adopt_cost,search_cost),
                                   gains_strong=search_gains_std(theta_draw,emit_imp_strong,adopt_cost,search_cost),
                                   tcp_weak=tcp(emit_imp_weak,theta_sent = theta_draw,tau_sent = 0),
                                   tcp_strong=tcp(emit_imp_strong,theta_sent = theta_draw,tau_sent = 0),
                                    e_weak=tcp_weak-(gains_weak>=0)*(gains_weak+search_cost),
                                    e_strong=tcp_strong-(gains_strong>=0)*(gains_strong+search_cost),
            
    )
  
  theta_grid<-theta_grid %>% mutate(search_cond=((gains_weak<0)*(gains_strong>0)==1),
                                    one_search=((gains_weak>=0)+(gains_strong>=0)>=1),
                                    no_search=((gains_weak<=0)+(gains_strong<=0)>1),
                                    both_search=((gains_weak>=0)+(gains_strong>=0)>1))
  
  theta_grid<-theta_grid %>% group_by(theta_draw)%>%
    mutate(test2=first(which(abs(theta_grid$tcp_strong-tcp_weak)==min(abs(theta_grid$tcp_strong-tcp_weak)))))%>%
    ungroup()%>% mutate(theta_min=theta_draw[test2])
  
  theta_grid<-theta_grid %>% mutate(porter_cond=(theta_min<res_tech_strong),#I would adopt a tech in POrter range
                                    porter_likely=search_prob(theta_min)) #likelihood of getting that draw
  

  
  min_high_search<-theta_grid %>% filter(gains_strong>0) %>% summarize(value=min(theta_draw)) %>% as.numeric()
  min_low_search<-theta_grid %>% filter(gains_weak>0) %>% summarize(value=min(theta_draw)) %>% as.numeric()
  
  theta_std<<-theta_grid
    col_vec<-c("grey99","grey80","grey50","black")
  #col_vec<-colors_tableau10()
  tcp_plot<-
    ggplot(theta_grid)+
    #geom_ribbon(data=theta_grid%>%filter(search_cond,porter_cond),aes(ymax=theta_min,ymin=0,x=theta_draw),alpha=0.5, colour="blue",fill="blue")+
    geom_line(aes(x=theta_draw,y=theta_draw))+
    geom_ribbon(aes(x=theta_draw,ymax=res_tech_strong,ymin=0,fill="A0"),alpha=1)+
    geom_ribbon(data=theta_grid%>%filter(!no_search),aes(ymax=res_tech_strong,ymin=0,x=theta_draw,fill="A1"),alpha=1,color="black")+
    geom_text_repel(data=tibble(1),aes(y=porter_max-12,x=porter_max-12),nudge_x = -10,label="45\u00B0 line",size=rel(5))+
    
    #geom_ribbon(aes(x=theta_draw,ymax=theta_min,ymin=0,fill="A2"),alpha=1, colour="black")+
    
    geom_ribbon(data=theta_grid%>%filter(one_search),aes(ymax=theta_min,ymin=0,x=theta_draw,fill="C"),alpha=1)+
    geom_ribbon(data=theta_grid%>%filter(both_search),aes(x=theta_draw,ymax=res_tech_strong,ymin=0,fill="A2"),alpha=1)+
    
    geom_segment(aes(x=min_high_search,xend=min_high_search,y=0,yend=min_high_search),lty="11")+
    geom_segment(aes(x=min_low_search,xend=min_low_search,y=0,yend=min_low_search),lty="11")+
    geom_line(aes(x=theta_draw,y=res_tech_strong),lty="11")+
    
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    expand_limits(y=c(0,porter_max),x=c(0,porter_max))+
    scale_fill_manual("",values=col_vec,labels=c("Adoptible ex-post values of \u03B8, but search not optimal with either standard",                                
                                                 "Search and adoption optimal under stringent standard, but ex-post values of adopted \u03B8 not consistent with Porter hypothesis",
                                                "Search condition not satisfied (search optimal with weak and strong standards)",
                                                "Search optimal only under stringent standard, ex-post values of adopted \u03B8 consistent with Porter hypothesis"))+
    scale_color_manual("",values=c("black","grey50","grey70"))+
    scale_linetype_manual("",values=c("11"))+
    labs(x="Initial value of \u03B8",y="Ex-post value of \u03B8")+
    theme_classic()+
    guides(fill=guide_legend(nrow=4,byrow=FALSE,keywidth = 1),lty="none")+
    theme(
      #legend.position=c(.25,.9),
      text=element_text(size=20),legend.text = element_text(size=18),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      
      legend.position = "bottom",
    )
  tcp_plot
  save(file="porter_plot_std.rdata",tcp_plot)
  
  if(file_exp==1)
    ggsave(tcp_plot,file=file_sent, width = width_sent, height = height_sent,dpi=300)
}

fig_porter_std(file_exp =1,inc=0.10,porter_max=300,width_sent = 16,height_sent = 17,theta_sent=180,adopt_cost =c_a,search_cost = 500)





e_porter_plot_std<-function(theta_data_sent){
  #testing theta_data_sent<-theta_std
  min_high_search<-theta_data_sent %>% filter(gains_strong>0) %>% summarize(value=min(theta_draw)) %>% as.numeric()
  min_low_search<-theta_data_sent %>% filter(gains_weak>0) %>% summarize(value=min(theta_draw)) %>% as.numeric()
  mid<-(min_high_search+min_low_search)/2
  e_plot<-ggplot(theta_data_sent)+
    geom_line(aes(x=theta_draw,y=e_weak,linetype="Less stringent standard"),size=1.5)+
    geom_line(aes(x=theta_draw,y=e_strong,linetype="More stringent standard"),size=1.5)+
    geom_vline(aes(xintercept=min_high_search),lty="11")+
    geom_vline(aes(xintercept=min_low_search),lty="11")+
    annotate("text",x=mid,y=600,
             label=str_wrap("Firm would only search when facing the stringent standard and ex post compliance costs are lower with increased policy stringency",12))+
    annotate("text",x=180,y=600,
             label=str_wrap("Firm would not search under either standard. Total compliance costs are a deterministic function of the initial value of \u03B8",40))+
    annotate("text",x=300,y=600,
             label=str_wrap("Firm would search under either standard, and expected total compliance costs are lower with the less stringent standard.",20))+
    scale_color_manual("",values=c("black","grey30","grey50","grey70"))+
    
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    
    expand_limits(x=325)+
    scale_linetype_manual("",values=c("solid","32"))+
    labs(x="Initial value of \u03B8",y="Expected total cost of pollution control")+
    theme_classic()+
    guides(color=guide_legend(nrow=1,byrow=TRUE,keywidth = 5),lty=guide_legend(nrow=1,byrow=TRUE,keywidth = 5))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
    )
  e_plot
  ggsave(e_plot,file="expected_cost_std.png", width = 12, height = 6,dpi=300)
  
  
}
e_porter_plot_std(theta_std)



#repeat for standards

#equil_emissions are equal to the standard

#set allocations in tcp to account for changes in allocation of credits
e_std<-7
theta_sent<-240
tcp(e_sent = e_std,theta_sent = theta_sent,tau=0,alloc = 0)
#reservation technology under standards changes because you only get TAX savings




e_solve<-equil_emissions(tau_sent=sv,theta_sent=theta)
theta_res<-res_tech_std(theta_sent,adopt_cost=c_a,e_s=e_solve)
tcp(e_sent = e_solve,theta_sent = theta_sent,tau=0,alloc = 0)
tcp(e_sent = e_solve,theta_sent = theta_res,tau=0,alloc = 0)





new_search_gains_std<-function(theta_sent,sv_sent,adopt_cost=c_a,search_cost=c_s,inc=10^(-1)){ #uses global search and adoption costs by default, but you can send your own
  #inc is the step in the numerical expected value of search
  #theta_sent<-240 # for testing
  #sv_sent<-50  # send shadow value of standard
  #reservation technology level
  #will adopt a tech with theta less than
  #emissions under the standard
  e_solve<-equil_emissions(tau_sent=sv_sent,theta_sent=theta_sent)
  res_theta<-max(10^(-6),res_tech_std(theta_sent,adopt_cost=adopt_cost,e_s=e_solve))
  #fixed it to 10^(-5) for cases where there is a negative
  #print(paste("will adopt any tech less than",res_theta))
  #likelihood of that draw is
  #print(paste("likelihood of that draw is",search_prob(res_theta)))
  #define existing cost given optimal choices 
  old_cost<-tcp(e_solve,theta_sent = theta_sent,tau_sent = 0)
  #print(paste("Old cost is", old_cost))
  #numerical expected value of savings of draws less than res tech
  #cost_at_res<-tcp(equil_emissions(tau_sent = tau_sent, theta_sent = res_theta),theta_sent = res_theta,tau_sent = tau_sent)
  #cost_at_res
  num_test<-tibble(theta_draw=seq(10^(-6),res_theta,inc))%>%mutate(prob=d_search_prob(theta_draw))%>%
    group_by(theta_draw) %>% mutate(new_cost=tcp(e_solve,theta_sent = theta_draw,tau_sent = 0)+adopt_cost,
                                    savings=old_cost-new_cost)%>%ungroup()%>%
    summarize(e_savings=sum(prob*savings)/sum(prob))%>% as.numeric()
  
  #print(paste("expected savings conditional on draw below res tech (",res_theta,") is ",as.numeric(num_test)))
  
  search_gains<-search_prob(res_theta)*num_test-search_cost
  #print(paste("expected savings conditional on searching is",search_gains))
  search_gains
  
}

#new_search_gains_std(150,150)


new_e_gains_plot_std<-function(file_sent="new_e_gains_plot_std.png",width_sent=8,height_sent=5,file_exp=1,sv_high=150,sv_low=50){
  #solve for optimal e at each theta_sent
  
  test<-tibble(theta_test=seq(10^(-3),sv_high*5,10))%>%group_by(theta_test)%>%mutate(gains_low=new_search_gains_std(theta_test,sv_sent = sv_low,search_cost = 0,adopt_cost = 0),
                                                                          gains_high=new_search_gains_std(theta_test,sv_sent = sv_high,search_cost = 0,adopt_cost = 0),
                                                                          gains_low_high_a =new_search_gains_std(theta_test,sv_sent = sv_low,search_cost = 0,adopt_cost = 500),
                                                                          gains_high_high_a=new_search_gains_std(theta_test,sv_sent = sv_high,search_cost = 0,adopt_cost =500),
                                                                          emit_high=equil_emissions(tau_sent=sv_high,theta_sent=theta_test),
                                                                          emit_low=equil_emissions(tau_sent=sv_low,theta_sent=theta_test),
                                                                          cost_high=tcp(emit_high,theta_sent = theta_test,tau_sent = 0),
                                                                          cost_low=tcp(emit_low,theta_sent = theta_test,tau_sent = 0)                                       
                                                                          )
  
  
  
  
  implied_standard_plot<-
    ggplot(test)+
    geom_vline(xintercept = sv_low)+
    geom_vline(xintercept = sv_high)+
    geom_line(aes(color="Stringent Standard Emissions",lty="Stringent Standard Emissions",x=theta_test,y=emit_high),size=rel(1.25))+
    geom_line(aes(color="Weak Standard Emissions",lty="Weak Standard Emissions",x=theta_test,y=emit_low),size=rel(1.25))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    #expand_limits(y=c(0,290))+
    scale_color_manual("",values=c("black","black","grey50","grey50"))+
    scale_linetype_manual("",values=c("solid","32","solid","3211"))+
    labs(x="Initial value of \u03B8",y="Implied Emissions Limit")+
    theme_classic()+
    guides(color=guide_legend(nrow=2,byrow=TRUE,keywidth = 4))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      #axis.text = element_blank(),
      #axis.ticks = element_blank(),
      legend.position = "bottom",
    )
  implied_standard_plot
  
  init_cost_plot<-
    ggplot(test)+
    geom_vline(xintercept = sv_low)+
    geom_vline(xintercept = sv_high)+
    geom_line(aes(color="Stringent Standard Cost",lty="Stringent Standard Cost",x=theta_test,y=cost_low),size=rel(1.25))+
    geom_line(aes(color="Weak Standard Cost",lty="Weak Standard Cost",x=theta_test,y=cost_high),size=rel(1.25))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    expand_limits(y=c(0,290))+
    scale_color_manual("",values=c("black","black","grey50","grey50"))+
    scale_linetype_manual("",values=c("solid","32","solid","3211"))+
    labs(x="Initial value of \u03B8",y="Total compliance cost before search ($)")+
    theme_classic()+
    guides(color=guide_legend(nrow=2,byrow=TRUE,keywidth = 4))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      #axis.text = element_blank(),
      #axis.ticks = element_blank(),
      legend.position = "bottom",
    )
  init_cost_plot
  
    e_gains_plot<-
    ggplot(test)+
    geom_vline(xintercept = sv_low)+
    geom_vline(xintercept = sv_high)+
    geom_line(aes(color="Stringent Standard",lty="Stringent Standard",x=theta_test,y=gains_high),size=rel(1.25))+
    #geom_line(aes(color="Stringent Standard, High Adoption Cost",lty="Stringent Standard, High Adoption Cost",x=theta_test,y=gains_low_high_a),size=rel(1.25))+
    geom_line(aes(color="Weak Standard",lty="Weak Standard",x=theta_test,y=gains_low),size=rel(1.25))+
    #geom_line(aes(color="Weak Standard, High Adoption Cost",lty="Weak Standard, High Adoption Cost",x=theta_test,y=gains_high_high_a),size=rel(1.25))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    expand_limits(y=c(0,290))+
    scale_color_manual("",values=c("black","black","grey50","grey50"))+
    scale_linetype_manual("",values=c("solid","32","solid","3211"))+
    labs(x="Initial value of \u03B8",y="Critical value of search costs ($)")+
    theme_classic()+
    guides(color=guide_legend(nrow=2,byrow=TRUE,keywidth = 4))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      #axis.text = element_blank(),
      #axis.ticks = element_blank(),
      legend.position = "bottom",
    )
  e_gains_plot
  
  grid.arrange(implied_standard_plot,init_cost_plot,e_gains_plot)
  
  if(file_exp==1)
    ggsave(e_gains_plot,file=file_sent, width = width_sent, height = height_sent,dpi=300)
}


new_e_gains_plot_std()




new_e_gains_plot_tax<-function(file_sent="new_e_gains_plot_std.png",width_sent=8,height_sent=5,file_exp=1,tau_high=150,tau_low=50){
  #solve for optimal e at each theta_sent
  
  test<-tibble(theta_test=seq(10^(-3),sv_high*5,10))%>%group_by(theta_test)%>%mutate(gains_low=search_gains(theta_test,tau_sent = tau_low,search_cost = 0,adopt_cost = 0),
                                                                                     gains_high=search_gains(theta_test,tau_sent = tau_high,search_cost = 0,adopt_cost = 0),
                                                                                     emit_high=equil_emissions(tau_sent=tau_high,theta_sent=theta_test),
                                                                                     emit_low=equil_emissions(tau_sent=tau_low,theta_sent=theta_test),
                                                                                     cost_high=tcp(emit_high,theta_sent = theta_test,tau_sent = tau_high),
                                                                                     cost_low=tcp(emit_low,theta_sent = theta_test,tau_sent = tau_low)                                       
  )
  
  
  
  
  emit_tax_plot<-
    ggplot(test)+
    geom_vline(xintercept = tau_low)+
    geom_vline(xintercept = tau_high)+
    geom_line(aes(color="Stringent Tax Emissions",lty="Stringent Tax Emissions",x=theta_test,y=emit_high),size=rel(1.25))+
    geom_line(aes(color="Weak Tax Emissions",lty="Weak Tax Emissions",x=theta_test,y=emit_low),size=rel(1.25))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    #expand_limits(y=c(0,290))+
    scale_color_manual("",values=c("black","black","grey50","grey50"))+
    scale_linetype_manual("",values=c("solid","32","solid","3211"))+
    labs(x="Initial value of \u03B8",y="Pre-Search Emissions")+
    theme_classic()+
    guides(color=guide_legend(nrow=2,byrow=TRUE,keywidth = 4))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      #axis.text = element_blank(),
      #axis.ticks = element_blank(),
      legend.position = "bottom",
    )
  emit_tax_plot
  
  tax_cost_plot<-
    ggplot(test)+
    geom_vline(xintercept = tau_low)+
    geom_vline(xintercept = tau_high)+
    geom_line(aes(color="Stringent Tax Cost",lty="Stringent Tax Cost",x=theta_test,y=cost_low),size=rel(1.25))+
    geom_line(aes(color="Weak Tax Cost",lty="Weak Tax Cost",x=theta_test,y=cost_high),size=rel(1.25))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    expand_limits(y=c(0,290))+
    scale_color_manual("",values=c("black","black","grey50","grey50"))+
    scale_linetype_manual("",values=c("solid","32","solid","3211"))+
    labs(x="Initial value of \u03B8",y="Total compliance cost before search ($)")+
    theme_classic()+
    guides(color=guide_legend(nrow=2,byrow=TRUE,keywidth = 4))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      #axis.text = element_blank(),
      #axis.ticks = element_blank(),
      legend.position = "bottom",
    )
  tax_cost_plot
  
  tax_gains_plot<-
    ggplot(test)+
    geom_vline(xintercept = tau_low)+
    geom_vline(xintercept = tau_high)+
    geom_line(aes(color="Stringent Tax",lty="Stringent Tax",x=theta_test,y=gains_high),size=rel(1.25))+
    geom_line(aes(color="Weak Tax",lty="Weak Tax",x=theta_test,y=gains_low),size=rel(1.25))+
    scale_x_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    scale_y_continuous(breaks=pretty_breaks(),expand=c(0,0))+
    expand_limits(y=c(0,290))+
    scale_color_manual("",values=c("black","black","grey50","grey50"))+
    scale_linetype_manual("",values=c("solid","32","solid","3211"))+
    labs(x="Initial value of \u03B8",y="Critical value of search costs ($)")+
    theme_classic()+
    guides(color=guide_legend(nrow=2,byrow=TRUE,keywidth = 4))+
    theme(
      #legend.position=c(.25,.9),
      text = element_text(size = 12),
      #axis.text = element_blank(),
      #axis.ticks = element_blank(),
      legend.position = "bottom",
    )
  tax_gains_plot
  
  grid.arrange(emit_tax_plot,tax_cost_plot,tax_gains_plot)
  
  if(file_exp==1)
    ggsave(e_gains_plot,file=file_sent, width = width_sent, height = height_sent,dpi=300)
}




