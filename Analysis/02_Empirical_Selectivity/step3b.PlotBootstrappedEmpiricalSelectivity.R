#--plot SBS empirical selectivity bootstrapping results

dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"Analysis/02_Empirical_Selectivity");
setwd(dirThs);

load(file=file.path(dirThs,"rda_Step3a_EmpiricalSelectivityFromBootstrapping.RData"));

plotEmpSel<-function(dfr,colour,z_lim=c(22.5,187.5),n_min=0,wrap=TRUE){
  p = ggplot(dfr |> dplyr::filter(n_BSFRF+n_NMFS >= n_min),
             aes(x=z,y=emp_sel,colour={{colour}},fill={{colour}},group=paste(y,iB)));
  p = p + geom_line(alpha=0.1,size=0.1) + 
          geom_point(alpha=0.1,size=0.1,position=position_jitter(0.5)) + 
          geom_line(data=dfr |> dplyr::filter(iB==1),alpha=1,size=1);
  if (wrap){
    p = p + geom_smooth(mapping=aes(group=NULL),colour=NA,fill="gray25",
                        method="gam",formula=y~s(x,bs="cs",k=5));
  } else {
    p = p + geom_smooth(mapping=aes(x=z,y=emp_sel),colour=NA,fill="gray25",inherit.aes=FALSE,
                        method="gam",formula=y~s(x,bs="cs",k=5));
  }
  p = p + geom_hline(yintercept=c(0,1),linetype=2) + 
          labs(x="size (mm CW)",y="empirical selectivity",colour="year",fill="year") + 
          scale_y_continuous(limits=c(0,3),oob=scales::squish) + 
          scale_x_continuous(limits=z_lim,oob=scales::squish) + 
          scale_colour_discrete(aesthetics=c("colour","fill")) + 
          guides(colour=guide_legend(override.aes=list(alpha=1))) + 
          wtsPlots::getStdTheme() + 
          theme(panel.grid.major=element_line(colour="gray75"),
                panel.grid.minor=element_line(colour="white"));
  if (wrap) p = p + facet_wrap(~y,nrow=2);
  return(p);
}

#--combined sexes
dfrZCsp<-dfrZCs[dplyr::between(dfrZCs$z,50,190),];
dfrESsp<-dfrESs[dplyr::between(dfrZCs$z,50,190),];
plotEmpSel(dfrESsp,colour=factor(y),z_lim=c(47.5,192.5),n_min=0,wrap=TRUE);
plotEmpSel(dfrESsp,colour=factor(y),z_lim=c(47.5,192.5),n_min=0,wrap=FALSE);
