#--use Thygesen et al approach to estimate relative selectivity ratio
##--using log-Gaussian Cox processes
require(TMB);
require(ggplot2);
dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"Analysis/04a_ThygesenIntercalibration");
dirICl = file.path("~/Work/Projects/SBS_SelectivityAnalysis_TannerCrab/Analysis/04a_ThygesenIntercalibration/Intercalibration");
dirCPP = file.path(dirICl,"gearcalib/src");

#--compile dynamic library
setwd(dirCPP);
compile("gearcalib.cpp",CXX="clang++")
if(is.loaded("EvalDoubleFunObject")) dyn.unload("gearcalib.so") #<-may error out
dyn.load("gearcalib.so")
setwd(dirThs);

#--get data
##--dfrCPUE: fleet, y, gis_station, sampling_factor, area_swept_variable, x, m, s, z, n, val, type
dfrCPUE = wtsUtilities::getObj(file.path(dirPrj,"Analysis/01_SBS_Data/rda_Step3_SBS_CrabAbundance.RData"))$dfrCPUE;
dfrCPUE = dfrCPUE |> dplyr::mutate(gear=factor(fleet,levels=c("BSFRF","NMFS")),
                                   group=paste0(y,"_",gis_station),
                                   ASp=area_swept_variable/sampling_factor);
dfrGr = dfrCPUE |> dplyr::select(group,x,n,ASp) |> 
          dplyr::group_by(group,x) |> 
          dplyr::summarize(totN=sum(n,na.rm=FALSE),
                           totA=sum(ASp,na.rm=FALSE)) |> 
          dplyr::ungroup() |>
          dplyr::arrange(group,x) |> 
          dplyr::filter(totN>0,is.finite(totA)) |> dplyr::select(group,x);
##--NOTE: sexes are combined as "all" (filter kept here for parallel w/ Tanner crab)
dfrA.N = dfrGr |> dplyr::filter(x=="all") |> 
          dplyr::inner_join(dfrCPUE |> dplyr::filter(dplyr::between(z,24,179)),by=c("group","x")) |> 
          dplyr::select(group,gear,z,n) |> 
          tidyr::pivot_wider(names_from=z,values_from=n) |> dplyr::arrange(group,gear);
dfrA.A = dfrGr |> dplyr::filter(x=="all") |> 
          dplyr::inner_join(dfrCPUE |> dplyr::filter(dplyr::between(z,24,179)),by=c("group","x")) |> 
          dplyr::select(group,gear,z,ASp) |> 
          tidyr::pivot_wider(names_from=z,values_from=ASp) |> dplyr::arrange(group,gear);

#--source the R code
source(file.path(dirICl,"gearcalib/R/gearcalib.R"));

#--convert to NEW gearcalib "d" format:
##--d$N: matrix with N in haul x size
##--d$SweptArea: matrix with area_swept/sampling_factor in haul x size
##--d$group: factor encoding hauls (size = nHauls)
##--d$Gear:  factor encoding gear  (size = nHauls)
###--sexes combined: x=="all"
zA = as.double(colnames(dfrA.N)[3:ncol(dfrA.N)]) + 3;
dA = list();
dA$L         = zA;
dA$N         = as.matrix(dfrA.N[,3:ncol(dfrA.N)]);
dA$SweptArea = as.matrix(dfrA.A[,3:ncol(dfrA.A)]);
dA$group     = factor(dfrA.N$group);
dA$Gear      = dfrA.N$gear; 

#--fit the model
fitA = gearcalibFit(dA,fit0=TRUE);
botA = boot(dA);
pltA = plot.gearcalibFit(fitA,boot=botA,Lvec=zA,ymax=NULL)


