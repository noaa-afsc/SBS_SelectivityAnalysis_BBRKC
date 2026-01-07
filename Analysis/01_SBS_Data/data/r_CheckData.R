#--compare SBS haul data characteristics
require(ggplot2)
dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"Analysis/01_SBS_Data/data");

dfr_haul = readr::read_csv(file.path(dirThs,"final_haul.csv"));
dfr_rkc  = readr::read_csv(file.path(dirThs,"final_rkc_sr_melt_acoustic.csv"));

match_hauls<-function(colnm){
  var = rlang::sym(colnm);
  dfr = dfr_haul |> 
           dplyr::select(year,stationid,start_date,charter,!!var) |> 
           tidyr::pivot_wider(id_cols=c("year","stationid","start_date"),
                              names_from="charter",values_from=colnm);
  return(dfr);
}

compareHD<-function(colnm){
  dfr = match_hauls(colnm);
  p = ggplot(dfr,aes(x=NMFS,y=BSFRF,colour=factor(year))) + 
        geom_abline(slope=1) + 
        geom_point() + 
        coord_fixed(ratio=1) +
        labs(colour="year",subtitle=colnm) +
        wtsPlots::getStdTheme();
  return(p)
}

#--check for unique matches 
#--the number of distinct hauls by charter should be identical to the number of matched hauls
dfr_match = match_hauls("haul");
nr_match = nrow(dfr_match);
nr_NMFS  = nrow(dfr_match |> dplyr::distinct(year,stationid,start_date,NMFS));
nr_BSFRF = nrow(dfr_match |> dplyr::distinct(year,stationid,start_date,BSFRF));

print(compareHD("gear_depth"))
print(compareHD("gear_temp"))
print(compareHD("ebssed2_phi"))

#--summarize dfr_rkc data. 
#----mean size here is count-weighted mean of means (not same as Bob's writeup)
dfr = dfr_rkc |> dplyr::select(year,stationid,pairnum,size=length,
                               BSFRF=freq.BSFRF,
                               NMFS=freq.NMFS) |> 
        tidyr::pivot_longer(c("BSFRF","NMFS"),names_to="gear",values_to="num") |> 
        dplyr::group_by(gear,year,stationid,pairnum) |> 
        dplyr::summarize(num_cnt=sum(num),
                         mn_size=sum(size*num)/sum(num>0)) |>
        dplyr::group_by(gear,year) |> 
        dplyr::summarize(num_stns = dplyr::n(),
                         num_non0 = sum(num_cnt>0),
                         num_cnt  = sum(num_cnt),
                         mn_size  = mean(mn_size,na.rm=TRUE)) |> 
        dplyr::ungroup();
#----mean size here is cpue-weighted mean of means
#----mean cpue here will not match Bob's: his is mean over lengths(?)
#------NMFS same as Bob's writeup, BSFRF off by 1 mm
dfr = dfr_rkc |> dplyr::select(year,stationid,pairnum,size=length,
                               BSFRF=cpuenum_km2_new.BSFRF,
                               NMFS=cpuenum_km2_new.NMFS) |> 
        tidyr::pivot_longer(c("BSFRF","NMFS"),names_to="gear",values_to="cpue") |> 
        dplyr::group_by(gear,year,stationid,pairnum) |> 
        dplyr::summarize(tot_cpue=sum(cpue),
                         mn_size=sum(size*cpue)/sum(cpue)) |>
        dplyr::group_by(gear,year) |> 
        dplyr::summarize(num_stns = dplyr::n(),
                         num_non0 = sum(tot_cpue>0),
                         mn_cpue  = mean(tot_cpue),
                         mn_size  = round(sum(mn_size*tot_cpue,na.rm=TRUE)/sum(tot_cpue))) |> 
        dplyr::ungroup();
#----mean size here is count-weighted mean over all (not same as Bob's writeup)
dfr = dfr_rkc |> dplyr::select(year,stationid,pairnum,size=length,
                               BSFRF=freq.BSFRF,
                               NMFS=freq.NMFS) |> 
        tidyr::pivot_longer(c("BSFRF","NMFS"),names_to="gear",values_to="num") |> 
        dplyr::group_by(gear,year) |> 
        dplyr::summarize(num_cnt=sum(num),
                         mn_size=sum(size*num)/sum(num)) |> 
        dplyr::ungroup();
#----mean size here is CPUE-weighted mean over all,
#------NMFS same as Bob's writeup, BSFRF off by 1 mm
dfr = dfr_rkc |> dplyr::select(year,stationid,pairnum,size=length,
                               BSFRF=cpuenum_km2_new.BSFRF,
                               NMFS=cpuenum_km2_new.NMFS) |> 
        tidyr::pivot_longer(c("BSFRF","NMFS"),names_to="gear",values_to="cpue") |> 
        dplyr::group_by(gear,year) |> 
        dplyr::summarize(mn_size=round(sum(size*cpue)/sum(cpue))) |> 
        dplyr::ungroup();

#--plot counts by 5-mm size bins
dfr = dfr_rkc |> dplyr::select(year,stationid,pairnum,size=length,
                               BSFRF=freq.BSFRF,
                               NMFS=freq.NMFS) |> 
        tidyr::pivot_longer(c("BSFRF","NMFS"),names_to="gear",values_to="num") |> 
        dplyr::mutate(size=5*floor(size/5)) |>
        dplyr::group_by(year,size,gear) |> 
        dplyr::summarize(tot=sum(num)) |> 
        dplyr::ungroup();
ggplot(dfr,aes(x=size,y=tot,colour=gear)) + 
  geom_errorbar(aes(ymax=tot),ymin=0,position=position_identity()) + 
  facet_grid(gear~year);

#--plot histograms of counts by 5-mm size bins
dfr = dfr_rkc |> dplyr::select(year,stationid,pairnum,size=length,
                               BSFRF=freq.BSFRF,
                               NMFS=freq.NMFS) |> 
        tidyr::pivot_longer(c("BSFRF","NMFS"),names_to="gear",values_to="num"); 
ggplot(dfr,aes(x=size,colour=gear,weight=num)) + 
  geom_histogram(breaks=seq(0,200,5)) + 
  facet_grid(gear~year);

#--plot mean cpue by 1-mm size bin, year
dfr = dfr_rkc |> dplyr::select(year,stationid,pairnum,size=length,
                               BSFRF=cpuenum_km2_new.BSFRF,
                               NMFS=cpuenum_km2_new.NMFS) |> 
        tidyr::pivot_longer(c("BSFRF","NMFS"),names_to="gear",values_to="cpue");
dfrm = dfr |>
        dplyr::group_by(year,size,gear) |> 
        dplyr::summarize(mean_cpue=mean(cpue)) |> 
        dplyr::ungroup();
ggplot(dfr,aes(x=size,y=cpue,colour=gear)) + 
  geom_errorbar(aes(ymax=cpue),ymin=0,position=position_jitter(1)) + 
  facet_grid(gear~year)
ggplot(dfrm,aes(x=size,y=mean_cpue,colour=gear)) + 
  geom_errorbar(aes(ymax=mean_cpue),ymin=0) + 
  facet_grid(gear~year)



