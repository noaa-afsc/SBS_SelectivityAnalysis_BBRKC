#--extract data for selectivity analysis
#--Need:                                        available from:
#----year, start_date, stationid, pairnum,       dfr_rkc 
#----gear_depth and gear_temp from NMFS hauls    dfr_hauls
#----aswept_km2_new.BSFRF,
#----freq.BSFRF, aswept_km2_new.BSFRF,           dfr_rkc
#----freq.NMFS,  aswept_km2_new.NMFS,            dfr_rkc
#----ebssed2_phi, ebssed2_sort, ebssed2_skew     dfr_rkc
require(ggplot2)
dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"Analysis/01_SBS_Data/data");

#--read csv files
dfr_haul = readr::read_csv(file.path(dirThs,"final_haul.csv"));
dfr_rkc  = readr::read_csv(file.path(dirThs,"final_rkc_sr_melt_acoustic.csv"));

#--extract catch abundance and "melt"
dfr_ext = dfr_rkc |>
            dplyr::select(year,stationid,pairnum,
                          length,freq.BSFRF,freq.NMFS) |>
            tidyr::pivot_longer(cols=c("freq.BSFRF","freq.NMFS"),
                                values_to="N") |> 
            dplyr::mutate(charter=stringr::str_sub(name,6,-1)) |> 
            dplyr::select(year,stationid,pairnum,length,charter,N);

#--bin catch abundance-at-length to 5mm CL bins, calc total abundance and proportion NMFS
cutpts = seq(0,250,5);
dfr_ext$bin5 = cutpts[cut(dfr_ext$length,breaks=cutpts)]+2.5;
dfrProps = dfr_ext |>
            dplyr::mutate(bin5 = cutpts[cut(length,breaks=cutpts)]+2.5) |>
            dplyr::group_by(year,stationid,pairnum,charter,bin5) |> 
            dplyr::summarize(N=wtsUtilities::Sum(N)) |> 
            dplyr::ungroup() |>  
            tidyr::pivot_wider(id_cols=c("year","stationid","pairnum","bin5"),
                               names_from="charter",names_prefix="N.",values_from="N") |> 
            dplyr::mutate(N.tot=N.NMFS+N.BSFRF,
                          propNMFS=N.NMFS/(N.NMFS+N.BSFRF));

#--extract area swept estimates
dfr_AS = dfr_haul |> 
           tidyr::pivot_wider(id_cols=c("year","stationid","pairnum"),
                              names_from="charter",names_prefix="as.",values_from="aswept_km2_new") |> 
           dplyr::mutate(q=(1/as.BSFRF)/(1/as.NMFS)); # q = expfBSFRF/expfNMFS such that logit(propNMFS) = ln(r(z))+ln(q)
           
#--join area swept estimates 
dfrProps = dfrProps |> 
             dplyr::inner_join(dfr_AS,by=c("year","stationid","pairnum"))
             
#--extract haul and environmental data from NMFS hauls
dfrNMFS = dfr_haul |> 
            dplyr::filter(charter=="NMFS") |>
            dplyr::select(year,stationid,pairnum,
                          gear_depth,gear_temp,ebssed2_phi,ebssed2_sort,ebssed2_skew);

#--merge environmental data
dfrAllSBS = dfrProps |>
              dplyr::inner_join(dfrNMFS,by=c("year","stationid","pairnum")) %>% 
              dplyr::arrange(year,stationid,pairnum,gear_depth,gear_temp,ebssed2_phi,ebssed2_sort,ebssed2_skew,
                             as.NMFS,as.BSFRF,q,bin5,N.NMFS,N.BSFRF,N.tot,propNMFS) |> 
              dplyr::mutate(lnR=log(propNMFS/(1-propNMFS))-log(q));
wtsUtilities::saveObj(dfrAllSBS,"rda_dfrAllSBS.RData");



