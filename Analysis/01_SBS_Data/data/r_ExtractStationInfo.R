#--extract trawl station info for BBRKC SBS study
require(ggplot2)
dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"Analysis/01_SBS_Data/data");

#--processing info
minYr = 2013;
maxYr = 2016;
verbosity = 0;

#--read NMFS strata definitions for RKC in standard format
dirData_NMFS   <-"~/Work/StockAssessments-Crab/Data/Survey.NMFS.EBS/Current";
fnStrata       <-file.path(dirData_NMFS,"AllCrab_SurveyStrata.csv");
dfrStrata<-readr::read_csv(fnStrata) |>
             dplyr::select(STATION_ID,DISTRICT,TOWS,SURVEY_YEAR,
                           LATITUDE,LONGITUDE,STRATUM,TOTAL_AREA_SQ_NM) |> 
             dplyr::filter(dplyr::between(SURVEY_YEAR,minYr,maxYr));
dfrSD<-selectStrata.TrawlSurvey(dfrStrata,
                                species="RKC",
                                strataType="2015",
                                export=FALSE,
                                verbosity=verbosity) |> 
         dplyr::filter(dplyr::between(YEAR,minYr,maxYr));
rm(dfrStrata);

#--read McConnaughey haul csv files
dfr_haul = readr::read_csv(file.path(dirThs,"final_haul.csv")) |> 
             dplyr::mutate(MID_LATITUDE=0.5*(start_latitude+end_latitude),
                           MID_LONGITUDE=0.5*(start_longitude+end_latitude)) |> 
             dplyr::select(type=charter,
                           YEAR=year,
                           GIS_STATION=stationid,
                           HAULJOIN=pairnum,
                           HAUL_TYPE=3,
                           START_DATE=start_date,
                           START_HOUR=start_time,
                           MID_LATITUDE,MID_LONGITUDE,
                           AREA_SWEPT_VARIALE=aswept_nm2_new, #--station area in nm^2
                           GEAR_TEMPERATURE=gear_temp,
                           BOTTOM_DEPTH=gear_depth);

#--extract distinct SBS station data in standard format, 
#----dropping 19 stations not in the Bristol Bay stratum
dfrSD_SBS = dfrSD |> dplyr::right_join(dfr_haul |> dplyr::distinct(YEAR, GIS_STATION),
                                       by=c("YEAR", "GIS_STATION")) |> 
                     dplyr::filter(STRATUM=="Bristol Bay");

#--extract NMFS haul data
dfrHD_NMFS = dfr_haul |> dplyr::filter(type=="NMFS") |> dplyr::select(!type);
#--extract BSFRF haul data
dfrHD_BSFRF = dfr_haul |> dplyr::filter(type=="BSFRF") |> dplyr::select(!type);

#--read McConnaughey combined data csv files
dfr_cd = readr::read_csv(file.path(dirThs,"final_rkc_sr_melt_acoustic.csv")) |> 
           dplyr::select(HAULJOIN=pairnum,SIZE=length,freq.BSFRF,freq.NMFS);
#--extract NMFS ID data
dfrID_NMFS = dfr_cd |> dplyr::select(HAULJOIN,
                                     numIndivs=freq.NMFS,
                                     SIZE) |> 
               dplyr::mutate(SEX="UNDETERMINED",
                             MATURITY="UNDETERMINED",
                             SHELL_CONDITION="UNDETERMINED",
                             SEX_CODE=1,
                             SHELL_CONDITION_CODE=1,
                             EGG_COLOR=NA,
                             EGG_CONDITION=NA,
                             CLUTCH_SIZE=NA,
                             CHELA_HEIGHT=NA,
                             SAMPLING_FACTOR=1,
                             WEIGHT=0,
                             CALCULATED_WEIGHT=0) |> 
                dplyr::filter(numIndivs>0);
dfrID_BSFRF = dfr_cd |> dplyr::select(HAULJOIN,
                                     numIndivs=freq.BSFRF,
                                     SIZE) |> 
               dplyr::mutate(SEX="UNDETERMINED",
                             MATURITY="UNDETERMINED",
                             SHELL_CONDITION="UNDETERMINED",
                             SEX_CODE=1,
                             SHELL_CONDITION_CODE=1,
                             EGG_COLOR=NA,
                             EGG_CONDITION=NA,
                             CLUTCH_SIZE=NA,
                             CHELA_HEIGHT=NA,
                             SAMPLING_FACTOR=1,
                             WEIGHT=0,
                             CALCULATED_WEIGHT=0) |> 
                dplyr::filter(numIndivs>0);


