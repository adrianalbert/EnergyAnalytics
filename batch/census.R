clockMean = function(v) {
  v[v > 12] = -1 * (24 - v[v > 12])
  mean(v)
}

#clockMean(c(13,11))

#############################
# ACS data summary
#
# the data files that contain demographic, housing, income, etc. data by ZCTA are from the 
# American Community Survey (ACS) 2011. The 5yr estimate data spans 2007 to 2012.
# Go here: http://factfinder2.census.gov/faces/nav/jsf/pages/searchresults.xhtml
# choose 5-digit Zip Code Tabulation Areas under "geographies" search option, then choose CA, then "all 5-digit...within CA"
# to locate ACS summary data, choose "Data Profile" under topic : "Product Type"
# This will produce a very manageable list, but to be very specific, you can also choose
# "2011 ACS 5-year estimates" under topic: "Dataset" as well
# the resulting DP* data files will be aggregated by ZCTA and span the breadth of topics covered by the ACS
# Check boxes next to all the data you want to download and wait for the zip file to be generated and download
# 
# The tables of greates interest are:
# DP02: SELECTED SOCIAL CHARACTERISTICS IN THE UNITED STATES (household types, marital status, education, citizenship, ancestry)
# DP03: SELECTED ECONOMIC CHARACTERISTICS (employment status, job, income, poverty)
# DP04: SELECTED HOUSING CHARACTERISTICS (occupancy, units, vintage, tenure, value, costs, fuels)
# DP05: ACS DEMOGRAPHIC AND HOUSING ESTIMATES (sex, age, race)
# file naming convention is ACS_<yr>_<span>_<table>_with_ann.csv
# for example: ACS_11_5YR_DP02_with_ann.csv is the table with social characteristics
# columns are named HC01_VC03, HC02_VC03, HC03_VC03, HC04_VC03, etc.
# HC01 is for the estimated value
# HC02 is for the margin of error for the estimated value
# HC03 is for the estimated percentage of the value
# HC04 is for the percentage margin of error
#
# ACS_<yr>_<span>_<table>_metadata.csv provides human readable names for the coded columns in the data files
# heirarchical levels of categorical breakdowns are separate by '-', 
# with the top level also containing the type of observation followed by a ';'
# For example: HC04_VC12 is Percent Margin of Error; HOUSEHOLDS BY TYPE - Family households (families) - Female householder, no husband present, family - With own children under 18 years
# So this field provides the % margin of error for:
# HOUSEHOLDS BY TYPE
#   Family households (families)
#     Female householder, no husband present, family
#       With own children under 18 years
#
# for entries in the data file:
# (X) means NA
# - means too few samples to compute an estimate
# ** means too few samples to calc margin of error (error cols only)
# a - following a median means the median is in the lowest interval of the distro
# a + following a median means the median is in the highest interval of the distro
# a *** in the margin of error column means that the estimate is in the lowest or highest interval (i.e. stat test inappropriate)
# a ***** in the margin of error means that the estimate is controlled, so a statistical test is inappropriate
# an 'N' in the margin of error indicates that the data cannot be displayed because the sample size is too small
# margin of error means that prob that the true value lies in the range defined by the estimate +- the margin is 90%
# the upshot is essentailly that values that cannot be parsed as numerical can be treated as missing data for the puroses of analysis
# these will be turned into NA values so lm, and functions with rm.na, etc. will ignore them.


loadACS = function(filterErr=T) {
  ACS = acsSocial(filterErr=filterErr)
  ACS = merge(ACS,acsEcon(filterErr=filterErr))     # will merge on both GEO.id and ZCTA
  ACS = merge(ACS,acsHousing(filterErr=filterErr))  # will merge on both GEO.id and ZCTA
  ACS = merge(ACS,acsDemo(filterErr=filterErr))     # will merge on both GEO.id and ZCTA
  return(ACS)
}

acsSocial = function(filterErr=T) {
  fields = c('VC20','VC17','VC18','VC93','VC94','VC118','VC119')
  names  = c('avg_hh_size','hh_under_18','hh_over_65','pop_past_highschool','pop_past_bachelors','res_same_1yr','res_diff_1yr')
  
  return (loadACSTable(table='DP02',colList=fields,colNames=names,filterErr=filterErr))
}
acsEcon = function(filterErr=T) {
  fields = c('VC85','VC86','VC112','VC113','VC115','VC156')
  names  = c('median_hh_income','mean_hh_income','median_fam_income','mean_fam_income','per_cap_income','pct_below_poverty')
  
  return (loadACSTable(table='DP03',colList=fields,colNames=names,filterErr=filterErr))
}
acsHousing = function(filterErr=T) {
  fields = c('VC04','VC63','VC64','VC66','VC67','VC125')
  names  = c('occupied_units','owner_occupied','renter_occupied','owner_hh_size','renter_hh_size','median_home_value')
  
  return (loadACSTable(table='DP04',colList=fields,colNames=names,filterErr=filterErr))
}
acsDemo = function(filterErr=T) {
  fields = c('VC21','VC23','VC29','VC30','VC26','VC33','VC34')
  names  = c('median_pop_age','pop_above_18','pop_above_18_M','pop_above_18_F','pop_above_65','pop_above_65_M','pop_above_65_F')
  
  return (loadACSTable(table='DP05',colList=fields,colNames=names,filterErr=filterErr))
}
# sapply(paste('a_',b,sep=''),FUN=function(x) grep(x,a))

loadACSTable = function(table='DP02',colList=NULL,colNames=NULL,filterErr=T){
  fileName = paste('census/ACS_11/ACS_11_5yr_',table,'_with_ann.csv',sep='')
  ACS <-read.csv(fileName, header=TRUE,as.is=T,na.strings=c('(X)','-','**','***','*****','N'))
  # GEO.id2 is the ZCTA
  colnames(ACS)[2] <- 'ZCTA'
  if(filterErr) {
    valCols = grep('^HC01',colnames(ACS))
    pctCols = grep('^HC03',colnames(ACS))
    ACS = ACS[,c(1,2,valCols,pctCols)]
  }
  if(!is.null(colList)) {
    newCols = c()
    newNames = c()
    suffix = c('val','valerr','pct','pcterr')
    for(i in 1:4) {
      prefix = paste('HC0',i,sep='')
      end    = suffix[i]
      #print(prefix)
      subCols = as.vector(sapply(colList,FUN=function(x) grep(paste(prefix,'_',x,sep=''),colnames(ACS))))
      if(is.numeric(subCols)) { # sapply with no grep matches returns a list, which should be ignored
        #print(class(subCols))
        newCols = c(newCols,subCols)
        #print(newCols)
        if(! is.null(colNames) ) {
          newNames = c(newNames,paste(colNames,'_',end,sep=''))
          #print(newNames)
        }
      }
    }
    ACS = ACS[,c(1,2,newCols)]
    colnames(ACS)[c(-1,-2)] <- newNames
  }
  #print(names(ACS))
  #ACS[ACS == '(X)']        = NA
  #ACS[ACS == '-']          = NA
  #ACS[ACS == '**']         = NA
  #ACS[ACS == '***']        = NA
  #ACS[ACS == '*****']      = NA
  #ACS[ACS == 'N']          = NA
  for (colName in colnames(ACS)[c(-1,-2)]) {
    data = ACS[,colName]
    if(is.character(data)) {
      data = gsub('[^0-9]','',data)
    }
    ACS[,colName] = as.numeric(data)
  }
  return(ACS)
}


# Basic idea behind merging files: Thanks to https://www.braintreepayments.com/braintrust/vaulted-credit-card-maps-with-R
# and https://gist.github.com/braintreeps/5006126#file-map-r
# crosswalk file Zip_to_ZCTA_Crosswalk_2011_JSI.csv from http://udsmapper.org/ziptozctacrosswalk.cfm
# census gazeteer file from http://www.census.gov/geo/maps-data/data/gazetteer2010.html

# merga an arbitrary data frame with zip code information with the associated 2010 ZCTAs
# under column name "ZCTA"
ZIP_TO_ZCTA = NULL
zipToZCTA = function() {
  if(is.null(ZIP_TO_ZCTA)) {
    #read in the zip to zcta mapping (crosswalk) file
    ZIP_TO_ZCTA <-read.csv("census/Zip_to_ZCTA_Crosswalk_2011_JSI.csv", header=TRUE, colClasses=c(rep("character", 5)))
    colnames(ZIP_TO_ZCTA)[5] = 'ZCTA'
  }
  return(ZIP_TO_ZCTA)
}

# merge an arbitrary data frame with ZCTA information with ZCTA data from the associated 2010 gazeteer file
CENSUS_GAZ = NULL
censusGaz = function() {
  if(is.null(CENSUS_GAZ)) {
    #read in the "gazeteer" file that has population, latitude, and longitude for the zcta's
    CENSUS_GAZ <- read.table("census/Gaz_zcta_national.txt", header=TRUE, colClasses=c("character", rep("numeric", 8)))
    # fun fact: the geoid in the gazeteer file IS the ZCTA
    colnames(CENSUS_GAZ)[1] = c('ZCTA')
  }
  return(CENSUS_GAZ)
}

addZCTA = function(df,zipCol='zip5') {
  return( merge(df,zipToZCTA()[,c('ZIP','ZCTA')],by.x=zipCol,by.y='ZIP'))
}

mergeCensus = function(df,zctaCol='ZCTA',zipCol='zip5',censusStats=NULL) {
  # if one isw not already present, use the zip code column to add a ZCTA column
  if(! 'ZCTA' %in% colnames(df)) {
    df = addZCTA(df,zipCol)
  }
  if(is.null(censusStats)) { censusStats = loadACS() }
  # merge the dataframe passed in with ZCTA distinguished rows from the gazeteer file.
  return( merge(df,censusStats,by.x='ZCTA',by.y=zctaCol))
}
  
  
  