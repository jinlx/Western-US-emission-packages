# This script is used to process EPA hourly data in the western US. We select the data of interest based on start/end date, MDL setting, state of interest. 
# In order to make use of this script, we need to 
# 1) change the default directory; 2) default name; 3) default year; 4) default date (the original setting is for summer)
# 5) default MDL setting; 6) default state setting (default is for western US sites) 
# EPA data: https://aqs.epa.gov/aqsweb/airdata/download_files.html#Raw

# Further anlaysis of playing with missing hours is in process_missing_hour.R
# ==================
# customized setting
# ==================
# set the working directory
setwd("./")

#name = "CO"
name = "O3"
#name ="PM25"
year ="2018"

# standard unit is ppm for CO and ozone 
# standard unit is ug/m3 for PM25

# Remove data where MDL is <= 0.020 (0.50000, 0.16000, 0.04000, 0.02000, 0.01100, and 0.00005) for CO
# Remove data where MDL is <= 0.005 (8.0000, 5.0000, 0.0050, 0.0030, 0.0015, and 0.0006) for O3
# Remove data where MDL is <= 0.1 (5.0, 2.0, and 0.1) for PM25


#MDL_value=0.04
MDL_value=0.0050
#MDL_value=5

# customized dates 
start_date = as.Date('2018-06-01')
end_date = as.Date('2018-10-01')

hourly_data <- read.csv(file = paste("hourly_",name,"_",year,".csv",sep = ""))
unique(hourly_data$State.Name)

unique(hourly_data$MDL)

# look at dataframe
#head(hourly_data)
# MDL: method detection limit
#unique(hourly_data$MDL)

# remove unnecessary columns
hourly_data$State.Code <- NULL
hourly_data$County.Code <- NULL
hourly_data$Parameter.Code <- NULL
hourly_data$POC <- NULL
hourly_data$Datum <- NULL
hourly_data$Parameter.Name <- NULL
hourly_data$Units.of.Measure <- NULL
hourly_data$Uncertainty <- NULL
hourly_data$Qualifier <- NULL
hourly_data$Method.Type <- NULL
hourly_data$Method.Code <- NULL
hourly_data$Method.Name <- NULL
hourly_data$Date.of.Last.Change <- NULL
#hourly_data$Site.Num <- NULL
hourly_data$Date.Local <- NULL
hourly_data$Time.Local <- NULL
#head(hourly_data)
# remove GMT time format in to a numeric number.
hourly_data$hours = substr(hourly_data$Time.GMT, start = 1, stop = 2)
hourly_data$Time.GMT <- NULL
#head(hourly_data)

# select months
hourly_data$month <- months(as.Date(hourly_data$Date.GMT))

# select states for western US, Alaska is not included. Customized
hourly_data_states = hourly_data[hourly_data$State.Name == "Arizona" | hourly_data$State.Name == "California" |
                                   hourly_data$State.Name == "Colorado" | hourly_data$State.Name == "Idaho" | hourly_data$State.Name == "Iowa" | 
                                                                      hourly_data$State.Name == "Montana" | hourly_data$State.Name == "Nevada" | hourly_data$State.Name == "New Mexico" |
                                                                                                         hourly_data$State.Name == "Oregon" | hourly_data$State.Name == "Utah" | hourly_data$State.Name == "Washington" |
                                                                                                                                            hourly_data$State.Name == "Wyoming",]

                                                                                                                                            head(hourly_data_states)
                                                                                                                                            unique(hourly_data_states$State.Name)
# select cities in Laing's paper
#hourly_data_JAS_states_Laing = hourly_data_JAS_states[hourly_data_JAS_states$County.Name == 'King' | hourly_data_JAS_states$County.Name == 'Multnomah' | 
#                                                       hourly_data_JAS_states$County.Name == 'Ada' | hourly_data_JAS_states$County.Name == 'Boise' | hourly_data_JAS_states$County.Name == 'Denver' |
#                                                       hourly_data_JAS_states$County.Name == 'San Joaquin' | hourly_data_JAS_states$County.Name == 'Fresno' |
#                                                       hourly_data_JAS_states$County.Name == 'Washoe' | hourly_data_JAS_states$County.Name == 'Butte' | hourly_data_JAS_states$County.Name == 'Missoula',]
#head(hourly_data_JAS_states)

#unique(hourly_data_states$County.Name)
#site_num = unique(hourly_data_states$Site.Num)
#length(site_num)
#unique(hourly_data_states$MDL)

# select data with MDL limit
hourly_data_states_MDL =  hourly_data_states[hourly_data_states$MDL <= MDL_value,]
#head(hourly_data_states_MDL)

#clock = 0
#for (nums in site_num) {
#  clock = clock+1
#  test = hourly_data_states_MDL
#  #lat_tmp = test$Latitude
#  #lon_tmp  =test$Longitude
#  state_names_tmp = test$State.Name
#  county_names_tmp = test$County.Name
#  print('!!!!!!!!!!!!!!!')
#  print(clock)
#  print('!!!!!!!!!!!!!!!')
#  print(unique(state_names_tmp))
#  #print(unique(county_names_tmp))
#}

# summary of the data
#summary(as.factor(hourly_data_states_MDL$County.Name))
#head(hourly_data_states_MDL)

# select the date
hourly_data_final = hourly_data_states_MDL[  hourly_data_states_MDL$Date.GMT >= start_date &  hourly_data_states_MDL$Date.GMT <= end_date,]
dfout = hourly_data_final

# make measurement lt MDL as half of the MDL
ct <- 0
for(i in 1:length(dfout$MDL)){
      if (dfout$Sample.Measurement[i] < dfout$MDL[i])  {
              ct=ct+1
                  dfout$Sample.Measurement[i] = 0.5*dfout$MDL[i]
                    }
}
print(ct)

# check the if there is a missing hour for each day!!!!!!
summary(as.factor(dfout$hours))


# concatenate two string columns;
dfout$location <- paste0(dfout$State.Name, ', ', dfout$County.Name)
#unique(dfout$location)

# choose latitude and longitude
#unique(dfout$Site.Num)
#unique(dfout$Latitude)
#unique(dfout$Longitude)

# check the length of differnet site#, and locations.
length(unique(dfout$Site.Num))
length(unique(dfout$Latitude))
length(unique(dfout$Longitude))

# interesting... The length of differnet site number is different from length of latitude and longitude

# save out this dataframe
write.csv(dfout,paste("hourly_EPA_",name,"_",year,".csv",sep = ""),row.names = FALSE)

# print out done
print('Done!')

# check lat and longitude: 
# Differnet city has differenment number sites, in order to pick up lat and lon,
# we need to know which site we need
#unique(dfout$location)
#unique(dfout$Longitude)
#unique(dfout$Latitude)

#length(unique(dfout$location))
#length(unique(dfout$Longitude))
#length(unique(dfout$Latitude))

# check the number of element/length
#print(length(dfout$Sample.Measurement))


