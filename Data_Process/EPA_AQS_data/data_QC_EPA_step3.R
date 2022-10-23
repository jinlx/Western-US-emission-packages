
# ==========================================================================================
# This script is used to provide lat and lon and good lon and lat for available ground sites
# ==========================================================================================

# ==================
# customized setting
# ==================
# set the working directory
#setwd("C:/Users/lj152920/Dropbox/data/Ground/")
setwd("./")
name = "CO"
#name = "O3"
#name ="PM25"

year ="2018"

# standard unit is ppm for CO and ozone 
# standard unit is ug/m3 for PM25

# Remove data where MDL is <= 0.005 (8.0000, 5.0000, 0.0050, 0.0030, 0.0015, and 0.0006) for O3
# Remove data where MDL is <= 0.020 (0.50000, 0.16000, 0.04000, 0.02000, 0.01100, and 0.00005) for CO
# Remove data where MDL is <= 0.1 (5.0, 2.0, and 0.1) for PM25

hourly_data <- read.csv(file = paste("hourly_EPA_",name,"_",year,"_fixed.csv",sep = ""))

head(hourly_data)
length(unique(hourly_data$location))
length(unique(hourly_data$Latitude))
length(unique(hourly_data$Longitude))

Lats = unique(hourly_data$Latitude)
Lons = unique(hourly_data$Longitude)

dataoutput = data.frame(Lats, Lons)



# save out this dataframe
write.csv(dataoutput,paste("hourly_EPA_",name,"_",year,"_fixed_locations.csv",sep = ""),row.names = FALSE)
print('Done!')

# avaialbe data for each site
pct = c()
for (lat in Lats) {
  temp = hourly_data[hourly_data$Latitude == lat,]
  temp_value = temp$Sample.Measurement
  pct = c(pct, sum(is.na(temp_value))/length(temp_value))
}
pct = 100 - pct * 100


pct_each_day = c(0,0,0)
# avaiable data for each day and each site
# use this to get sites that they have >80% avaialbe data for each day
Days = unique(hourly_data$Date.GMT)
Sites = unique(hourly_data$Site.Num)

goodsite_lat = c()
ct = 0 

for (lat in Lats) {
  
  test = c()
  
  for (day in Days) {
  
    temp = hourly_data[hourly_data$Latitude == lat & hourly_data$Date.GMT == day,]
    
    temp_value = temp$Sample.Measurement
    
    tmp = 100 - sum(is.na(temp_value))/length(temp_value)*100
    
    test = c(test,tmp)
  }
  
  
  # 75 sites still satisfy the requirement
  max = 0
  for (tt in test) {
    if (tt<80){
      break
    }
    if (tt>=80){
      max = max + 1
    }

    if (max == length(test)){
      goodsite_lat = c(goodsite_lat,lat)
      ct = ct + 1
      print(ct)
      print(lat)
    }
  }
}

goodlat = c()
goodlon = c()
# Save out data with good sites
for (lat in goodsite_lat){
  gooddata = hourly_data[hourly_data$Latitude==lat, ]
  goodlat = c(goodlat,unique(gooddata$Latitude))
  goodlon = c(goodlon,unique(gooddata$Longitude))
}
  
length(goodsite_lat)
length(goodlat)
length(goodlon)

gooddataoutput = data.frame(goodlat, goodlon)

# save out this dataframe
write.csv(gooddataoutput,paste("hourly_EPA_",name,"_",year,"_fixed_locations_QC.csv",sep = ""),row.names = FALSE)
print('Done!')

# =========
# plotting
# =========
# histogram plot for QC
# # of sites need to be changed
h <- hist(pct,
          main="Data Coverage frequency over all sites",
          xlab="Data Coverage over WE-CAN period",
          xlim=c(0,100),
          col="darkmagenta",
          freq=TRUE
)


text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.4))


