# set the working directory
# This script to fill in missing hour
#setwd("C:/Users/jinlx/CSV_preprocessing/")
#df <- read.csv(file = 'data/hourly_EPA_PM25_2018.csv', check.names = F, stringsAsFactors = F)

setwd("./")
#name = "CO"
name ="O3"
#name="PM25"

year ="2018"
df <- read.csv(file = paste("hourly_EPA_",name,"_",year,".csv",sep = ""), check.names = F, stringsAsFactors = F)
#df <- read.csv(file = paste("hourly_EPA_",name,"_",year,"_fixed.csv",sep = ""), check.names = F, stringsAsFactors = F)

head(df)
dates <- unique(df$Date.GMT) # Gets all unique dates in a vector
sites <- unique(df$location) # Gets unique location names
summary(as.factor(sites))

# Alternatively, install just lubridate:
# install.packages("lubridate")
library(lubridate)
library(dplyr)
## get dates and sites, and create an empty dataframe
merged_data <- df[0,] # Creates an empty dataframe with the same column structure as your data
head(merged_data)


for(i in 1:length(dates)){
    print(dates[i])
    for(j in 1:length(sites)){
        # only return the statement satisfy the filter into the variable
        df.filt <- df %>% filter(location == sites[j])           
        # Gets unique site numbers for each location name
        site.number <- unique(df.filt$Site.Num)
        #print(site.number)
        #break                                
        # for data with multiple site numbers in one location
        for(h in 1:length(site.number)){
            # Extracts data for only one date and site into a dataframe
            single_date <- df %>% filter(Date.GMT == dates[i] & location == sites[j] & Site.Num == site.number[h]) 
            #print(dates[i])
            #print(sites[j])
            #print(site.number[h])
            #print(length(single_date$location))
            #print(single_date)
            #stop

            if(length(single_date$location) == 0){ # Works to NA an entire date if a site is missing an entire date.
                k = 1
                while(length(single_date$location) == 0){
                    single_date <- df %>% filter(Date.GMT == dates[k] & location == sites[j] & Site.Num == site.number[h])
                    k = k+1
                    #print(dates[k])
                    #print(sites[j])
                    #print(site.number[h])
                }
                hours_frame <- data.frame('hours'=seq(0,23),'Site.Num' = single_date$Site.Num[1],"Latitude" = single_date$Latitude[1]
                                          ,"Longitude" = single_date$Longitude[1],"Date.GMT" = dates[i],"MDL" = single_date$MDL[1]
                                          ,"State.Name" = single_date$State.Name[1],"County.Name" = single_date$County.Name[1]
                                          ,"month" = lubridate::month(dates[i], label = T, abbr = F),"location" = single_date$location[1])
                merged_data <- merge(merged_data,hours_frame,all=T)
            }else{
                hours_frame <- data.frame('hours'=seq(0,23),'Site.Num' = single_date$Site.Num[1],"Latitude" = single_date$Latitude[1]
                                          ,"Longitude" = single_date$Longitude[1],"Date.GMT" = single_date$Date.GMT[1],"MDL" = single_date$MDL[1]
                                          ,"State.Name" = single_date$State.Name[1],"County.Name" = single_date$County.Name[1]
                                          ,"month" = single_date$month[1],"location" =  single_date$location[1])
                semi_merged_data <- merge(hours_frame, single_date, all=T)
                #print(semi_merged_data)
                merged_data <- merge(merged_data,semi_merged_data,all=T)
            }
        }
    }
}
# save out this dataframe
write.csv(merged_data,paste("hourly_EPA_",name,"_",year,"_fixed.csv",sep = ""),row.names = FALSE)
print('Done!')


