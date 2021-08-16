library(tidyverse)
library(lubridate)
library(reshape)
library(readxl)

# load data
jb_data <- read_excel("./data/jb_project_resights.xlsx", sheet = "Sheet1")
banding_records <- read_excel("./data/banding_records.xlsx", sheet = "REKNcaptures")

# start cleaning data
# start filtering process
all_resights <- jb_data %>% 
  mutate(year = year(ResightDate)) %>% 
  filter(year %in% 2009:2018)

# remove uncertain resights
all_resights_clean <- all_resights %>% 
  mutate(ResightCertainty = str_replace(ResightCertainty, "^75% \\(COULD BE TV9\\)\\?$", "75")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^100%$", "100")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^0.95$", "95")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^0.75$", "75")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^NOT 100% SURE$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^POSSIBLY H$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^90%$", "90")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^1$", "100")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^0.8$", "80")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^75%$", "75")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^NOT 100% - COULD HAVE BEEN AHU\\?$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^N VERY FADED$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^0.5$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^80%$", "80")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^0.9$", "90")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^COULD BE N9V$", "50")) %>% 
  mutate(ResightCertainty = str_replace(ResightCertainty, "^Y$", "100")) %>% 
  mutate(ResightCertainty = str_replace(ResightCertainty, "^C$", "100")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^\\?$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^N$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^L$", "50"))

all_resights_clean$ResightCertainty <- as.numeric(all_resights_clean$ResightCertainty)  

all_resights_clean <- all_resights_clean %>% 
  filter(!str_detect(FlagCode, "Q")) %>% 
  filter(is.na(ResightCertainty) | ResightCertainty > 94 | ResightCertainty == 100 )

# remove/correct manually identified suspect detections
unique(all_resights_clean$FlagID)

suspect_resights <- all_resights_clean %>%
  filter(FlagID %in% c("FEY", "FEBK", "EY", "FEDP", "NA", "CB", "FEGY"))

all_resights_clean <- all_resights_clean %>% 
  filter(!FlagID %in% c("FEY", "FEBK", "EY", "FEDP", "CB", "FEGY"))

# pull out band numbers from banding records
band_numbers <- banding_records %>% 
  select(MetalID) %>% 
  distinct() %>% pull()

# remove resights that don't have a banding record
all_resights_clean <- all_resights_clean %>% 
  filter(MetalID %in% band_numbers)

# create encounter histories for CJS model
resights_enchist_cjs <- all_resights_clean %>% 
  arrange(ResightDate) %>% 
  select(year, BirdID) %>% 
  distinct() %>% 
  arrange(year) %>% 
  mutate(Seen = 1)

length(unique(resights_enchist_cjs$BirdID))

# create encounter histories
enchist_cjs <- resights_enchist_cjs %>% 
  distinct() %>% 
  arrange(year) %>% 
  tidyr::pivot_wider(id_cols = BirdID, names_from = year, values_from = Seen) %>% 
  replace(is.na(.), 0)

# get number of rekn resighted per period
apply(enchist_cjs, 2, function(x) length(which(!is.na(x))))

# export encounter histories
saveRDS(enchist_cjs, file = "./data/rekn_enchist.rds")

# create INP file to check assumptions
# now do formatting of encounter history file for .inp file for MARK
enchist <- enchist_cjs

enchist$eh <- apply(enchist[2:ncol(enchist)],1,paste,collapse="") # concatenates encounter columns into eh
enchist[2:(ncol(enchist)-1)] <- NULL # drops individual encounter columns

# create commented tag
enchist$BirdID <- paste("/*", enchist$BirdID, "*/", sep=" ")

# sort by descending encounter histories
enchist <- enchist[order(enchist$eh,decreasing=TRUE),]

# tack on the frequency for the individual
enchist$end <- "1;"

write.table(enchist, file = "./data/rekn_enchist.inp", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
  
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  nind <- dim(CH)[1]
  
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

enchist_cjs <- enchist_cjs[,2:39]

# # remove capture histories of individuals marked in last occasion
# get.first <- function(x) min(which(x!=0))
# first <- apply(enchist_cjs, 1, get.first)
# last.only <- which(first==dim(enchist_cjs)[1])
# if (length(last.only) > 0) enchist_cjs <- enchist_cjs[-last.only,]

marr <- marray(enchist_cjs)

write.csv(marr, file = "./processed-data/rekn-cjs-marray.csv", row.names = FALSE)
saveRDS(marr, file = "./processed-data/rekn-cjs-marray.rds")

