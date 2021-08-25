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

# separate out brazil and argentina - don't have the banding records for them so keep all
brar_resights <- all_resights_clean %>% 
  filter(FlagID %in% c("FDB", "FEDB", "FO", "FEO"))

other_resights <- all_resights_clean %>% 
  filter(!FlagID %in% c("FDB", "FEDB", "FO", "FEO"))

# remove resights that don't have a banding record
other_resights <- other_resights %>% 
  filter(MetalID %in% band_numbers)

# add brazil and argentina resights back in
all_resights_clean <- bind_rows(other_resights, brar_resights)

# create encounter histories for CJS model
resights_enchist_cjs <- all_resights_clean %>% 
  arrange(ResightDate) %>% 
  dplyr::select(year, BirdID) %>% 
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

################################################################################

# summarize number of observations per year by flag colour

bb_resights <- all_resights_clean %>% 
  dplyr::select(year, FlagID, BirdID) %>% 
  mutate(FlagID = replace(FlagID, str_detect(FlagID, "FELG"), "FLG")) %>% 
  mutate(FlagID = replace(FlagID, str_detect(FlagID, "FEDG"), "FDG")) %>% 
  mutate(FlagID = replace(FlagID, str_detect(FlagID, "FEDB"), "FDB")) %>% 
  mutate(FlagID = replace(FlagID, str_detect(FlagID, "FEO"), "FO")) %>% 
  mutate(FlagID = replace(FlagID, str_detect(FlagID, "FEW"), "FW")) %>% 
  mutate(FlagID = replace(FlagID, str_detect(FlagID, "FER"), "FR")) %>% 
  distinct()

resights_summary <- bb_resights %>% 
  group_by(year, FlagID) %>% 
  tally()

year_summary <- resights_summary %>% 
  group_by(year) %>% 
  summarize(n = sum(n))

flag_summary <- resights_summary %>% 
  group_by(FlagID) %>% 
  summarize(n = sum(n))
