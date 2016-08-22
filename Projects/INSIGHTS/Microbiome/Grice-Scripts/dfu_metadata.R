# packages
library(loescher)
library(qiimer)
library(dplyr)
library(tidyr)

# Set working directory
setwd('~/Coding/Club_Grice/scripts/loesche/dfu_uclust/metadata/')

# Read in mapping file and key ----------------------------------------------------------------

## Give mapping file key to subject visit pairs
## read in mapping file
map <- read_qiime_mapping_file('./dfu_map.txt') %>% 
  mutate(Plate = as.factor(Plate), WellPosition = as.character(WellPosition)) %>%
  select(SampleID,WellPosition,Plate)

## read in sample key
map <- read.csv('./v1v3_sample_key.csv') %>% 
  separate(SubjectID, c('SubjectID','Visit'),sep = '-') %>%
  select(-ProjectID) %>%
  mutate(Visit = as.integer(Visit), 
         SampleID = as.character(SampleID),
         Notes = as.character(Notes)) %>%
  full_join(x = map, y = .,by = 'SampleID') %>%
  filter(!SubjectID %in% c('135','150'))

# Add number of raw samples and the last visit
map <- map %>% group_by(SubjectID) %>% 
  summarise(num_raw_samples = length(Visit)) %>%
  left_join(x = map, y = ., by = 'SubjectID')


# Read in clinical metadata -------------------------------------------------------------------

# Read in table
dfu <- read.table('complete_metadata_wide.txt',sep='\t',header=T) %>%
  # Make into long format
  reshape(timevar = 'visit', idvar = 'studyid', varying = 20:257, 
          sep = '', direction = 'long') %>% 
  # rename variables
  rename(SubjectID = studyid, Visit = visit) %>%
  # format variables
  mutate(SubjectID = as.character(SubjectID)) %>%
  # arrange entries by subject and visit
  arrange(SubjectID, Visit) %>%
  # remove subjects that we have no samples for
  filter(SubjectID %in% map$SubjectID) %>% droplevels

# Remove blank rows, ie visits that didn't happen
dfu <- dfu[rowSums(is.na(dfu[,21:37])) != 17,]

# Merge with mapping file
dfu <- full_join(dfu, map, by = c('SubjectID','Visit')) %>%
  # temporarily remove control samples, will be added back later. For the sake of reducing NA's
  filter(SubjectID != 'Control') %>% droplevels

# Make 0/1 variables into booleans
dfu$enrolled_new_criteria <- dfu$enrolled_new_criteria == 1
dfu$antbiostart <- dfu$antbiostart == 1
dfu$topical_tx_antibiotic <- dfu$topical_tx_antibiotic == 1
dfu$topical_tx_other <- dfu$topical_tx_other == 1
dfu$wd_clinical <- dfu$wd_clinical == 1
dfu$wd_size <- dfu$wd_size == 1
dfu$osteomye <- dfu$osteomye == 1

# Turn coded data into factors
dfu$DiabType <- factor(dfu$DiabType, levels = c(1,2), labels = c('Type 1', 'Type 2'))
dfu$Sex <- factor(dfu$Sex, levels = c(2,1), labels = c('Female', 'Male'))
dfu$race_cat <- factor(dfu$race_cat, levels = c(1,2,9), labels = c('White', 'Other', 'Unreported'))
dfu$EOSReas <- factor(dfu$EOSReas, levels = c(1,2,4,5,6,8,9), labels = c('Unhealed','Dropped-Subject','Amputation','Other Infection','Healed','Dropped-Study','Other'))
dfu$ulcerloc_cat <- factor(dfu$ulcerloc_cat, levels = c(1:3), labels = c('Forefoot','Midfoot','Hindfoot'))
dfu$necrotisprior00 <- NULL
dfu$SampleID <- as.character(dfu$SampleID)

## Make new variables
dfu <- dfu %>% group_by(SubjectID) %>% 
  mutate(idx.last.visit = max(Visit),
         last.visit = idx.last.visit == Visit,
         days = round(weeks * 7),
         UlcerDur = UlcerDur + (Visit * 2),
         healed = EOSReas == 'Healed',
         idx.b4.healing = visit_healed - 1)

### Code visit types
# Did this visit precede healing?
dfu$v.b4.healing <- unlist(
  sapply(unique(dfu$SubjectID), simplify = F, USE.NAMES = T, function(x) {
    tmp <- subset(dfu, SubjectID == x)
    if(any(is.na(tmp$healed))) {
      return(rep(NA,nrow(tmp)))
    } else if(any(tmp$healed)) {
      v.b4.healing <- tmp$idx.b4.healing == tmp$Visit
    } else {
      v.b4.healing <- rep(F, nrow(tmp))
    }
    return(v.b4.healing)
  })
)

## Complication
# Condense wound deterioration types
dfu$v.wd <- dfu$wd_clinical | dfu$wd_size
# Boolean visit before wd, egrice style - does not distinguish between runs of complications
dfu$v.b4.wd.eg <- unlist(
  sapply(unique(dfu$SubjectID), simplify = F, USE.NAMES = T, function(x) {
    tmp <- subset(dfu, SubjectID == x)
    result <- rep(F, nrow(tmp))
    if(length(result) > 0) result[which(tmp$v.wd) - 1] <- T
    return(result)
  }))
# Boolean visit before wd, excluding visits that are currently in wd
dfu$v.b4.wd <- 
  sapply(1:nrow(dfu), function(x) {
    if(!dfu$v.b4.wd.eg[x]) return(F)
    if(is.na(dfu$v.wd[x])) {
      return(NA)
    } else if(dfu$v.wd[x]) {
      return(F)
    } else {
      return(T)
    }})

## Amputation
# Subject get amputation?
dfu$amputation <- dfu$EOSReas == "Amputation"
# Is this the visit before amputation?
dfu$v.b4.amp <- F
dfu$v.b4.amp[dfu$amputation] <- (dfu$Visit == dfu$idx.last.visit)[dfu$amputation]

## Osteo
# Does subject have osteomyelitis currently?
dfu$v.osteo <- dfu$osteomye
# Is this the visit before developing osteo?
dfu$v.b4.osteo <- unlist(
  sapply(unique(dfu$SubjectID), function(x){
    tmp <- droplevels(subset(dfu, SubjectID == x))
    v <- which(tmp$v.osteo)
    result <- rep(F, nrow(tmp))
    if(length(v) == 0) return(result)
    else {
      result[(v-1)] <- T
      return(result)
    }
  }) )

## Make combined complication column: wd, osteo, amputation
# Visit of complication, no samples for amputation
dfu$v.comp <- apply(dfu[,c('v.wd','v.osteo')],1,any,na.rm=T)
# Visit before complication excluding the runs of complications
dfu$v.b4.comp <- apply(dfu[,c('v.b4.wd','v.b4.amp','v.b4.osteo')],1,any,na.rm=T)
# Visit before complication excluding runs of complications, and making WD visits NA's
dfu$v.b4.comp.na <- dfu$v.b4.comp
dfu$v.b4.comp.na[dfu$v.comp] <- NA
# Visit before complication including the runs of complications
dfu$v.b4.comp.eg <- apply(dfu[,c('v.b4.wd.eg','v.b4.amp','v.b4.osteo')],1,any,na.rm=T)


## Outcomes Variables
dfu <- dfu %>% 
  # Heals in 12 wks from baseline
  mutate(heals.12wks = this2that(visit_healed <= 6, NA, F),
         # Heals in 4 wks from baseline
         heals.4wks = this2that(visit_healed <= 2, NA, F),
         # Overall complication: does subject develop complication at some point
         comp = any(v.comp,na.rm=T),
         # Complication 12 wks from baseline
         comp.12wks = any(Visit[v.comp] <= 6,na.rm=T),
         # Complication 4 wks from baseline
         comp.4wks = any(Visit[v.comp] <= 2,na.rm=T))
# Make a multi-level outcome variable that combines multiple outcomes
dfu <- dfu %>% group_by(SubjectID, Visit) %>% 
  summarise_each(funs(first)) %>%
  mutate(unhealed = EOSReas == 'Unhealed') %>% 
  select(SubjectID,Visit,amputation,unhealed,healed,heals.12wks,heals.4wks) %>%
  mutate(healed = healed & !heals.12wks,
         heals.12wks = heals.12wks & !heals.4wks) %>% 
  gather('Outcome','bool',amputation:heals.4wks) %>% 
  filter(bool) %>%
  select(-bool) %>% 
  left_join(x = dfu, y = ., by = c('SubjectID','Visit'))

# Remove variables
dfu <- select(dfu, -wks_to_heal,-weeks,-study_week_healed,-osteomye) %>% rename(hgba1c = hgb)

# Merge with mapping file to make DFU table
dfu <- full_join(dfu,map[map$SubjectID == 'Control',], by = c('SampleID','WellPosition','Plate','SubjectID','Visit','Notes','num_raw_samples')) %>%
  mutate(has_raw_sample = !is.na(SampleID))

# Make new subject id column that is ordered by number of visits
tmp <- dfu %>% group_by(SubjectID) %>% summarise(n = first(idx.last.visit)) %>% arrange(desc(n))
dfu$SubjOrder <- factor(dfu$SubjectID, levels = tmp$SubjectID)


# Combine Multiple Data Files -----------------------------------------------------------------

## Total qPCR bacterial burden
dfu <- read.csv('./bacterial_burden.csv',header=T) %>% 
  mutate(SubjectID = as.character(SubjectID), Notes = as.character(Notes)) %>%
  left_join(dfu, ., by = c('SubjectID','Visit','Notes'))

## Read in new data and merge with dfu data frame
dfu <- read.table('./strep_qpcr.txt',sep='\t',header=T) %>% 
  separate(SwabID, c('patient_id','visit'),'_',remove = T, convert = T) %>%
  rename(SubjectID = patient_id, Visit = visit) %>% 
  mutate(SubjectID = as.character(SubjectID), log_strep = inf2na(log10(strep.qpcr))) %>%
  select(-strep.qpcr) %>%
  left_join(dfu, ., by = c('SubjectID','Visit'))

## Add column for offloading therapy type
dfu <- read.csv('./offloading_therapy.csv',header=T) %>%
  mutate(SubjectID = as.character(SubjectID)) %>%
  left_join(dfu, ., by = 'SubjectID')

## Add necrotic tissue
dfu <- read.csv('./necrotic_tissue.csv',header=T) %>% 
  mutate(SubjectID = as.character(SubjectID)) %>%
  left_join(dfu, ., by = c('SubjectID','Visit'))


# Read in Antibiotics Data --------------------------------------------------------------------

# Read in antiobiotic tables
abx <- read.csv('./abx_data.csv', header = T, colClasses = c('character','integer', rep('factor',3), rep('numeric',3),'character'))

# Format factor levels
abx$abx.type <- factor(abx$abx.type, levels = c(1:8,'Missing'), labels = c('Aminoglycosides','Cephalosporins','Fluoroquinolones','Ketolides','Penicillin','Sulfonamides','Tetracycline','Misc','Missing'))
abx$abx.type[is.na(abx$abx.type)] <- 'Missing'
abx$abx.route <- factor(abx$abx.route, labels = c('IV','Oral'))
abx$abx.reason <- factor(abx$abx.reason, labels = c('Study Ulcer','Non-Study Ulcer','UTI','URI','Sinus Infection','Other'))

# Spread data
abx.w <- abx %>% select(SubjectID:abx.type) %>% 
  arrange(SubjectID, Visit, abx.reason) %>% 
  group_by(SubjectID,Visit,abx.type) %>% 
  summarise(n = length(abx.reason), abx.reason = first(abx.reason)) %>% 
  spread(key = abx.type, value = n, fill = 0)

# Condense patients with more than one visit. In practice this means these patients have antibiotics
# for more than one indication. The antibiotic reason levels are in order of precedence, so taking
# the first entry makes the most sense.
tmp <- abx.w %>% group_by(SubjectID,Visit) %>% 
  select(Aminoglycosides:Missing) %>% 
  summarise_each(funs(sum))
abx.w <- abx.w %>% group_by(SubjectID,Visit) %>% 
  summarise(abx.reason = first(abx.reason)) %>% 
  full_join(tmp, by = c('SubjectID','Visit'))

# Visit that antibiotics were administered
dfu$v.abx <- FALSE
dfu$v.abx[match(paste(abx.w$SubjectID,abx.w$Visit), paste(dfu$SubjectID,dfu$Visit))] <- TRUE
dfu$antbiostart <- NULL
# Create visit before and post antibiotic treatment
dfu <- sapply(unique(dfu$SubjectID), simplify = F, function(p) {
  # Subset and initialize new columns
  df <- subset(dfu, SubjectID == p)
  df$v.b4.abx <- FALSE
  df$v.post.abx <- FALSE
  # Identify visit number that abx were given 
  idx <- df$Visit[which(df$v.abx)]
  # If no antibiotics given, just return default values
  if(length(idx) == 0) return(df)
  # Loop through the antibiotic visits to identify the visit before and after
  for(i in idx) {
    # Mark visit before antibiotic column
    df$v.b4.abx[df$Visit == (i - 1)] <- TRUE
    # Mark visit after antibiotic column
    df$v.post.abx[df$Visit == (i + 1)] <- TRUE
  }
  return(df)
}) %>% do.call(what = 'rbind') %>% ungroup

# Make column to identify subjects that received antibiotics
dfu <- dfu %>% group_by(SubjectID) %>% mutate(antibiotics = any(v.abx,na.rm = T))

# Cleanup
rm(tmp,map)
detach("package:dplyr", unload=TRUE)
detach("package:loescher", unload=TRUE)
detach("package:qiimer", unload=TRUE)
detach("package:tidyr", unload=TRUE)
