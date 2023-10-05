library(here)
library(tidyverse)
library(lubridate)
library(readxl)
library(tidylog)
set.seed(2023)


## files to load
sheets <- list.files(here::here("data/annual_data_updates/data"), recursive = TRUE, full.names = TRUE)


##### prep survival data ####

## function to clean up survival data
prep.surv <- function(path) {
  df <- read_excel(path, sheet = "Survival")%>%
    filter_all(any_vars(!is.na(.)))
  herd.name <- read_excel(path, sheet = 1)[1, 2] %>%
    as.character()
  print.path <- basename(path)

  if (nrow(df) > 0) {
    ## checks
    year <- str_sub(path, -9, -6) %>%
      as.numeric()

    if (all((df$`Date exit` - df$`Date entry`) > 365)) { ## all within 1 year
      stop(paste0("check dates, >365 days ",   print.path))
    }

    if (!all(df$`Date entry` >= ymd(paste0(year - 1, "-04-01")))) { ## check start date
      stop(paste0("start date too early ",   print.path))
    }

    if (!all(df$`Date exit` <= ymd(paste0(year, "-03-31")))) { ## check end date
      stop(paste0("end date too late ",   print.path))
    }
    
    if (any(!df$Outcome%in%c("Alive", "Censored","Mortality"))) { 
      stop(paste0("Outcome has an invalid value ",   print.path))
    }
    
    if (any(!df$Sex%in%c("M","F"))) { 
      stop(paste0("Sex has an invalid value ",   print.path))
    }
    
    if (any(!df$Ageclass%in%c("Adult", "Juvenile", "Calf"))) { 
      stop(paste0("Ageclass has an invalid value ",   print.path))
    }
    
    if (any(duplicated(df$`Animal ID`,incomparables=NA))) { 
      stop(paste0("Animal ID duplicates ",   print.path))
    }

    if (any(duplicated(df$WLHID,incomparables=NA))) { 
      stop(paste0("WLHID duplicates ",   print.path))
    }
    
    if(is.na(herd.name)){
      stop(paste0("herd name missing ",   print.path))
    }

    df <- df %>%
      mutate(
        herd = herd.name,
        `Animal ID` = as.character(`Animal ID`)
      ) %>%
      select(herd, `Animal ID`, WLHID, Sex, Ageclass, `Date entry`, `Date exit`, Outcome)

    return(df)
  }
}

## loop over all data
path <- sheets
surv.list <- list()

for (i in 1:length(sheets)) {
  surv.list[[i]] <- prep.surv(path[i])
}

surv.add <- bind_rows(surv.list)

write_csv(surv.add, here::here("data/raw/surv.2022onwards.csv"))





##### prep count data ####


## function to clean up survival data
prep.count <- function(path) {
  df <- read_excel(path, sheet = "Counts", skip = 1, col_names = TRUE, na = c("", " ", NA, "N/A"))%>%
    filter_all(any_vars(!is.na(.)))
  herd.name <- read_excel(path, sheet = 1)[1, 2] %>%
    as.character()
  print.path <- basename(path)

  if (nrow(df) > 0) {
    ## checks
    year <- str_sub(path, -9, -6) %>%
      as.numeric()

    if (any(is.na(df$`Official Survey Timing`))) {
      stop(paste0("Missing month ",   print.path))
    }

    if (any(!is.na(df$`N collar available`))) {
    if (any(df$`N collars detected` > df$`N collar available`)) {
      stop(paste0("Too many collars detected ",   print.path))
    }
    }

    if (sum(df$`Total Counted (OTC)`,
            na.rm=TRUE)!= sum(df$`N. Calves (OTC)`,
      df$`N. Adult F (OTC)`,
      df$`N. Adult M (OTC)`,
      df$`N. Adult (Sex Unclassified) (OTC)`,
      df$`Unclassified Life Stage and Sex (OTC)`,
      df$`Yrling - Unclas Sex (OTC)`,
      na.rm = TRUE
    )) {
      stop(paste0("OTC parts don't add up to total ",   print.path))
    }
    
    if (sum(df$`Total Counted (OSC)`,
            na.rm=TRUE)!= sum(df$`N. Calves (OSC)`,
                              df$`N. Adult F (OSC)`,
                              df$`N. Adult M (OSC)`,
                              df$`N. Adult (Sex Unclassified) (OSC)`,
                              df$`Unclassified Life Stage and Sex (OSC)`,
                              df$`Yrling - Unclas Sex (OSC)`,
                              na.rm = TRUE
            )) {
      stop(paste0("OSC parts don't add up to total ",   print.path))
    }
    
    
    if (sum(df$`Total Counted (MNKA)`,
            na.rm=TRUE)!= sum(df$`N. Calves (MNKA)`,
                              df$`N. Adult F (MNKA)`,
                              df$`N. Adult M (MNKA)`,
                              df$`N. Adult (Sex Unclassified) (MNKA)`,
                              df$`Unclassified Life Stage and Sex (MNKA)`,
                              df$`Yrling - Unclas Sex (MNKA)`,
                              na.rm = TRUE
            )) {
      stop(paste0("MNKA parts don't add up to total ",   print.path))
    }



    df <- df %>%
      mutate(
        herd = herd.name,
        Year = year,
        NFG = NA,
        `NFG-Reason` = NA,
        `Estimate (MC, eg JHE)`=NA
      ) %>%
      select(herd,
        Year,
        `Official Survey Timing`,
        NFG,
        `NFG-Reason`,
        `N collar available`,
        `N collars detected`,
        `Survey Count (KMB = SO, Total count) OTC` = `Total Counted (OTC)`,
        `KMB = Observed sample count, OSC; Min count (survey for recruitment - cannot be used for trend analyses in the IPM or otherwise)`=`Total Counted (OSC)`,
        `MinCount (KMB = Minimum # known alive)`=`Total Counted (MNKA)`,
        `N. Calves (OTC)`,
        `Yrling - Unclas Sex (OTC)`,
        `N. Adult F (OTC)`,
        `N. Adult M (OTC)`,
        `N. Adult (Sex Unclassified) (OTC)`,
        `Unclassified Life Stage and Sex (OTC)`,
        `N. Calves (OSC)`,
        `Yrling - Unclas Sex (OSC)`,
        `N. Adult F (OSC)`,
        `N. Adult M (OSC)`,
        `N. Adult (Sex Unclassified) (OSC)`,
        `Unclassified Life Stage and Sex (OSC)`,
        `N. Calves (Min # Known Alive)` = `N. Calves (MNKA)`,
        `Yrling - Unclas Sex (Min # Known Alive)` = `Yrling - Unclas Sex (MNKA)`,
        `N. Adult F (Min # Known Alive)` = `N. Adult F (MNKA)`,
        `N. Adult M (Min # Known Alive)` = `N. Adult M (MNKA)`,
        `N. Adult (Sex Unclassified) (Min # Known Alive)` = `N. Adult (Sex Unclassified) (MNKA)`,
        `Unclassified Life Stage and Sex (Min # Known Alive)` = `Unclassified Life Stage and Sex (MNKA)`,
        `Estimate (MC, eg JHE)`
      )

    return(df)
  }
}

## loop over all data
path <- sheets
count.list <- list()

for (i in 1:length(sheets)) {
  count.list[[i]] <- prep.count(path[i])
}

count.add <- bind_rows(count.list)

write_csv(count.add, here::here("data/raw/count.2022onwards.csv"))



##### prep treatment data ####
prep.trt <- function(path) {
  df <- read_excel(path, sheet = "Treatment")%>%
    filter_all(any_vars(!is.na(.)))
  herd.name <- read_excel(path, sheet = 1)[1, 2] %>%
    as.character()
  print.path <- basename(path)


  ## checks
  year <- str_sub(path, -9, -6) %>%
    as.numeric()

  if (nrow(df) == 0) {
    stop(paste0("No treatment data ",   print.path))
  }

  df <- df %>%
    mutate(
      herd = herd.name,
      start.year = year,
      end.year = year
    ) %>%
    select(herd, treatment = Treatment, start.year, end.year, intensity = Intensity)

  return(df)
}

## loop over all data
path <- sheets
trt.list <- list()

for (i in 1:length(sheets)) {
  trt.list[[i]] <- prep.trt(path[i])
}

trt.add <- bind_rows(trt.list)

write_csv(trt.add, here::here("data/raw/trt.2022onwards.csv"))


## function for renaming
prep.trt <- function(path) {
  df <- read_excel(path, sheet = "Treatment")
  herd.name <- read_excel(path, sheet = 1)[1, 2] %>%
    as.character()

  path[i]
  ## checks
  year <- str_sub(path, -9, -6) %>%
    as.numeric()

  if (nrow(df) == 0) {
    stop("No treatment data")
  }

  df <- df %>%
    mutate(
      herd = herd.name,
      start.year = year,
      end.year = year
    ) %>%
    select(herd, treatment = Treatment, start.year, end.year, intensity = Intensity)%>%
    mutate(treatment=case_when(treatment%in%c("wolf reduction")~"reduce wolves",
                               treatment%in%c("moose reduction")~"reduce moose",
                               treatment%in%c("maternal penning", "mat pen")~"pen",
                               treatment%in%c("supplementary feeding", "feeding")~"feed",
                               TRUE~treatment),
           Exclude=NA)
  
  if (any(!df$treatment%in%c("reduce wolves","reduce moose","pen","feed","transplant", "none"))) { 
    stop(paste0("Treatment has an invalid value ",   print.path))
  }

  return(df)
}

## loop over all data
path <- sheets
trt.list <- list()

for (i in 1:length(sheets)) {
  trt.list[[i]] <- prep.trt(path[i])
}

trt.add <- bind_rows(trt.list)

write_csv(trt.add, here::here("data/raw/trt.2022onwards.csv"))
