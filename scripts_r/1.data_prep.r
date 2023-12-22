## ----render, eval=FALSE,include=FALSE----------------------------------------------------------------------------------------------------------------------------------
## rmarkdown::render(here::here("data", "dataprep.Rmd"),
##   output_file = "README.md"
## )
## 
## knitr::purl(
##   input = here::here("data", "dataprep.Rmd"),
##   output = here::here("scripts_r", "1.data_prep.r")
## )


## ----Load packages and data, results='hide', message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
library(renv)
## to pull packages
# restore(repos="https://cloud.r-project.org")
library(styler)
library(here)
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(scales)
library(hrbrthemes)
library(gt)
library(knitr)
library(tidylog)
set.seed(2023)



#### Year cut off, for trimming of data that is being updated
year.cut <- 2023

count.raw <- read_csv(here::here("data", "raw", "count.csv"))
surv.raw <- read_csv(here::here("data", "raw", "surv.csv"))
treat.raw <- read_csv(here::here("data", "raw", "treatment.csv"))

## AB
ab <- read_csv(here::here("data", "raw", "ab", "Caribou_Demographic_Vital_Rates_22July2021_forSerrouya.csv")) %>%
  mutate(
    Reporting_Pop = case_when(Reporting_Pop %in% "Narraway" ~ "Narraway AB", TRUE ~ Reporting_Pop), ## break out AB data
    Year_START = Year_END - 1
  ) ## switch reporting year to start year


## TWEEDS
tweed <- read_csv(here::here("data", "raw", "TweedsmuirSummary.csv"))

## Tonquin
tonq.surv <- read_csv(here::here("data", "raw", "TonquinSurvivalAnnualMedian(nodatafor2022).csv")) ## posteriors by age-sex in appendix
tonq.abund <- read_csv(here::here("data", "raw", "TonquinAbundanceAnnualEstimates.csv")) ## overall totalN from Layla Nov 16, 2023


## ----filter, results='hide', message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------
herds <- treat.raw %>%
  filter(!Exclude %in% "Y") %>% ## remove herds that don't have enough data/years
  filter(!Herd %in% "Central Selkirks") %>% ## now split into Nakusp/Duncan
  filter(!Herd %in% "Quintette Full") %>% ## pull from here and add in below to QT Full is end of list for easy filtering w/o messing up indexing numbers
  distinct(Herd) %>%
  arrange(Herd) %>%
  pull(Herd)

herds <- c(herds, "Quintette Full") ## add to end of list

treat.raw <- treat.raw %>%
  filter(Herd %in% herds)

count.raw <- count.raw %>%
  rename(herd = Herd) %>%
  filter(herd %in% herds) %>%
  filter(!(herd == "Tweedsmuir" & Year == 1987)) ## only keep one count

surv.raw <- surv.raw %>%
  dplyr::select(-Herd) %>%
  rename(herd = `Herd (RS, CL)`) %>%
  filter(herd %in% herds)


## ----add new data, results='hide', message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------
## update the data
source("/Users/claytonlamb/Dropbox/Documents/University/Work/WSC/CaribouIPM_BCAB/data/prepupdates.R")

unique(surv.raw$Outcome)
unique(surv.raw$`Sex (F, M)`)
unique(surv.raw$Age.when.collard)

## SURVIVAL: format dates, calculate monitoring days, and trim to only adult females
surv <- surv.raw %>%
  mutate(
    DateEntry = ymd(`Date entry`), # get dates dialed
    DateExit = ymd(`Date exit`)
  ) %>%
  select(herd, `Animal ID` = `Animal ID (RS, CL)`, WLHID, Sex = `Sex (F, M)`, Ageclass = Age.when.collard, DateEntry, DateExit, Outcome) %>% ## slim down to only needed columns
  # filter(!(herd%in%"Columbia North" & Comment%in%c("wild","release")))%>% ##to remove penned animals in CN
  rbind(read_csv(here::here("data/raw/surv.2022onwards.csv")) %>%
    mutate(
      DateEntry = ymd(`Date entry`), # get dates dialed
      DateExit = ymd(`Date exit`)
    ) %>%
    select(herd, `Animal ID`, WLHID, Sex, Ageclass, DateEntry, DateExit, Outcome)) %>% ## add in 2022 onwards data
  mutate(
    id = case_when(
      is.na(WLHID) ~ `Animal ID`, ## fix ID's
      TRUE ~ WLHID
    ),
    MonitoringDays = (DateExit - DateEntry) %>% as.numeric()
  ) %>%
  mutate(event = case_when(
    Outcome %in% c("Mortality") ~ 1,
    TRUE ~ 0
  )) %>%
  filter(
    !`Sex` %in% c("M"),
    !Ageclass %in% c("Calf"),
    MonitoringDays > 0
  ) %>%
  dplyr::select(id, herd, DateEntry, DateExit, event)


## COUNTS
count.raw <- count.raw %>%
  select(
    herd,
    Year,
    `Official Survey Timing`,
    NFG,
    `NFG-Reason`,
    `N collar available`,
    `N collars detected`,
    `Survey Count (KMB = SO, Total count) OTC`,
    `KMB = Observed sample count, OSC; Min count (survey for recruitment - cannot be used for trend analyses in the IPM or otherwise)`,
    `MinCount (KMB = Minimum # known alive)`,
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
    `N. Calves (Min # Known Alive)`,
    `Yrling - Unclas Sex (Min # Known Alive)`,
    `N. Adult F (Min # Known Alive)`,
    `N. Adult M (Min # Known Alive)`,
    `N. Adult (Sex Unclassified) (Min # Known Alive)`,
    `Unclassified Life Stage and Sex (Min # Known Alive)`,
    `Estimate (MC, eg JHE)`
  ) %>%
  rbind(read_csv(here::here("data/raw/count.2022onwards.csv")))


## TREATMENTS
# clean up names of treatments
treat.raw <- treat.raw %>%
  mutate(treatment = case_when(
    treatment %in% "sterilization" ~ "sterilize wolves",
    treatment %in% "pen" ~ "pen",
    treatment %in% "feeding" ~ "feed",
    treatment %in% "wolf reduction" ~ "reduce wolves",
    treatment %in% "moose reduction" ~ "reduce moose",
    treatment %in% "moose reduction" ~ "reduce moose",
    TRUE ~ treatment
  )) %>%
  select(herd = Herd, treatment, start.year, end.year, intensity, Exclude) %>%
  rbind(read_csv(here::here("data/raw/trt.2022onwards.csv")))


## ----survival prep, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------------------
## remove some duplicated columns
surv <- surv %>%
  distinct() # none

## fix trailing letter in Chase and Wolverine animal names
surv <- surv %>%
  mutate(id = case_when(
    str_starts(id, "CN") & herd %in% c("Klinse-Za", "Chase", "Wolverine") ~ str_sub(id, 1, -2),
    TRUE ~ id
  ))


## offset dates to reflect census year
## most herds have march census, but Chilco herds (Ithcas, CA, Rainb) has June. So, offset. Parks AB has fall

offset.bioyr <- ymd("2020-April-01") - ymd("2020-Jan-01")
offset.bioyr_spring <- ymd("2020-June-01") - ymd("2021-Jan-01")
spring.herds <- c("Itcha-Ilgachuz", "Charlotte Alplands", "Rainbows")

offset.bioyr_fall <- ymd("2020-Oct-01") - ymd("2021-Jan-01")
fall.herds <- c("Brazeau", "Maligne", "Tonquin")

surv <- surv %>%
  mutate(
    DateEntry.bio = case_when(
      !herd %in% c(spring.herds, fall.herds) ~ DateEntry - offset.bioyr,
      herd %in% c(spring.herds) ~ DateEntry - offset.bioyr_spring,
      herd %in% c(fall.herds) ~ DateEntry - offset.bioyr_fall
    ),
    DateExit.bio = case_when(
      !herd %in% c(spring.herds, fall.herds) ~ DateExit - offset.bioyr,
      herd %in% c(spring.herds) ~ DateExit - offset.bioyr_spring,
      herd %in% c(fall.herds) ~ DateExit - offset.bioyr_fall
    )
  )


## transform into annual survival (first to daily)
source(here::here("other", "seasonalsurvival_fn.r"))
surv.day <- surv %>%
  dplyr::select(id, herd, DateEntry.bio, DateExit.bio, event) %>%
  mutate(id_period = paste0(id, DateEntry.bio)) %>%
  ungroup() %>%
  stretch_survival_data("1 day") %>%
  mutate(year = year(DateEntry.bio))

surv.yr <- surv.day %>%
  group_by(id, herd, year) %>%
  summarise(
    enter.date = min(DateEntry.bio),
    exit.date = max(DateExit.bio) - 1,
    event = max(dead)
  ) %>%
  mutate(
    enter = month(enter.date) - 1,
    exit = month(exit.date),
    time = (exit.date - enter.date) %>% as.numeric()
  )

## save
write_csv(surv.day, here::here("data", "raw", "survival_day_noCNpen.csv"))
write_csv(surv.yr, here::here("data", "raw", "survival_yrly.csv"))


## ----check surv, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------
fit <- survfit(Surv(time, event) ~ herd, data = surv.yr)
ggsurvplot(fit,
  pval = TRUE,
  # risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  risk.table.height = 0.5,
  ggtheme = theme_ipsum(), # Change ggplot2 theme
  ylim = c(0.5, 1),
  xlim = c(0, 360)
)


# summary(survfit(Surv(time, event) ~ herd, data = surv.yr), times = 360)


## ----surv est, fig.height=10, fig.width=10, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------
surv.yr$herd <- as.character(surv.yr$herd)
surv.herds <- unique(surv.yr$herd) %>% as.character()

## extract survival estimates from each surv.fit model
surv.yr.est <- data.frame()
for (i in 1:length(surv.herds)) {
  out <- summary(
    survfit(
      Surv(enter, exit, event) ~ year,
      conf.type = "log-log",
      data = surv.yr %>% dplyr::filter(herd %in% !!surv.herds[i]) %>% dplyr::mutate(year = as.factor(year))
    ),
    times = 12,
    extend = TRUE
  )

  surv.yr.est <- rbind(
    surv.yr.est,
    data.frame(
      herd = surv.herds[i],
      year = str_sub(out$strata, start = 6, end = -1),
      est = out$surv,
      se = out$std.err,
      lower = out$surv - (1.96 * out$std.err),
      upper = out$surv + (1.96 * out$std.err),
      n = out$n,
      # sd=out$std.err*sqrt(out$n.risk+out$n.event+out$n.censor),
      sd = out$std.err, ## consistent with Eacker App
      type = "Frequentist"
    )
  )
}



## add in AB survival
surv.yr.est <- rbind(
  surv.yr.est,
  data.frame(
    herd = ab$Reporting_Pop,
    year = ab$Year_START,
    est = ab$Survival_Mean,
    se = ab$Survival_SD,
    lower = ab$Survival_LCL95,
    upper = ab$Survival_UCL95,
    n = NA,
    sd = ab$Survival_SD,
    type = "Frequentist"
  )
)



## add in Tweeds survival
surv.yr.est <- rbind(
  surv.yr.est,
  data.frame(
    herd = "Tweedsmuir",
    year = tweed$Year,
    est = tweed$`Survival Rate`,
    se = sqrt((tweed$`Survival Rate` * (1 - tweed$`Survival Rate`)) / tweed$N2),
    lower = NA,
    upper = NA,
    n = tweed$N2,
    sd = sqrt((tweed$`Survival Rate` * (1 - tweed$`Survival Rate`)) / tweed$N2),
    type = "Frequentist"
  ) %>%
    drop_na(est)
)

## add in Tonquin
surv.yr.est <- rbind(
  surv.yr.est,
  data.frame(
    herd = "Tonquin",
    year = tonq.surv$Year,
    est = tonq.surv$Mean,
    se = tonq.surv$SD,
    lower = tonq.surv$LCL,
    upper = tonq.surv$UCL,
    n = NA,
    sd = tonq.surv$SD,
    type = "Frequentist"
  ) %>%
    drop_na(est)
)



ggplot(
  surv.yr.est,
  aes(x = as.numeric(year), y = est)
) +
  geom_pointrange(aes(ymin = est - sd, ymax = est + sd), shape = 21) +
  theme_ipsum() +
  ylab("Survival") +
  xlab("Year") +
  facet_wrap(vars(herd), scales = "free_x") +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  labs(x = "Year", title = "KM Survival") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  )

ggsave(here::here("data", "plots", "input_survival1.png"), width = 11, height = 9, bg = "white")



mean(surv.yr.est$est)
mean(surv.yr.est$se, na.rm = TRUE)


## ----surv est low samp, fig.height=10, fig.width=10, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------
## remove herd-years with <=2 animals (means survival est could only be 0,0.5,or 1)
surv.yr.est <- surv.yr.est %>%
  mutate(
    n = case_when(is.na(n) ~ 15, TRUE ~ n), # fill in sample size for AB data, this is approximate, but doesn't matter much, rarely surv==1 or 0 w/ no error
    est = case_when(
      n <= 2 ~ NA_real_,
      TRUE ~ est
    )
  ) %>%
  drop_na(est)


## add in bootstrapped error to deal with no error when survival ==1 or 0
surv.rep <- surv.yr.est %>%
  filter(est %in% c(1, 0)) %>%
  group_by(herd, year) %>%
  slice(rep(1:n(), each = n)) %>%
  mutate(surv = case_when(row_number() == 1 ~ 0, TRUE ~ 1))

## BOOT
surv.boot <- tibble()
for (i in 1:1000) {
  surv.boot.i <- surv.rep %>%
    dplyr::group_by(herd, year) %>%
    dplyr::sample_frac(1, replace = TRUE) %>%
    dplyr::summarise(
      mean = mean(surv),
      n = max(n),
      .groups = "keep"
    ) %>%
    dplyr::mutate(iter = i)

  surv.boot <- dplyr::bind_rows(surv.boot, surv.boot.i)
}

## summarise
surv.boot.summary <- surv.boot %>%
  group_by(herd, year) %>%
  summarise(
    av = mean(mean),
    sd.boot = sd(mean),
    n = mean(n)
  )

## add in to data
surv.yr.est <- surv.yr.est %>%
  left_join(surv.boot.summary %>% dplyr::select(-n) %>% ungroup(), by = c("herd", "year")) %>%
  mutate(
    se = case_when(
      est %in% c(1, 0) ~ sd.boot,
      TRUE ~ se
    ),
    est = case_when(
      est %in% c(0) ~ 0.5,
      TRUE ~ est
    )
  ) %>%
  mutate(
    lower = est - (1.96 * se),
    upper = est + (1.96 * se)
  ) %>%
  mutate(
    sd = se,
    year = as.numeric(year)
  ) %>%
  ungroup() %>%
  dplyr::select(-av, -sd.boot)


## force upper and lowers between 0-1 (has to be true for survival)
surv.yr.est <- surv.yr.est %>%
  mutate(
    lower = case_when(
      lower < 0 ~ 0,
      TRUE ~ lower
    ),
    upper = case_when(
      upper > 1 ~ 1,
      TRUE ~ upper
    )
  )

## plot
ggplot(
  surv.yr.est,
  aes(x = year, y = est)
) +
  geom_pointrange(aes(ymin = est - sd, ymax = est + sd), shape = 21) +
  theme_ipsum() +
  ylab("Survival") +
  xlab("Year") +
  facet_wrap(vars(herd), scales = "free_x") +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  labs(x = "Year", title = "KM Survival") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  )

ggsave(here::here("data", "plots", "input_survival2_fixsmallSS.png"), width = 11, height = 9, bg = "white")



## save
write_csv(surv.yr.est %>% dplyr::select(herd, year, est, sd), here::here("data", "clean", "survival.csv"))


## ----sex ratio, fig.height=10, fig.width=12, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
## extract sex ratio data from counts
sr <- count.raw %>%
  filter(!NFG %in% c("Y")) %>%
  dplyr::select(herd,
    year = Year,
    adF = `N. Adult F (OTC)`,
    adM = `N. Adult M (OTC)`,
    calf = `N. Calves (OTC)`,
    total = `Survey Count (KMB = SO, Total count) OTC`,
    `Official Survey Timing`
  ) %>%
  mutate(type = "OTC") %>%
  rbind(count.raw %>%
    filter(!NFG %in% c("Y")) %>%
    dplyr::select(herd,
      year = Year,
      adF = `N. Adult F (OSC)`,
      adM = `N. Adult M (OSC)`,
      calf = `N. Calves (OSC)`,
      total = `KMB = Observed sample count, OSC; Min count (survey for recruitment - cannot be used for trend analyses in the IPM or otherwise)`,
      `Official Survey Timing`
    ) %>%
    mutate(type = "OSC")) %>%
  rbind(count.raw %>%
    filter(!NFG %in% c("Y")) %>%
    dplyr::select(herd,
      year = Year,
      adF = `N. Adult F (Min # Known Alive)`,
      adM = `N. Adult M (Min # Known Alive)`,
      calf = `N. Calves (Min # Known Alive)`,
      total = `MinCount (KMB = Minimum # known alive)`,
      `Official Survey Timing`
    ) %>%
    mutate(type = "MNKA"))


## denote seasons
sr <- sr %>%
  mutate(season = case_when(
    `Official Survey Timing` %in% c("Winter", "March", "April", "February") ~ "winter",
    `Official Survey Timing` %in% c("August", "September", "October", "November") ~ "fall",
    `Official Survey Timing` %in% c("June", "July") ~ "spring",
    is.na(`Official Survey Timing`) ~ NA_character_,
    TRUE ~ "needs class"
  )) %>%
  filter(!season %in% c("needs class")) %>%
  dplyr::select(-`Official Survey Timing`)




### Filter to only good  counts for identifying sex ratio
sr <- sr %>%
  mutate(
    ad.total = (total - calf),
    ad.id = (adF + adM) / ad.total
  ) %>%
  filter(ad.id > 0.8, (adF + adM) > 20) ## keep only data where sex was identified in >80% of adults during survey, and at least 20 animals seen.


## get average sex ratio
sr <- sr %>%
  mutate(
    sratio = adF / (adF + adM),
    sratio.sd = sqrt(sratio * (1 - sratio) / (adF + adM))
  )

## summarize
sr.summary <- sr %>%
  filter(season %in% c("fall", "winter")) %>% ## keep only fall and winter surveys when sexes are generally better mixed
  group_by(type) %>%
  summarise(
    mean = mean(sratio) %>% round(2),
    sd = sd(sratio) %>% round(2),
    n = n()
  )

ggplot(data = sr, aes(x = sratio, fill = type)) +
  geom_density(alpha = 0.3) +
  theme_ipsum()


ggplot(data = sr, aes(x = sratio, color = type, y = herd)) +
  geom_point(alpha = 0.5)


ggplot(data = sr %>% group_by(herd, type) %>% summarise(sratio = mean(sratio)), aes(x = sratio, color = type, y = herd, group = herd)) +
  geom_path(color = "grey") +
  geom_point(alpha = 0.8) +
  theme_ipsum()

### sr for OTC only (OSC and MNKA female-biased due to flying to Female-only collars)
sr.otc <- sr %>%
  filter(season %in% c("fall", "winter")) %>% ## keep only fall and winter surveys when sexes are generally better mixed
  filter(type %in% "OTC")

## add missing KZ Sex Ratio
kz.sr <- tibble(
  herd = "Klinse-Za",
  year = c(2014:2019),
  adF = NA,
  adM = NA,
  calf = NA,
  total = NA,
  type = "OTC",
  season = "winter",
  ad.total = NA,
  ad.id = NA,
  sratio = c(0.68, 0.71, 0.64, 0.65, 0.51, 0.55),
  sratio.sd = 0.05
)

## export raw sr data
write_csv(sr.otc %>% rbind(kz.sr), here::here("data", "clean", "sexratios.csv"))


## raw distribution to characterize priors
sr.median <- median(sr.otc %>% rbind(kz.sr) %>% pull(sratio))
sr.sd <- sd(sr.otc %>% rbind(kz.sr) %>% pull(sratio))

## summarize average for priors
write_csv(tibble(sr = sr.median, sr.sd = sr.sd), here::here("data", "clean", "sexratio_summary.csv"))


## ----recruitment, fig.height=10, fig.width=10, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------
## extract all 3 types of counts and remove unknown age+sex from the totals (assumes these animals have same calf/cow ratio as pop. Likely if these are just whole groups in the trees)
##* sorry that some of these column names are terribly long, I shorten them here**
recruitment <- count.raw %>%
  mutate(
    `Unclassified Life Stage and Sex (OTC)` = case_when(
      is.na(`Unclassified Life Stage and Sex (OTC)`) ~ 0, ### Make the NA's in unclassified ==0
      TRUE ~ `Unclassified Life Stage and Sex (OTC)`
    ),
    `Unclassified Life Stage and Sex (OSC)` = case_when(
      is.na(`Unclassified Life Stage and Sex (OSC)`) ~ 0,
      TRUE ~ `Unclassified Life Stage and Sex (OSC)`
    ),
    `Unclassified Life Stage and Sex (Min # Known Alive)` = case_when(
      is.na(`Unclassified Life Stage and Sex (Min # Known Alive)`) ~ 0,
      TRUE ~ `Unclassified Life Stage and Sex (Min # Known Alive)`
    )
  ) %>%
  mutate(
    `N. Adults (OTC)` = `Survey Count (KMB = SO, Total count) OTC` - `N. Calves (OTC)` - `Unclassified Life Stage and Sex (OTC)`, ### get to adults only
    `N. Adults (OSC)` = `KMB = Observed sample count, OSC; Min count (survey for recruitment - cannot be used for trend analyses in the IPM or otherwise)` - `N. Calves (OSC)` - `Unclassified Life Stage and Sex (OSC)`,
    `N. Adults (MNA)` = `MinCount (KMB = Minimum # known alive)` - `N. Calves (Min # Known Alive)` - `Unclassified Life Stage and Sex (Min # Known Alive)`
  ) %>%
  mutate(
    r.OTC = `N. Calves (OTC)` / `N. Adults (OTC)`,
    r.OSC = `N. Calves (OSC)` / `N. Adults (OSC)`,
    r.MNA = `N. Calves (Min # Known Alive)` / `N. Adults (MNA)`
  ) %>%
  ungroup() %>%
  mutate(season = case_when(
    `Official Survey Timing` %in% c("Winter", "March", "April", "February", "January") ~ "winter",
    `Official Survey Timing` %in% c("August", "September", "October", "November") ~ "fall",
    `Official Survey Timing` %in% c("June", "July", "May") ~ "spring",
    is.na(`Official Survey Timing`) ~ NA_character_,
    TRUE ~ "needs class"
  )) %>%
  dplyr::select(herd,
    year = Year, season, r.OTC, `N. Calves (OTC)`, `N. Adults (OTC)`,
    r.OSC, `N. Calves (OSC)` = `N. Calves (OSC)`, `N. Adults (OSC)`,
    r.MNA, `N. Calves (MNA)` = `N. Calves (Min # Known Alive)`, `N. Adults (MNA)`
  ) %>%
  drop_na(season)

if (any(recruitment$season == "needs class")) {
  stop("Missing seasons")
}

## remove one record, Itcha's 1978, May survey
recruitment <- recruitment %>% filter(!(herd == "Itcha-Ilgachuz" & year == 1978 & `N. Adults (OSC)` == 1 & is.na(r.OTC)))

## have a look at OTC vs MNA when they are both present
recruitment %>%
  filter(season == "winter") %>%
  group_by(herd, year) %>%
  mutate(count = sum(!is.na(r.OTC), !is.na(r.MNA))) %>%
  filter(count > 1, r.MNA <= 1) %>%
  ggplot() +
  geom_point(aes(x = r.OTC, y = r.MNA, size = `N. Adults (OTC)`)) +
  theme_ipsum() +
  geom_abline(slope = 1, intercept = 0)


## have a look at OTC vs OSC when they are both present
recruitment %>%
  filter(season == "winter") %>%
  group_by(herd, year) %>%
  mutate(count = sum(!is.na(r.OTC), !is.na(r.OSC))) %>%
  filter(count > 1) %>%
  ggplot() +
  geom_point(aes(x = r.OTC, y = r.OSC, size = `N. Adults (OTC)`)) +
  theme_ipsum() +
  geom_abline(slope = 1, intercept = 0)



## reorganize recruitment dataframe
recruitment <- tibble(
  herd = recruitment$herd,
  year = recruitment$year,
  season = recruitment$season,
  type = "OTC",
  calves = recruitment$`N. Calves (OTC)`,
  adults = recruitment$`N. Adults (OTC)`,
  recruitment = calves / adults
) %>%
  rbind(
    tibble(
      herd = recruitment$herd,
      year = recruitment$year,
      season = recruitment$season,
      type = "OSC",
      calves = recruitment$`N. Calves (OSC)`,
      adults = recruitment$`N. Adults (OSC)`,
      recruitment = calves / adults
    )
  ) %>%
  rbind(
    tibble(
      herd = recruitment$herd,
      year = recruitment$year,
      season = recruitment$season,
      type = "MNKA",
      calves = recruitment$`N. Calves (MNA)`,
      adults = recruitment$`N. Adults (MNA)`,
      recruitment = calves / adults
    )
  )



## choose which recruitment to use, hierarchy is MNKA>OTC>OSC

## make sure this works, test an example
# tibble(calves_MNKA=c(0,1,2,3),
#        calves_OTC=c(10,11,12,13),
#        calves_OSC=c(110,111,112,113),
#        MNKA=c(1,2,NA,NA),
#        OTC=c(NA,2,2,NA),
#        OSC=c(1,NA,NA,1))%>%
#    mutate(calves2=case_when(!is.na(MNKA)~calves_MNKA,
#                           !is.na(OTC)~calves_OTC,
#                           !is.na(OSC)~calves_OSC))

## find duplicates
# recruitment%>%dplyr::group_by(herd, year, season, type) %>%
#     dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#     dplyr::filter(n > 1L)

recruitment <- recruitment %>%
  ungroup() %>%
  tidyr::pivot_wider(id_cols = c(herd, year, season), names_from = type, values_from = c(calves, adults, recruitment)) %>%
  mutate(
    calves = case_when(
      !is.na(recruitment_OSC) ~ calves_OSC,
      !is.na(recruitment_MNKA) ~ calves_MNKA,
      !is.na(recruitment_OTC) ~ calves_OTC
    ),
    adults = case_when(
      !is.na(recruitment_OSC) ~ adults_OSC,
      !is.na(recruitment_MNKA) ~ adults_MNKA,
      !is.na(recruitment_OTC) ~ adults_OTC
    ),
    recruitment = case_when(
      !is.na(recruitment_OSC) ~ recruitment_OSC,
      !is.na(recruitment_MNKA) ~ recruitment_MNKA,
      !is.na(recruitment_OTC) ~ recruitment_OTC
    ),
    count.used = case_when(
      !is.na(recruitment_OSC) ~ "OSC",
      !is.na(recruitment_MNKA) ~ "MNKA",
      !is.na(recruitment_OTC) ~ "OTC"
    )
  ) %>%
  dplyr::select(herd, year, season, calves, adults, recruitment, count.used) %>%
  drop_na(recruitment)

## summary stats by type
recruitment %>%
  filter(!recruitment %in% c(NA, Inf)) %>%
  group_by(count.used) %>%
  summarise(
    mean = mean(recruitment) %>% round(2),
    sd = sd(recruitment) %>% round(2),
    max = max(recruitment) %>% round(2),
    n = n()
  ) %>%
  kable()



## remove herd-years with <=2 adults seen (means r est could only be 0,0.5,or 1)
recruitment <- recruitment %>%
  mutate(adults = case_when(
    adults <= 2 ~ NA_real_,
    TRUE ~ adults
  )) %>%
  drop_na(adults)



## estimate R error
recruitment <- recruitment %>%
  group_by(herd, year) %>%
  mutate(
    se = sqrt((recruitment * (1 - recruitment)) / sum(calves + adults)),
    # sd=se*sqrt(sum(calves+adults)),
    sd = se,
    lower = recruitment - (se * 1.96),
    upper = recruitment + (se * 1.96)
  ) %>%
  drop_na(recruitment) %>%
  dplyr::select(herd, year, season, est = recruitment, se, lower, upper, sd, calves, adults, count.used)





## add in Tweeds recruitment
recruitment <- rbind(
  recruitment,
  data.frame(
    herd = "Tweedsmuir",
    year = tweed$Year + 1,
    season = "winter",
    est = tweed$R,
    se = sqrt((tweed$R * (1 - tweed$R)) / tweed$N4),
    lower = NA,
    upper = NA,
    sd = sqrt((tweed$R * (1 - tweed$R)) / tweed$N4),
    calves = round(tweed$N4 * tweed$R, 0),
    adults = tweed$N4,
    count.used = "OSC"
  ) %>%
    drop_na(est)
) %>%
  drop_na(est)


## Shift years for surveys not on march timing, these refer to a different cohort of calves
recruitment <- recruitment %>%
  mutate(year = case_when(
    season %in% c("spring", "fall") ~ year + 1,
    TRUE ~ year
  ))

## Prioritize surveys in a given year based on winter>spring>fall
recruitment <- recruitment %>%
  mutate(season_order = case_when(
    season %in% "winter" ~ 1,
    season %in% "fall" ~ 3,
    season %in% "spring" ~ 2
  )) %>%
  group_by(herd, year) %>%
  filter(season_order <= min(season_order)) %>%
  mutate(season_int = case_when(
    season %in% "winter" ~ 1,
    season %in% "fall" ~ 2,
    season %in% "spring" ~ 3
  ))


## offset OSC and MNKA b/c sex ratios are higher. Use observed sex ratios when possible.
otc.sr <- sr.summary %>%
  filter(type %in% "OTC") %>%
  pull(mean)
osc.sr <- sr.summary %>%
  filter(type %in% "OSC") %>%
  pull(mean)
mnka.sr <- sr.summary %>%
  filter(type %in% "MNKA") %>%
  pull(mean)

recruitment <- recruitment %>%
  left_join(sr %>% filter(type %in% c("OSC", "MNKA")) %>% dplyr::select(herd, year, season, sratio, type), by = c("herd", "season", "year", "count.used" = "type")) %>%
  mutate(
    sratio = case_when(
      count.used %in% "OSC" & is.na(sratio) ~ osc.sr,
      count.used %in% "MNKA" & is.na(sratio) ~ mnka.sr,
      TRUE ~ sratio
    ),
    est.adj = case_when(
      count.used %in% c("OSC", "MNKA") ~ est * (otc.sr / sratio),
      TRUE ~ est
    )
  )



write_csv(recruitment %>% dplyr::select(herd, year, est = est.adj, sd, season_int), here::here("data", "clean", "recruitment.csv"))
write_csv(recruitment %>% dplyr::select(herd, year, season, type = count.used, calves, adults, recruitment = est), here::here("data", "clean", "recruitment_counts.csv"))


ggplot(recruitment, aes(x = year, y = est, color = season)) +
  geom_pointrange(aes(ymin = est - sd, ymax = est + sd), shape = 21) +
  theme_ipsum() +
  theme(legend.position = "bottom") +
  ylab("Calves (proportion)") +
  xlab("Year") +
  facet_wrap(vars(herd), scales = "free_x") +
  labs(x = "Year", title = "Recruitment") +
  expand_limits(y = 0) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  )

ggsave(here::here("data", "plots", "input_recruitment.png"), width = 15, height = 12, bg = "white")


ggplot(recruitment, aes(x = est, y = est.adj, color = count.used)) +
  geom_point(alpha = 0.5) +
  theme_ipsum() +
  theme(legend.position = "bottom") +
  labs(title = "Recruitment") +
  expand_limits(y = 0) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  )


## ----counts, fig.height=10, fig.width=12, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------
## clean up counts
counts <- count.raw %>%
  filter(!NFG %in% c("Y")) %>%
  mutate(
    `N collars detected` = case_when(is.na(`Survey Count (KMB = SO, Total count) OTC`) ~ NA_real_, `N collar available` == 0 ~ NA_real_, TRUE ~ `N collars detected`),
    `N collar available` = case_when(is.na(`Survey Count (KMB = SO, Total count) OTC`) ~ NA_real_, `N collar available` == 0 ~ NA_real_, TRUE ~ `N collar available`)
  ) %>%
  mutate(`Survey Count (KMB = SO, Total count) OTC` = case_when(
    herd %in% c("Chase", "Wolverine") & !is.na(`Estimate (MC, eg JHE)`) ~ (`Estimate (MC, eg JHE)` * (`N collars detected` / `N collar available`)),
    TRUE ~ `Survey Count (KMB = SO, Total count) OTC`
  )) %>% ## back transform Chase and Wolverine
  dplyr::select(herd,
    year = Year,
    timing = `Official Survey Timing`,
    SurveyCount = `Survey Count (KMB = SO, Total count) OTC`,
    MinCount = `MinCount (KMB = Minimum # known alive)`,
    CollarsSeen = `N collars detected`,
    CollarsAvailable = `N collar available`
  ) %>%
  mutate(
    Sightability = CollarsSeen / CollarsAvailable,
    se = sqrt(Sightability * (1 - Sightability) / CollarsAvailable),
    sd = se
  )

## sightability by herd
herd.sight <- counts %>%
  group_by(herd) %>%
  summarise(
    est = sum(CollarsSeen, na.rm = TRUE) / sum(CollarsAvailable, na.rm = TRUE),
    se = sqrt(est * (1 - est) / mean(CollarsAvailable, na.rm = TRUE)),
    sd = se,
    n = sum(!is.na(CollarsSeen))
  )


## to get rough sight-corrected sightability (Est_CL) add in mean sightability where there is none, or if 0,
counts <- counts %>%
  group_by(herd) %>% ## by herd
  mutate(
    Sightability.pooled = case_when(
      is.na(Sightability) | Sightability %in% 0 ~ mean(Sightability, na.rm = TRUE),
      TRUE ~ Sightability
    ),
    Est_CL = case_when(
      Sightability.pooled > 0 ~ SurveyCount / Sightability.pooled,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup() %>% ## across herds if no sightability at all for that herd
  mutate(
    Sightability.pooled = case_when(
      is.na(Sightability) | Sightability %in% 0 ~ mean(Sightability, na.rm = TRUE),
      TRUE ~ Sightability
    ),
    Est_CL = case_when(
      Sightability.pooled > 0 ~ SurveyCount / Sightability.pooled,
      TRUE ~ NA_real_
    )
  )



## if estimate is <min count, use min count but force sightability to 0.99
counts <- counts %>%
  mutate(
    count = case_when(
      Est_CL >= MinCount ~ SurveyCount,
      Est_CL < MinCount ~ MinCount,
      is.na(MinCount) ~ SurveyCount,
      is.na(SurveyCount) ~ MinCount
    ),
    Sightability = case_when(
      is.na(MinCount) ~ Sightability,
      Est_CL >= MinCount ~ Sightability,
      Est_CL < MinCount ~ 0.99,
      is.na(SurveyCount) ~ 0.99
    ),
    MinUsed = case_when(
      Est_CL < MinCount ~ 1,
      is.na(SurveyCount) ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(sd = case_when(
    Sightability == 0.99 ~ 0.15, # add error around min counts. had mean(sd, na.rm = TRUE) which was ~0.08, but more error so counts are less precise is good.
    TRUE ~ sd
  ))


## add in bootstrapped error to deal with when sightability ==1 or 0 (no error)
sight.rep <- counts %>%
  filter(Sightability %in% c(1, 0)) %>%
  dplyr::select(herd, year, CollarsSeen, CollarsAvailable) %>%
  group_by(herd, year) %>%
  slice(rep(1:n(), each = CollarsAvailable)) %>%
  mutate(sight = case_when(row_number() == 1 ~ 0, TRUE ~ 1))

## BOOT
sight.boot <- tibble()
for (i in 1:1000) {
  sight.boot.i <- sight.rep %>%
    dplyr::group_by(herd, year) %>%
    dplyr::sample_frac(1, replace = TRUE) %>%
    dplyr::summarise(mean = mean(sight), collars = max(CollarsAvailable), .groups = "keep") %>%
    dplyr::mutate(iter = i)

  sight.boot <- bind_rows(sight.boot, sight.boot.i)
}

## summarise
sight.boot.summary <- sight.boot %>%
  group_by(herd, year) %>%
  summarise(
    av = mean(mean),
    sd.boot = sd(mean),
    n = mean(collars)
  )


## add in
counts <- counts %>%
  left_join(sight.boot.summary %>% ungroup(), by = c("herd", "year")) %>%
  mutate(sd = case_when(
    Sightability %in% c(1, 0) ~ sd.boot,
    TRUE ~ sd
  )) %>%
  ungroup() %>%
  dplyr::select(-av, -sd.boot, -n)


## create final count data
counts <- counts %>%
  ungroup() %>%
  mutate(
    Est_CL = case_when(
      is.na(Sightability) | Sightability %in% 0 ~ count / Sightability.pooled,
      !is.na(Sightability) ~ count / Sightability
    ),
    Est_CL.min = case_when(
      is.na(Sightability) | Sightability %in% 0 ~ count / (Sightability.pooled + (1.64485 * mean(se, na.rm = TRUE))),
      !is.na(Sightability) & (Sightability + (1.64485 * sd)) < 1 & Sightability != 0 ~ count / (Sightability + (1.64485 * sd)),
      !is.na(Sightability) & (Sightability + (1.64485 * sd)) >= 1 & Sightability != 0 ~ count / 1
    ),
    Est_CL.max = case_when(
      is.na(Sightability) | Sightability %in% 0 ~ count / (Sightability.pooled - (1.64485 * mean(se, na.rm = TRUE))),
      !is.na(Sightability) & (Sightability - (1.64485 * sd)) > 0 & Sightability != 0 ~ count / (Sightability - (1.64485 * sd)),
      !is.na(Sightability) & (Sightability - (1.64485 * sd)) <= 0 & Sightability != 0 ~ count / (Sightability * 0.75)
    )
  ) %>%
  drop_na(count)

## Add AB DNA counts, add error in via imperfect sightability. These are single counts for each herd, so dont impact trend, just anchor it.
counts <- counts %>%
  rbind(
    tibble(
      herd = c("A La Peche", "Narraway AB", "Redrock/Prairie Creek"),
      year = c(2018, 2019, 2019),
      timing = "Winter",
      SurveyCount = c(152, 56, 153) * 0.9,
      MinCount = NA,
      CollarsSeen = NA,
      CollarsAvailable = NA,
      Sightability = rep(0.9, times = 3),
      se = NA,
      sd = rep(0.05, times = 3),
      Sightability.pooled = NA,
      Est_CL = c(152, 56, 153),
      count = c(152, 56, 153),
      MinUsed = 0,
      Est_CL.min = c(152, 56, 153) * .95,
      Est_CL.max = c(152, 56, 153) * 1.05
    ) %>% ungroup()
  )

## add AB sight
herd.sight <- herd.sight %>%
  filter(!herd %in% c("A La Peche", "Narraway AB", "Redrock/Prairie Creek")) %>%
  rbind(
    tibble(
      herd = c("A La Peche", "Narraway AB", "Redrock/Prairie Creek"),
      est = rep(0.9, times = 3),
      se = NA,
      sd = rep(0.05, times = 3),
      n = 1
    )
  )



# Add Tonquin abundance
counts <- counts %>%
  rbind(tonq.abund %>%
    mutate(
      SurveyCount = Mean * 0.9,
      Sightability = 0.9
    ) %>%
    mutate(
      herd = c("Tonquin"),
      year = Year,
      timing = "October",
      MinCount = NA,
      CollarsSeen = NA,
      CollarsAvailable = NA,
      Sightability = Sightability,
      se = NA,
      sd = (SD * Sightability^2) / (SD * Sightability + SurveyCount),
      Sightability.pooled = NA,
      Est_CL = SurveyCount / Sightability,
      count = SurveyCount / Sightability,
      MinUsed = 0,
      Est_CL.min = (SurveyCount / (Sightability + sd)),
      Est_CL.max = (SurveyCount / (Sightability - sd))
    ) %>%
    dplyr::select(colnames(counts)))


## plot
ggplot(
  counts %>%
    mutate(year = as.integer(year)) %>%
    pivot_longer(SurveyCount:Est_CL) %>%
    filter(name %in% c("SurveyCount", "MinCount", "Estimate", "Est_CL")),
  aes(x = year, y = value, color = name)
) +
  geom_point() +
  theme_ipsum() +
  ylab("Abundance") +
  xlab("Year") +
  facet_wrap(vars(herd), scales = "free") +
  labs(x = "Year", title = "Population Trajectory") +
  expand_limits(y = 0) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  )

ggplot(
  counts,
  aes(x = year, y = Est_CL)
) +
  geom_point() +
  theme_ipsum() +
  ylab("Abundance") +
  xlab("Year") +
  facet_wrap(vars(herd), scales = "free") +
  labs(x = "Year", title = "Population Trajectory") +
  expand_limits(y = 0) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  )
ggsave(here::here("data", "plots", "abundance1.png"), width = 12, height = 9, bg = "white")


ggplot(
  counts,
  aes(x = year, y = Est_CL, shape = as.factor(MinUsed))
) +
  geom_point() +
  theme_ipsum() +
  ylab("Abundance") +
  xlab("Year") +
  facet_wrap(vars(herd), scales = "free") +
  labs(x = "Year", title = "Population Trajectory", shape = "Minimum count?") +
  expand_limits(y = 0) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  )

ggsave(here::here("data", "plots", "abundance2.png"), width = 12, height = 10, bg = "white")



write_csv(counts %>% dplyr::select(herd, year, count, Sightability, sd, MinUsed, Est_CL, Est_CL.min, Est_CL.max), here::here("data", "clean", "counts.csv"))

write_csv(herd.sight %>% drop_na(est), here::here("data", "clean", "herd_sightability.csv"))


counts %>%
  filter(CollarsAvailable > 0) %>%
  group_by(herd) %>%
  summarise(
    available = sum(CollarsAvailable),
    seen = sum(CollarsSeen)
  ) %>%
  mutate(sightability = seen / available) %>%
  write_csv(here::here("data", "clean", "herd_sightability_summary.csv"))


## ----treatment, fig.height=8, fig.width=12, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------
trt_raw <- treat.raw %>%
  filter(
    herd %in% c(counts$herd, surv$herd, recruitment$herd),
    !treatment %in% "all"
  ) %>%
  group_by(herd, treatment) %>%
  mutate(
    start.year = min(start.year),
    end.year = max(end.year)
  ) %>%
  slice(1) %>%
  as.data.frame()


# shift start year of treatments +1 to reflect lag in detecting change
trt_raw <- trt_raw %>%
  group_by(herd) %>%
  mutate(trt = n_distinct(treatment)) %>%
  ungroup() %>%
  mutate(
    start.year = case_when(
      !treatment %in% "none" ~ start.year + 1,
      TRUE ~ start.year
    ),
    end.year = case_when(
      !treatment %in% "none" ~ end.year + 1,
      TRUE ~ end.year
    )
  )


# Herds and herd number
hn <- tibble(
  herd = herds,
  herd_num = seq(1:length(herds))
)
nherd <- nrow(hn)

#  Years of study
yrs <- seq(from = min(trt_raw$start.year, na.rm = TRUE), to = max(trt_raw$end.year, na.rm = TRUE), by = 1)
nyr <- length(yrs)

# Possible treatments (excluding "none")
trt_only <- trt_raw %>%
  dplyr::filter(treatment != "none") %>%
  left_join(hn, by = "herd")
u_trts <- unique(trt_only$treatment)
ntrts <- length(u_trts)
possible_trts <- data.frame(
  year = rep(yrs, times = ntrts + 1),
  treatment = rep(c("none", u_trts), each = nyr),
  applied = 0
)
# Generate treatment data per herd
herd_trt_ls <- list()
for (i in 1:nherd) {
  possible_trts_h <- data.frame(herd = rep(herds[i],
    times = nrow(possible_trts)
  )) %>%
    cbind(possible_trts)
  tmp_dat <- trt_raw %>%
    dplyr::filter(herd == herds[i])
  u_trts_h <- unique(tmp_dat$treatment)

  none_st <- trt_raw %>%
    dplyr::filter(herd == herds[i]) %>%
    dplyr::filter(treatment == "none")
  none_st <- none_st$start.year

  first_st <- tmp_dat %>%
    arrange(start.year) %>%
    slice(1)
  first_st <- first_st$start.year

  if (sum(none_st) == 0) {
    none_st <- first_st ## for QT full, which had no "none"
  }

  trt_ls <- list()
  for (j in 1:length(u_trts_h)) {
    trt_h <- data.frame(
      herd = rep(herds[i], times = nyr),
      year = yrs,
      treatment = rep(u_trts_h[j], times = nyr),
      applied = NA
    )

    trt_dur <- tmp_dat %>%
      dplyr::filter(treatment == u_trts_h[j])
    trt_st <- trt_dur$start.year
    trt_ed <- trt_dur$end.year

    trt_h <- trt_h %>%
      dplyr::mutate(applied = ifelse(year >= first_st, 0, NA)) %>%
      dplyr::mutate(applied = ifelse(year >= none_st, 0, applied)) %>%
      dplyr::mutate(applied = ifelse(year >= trt_st & year <= trt_ed, 1, applied))

    trt_ls[[j]] <- trt_h
  }

  all_trt_h <- do.call(rbind, trt_ls)
  missing_trts <- possible_trts_h %>%
    dplyr::anti_join(all_trt_h, by = c("herd", "year", "treatment")) %>%
    dplyr::mutate(applied = ifelse(year < none_st, NA, applied))
  all_trt_h <- all_trt_h %>%
    rbind(missing_trts)
  herd_trt_ls[[i]] <- all_trt_h
}
trt_long <- do.call(rbind, herd_trt_ls) %>%
  dplyr::filter(treatment != "none")
# 0 = treatment did not occur for that herd in that year
# 1 = treatment did occur
# NA = herd not monitored

## add intensity
trt_int <- treat.raw %>%
  filter(intensity == "low")

trt_intensity.df <- tibble()
for (i in 1:nrow(trt_int)) {
  a <- tibble(
    herd = trt_int[i, ]$herd,
    year = c((trt_int[i, ]$start.year + 1):(trt_int[i, ]$end.year + 1)),
    treatment = trt_int[i, ]$treatment,
    intensity = trt_int[i, ]$intensity
  )

  trt_intensity.df <- bind_rows(trt_intensity.df, a)
}


trt_long <- trt_long %>%
  left_join(trt_intensity.df, by = c("herd", "year", "treatment"))

## remove years that exceed available demographic data (happens due to +1 above)
trt_long <- trt_long %>%
  filter(year <= year.cut)

write_csv(trt_long, here::here("data", "clean", "treatments.csv"))

trt_long %>%
  mutate(applied = as.factor(applied)) %>%
  drop_na(applied) %>%
  ggplot(aes(x = year, y = herd, color = applied)) +
  geom_point() +
  facet_wrap(vars(treatment)) +
  labs(x = "Year", y = "Herd", title = "Treatments Applied") +
  expand_limits(y = 0) +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  )


## ----assess herds, results='hide', message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------
herds.present <- treat.raw %>%
  filter(!Exclude %in% "Y") %>%
  distinct(herd) %>%
  left_join(counts %>%
    distinct(herd) %>%
    dplyr::select(herd) %>%
    mutate(count = 1)) %>%
  left_join(recruitment %>%
    ungroup() %>%
    distinct(herd) %>%
    dplyr::select(herd) %>%
    mutate(recruitment = 1)) %>%
  left_join(surv.yr.est %>%
    distinct(herd) %>%
    dplyr::select(herd) %>%
    mutate(survival = 1)) %>%
  left_join(trt_long %>%
    distinct(herd) %>%
    dplyr::select(herd) %>%
    mutate(treatment = 1))

write_csv(herds.present, here::here("data", "herds.present.csv"))


herds.present %>%
  pivot_longer(count:treatment) %>%
  replace(is.na(.), 0) %>%
  ggplot(aes(x = name, y = herd, color = as.factor(value))) +
  geom_point() +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13)
  ) +
  labs(color = "Present?")


## ----bp, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------------------------------------
# sight_group <- tibble(herd=c("Columbia North", "Columbia South", "Central Rockies", "Frisby-Boulder", "Nakusp", "Duncan",
#                               "Burnt Pine",  "Kennedy Siding", "Quintette", "Klinse-Za","Narraway BC",
#                              "Purcells Central", "Purcells South", "South Selkirks","Central Selkirks","Monashee",
#                              "Narraway AB", "Redrock/Prairie Creek", "A La Peche",
#                               "Graham",
#                              "Barkerville","Chase","George Mountain","Hart North","Hart South","Narrow Lake","North Cariboo","Wells Gray North","Wells Gray South",
#                              "Charlotte Alplands","Itcha-Ilgachuz","Rainbows","Tweedsmuir"),
#                       sight_grp=c(rep(1,times=6),
#                                   rep(2,times=5),
#                                   rep(3,times=5),
#                                   rep(4,times=3),
#                                   rep(5,times=1),
#                                   rep(6,times=9),
#                                   rep(7,times=4)))

sight_group <- herd.sight %>%
  mutate(sight_grp = case_when(
    est <= 0.5 ~ 1,
    est > 0.5 & est <= 0.85 ~ 2,
    est > 0.85 & est <= 1 ~ 3,
    is.na(est) ~ 2
  )) %>%
  dplyr::select(herd, sight_grp) %>%
  distinct()

bp <- hn %>%
  left_join(sight_group)

write_csv(bp, here::here("data", "clean", "blueprint.csv"))


## ----data check, eval=FALSE, message=FALSE, warning=FALSE, include=FALSE-----------------------------------------------------------------------------------------------
## for (i in 1:length(herds)) {
##   if (herds[i] != "Redrock/Prairie Creek") {
##     ## make folder structure
##     dir.create(here::here("data", "herd_checks", herds[i]))
##     dir.create(here::here("data", "herd_checks", herds[i], "cleaned_data"))
##     dir.create(here::here("data", "herd_checks", herds[i], "raw_data"))
##     dir.create(here::here("data", "herd_checks", herds[i], "plots"))
## 
##     ## export raw data
##     write_csv(trt_raw %>% filter(herd %in% herds[i]) %>% dplyr::select(herd, treatment, start.year, end.year), here::here("data", "herd_checks", herds[i], "raw_data", "treatments_raw.csv"))
##     write_csv(count.raw %>% filter(herd %in% herds[i]), here::here("data", "herd_checks", herds[i], "raw_data", "counts_raw.csv"))
##     write_csv(recruitment %>% dplyr::select(herd, year, calves, adults) %>% filter(herd %in% herds[i]), here::here("data", "herd_checks", herds[i], "raw_data", "recruitment_raw.csv"))
##     write_csv(surv.raw %>% filter(herd %in% herds[i]), here::here("data", "herd_checks", herds[i], "raw_data", "survival_raw.csv"))
## 
##     ## export cleaned data
##     write_csv(trt_long %>% filter(herd %in% herds[i]), here::here("data", "herd_checks", herds[i], "cleaned_data", "treatments.csv"))
##     write_csv(counts %>% dplyr::select(herd, year, count, Sightability, sd, MinUsed, Est_CL) %>% filter(herd %in% herds[i]), here::here("data", "herd_checks", herds[i], "cleaned_data", "counts.csv"))
##     write_csv(recruitment %>% dplyr::select(herd, year, est, sd) %>% filter(herd %in% herds[i]), here::here("data", "herd_checks", herds[i], "cleaned_data", "recruitment.csv"))
##     write_csv(surv.yr.est %>% dplyr::select(herd, year, est, sd) %>% filter(herd %in% herds[i]), here::here("data", "herd_checks", herds[i], "cleaned_data", "survival.csv"))
## 
##     ## plots
##     recruitment %>%
##       dplyr::select(herd, year, est, sd) %>%
##       filter(herd %in% herds[i]) %>%
##       ggplot(aes(x = year, y = est)) +
##       geom_pointrange(aes(ymin = est - sd, ymax = est + sd), shape = 21) +
##       theme_ipsum() +
##       theme(legend.position = "none") +
##       ylab("Calves (proportion)") +
##       xlab("Year") +
##       labs(x = "Year", title = "Recruitment") +
##       scale_x_continuous(breaks = pretty_breaks(n = 3)) +
##       theme(
##         axis.title.x = element_text(size = 15),
##         axis.title.y = element_text(size = 15),
##         strip.text.x = element_text(size = 15),
##         strip.text.y = element_text(size = 15),
##         axis.text = element_text(size = 10),
##         legend.text = element_text(size = 13),
##         legend.title = element_text(size = 15)
##       )
##     ggsave(here::here("data", "herd_checks", herds[i], "plots", "recruitment_clean.png"), width = 7, height = 5, bg = "white")
## 
## 
##     surv.yr.est %>%
##       dplyr::select(herd, year, est, sd) %>%
##       filter(herd %in% herds[i]) %>%
##       ggplot(aes(x = year, y = est)) +
##       geom_pointrange(aes(ymin = est - sd, ymax = est + sd), shape = 21) +
##       theme_ipsum() +
##       theme(legend.position = "none") +
##       ylab("Annual survial") +
##       xlab("Year") +
##       labs(x = "Year", title = "Survival") +
##       scale_x_continuous(breaks = pretty_breaks(n = 3)) +
##       theme(
##         axis.title.x = element_text(size = 15),
##         axis.title.y = element_text(size = 15),
##         strip.text.x = element_text(size = 15),
##         strip.text.y = element_text(size = 15),
##         axis.text = element_text(size = 10),
##         legend.text = element_text(size = 13),
##         legend.title = element_text(size = 15)
##       )
## 
##     ggsave(here::here("data", "herd_checks", herds[i], "plots", "survival_clean.png"), width = 7, height = 5, bg = "white")
## 
## 
##     recruitment %>%
##       dplyr::select(herd, year, calves, adults) %>%
##       filter(herd %in% herds[i]) %>%
##       pivot_longer(calves:adults) %>%
##       ggplot(aes(x = year, y = value, color = name, group = name)) +
##       geom_point() +
##       theme_ipsum() +
##       ylab("Animals (#)") +
##       xlab("Year") +
##       labs(x = "Year", title = "Recruitment (raw)", color = "Class") +
##       scale_x_continuous(breaks = pretty_breaks(n = 3)) +
##       theme(
##         axis.title.x = element_text(size = 15),
##         axis.title.y = element_text(size = 15),
##         strip.text.x = element_text(size = 15),
##         strip.text.y = element_text(size = 15),
##         axis.text = element_text(size = 10),
##         legend.text = element_text(size = 13),
##         legend.title = element_text(size = 15)
##       )
## 
##     ggsave(here::here("data", "herd_checks", herds[i], "plots", "recruitment_raw.png"), width = 7, height = 5, bg = "white")
## 
## 
## 
##     surv.raw.plot <- surv.raw %>%
##       filter(herd %in% herds[i]) %>%
##       dplyr::select(`Animal ID (RS, CL)`, `Date entry`, `Date exit`, Outcome) %>%
##       group_by(`Animal ID (RS, CL)`) %>%
##       mutate(
##         `Date exit` = ymd(`Date exit`),
##         uID = paste(`Animal ID (RS, CL)`, row_number(), sep = "_")
##       ) %>%
##       pivot_longer(`Date entry`:`Date exit`)
##     ggplot() +
##       geom_line(data = surv.raw.plot, aes(x = value, y = `Animal ID (RS, CL)`, group = uID)) +
##       geom_point(data = surv.raw.plot %>% filter(name %in% "Date exit"), aes(x = value, y = `Animal ID (RS, CL)`, group = `Animal ID (RS, CL)`, color = Outcome)) +
##       theme_ipsum() +
##       ylab("Animal ID") +
##       xlab("Date") +
##       labs(x = "Date", title = "Survival (raw)", color = "Outcome") +
##       scale_x_date(labels = date_format("%m-%Y")) +
##       theme(
##         axis.title.x = element_text(size = 15),
##         axis.title.y = element_text(size = 15),
##         strip.text.x = element_text(size = 15),
##         strip.text.y = element_text(size = 15),
##         axis.text = element_text(size = 10),
##         legend.text = element_text(size = 13),
##         legend.title = element_text(size = 15)
##       )
## 
##     ggsave(here::here("data", "herd_checks", herds[i], "plots", "survival_raw.png"), width = 7, height = 15, bg = "white")
## 
## 
## 
##     ggplot(
##       counts %>% filter(herd %in% herds[i]) %>%
##         mutate(year = as.integer(year)) %>%
##         pivot_longer(SurveyCount:Est_CL) %>%
##         filter(name %in% c("SurveyCount", "MinCount", "Estimate", "Est_CL")),
##       aes(x = year, y = value, color = name)
##     ) +
##       geom_point() +
##       theme_ipsum() +
##       ylab("Abundance") +
##       xlab("Year") +
##       labs(x = "Year", title = "Population Trajectory") +
##       expand_limits(y = 0) +
##       scale_x_continuous(breaks = pretty_breaks(n = 3)) +
##       theme(
##         axis.title.x = element_text(size = 15),
##         axis.title.y = element_text(size = 15),
##         strip.text.x = element_text(size = 15),
##         strip.text.y = element_text(size = 15),
##         axis.text = element_text(size = 10),
##         legend.text = element_text(size = 13),
##         legend.title = element_text(size = 15)
##       )
## 
##     ggsave(here::here("data", "herd_checks", herds[i], "plots", "counts.png"), width = 7, height = 5, bg = "white")
##   }
## }


## ----ss, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------------------------------------
## herds
nherd

## n abundance estimates
nrow(counts)
length(unique(counts$herd)) ## n herds with counts
sum(counts$MinUsed == 1)
sum(counts$MinUsed != 1)

## n sightability estimates
sight.count <- counts %>% filter(!is.na(Sightability) & MinUsed != 1 & !is.na(CollarsSeen))
nrow(sight.count)
length(unique(sight.count$herd)) ## n herds with sightability

## n recruitment estimates
nrow(recruitment)
length(unique(recruitment$herd)) ## n herds with recruit

## n survival estimates
nrow(surv.yr.est)
length(unique(surv.yr.est$herd)) ## n herds with recruit
sum(surv.yr$time) / 365 # n animal-years for survival
sum(surv.yr$event) # n dead
length(surv.yr$id %>% unique()) # n animals for survival

## timeframe
min(counts$year, recruitment$year, surv.yr$year)

## treatments
trt_long %>%
  filter(applied == 1) %>%
  summarise(
    n = n(),
    n.herds = n_distinct(herd),
    min = min(year),
    max = max(year)
  )

trt_long %>%
  filter(applied == 1, year > 2004) %>%
  summarise(
    n = n(),
    n.herds = n_distinct(herd),
    min = min(year),
    max = max(year)
  )

trt_long %>%
  filter(applied == 1) %>%
  group_by(herd, year) %>%
  count() %>%
  filter(n > 1)

trt_long %>%
  filter(intensity == "low")

trt_long %>%
  filter(applied == 1) %>%
  group_by(treatment) %>%
  summarise(
    n = n(),
    n.herds = n_distinct(herd),
    min = min(year),
    max = max(year)
  )






trt_long %>%
  left_join(read_csv(here::here("tables", "demog.csv")) %>% dplyr::select(year = yrs, herd, totNMF)) %>%
  filter(applied == 1) %>%
  mutate(application = case_when(intensity == "low" | totNMF < 30 ~ "low", TRUE ~ "standard")) %>%
  group_by(herd, treatment, application) %>%
  summarise(
    start = min(year),
    end = max(year)
  ) %>%
  ungroup() %>%
  arrange(treatment, herd, application) %>%
  write_csv(here::here("tables", "appendix", "treatment.table.csv"))




## check timing
count.timing <- counts %>%
  group_by(herd, timing) %>%
  count() %>%
  pivot_wider(values_from = n, names_from = timing)

## only herds that matter are fall.herds and spring.herds
a <- counts %>%
  mutate(season = case_when(
    timing %in% c("Winter", "March", "April", "February") ~ "winter",
    timing %in% c("August", "September", "October", "November") ~ "fall",
    timing %in% c("June", "July") ~ "spring"
  )) %>%
  dplyr::select(herd, year, season, est = Est_CL) %>%
  mutate(param = "count") %>%
  rbind(recruitment %>% dplyr::select(herd, year, season, est) %>% mutate(param = "R")) %>%
  filter(herd %in% c(fall.herds, spring.herds)) %>%
  group_by(herd, season, param) %>%
  count()

