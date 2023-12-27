BC AB Caribou IPM run
================
Clayton T. Lamb
27 December, 2023

## Load Data

``` r
library(renv)
##to pull packages
#restore(repos="https://cloud.r-project.org")
library(tidybayes)
library(ggmcmc)
library(hrbrthemes)
library(tidylog)
library(tidyverse)

# Load data ---------------------------------------------------------------


hd <- read.csv(here::here("data/clean/blueprint.csv"))
hn <- hd %>%
    dplyr::select(herd, herd_num)
sgt <- read.csv(here::here("data/clean/herd_sightability.csv")) %>%
    arrange(herd) %>%
    left_join(hn, by = "herd")
afs <- read.csv(here::here("data/clean/survival.csv")) %>%
    arrange(herd) %>%
    left_join(hn, by = "herd")
afr <- read.csv(here::here("data/clean/recruitment.csv")) %>%
    arrange(herd) %>%
    left_join(hn, by = "herd")
counts <- read.csv(here::here("data/clean/counts.csv")) %>%
    arrange(herd) %>%
    left_join(hn, by = "herd")
trt <- read.csv(here::here("data/clean/treatments.csv")) %>%
    arrange(herd) %>%
    left_join(hn, by = "herd")
ecotype <- read.csv(here::here("data/raw/treatment.csv")) %>%
  dplyr::select(herd=Herd,ECCC=ECCC_Recov_Grp, COSEWIC=COSEWIC_Grp, Heard_Vagt1998=Heard.and.Vagt.1998.grouping)%>%
  distinct()

sexratios <-  read.csv(here::here("data/clean/sexratios.csv")) %>%
  arrange(herd) %>%
  left_join(hn, by = "herd")
sexratio_summary <-  read.csv(here::here("data/clean/sexratio_summary.csv"))

# 0 = treatment did not occur for that herd in that year
# 1 = treatment did occur
# NA = herd not monitored? Possibly change to 0s
#trt$applied[is.na(trt$applied)] <- 0
```

## Prep data for IPM

``` r
# Herds and herd number
herds <- unique(hd$herd)
nherd <- nrow(hn)
hd <- hd%>%
  left_join(ecotype%>%
              add_row(herd="Quintette Full",
                      ECCC="Central Group",
                      COSEWIC="DU8",
                      Heard_Vagt1998="Northern"), by="herd")
hd$demog_grp <- hd$ECCC%>%factor%>%as.numeric() 
hd[hd$herd%in%c("Itcha-Ilgachuz"),"demog_grp"]<-4
nsight_grp <- length(unique(hd$sight_grp))
ndemog_grp <- length(unique(hd$demog_grp))

# Index for herd per demog group
hd <- hd %>%
    group_by(demog_grp) %>%
    mutate(hd_per_dg = 1:n()) %>%
    as.data.frame()

nhd_per_demog_grp <- hd %>%
    dplyr::select(demog_grp, hd_per_dg) %>%
    group_by(demog_grp) %>%
    mutate(max_hd_num = max(hd_per_dg)) %>%
    slice(1) %>%
    as.data.frame()
    

#  Years of study
yrs <-  seq(from = min(trt$year), to = max(trt$year), by = 1)
nyr <- length(yrs)
yr_idx <- seq(from = 1, to = nyr, by = 1)
yr_df <- as.data.frame(cbind(yrs, yr_idx))

# Age classes
nage <- 2



# Restrict year range for trouble shooting
rest_yr <- nyr


ntrt <- length(unique(trt$treatment))

# Generate treatment matrices
mtrt <- trt %>%
    dplyr::filter(treatment == "reduce moose") %>%
    left_join(yr_df, by = c("year" = "yrs")) %>%
    dplyr::select(herd_num, yr_idx, applied) %>%
    tidyr::spread(yr_idx, applied)
mtrt <- mtrt[,-1]
mtrt[is.na(mtrt)] <- 0


wtrt <- trt %>%
    dplyr::filter(treatment %in% c("reduce wolves")) %>%
    left_join(yr_df, by = c("year" = "yrs")) %>%
    dplyr::select(herd_num, yr_idx, applied) %>%
    tidyr::spread(yr_idx, applied)
wtrt <- wtrt[,-1]
wtrt[is.na(wtrt)] <- 0


ptrt <- trt %>%
    dplyr::filter(treatment == "pen") %>%
    left_join(yr_df, by = c("year" = "yrs")) %>%
    dplyr::select(herd_num, yr_idx, applied) %>%
    tidyr::spread(yr_idx, applied)
ptrt <- ptrt[,-1]
ptrt[is.na(ptrt)] <- 0


ftrt <- trt %>%
    dplyr::filter(treatment == "feed") %>%
    left_join(yr_df, by = c("year" = "yrs")) %>%
    dplyr::select(herd_num, yr_idx, applied) %>%
    tidyr::spread(yr_idx, applied)
ftrt <- ftrt[,-1]
ftrt[is.na(ftrt)] <- 0


strt <- trt %>%
  dplyr::filter(treatment == "sterilize wolves") %>%
  left_join(yr_df, by = c("year" = "yrs")) %>%
  dplyr::select(herd_num, yr_idx, applied) %>%
  tidyr::spread(yr_idx, applied)
strt <- strt[,-1]
strt[is.na(strt)] <- 0


ttrt <- trt %>%
  dplyr::filter(treatment == "transplant") %>%
  left_join(yr_df, by = c("year" = "yrs")) %>%
  dplyr::select(herd_num, yr_idx, applied) %>%
  tidyr::spread(yr_idx, applied)
ttrt <- ttrt[,-1]
ttrt[is.na(ttrt)] <- 0


# Recruitment
rdat <- afr %>%
  #filter(!herd%in%"Tonquin")%>%
    mutate(herd = herd_num,
        age = NA, 
        sex = 1,
        #mu = recruitment/100,
        mu = est,
        #tau = 500) %>%
        tau = 1/(sd*sd)) %>%
    left_join(yr_df, by = c("year" = "yrs")) %>%
    dplyr::filter(!is.na(mu)) %>%
    dplyr::filter(mu != "Inf") %>%
    dplyr::select(herd, year = yr_idx, age, sex, mu, tau) %>%
    dplyr::filter(year <= rest_yr) %>%
    arrange(herd, year)
rdat$tau[rdat$tau == "Inf"] <- mean(rdat$tau[rdat$tau != "Inf"], na.rm = TRUE)
nr <- nrow(rdat)    

##mean for prior
afr %>%
  filter(season_int==1,
         year<=2010)%>%
  summarise(meanR=mean(est))

# Survival 
sdat <- afs %>%
  filter(est<0.99 & est>0.5)%>%
    mutate(herd = herd_num,
        age = 2, 
        sex = 1,
        mu = est,
        #tau = 100) %>%
        tau = 1/(sd*sd)) %>%
    left_join(yr_df, by = c("year" = "yrs")) %>%
    dplyr::filter(!is.na(mu)) %>%
    dplyr::filter(mu != "Inf") %>%
    dplyr::select(herd, year = yr_idx, age, sex, mu, tau) %>%
    dplyr::filter(year <= rest_yr) %>%
    arrange(herd, year)
sdat$tau[sdat$tau == "Inf"] <- mean(sdat$tau[sdat$tau != "Inf"])
ns <- nrow(sdat)    

#sdat$mu[sdat$mu == 1.0000000] <- 0.999



# Sightability
# Estimates of 1 are not usable here, not sure why they are left in or how this 
# even works given qlogis(1) equals Inf
# We need to transform the uncertainty along with the mean, changes made in JAGS
grp_p <- sgt %>% 
    left_join(hd, by = c("herd", "herd_num")) %>%
    group_by(sight_grp) %>%
    mutate(
      n=n(),
      mean_grp_p = mean(est, na.rm = TRUE),
      mean_grp_pvar = mean(sd^2, na.rm = TRUE), # If this works I am willing to bet it is only valid when performed on the variance and not the sd 
      mean_grp_ptau = 1/(mean_grp_pvar)
    ) %>% 
    slice(1) %>% 
    as.data.frame() %>%
    dplyr::select(sight_grp, mean_grp_p, mean_grp_pvar, mean_grp_ptau) %>%
    arrange(sight_grp)
    

mean_grp_p <- grp_p$mean_grp_p
mean_grp_ptau <- grp_p$mean_grp_ptau

sgt_grp_ind <- matrix(nrow = nherd, ncol = nyr)
for(sg in 1:nrow(sgt_grp_ind)){
    sgt_grp_ind[sg,] <- hd$sight_grp[hd$herd_num == sg]
    }


    
# Abundance -- separate out rows where no sightability data provided

# Survey counts (have sightability)
cdat <- counts %>% 
    dplyr::filter(!is.na(Sightability)) %>%
    mutate(herd = herd_num, 
        age = NA, 
        sex = NA,
        mu = as.integer(count),
        tau = case_when(
          count > 0 ~ 1/((count*0.15)^2), 
          count %in% 0 ~ 1
        )
    ) %>%
        #tau = 1/(sd*sd)) %>%
    left_join(yr_df, by = c("year" = "yrs")) %>%
    dplyr::select(herd, year = yr_idx, age, sex, mu, tau) %>%
    dplyr::filter(year <= rest_yr) %>%
    arrange(herd, year)

nc <- nrow(cdat)


# Sightability data 
pdat <- counts %>% 
    left_join(hd, by = c("herd", "herd_num")) %>%
    dplyr::filter(!is.na(Sightability)) %>%
    mutate(herd = herd_num, 
        age = NA, 
        sex = NA,
        mu = Sightability,
        tau = 1/(sd*sd)) %>%
    left_join(yr_df, by = c("year" = "yrs")) %>%
    dplyr::select(herd, year = yr_idx, age, sex, mu, tau) %>%
    dplyr::filter(year <= rest_yr) %>%
    arrange(herd, year)
#pdat$pdat[pdat$pdat == 1] <- 0.95  
pdat$tau[pdat$tau == "Inf"] <- mean(pdat$tau[pdat$tau != "Inf"])
np <- nrow(pdat)

# Minimum Estimates
edat <- counts %>%
    dplyr::filter(is.na(Sightability)) %>%
    left_join(hd, by = c("herd", "herd_num")) %>%
    mutate(herd = herd_num, 
        age = NA, 
        sex = NA,
        mu = as.integer(count),
        tau = case_when(
          count > 0 ~ 1/((count*0.2)^2), 
          count %in% 0 ~ 10
        ),
        #tau = 1/(sd*sd),
        sight_grp = sight_grp) %>%
    left_join(yr_df, by = c("year" = "yrs")) %>%
    dplyr::select(herd, year = yr_idx, age, sex, mu, tau, sight_grp) %>%
    dplyr::filter(year <= rest_yr) %>%
    arrange(herd, year)

edat <- edat%>%
  rbind(edat%>%
          filter((herd==19 & year==48))%>%
          mutate(year=49,mu=99, tau=0.01))

ne <- nrow(edat)


# Survey timing: most are in March. Those NOT in March:
# 1 = March (reference), 2 = Fall, 3 = Spring

survey_timing <- afr%>%
  left_join(yr_df, by = c("year" = "yrs")) %>%
  dplyr::select(herd_num,yr_idx,season_int)%>%
  rbind(expand.grid(herd_num=(1:nherd)[!(1:nherd)%in%unique(afr$herd_num)],
               yr_idx=(1:nyr)[!(1:nyr)%in%unique(afr%>%
                                                   left_join(yr_df, by = c("year" = "yrs"))%>%pull(yr_idx))],
               season_int=1))%>%
  tidyr::spread(yr_idx, season_int)
survey_timing <- survey_timing[,-1]
survey_timing[is.na(survey_timing)] <- 1



########## UPDATE 11/2/22
# Calves born in May
# May/Jun survey = Spring survey (age = 1 month) --- season_int = 3
# Oct survey = Fall survey (age = 5 months) ---- season_int = 2
# Mar survey = winter survey (age = 10 months) ---- season_int = 1

month_offset <- survey_timing
month_offset[month_offset == 1] <- 0 
month_offset[month_offset == 2] <-  5 # 5 units offset
month_offset[month_offset == 3] <-  9 # 9 units offset


# # Starting population size of each herd
#take average of all counts that occured within 3 years of first count per herd
#   Using the "Est_CL" column 
# Starting population size of each herd
meancount <- counts %>%
    arrange(herd_num, year) %>%
    group_by(herd_num) %>%
    mutate(keep = year - min(year)) %>%
    dplyr::filter(keep < 4) %>%
    mutate(n1_mean = mean(Est_CL, na.rm = TRUE)) %>%
    slice(1) %>%
    as.data.frame() %>%
    dplyr::select(herd, herd_num, n1_mean) %>%
    right_join(hd, by = c("herd", "herd_num")) %>%
    mutate(n1_mean = ifelse(is.na(n1_mean), mean(n1_mean, na.rm = TRUE), n1_mean)) %>%
    arrange(herd_num)
n1s <- meancount$n1_mean

##starting age class distribution (sensitivity tested this in McNay et al. 2022 and it doesnt have a meaningful effect even if substantially changed)
n1 <- matrix(NA, nrow = nherd, ncol = nage)
for(h in 1:nherd){
    n1[h,1] <- (n1s[h]*0.5)*0.15
    n1[h,2] <- (n1s[h]*0.64)*0.85 
    }


# Sex ratio
srdat <- sexratios%>%
  mutate(herd = herd_num,
         mu = sratio,
         tau = 1/(sratio.sd*sratio.sd)) %>%
  left_join(yr_df, by = c("year" = "yrs")) %>%
  dplyr::filter(!is.na(mu)) %>%
  dplyr::filter(mu != "Inf") %>%
  dplyr::select(herd, year = yr_idx, mu, tau) %>%
  dplyr::filter(year <= rest_yr) %>%
  arrange(herd, year)%>%
  distinct(herd,year, .keep_all=TRUE)

nsr <- nrow(srdat)  

meansr <- array(NA, c(1,2))
meansr[1,1] <- sexratio_summary[1,1]
meansr[1,2] <- 1/(sexratio_summary[1,2]^2)


#  First year of data for each herd
first_s <- sdat %>%
    group_by(herd) %>%
    arrange(year) %>%
    slice(1) %>%
    as.data.frame() %>%
    dplyr::select(herd, first_s_year = year) 
first_r <- rdat %>%
    group_by(herd) %>%
    arrange(year) %>%
    slice(1) %>%
    as.data.frame() %>%
    dplyr::select(herd, first_r_year = year)    
first_c <- cdat %>%
    group_by(herd) %>%
    arrange(year) %>%
    slice(1) %>%
    as.data.frame() %>%
    dplyr::select(herd, first_c_year = year)    
first_e <- edat %>%
    group_by(herd) %>%
    arrange(year) %>%
    slice(1) %>%
    as.data.frame() %>%
    dplyr::select(herd, first_e_year = year)

first_per_herd <- first_s %>%
    full_join(first_r, by = "herd") %>%
    full_join(first_c, by = "herd") %>%
    full_join(first_e, by = "herd") %>%
    rowwise(.) %>%
    mutate(first_yr = min(first_s_year, first_r_year, first_c_year, first_e_year, na.rm = TRUE)) %>%
    as.data.frame() %>%
    arrange(herd)%>%
  left_join(hd%>%dplyr::select(name=herd,herd=herd_num))%>%
  left_join(yr_df%>%dplyr::select(first_yr=yr_idx, yrs))
first <- first_per_herd$first_yr

first <- rep(1,nherd)


#  First year of treatment for each herd
treatment_start<- trt%>%
  filter(applied==1 & !treatment%in%c("transplant","reduce moose"))%>%
  group_by(herd_num)%>%
  summarize(yr=min(year)-1)%>% ##year before trt starts
  full_join(tibble(herd_num=1:nherd,
                    yr_fill=2021),
             by="herd_num")%>%
  left_join(hd%>%dplyr::select(name=herd,herd_num))%>%
  mutate(yrs=case_when(is.na(yr)~yr_fill,
                       name%in%"Columbia North"~2004, ##CN reduce moose did work so keep this one
                      TRUE~yr))%>%
  left_join(yr_df)%>%
  arrange(herd_num)%>%
  pull(yr_idx)

startmean <- first_per_herd$first_yr

    
nyr <- rest_yr
```

## Gather data inputs in a list

``` r
ipm_dat <- list(
    nherd = nherd,
    nyr = nyr,
    first = first,
    
    month_offset = month_offset,
    treatment_start=treatment_start,
    
    nc = nc, 
    ne = ne,
    ns = ns,
    nr = nr,
    nsr = nsr,
    
    nsight_grp = nsight_grp,
    ndemog_grp = ndemog_grp,
    sight_grp = hd$sight_grp,
    demog_grp = hd$demog_grp,
    
    mean_grp_p = mean_grp_p,
    mean_grp_ptau = mean_grp_ptau,
    meansr = meansr,
    
    n1 = n1, 
    cdat = cdat,
    edat = edat,
    sdat = sdat,
    srdat = srdat,
    pdat = pdat,
    rdat = rdat,
    
    mtrt = mtrt, 
    wtrt = wtrt,
    ptrt = ptrt, 
    ftrt = ftrt,
    strt = strt,
    ttrt = ttrt)
    
    
#  Initial values for N to avoid parent node erros
#  Indexing was off here
Nst <- array(NA_integer_, c(nherd, nyr, nage))

for(h in 1:nherd) {
  for(a in 1:nage) {
    Nst[h, first[h]:nyr, a] <- as.integer(
      max(
        c(n1[h, a], cdat$mu[cdat$herd == h], edat$mu[edat$herd == h])))  
  }
}

ipm_inits <- function(){ 
  list(
    N = Nst
    )
}


#  Model parameters to monitor
model_parms <- c(

    "adj_totNMF", "totN", "N", 
    "totCalves", "totAdults", "totAdultsMF", "totCalvesMF", "totNMF", 

    "lambda", "logla",
    
    "S", "R",  "R_adj",
    "p_mu", "sight_est", "p",
    
    "muS", "muR",
    
    "surv_yr_dg", "recruit_yr_dg",
    "surv_tau_yr_dg", "recruit_tau_yr_dg",
    "surv_hd", "recruit_hd",
    
    "surv_sd_yr_dg","sight_tau_yr_sg",
    
    "offset",
    
    "pred_totNMF",
    "pred_totCalvesMF",
    "pred_totAdultsMF",
    "mean_R",
    "mean_S",
    "sexratio",
    
    "mtrt_eff_s", "mtrt_eff_r",
    "wtrt_eff_s", "wtrt_eff_r",
    "ptrt_eff_s", "ptrt_eff_r",
    "ftrt_eff_s", "ftrt_eff_r",
    "strt_eff_s", "strt_eff_r",
    "ttrt_eff_s", "ttrt_eff_r",
    
    "moosewolf",
    "penwolf"
    )
```

## Run IPM

``` r
# nth <- 90
# nbu <- 30000
# nch <- 3
# nad <- 60000
# nit <- 400000
# 
# nth <- 1
# nbu <- 200
# nch <- 3
# nad <- 100
# nit <- 300

nth <- 90
nbu <- 0
nch <- 3
nad <- 60000
nit <- 400000


out <- jagsUI::jags(data=ipm_dat, 
    inits = ipm_inits,
    parameters.to.save=model_parms,
    model.file = here::here("jags/BCAB_IPM_20231211.txt"),
    n.chains = nch,
    n.cores = nch,
    n.iter = nit,
    n.burnin = nbu,
    n.thin = nth,
    n.adapt = nad)
```

## Save outputs

``` r
#mcmcplots::mcmcplot(out$samples, par = c("S", "R", "totNMF")) 
saveRDS(out, file = here::here("jags/output/BCAB_CaribouIPM_23update.rds"))
```
