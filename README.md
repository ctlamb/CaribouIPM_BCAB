BC AB Caribou IPM Results
================
Clayton T. Lamb
11 March, 2024

## Load Data

``` r
library(renv)
## to pull packages
# restore(repos="https://cloud.r-project.org")
library(ggmap)
library(MCMCvis)
library(tidybayes)
library(ggmcmc)
library(boot)
library(cowplot)
library(hrbrthemes)
library(lme4)
library(tidymodels)
library(broom.mixed)
library(RColorBrewer)
library(ggridges)
library(ggrepel)
library(tidylog)
library(gt)
library(patchwork)
library(sf)
library(basemaps)
library(ggtext)
library(knitr)
library(terra)
library(tidyverse)

# Load data ---------------------------------------------------------------

## IPM Output
out <- readRDS(file = here::here("jags/output/BCAB_CaribouIPM_23update.rds"))

## IPM input to compare results
hd <- read.csv("data/clean/blueprint.csv")
hn <- hd %>%
  dplyr::select(herd, herd_num)
afs <- read.csv("data/clean/survival.csv") %>%
  arrange(herd) %>%
  left_join(hn, by = "herd")
afr <- read.csv("data/clean/recruitment.csv") %>%
  arrange(herd) %>%
  left_join(hn, by = "herd")
counts <- read.csv("data/clean/counts.csv") %>%
  arrange(herd) %>%
  left_join(hn, by = "herd")
trt <- read.csv("data/clean/treatments.csv") %>%
  arrange(herd) %>%
  left_join(hn, by = "herd")
ecotype <- read.csv("data/raw/treatment.csv") %>%
  dplyr::select(herd = Herd, ECCC = ECCC_Recov_Grp, COSEWIC = COSEWIC_Grp, Heard_Vagt1998 = Heard.and.Vagt.1998.grouping) %>%
  distinct()
labels <- read.csv("data/clean/labels.csv")
label.lookup <- read.csv("data/clean/label_lookup.csv")

#  Years of study
yrs <- seq(from = min(trt$year), to = max(trt$year), by = 1)
nyr <- length(yrs)
yr_idx <- seq(from = 1, to = nyr, by = 1)
yr_df <- as.data.frame(cbind(yrs, yr_idx))

rm(yrs)
rm(nyr)
rm(yr_idx)
```

## Data housekeeping

``` r
## set up colors for plotting
display.brewer.pal(8, "Accent")
```

![](README_files/figure-gfm/Data%20housekeeping-1.png)<!-- -->

``` r
cols <- RColorBrewer::brewer.pal(8, "Accent")

### Change Narraway BC to Bearhole Redwillow
trt <- trt %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
hd <- hd %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
counts <- counts %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
afr <- afr %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
afs <- afs %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
ecotype <- ecotype %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))

# Pull out some values, summarize to demog data frame ---------------------

#### TREATMENT COMBOS
treatment.combos <- trt %>%
  filter(applied == 1) %>%
  group_by(herd, year) %>%
  summarize(
    trt = paste(
      gsub(" ", "", treatment), # Remove spaces
      collapse = "-" # Add - between treatments
    ),
    .groups = "drop"
  ) %>%
  dplyr::select(herd, yrs = year, trt)


## pull posterior draws, add in herd, year, and treatment to each herd-year
ndraws <- 10000
demog.raw <- out %>%
  spread_draws(
    c(totNMF, totN, totAdults, S, R_adj, lambda, SR)[i, j],
    ndraws = ndraws
  ) %>%
  median_qi(.width = 0.9) %>%
  left_join(yr_df %>% rename(j = yr_idx), by = "j") %>%
  left_join(hd %>% dplyr::select(herd, i = herd_num), by = "i") %>%
  left_join(treatment.combos, by = c("herd", "yrs")) %>%
  mutate(trt = replace_na(trt, "Reference")) %>%
  rename(R = R_adj, R.lower = R_adj.lower, R.upper = R_adj.upper) %>%
  ungroup()

## determine which herds to keep
herds.keep <- demog.raw %>%
  filter(herd != "Quintette Full") %>%
  distinct(herd) %>%
  pull()


## treatment sample sizes
trt.n <- demog.raw %>%
  filter(herd %in% herds.keep) %>%
  group_by(trt) %>%
  summarize(
    n_herds = n_distinct(herd),
    n_yrs = n()
  )

write_csv(trt.n, here::here("tables", "trt_sample_sizes.csv"))


## filter herds to a period where not fully extirpated
extirpated.yr <- demog.raw %>%
  filter(round(totN, 0) == 0) %>%
  group_by(herd) %>%
  filter(yrs == min(yrs)) %>%
  ungroup()


herds <- unique(demog.raw$herd)

demog.trim <- tibble()
for (i in 1:length(herds)) {
  if (herds[i] %in% extirpated.yr$herd) {
    yr <- extirpated.yr %>%
      dplyr::filter(herd == !!herds[i])
    a <- demog.raw %>%
      dplyr::filter(herd == !!herds[i] & yrs <= yr$yrs)
  }

  if (!herds[i] %in% extirpated.yr$herd) {
    a <- demog.raw %>%
      dplyr::filter(herd == !!herds[i])
  }

  demog.trim <- bind_rows(a, demog.trim)
}

demog.raw <- demog.trim
rm(demog.trim)
rm(a)
rm(yr)

write_csv(demog.raw, "tables/demog.csv")

demog <- demog.raw %>% filter(herd %in% herds.keep)

demog.mod <- demog %>% filter(totAdults > 10 & totNMF > 20) ## modelling data that doesn't include functionally extirpated herds, demography gets unstable

ggplot(demog, aes(y = log(lambda), x = totNMF)) +
  geom_point(alpha = 0.3) +
  geom_vline(xintercept = 30, color = "red") +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15)
  ) +
  labs(
    title = "Population growth unstable when subpopulation too small",
    subtitle = "remove data below functionally extirpated abundance (red line)",
    x = "Subopulation abundance",
    y = "Instanteous rate of increase (r)"
  )
```

![](README_files/figure-gfm/Data%20housekeeping-2.png)<!-- -->

``` r
## how many herds had sightability
counts %>%
  filter(MinUsed == 0, !is.na(Sightability)) %>%
  distinct(herd)

counts %>%
  mutate(sight = case_when(!is.na(Sightability) ~ 1, TRUE ~ 0)) %>%
  count(sight)
```

## First year of demographic data for each herd

``` r
raw.demog <- rbind(
  afr %>% dplyr::select(herd, year, est) %>% mutate(type = "Recruit"),
  afs %>% dplyr::select(herd, year, est) %>% mutate(type = "Surv"),
  counts %>% dplyr::select(herd, year, est = Est_CL) %>% mutate(type = "Count")
)

### Identify first year of demographic data for each herd
first.yr <- raw.demog %>%
  dplyr::select(herd, first.year = year) %>%
  distinct() %>%
  group_by(herd) %>%
  filter(first.year == min(first.year))

ggplot(first.yr, aes(x = first.year)) +
  stat_bin(aes(y = cumsum(..count..)), geom = "step") +
  theme_ipsum() +
  labs(
    x = "Year", y = "Cumulative count",
    title = "First year of demographic data for the 41 SMC herds",
    subtitle = "Cumulative count of herds being monitored through time"
  ) +
  theme(
    plot.title = element_text(size = 18),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.position = "bottom"
  )
```

![](README_files/figure-gfm/Firstyr-1.png)<!-- -->

## Plot Abundance by Herd

``` r
## Prep data and layout for plot

## treatments for plot,and where they go
trt.plot <- trt %>%
  filter(applied == 1) %>%
  left_join(demog %>% group_by(herd) %>% summarize(max = max(totNMF.upper))) %>%
  mutate(y = case_when(
    treatment %in% "transplant" & herd %in% "South Selkirks" ~ 250,
    treatment %in% "reduce wolves" & herd %in% "Charlotte Alplands" ~ 150,
    treatment %in% "transplant" & herd %in% "Telkwa" ~ 150,
    treatment %in% "sterilize wolves" & herd %in% c("Wells Gray North") ~ (max + 20) - (max * 0.25),
    treatment %in% "reduce wolves" & herd %in% c("Wells Gray North") ~ (max + 20) - (max * 0.10),
    treatment %in% "sterilize wolves" & herd %in% "Barkerville" ~ 160,
    treatment %in% "reduce wolves" & herd %in% "Barkerville" ~ 135,
    treatment %in% "reduce wolves" & herd %in% "Kennedy Siding" ~ max - (max * 0.05),
    treatment %in% "feed" & herd %in% "Kennedy Siding" ~ max - (max * 0.20),
    treatment %in% "pen" & herd %in% "Columbia North" ~ 50,
    treatment %in% "reduce moose" & herd %in% "Columbia North" ~ -5,
    treatment %in% "reduce wolves" & herd %in% "Columbia North" ~ 90,
    treatment %in% "reduce moose" & herd %in% "Groundhog" ~ 88,
    treatment %in% "reduce wolves" & herd %in% "Groundhog" ~ 70,
    treatment %in% "reduce moose" & herd %in% "Hart North" ~ 420,
    treatment %in% "reduce wolves" & herd %in% "Hart North" ~ 340,
    treatment %in% "reduce wolves" ~ max - (max * 0.1),
    treatment %in% "transplant" ~ max - (max * 0.15),
    treatment %in% "reduce moose" ~ max - (max * 0.2),
    treatment %in% "sterilize wolves" ~ max - (max * 0.25),
    treatment %in% "pen" ~ max - (max * 0.3),
    treatment %in% "feed" ~ max - (max * 0.35)
  ))

### a couple functions for the red boxes around herd names of extirpated
element_textbox_highlight <- function(..., hi.labels = NULL, hi.fill = NULL,
                                      hi.col = NULL, hi.box.col = NULL, hi.family = NULL) {
  structure(
    c(
      element_textbox(...),
      list(hi.labels = hi.labels, hi.fill = hi.fill, hi.col = hi.col, hi.box.col = hi.box.col, hi.family = hi.family)
    ),
    class = c("element_textbox_highlight", "element_textbox", "element_text", "element")
  )
}

element_grob.element_textbox_highlight <- function(element, label = "", ...) {
  if (label %in% element$hi.labels) {
    element$fill <- element$hi.fill %||% element$fill
    element$colour <- element$hi.col %||% element$colour
    element$box.colour <- element$hi.box.col %||% element$box.colour
    element$family <- element$hi.family %||% element$family
  }
  NextMethod()
}

## prep simulated counterfactuals

sims <- out %>%
  gather_draws(c(pred_totNMF, totNMF, totAdults)[i, j], ndraws = ndraws) %>%
  left_join(yr_df %>% rename(j = yr_idx), by = "j") %>%
  left_join(hd %>% dplyr::select(herd, i = herd_num), by = "i") %>%
  dplyr::left_join(treatment.combos, by = c("herd", "yrs")) |>
  dplyr::mutate(
    trt = tidyr::replace_na(trt, "Reference")
  ) %>%
  filter(herd %in% herds.keep)

sims %>% write_csv(here::here("tables", "draws", "sims.draws.csv"))

sims.plot <- sims %>%
  dplyr::group_by(yrs, .variable, .draw) %>%
  dplyr::summarize(
    mean = sum(.value),
    .groups = "drop"
  ) %>%
  dplyr::group_by(yrs, .variable) %>%
  dplyr::summarize(
    LCL = quantile(mean, 0.05),
    UCL = quantile(mean, 0.95),
    mean = mean(mean),
    .groups = "drop"
  )

ext.yr <- demog %>%
  dplyr::select(herd, yrs, totAdults, totNMF) %>%
  group_by(herd) %>%
  mutate(
    viable = case_when(yrs == 2021 & round(totAdults, 0) > 8 & round(totNMF, 0) > 20 ~ 1, TRUE ~ 0),
    viable = max(viable)
  ) %>%
  filter(round(totAdults, 0) <= 8 | round(totNMF, 0) <= 20, viable == 0) %>% ## pull if 8 or fewer females present, or 20 or fewer total
  filter(yrs == (min(yrs))) %>%
  left_join(labels, by = "herd") %>%
  mutate(herd2 = paste0(number_label, ".", herd, " (", human, "%)")) %>%
  left_join(
    sims.plot %>%
      filter(.variable %in% "totNMF") %>%
      dplyr::select(yrs, mean),
    by = "yrs"
  ) %>%
  ungroup()

sims.plot <- sims.plot %>%
  filter(.variable %in% c("totNMF", "pred_totNMF"))

## Herd Abundance
ggplot() +
  geom_ribbon(
    data = demog %>%
      ungroup() %>%
      left_join(labels, by = "herd") %>%
      mutate(herd = paste0(number_label, ".", herd, " (", human, "%)") %>%
        fct_reorder(.f = ., .x = number_label, .fun = max)),
    aes(x = yrs, y = totNMF, ymin = totNMF.lower, ymax = totNMF.upper), alpha = 0.4, color = NA, fill = cols[3]
  ) +
  geom_line(
    data = demog %>%
      ungroup() %>%
      left_join(labels, by = "herd") %>%
      mutate(herd = paste0(number_label, ".", herd, " (", human, "%)") %>%
        fct_reorder(number_label)),
    aes(x = yrs, y = totNMF, ymin = totNMF.lower, ymax = totNMF.upper), size = 1, color = cols[3]
  ) +
  geom_rug(
    data = raw.demog %>%
      filter(herd %in% herds.keep) %>%
      ungroup() %>%
      left_join(labels, by = "herd") %>%
      mutate(herd = paste0(number_label, ".", herd, " (", human, "%)") %>%
        fct_reorder(number_label)),
    aes(x = year), sides = "t", length = unit(0.05, "npc"), alpha = 0.5
  ) +
  theme_ipsum() +
  theme(legend.position = "none") +
  facet_wrap(vars(herd), scales = "free_y", ncol = 7) +
  labs(x = "", y = "", title = "Abundance") +
  expand_limits(y = 0) +
  scale_x_continuous(
    limits = c(1973, 2026),
    breaks = seq(1980, 2020, by = 20)
  ) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text = element_textbox_highlight(
      size = 14,
      fill = "white", box.color = "white",
      halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
      padding = margin(b = 2, t = 2, r = 2, l = 2), margin = margin(b = 8, t = 2),
      hi.labels = ext.yr$herd2,
      hi.fill = "firebrick", hi.box.col = "firebrick", hi.col = "white"
    ),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank()
    # axis.line = element_line(colour = "black")
  ) +
  geom_point(
    data = counts %>%
      filter(herd %in% herds.keep) %>%
      ungroup() %>%
      left_join(labels, by = "herd") %>%
      mutate(herd = paste0(number_label, ".", herd, " (", human, "%)") %>%
        fct_reorder(number_label)),
    aes(x = year, y = Est_CL), size = 0.5, alpha = 0.7
  ) +
  geom_linerange(
    data = counts %>%
      filter(herd %in% herds.keep) %>%
      ungroup() %>%
      mutate(Est_CL.max = case_when(Est_CL.max > 5000 ~ 5000, herd == "Rainbows" & Est_CL.max > 500 ~ 500, TRUE ~ Est_CL.max)) %>%
      left_join(labels, by = "herd") %>% mutate(herd = paste0(number_label, ".", herd, " (", human, "%)") %>%
        fct_reorder(number_label)),
    aes(x = year, ymin = Est_CL.min, ymax = Est_CL.max), alpha = 0.5
  ) +
  geom_point(
    data = trt.plot %>%
      filter(herd %in% herds.keep) %>%
      ungroup() %>%
      left_join(labels, by = "herd") %>%
      mutate(herd = paste0(number_label, ".", herd, " (", human, "%)") %>%
        fct_reorder(number_label)),
    aes(x = year, y = y, group = treatment), size = 0.5
  ) +
  scale_color_manual(values = cols[-4]) +
  geom_text(
    data = trt.plot %>%
      filter(herd %in% herds.keep) %>%
      ungroup() %>%
      distinct(herd, treatment, y) %>%
      mutate(t = str_remove(treatment, "reduce ") %>% str_sub(1, 1)) %>%
      left_join(labels, by = "herd") %>% mutate(herd = paste0(number_label, ".", herd, " (", human, "%)") %>%
        fct_reorder(number_label)),
    aes(label = t, x = 2026, y = y),
    direction = "y"
  ) +
  coord_cartesian(clip = "off")
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

    ## Warning in geom_line(data = demog %>% ungroup() %>% left_join(labels, by = "herd") %>% : Ignoring unknown aesthetics: ymin
    ## and ymax

    ## Warning in geom_text(data = trt.plot %>% filter(herd %in% herds.keep) %>% : Ignoring unknown parameters: `direction`

![](README_files/figure-gfm/Plot%20herd%20abundance-1.png)<!-- -->

``` r
ggsave(here::here("plots", "abundance.png"), width = 15, height = 11, bg = "white", dpi=600)

# #for Fuse
# library(Cairo)
# ggsave(here::here("plots", "abundance_forFUSE.svg"), width = 15, height = 11, bg = "transparent")
```

## Plot Abundance by Herd

## Plot Abundance across Herds

### Considered only 1991-2023, as most herds don’t have demographic data before 1991,pre-1991 trends for most SMC herds are informed solely by the priors and shared paramaters with other herds

``` r
#### Summarize bou pop in '91 vs 2023####
sims.summary <- sims.plot %>%
  group_by(.variable) %>%
  filter(yrs == 1991 | yrs == 2023) %>%
  arrange(yrs) %>%
  ungroup()

write_csv(sims.summary %>% mutate(.variable = case_when(.variable == "pred_totNMF" ~ "No action counterfactual", TRUE ~ "Actual")), here::here("tables/pop.size.csv"))

## what was the decline?
p.decline <- sims.summary %>%
  filter(.variable == "totNMF") %>%
  dplyr::select(yrs, mean) %>%
  pivot_wider(names_from = yrs, values_from = mean) %>%
  summarise(dif = ((`1991` - `2023`) / `1991`) * 100) %>%
  round(0)


## how many more caribou now than status quo w/ error
sims.draws <- out %>%
  gather_draws(
    pred_totNMF[i, j], totNMF[i, j], totCalvesMF[i, j], pred_totCalvesMF[i, j],
    ndraws = ndraws
  ) %>%
  left_join(hd %>% dplyr::select(herd, i = herd_num), by = "i") %>%
  filter(herd %in% herds.keep)

n.recovery.all <- sims.draws %>%
  ungroup() %>%
  filter(j == max(j), .variable %in% c("totNMF", "pred_totNMF")) %>%
  group_by(.draw, .variable) %>%
  summarise(across(.value, ~ sum(.x))) %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  mutate(dif = totNMF - pred_totNMF) %>%
  pull(dif)

## with error
median(n.recovery.all)
```

    ## [1] 1547.97

``` r
quantile(n.recovery.all, c(0.05, 0.5, 0.95)) %>% round(0)
```

    ##   5%  50%  95% 
    ## 1175 1548 1942

``` r
n.recovery <- median(n.recovery.all) %>% round(0)

write.csv(quantile(n.recovery.all, c(0.05, 0.5, 0.95)) %>% round(0) %>% as.data.frame(), here::here("tables/adults.recovered.csv"), row.names = TRUE)

## do again but just for calves
## unlike above, also keep all year so we can compare how many more calves were born over the entire period
calves.recovered.all <- sims.draws %>%
  filter(.variable %in% c("totCalvesMF", "pred_totCalvesMF")) %>%
  group_by(.draw, .variable, j) %>%
  summarise(across(.value, ~ sum(.x))) %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  mutate(dif = totCalvesMF - pred_totCalvesMF) %>%
  group_by(.draw) %>%
  summarise(dif = sum(dif)) %>%
  pull(dif)

quantile(calves.recovered.all, c(0.05, 0.5, 0.95)) %>% round(0)
```

    ##   5%  50%  95% 
    ## 1140 1552 1934

``` r
calves.recovered <- median(calves.recovered.all) %>% round(0)

sum(out$mean$totCalvesMF - out$mean$pred_totCalvesMF)
```

    ## [1] 1572.917

``` r
write.csv(quantile(calves.recovered.all, c(0.05, 0.5, 0.95)) %>% round(0) %>% as.data.frame(), here::here("tables/calves.recovered.csv"), row.names = TRUE)


#### Total Abundance####
abundance.all.plot <- ggplot(data = sims.plot %>%
  mutate(.variable = case_when(
    .variable == "totNMF" ~ "With recovery\nactions",
    TRUE ~ "Status quo"
  ))) +
  geom_ribbon(alpha = 0.3, aes(x = yrs, y = mean, ymin = LCL, ymax = UCL, fill = fct_relevel(.variable, "Status quo", "With recovery\nactions")), color = NA) +
  geom_line(size = 1, aes(x = yrs, y = mean, ymin = LCL, ymax = UCL, color = fct_relevel(.variable, "Status quo", "With recovery\nactions"), fill = fct_relevel(.variable, "Status quo", "With recovery\nactions"))) +
  geom_text(data = sims.plot %>% filter(yrs == 2023) %>%
    mutate(.variable = case_when(
      .variable == "totNMF" ~ "With recovery\nactions",
      TRUE ~ "Status quo"
    )), aes(label = fct_relevel(.variable, "Status quo", "With recovery\nactions"), colour = .variable, x = Inf, y = mean), hjust = 0, size = 4) +
  geom_jitter(data = ext.yr, size = 1, aes(x = yrs, y = mean), alpha = 0.5) +
  theme_ipsum() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("Year") +
  labs(
    x = "Year", y = "Abundance", title = "b) Population Trend",
    subtitle = paste0(p.decline, "% decline since 1991, +", n.recovery, " caribou from recovery actions")
  ) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 100), limits = c(0, 12000)) +
  scale_x_continuous(breaks = seq(1980, 2020, by = 10), limits = c(1991, 2023)) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 13),
    plot.subtitle = element_text(size = 15),
    legend.title = element_text(size = 15),
    plot.margin = unit(c(1, 5, 1, 1), "lines")
  ) +
  geom_text(data = trt %>% filter(applied %in% 1) %>% group_by(year) %>% summarise(n = n_distinct(herd)) %>% filter(year %% 2 == 1), aes(x = year, y = 200, label = n), size = 3.5) +
  scale_color_manual(values = cols[c(3, 1)]) +
  annotate(geom = "text", x = 1992, y = 2500, label = "Subpopulations w/\nrecovery actions", hjust = "left", size = 5) +
  annotate(
    geom = "curve", x = 1996, y = 1800, xend = 1998, yend = 500,
    curvature = 0, arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(geom = "text", x = 2010, y = 10000, label = "Subpopulations w/\ndemographic data", hjust = "left", size = 5) +
  annotate(
    geom = "curve", x = 2015, y = 10600, xend = 2013, yend = 11500,
    curvature = 0, arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(geom = "text", x = 1991, y = 7000, label = "Subpopulation\nextirpation event", hjust = "left", size = 5) +
  annotate(
    geom = "curve", x = 1998, y = 7800, xend = (ext.yr %>% ungroup() %>% filter(herd == "Banff") %>% pull(yrs)) + 0.4, yend = (ext.yr %>% ungroup() %>% filter(herd == "Banff") %>% pull(mean)) - 200,
    curvature = 0, arrow = arrow(length = unit(2, "mm"))
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  geom_text(data = raw.demog %>% group_by(year) %>% summarise(n = n_distinct(herd)) %>% filter(year %% 2 == 0), aes(x = year, y = 11800, label = n), size = 3.5)
abundance.all.plot
```

![](README_files/figure-gfm/Plot%20total%20abundance-1.png)<!-- -->

``` r
ggsave(plot = abundance.all.plot, here::here("plots", "abundance_all.png"), width = 6, height = 6, bg = "white", dpi=600)
```

# Summarize effects of treatments

## Population growth by treatment

``` r
## Gather draws
demog.draws <- out %>%
  gather_draws(logla[i, j], S[i, j], R_adj[i, j], totNMF[i, j], ndraws = ndraws) %>%
  mutate(
    .variable = case_when(
      .variable == "R_adj" ~ "R",
      .variable == "logla" ~ "r",
      TRUE ~ .variable
    )
  ) %>%
  ### pivot wider
  ungroup() %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  ## add in years and herd names
  left_join(yr_df %>% rename(j = yr_idx), by = "j") %>%
  left_join(trt %>% distinct(herd_num, herd) %>% dplyr::select(i = herd_num, herd), by = c("i")) %>%
  left_join(treatment.combos) %>%
  left_join(trt %>% distinct(herd, year, intensity) %>% dplyr::select(herd, yrs = year, intensity) %>% filter(intensity %in% "low"), by = c("herd", "yrs")) %>%
  replace_na(
    list(trt = "Reference")
  ) %>%
  ## remove first year lambda for each herd, as lambda==1
  dplyr::mutate(
    r = replace(r, yrs == 1973, NA_real_)
  ) %>%
  drop_na(r) %>%
  filter(herd %in% herds.keep)



## filter herds to a period where not extirpated
demog.draws.trim <- tibble()
for (i in 1:length(herds)) {
  if (herds[i] %in% ext.yr$herd) {
    yr <- ext.yr %>%
      dplyr::filter(herd == !!herds[i])
    a <- demog.draws %>%
      dplyr::filter(herd == !!herds[i] & yrs < yr$yrs)
  }

  if (!herds[i] %in% ext.yr$herd) {
    a <- demog.draws %>%
      dplyr::filter(herd == !!herds[i])
  }

  demog.draws.trim <- bind_rows(a, demog.draws.trim)
}

demog.draws <- demog.draws.trim
rm(demog.draws.trim)

## keep treatment years when applied to pops before functional extirpation
demog.draws <- demog.draws %>%
  filter(paste(i, j, sep = "_") %in% paste(demog.mod$i, demog.mod$j, sep = "_"))



### Remove years where transplant was being done (impossibly high lambda due to adding translocated individuals, not population response)
demog.draws <- demog.draws %>%
  filter(
    !(herd %in% "Charlotte Alplands" & yrs %in% c(1984:1991)),
    !(herd %in% "Telkwa" & yrs %in% c(1997:1999)),
    !(herd %in% "South Selkirks" & yrs %in% c(1987:1989)),
    !(herd %in% "Purcell South" & yrs %in% c(2012:2013))
  )

demog.draws %>% write_csv(here::here("tables", "draws", "demog.draws.csv"))

# ###how to calculate from draws
# lambda.compare <- demog%>%
#   ungroup%>%
#   dplyr::select(herd,yrs,lambda)%>%
#   left_join(demog.draws %>%
#               group_by(herd,yrs)%>%
#               summarise(median=median(lambda),
#                         mean=mean(lambda)))
#
# lm(lambda~median, data=lambda.compare)%>%summary
# lm(lambda~mean, data=lambda.compare)%>%summary
#
# ##definitely median when calculating for each posterior. Geo mean likely more appropriate for summarizing pooled distributions (i.e., overall means)

demog.draws.combotreat <- demog.draws %>%
  # filter(!(intensity%in%"low"|totNMF<20))%>%
  group_by(.draw, trt) %>%
  summarise(r = mean(r, na.rm = TRUE)) %>% ## mean per treatment-draw
  group_by(.draw, trt) %>%
  summarise(r = mean(r)) %>% ## mean lambda per treatment-draw
  mutate(trt = case_when(is.na(trt) ~ "Reference", TRUE ~ trt)) %>%
  mutate(group = case_when(trt %in% "Reference" ~ "Reference", TRUE ~ "Treatment"))

demog.draws.combotreat.rug <- demog.draws %>%
  # filter(!(intensity%in%"low"|totNMF<20))%>%
  group_by(trt, herd) %>%
  summarise(r = mean(r, na.rm = TRUE)) %>% ## mean per treatment-draw
  mutate(trt = case_when(is.na(trt) ~ "Reference", TRUE ~ trt)) %>%
  mutate(group = case_when(trt %in% "Reference" ~ "Reference", TRUE ~ "Treatment"))

order <- demog.draws.combotreat %>%
  group_by(trt) %>%
  summarise(med = median(r))

demog.draws.combotreat %>%
  left_join(order) %>%
  left_join(label.lookup, by = "trt") %>%
  filter(trt != "transplant") %>%
  ggplot(aes(x = r, y = fct_reorder(new, med), fill = group)) +
  geom_density_ridges( # scale = 1.5,
    # scale = 1.3,
    rel_min_height = .01,
    size = 0.15,
    alpha = 0.9
  ) +
  theme_ipsum() +
  geom_point(
    data = demog.draws.combotreat.rug %>%
      left_join(order) %>%
      left_join(label.lookup, by = "trt") %>%
      filter(trt != "transplant"),
    aes(y = fct_reorder(new, med), x = r),
    shape = "|",
    size=0.8
  ) +
  theme(
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    strip.text.x = element_text(size = 6),
    strip.text.y = element_text(size = 6),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.position = "none",
    plot.title = element_text(size = 7),
    panel.grid.minor=element_line(linewidth = 0.1),
    panel.grid.major=element_line(linewidth = 0.1),
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size=0.3, alpha=0.8) +
  labs(
    x = "Instantaneous rate of increase (r)",
    y = "Recovery action(s)",
    fill = "",
    title = "Population Growth by Recovery Action"
  ) +
  scale_fill_manual(values = cols[c(3, 1)]) +
  xlim(-0.25, 0.3)
```

![](README_files/figure-gfm/trt%20eff-%20r-1.png)<!-- -->

``` r
#ggsave(here::here("plots", "lambda_treatments.png"), width = 8, height = 7, bg = "white", dpi=600)
ggsave(here::here("plots", "final","Fig3.tiff"), width = 8.5, height = 8, bg = "white", units="cm", dpi=600)


lambda.table <- demog.draws.combotreat %>%
  group_by(trt) %>%
  summarise(
    r.med = median(r, na.rm = TRUE) %>% round(2),
    lower = quantile(r, 0.05, na.rm = TRUE) %>% round(2),
    upper = quantile(r, 0.95, na.rm = TRUE) %>% round(2)
  ) %>%
  arrange(-r.med) %>%
  mutate(r = paste0(r.med, " (", lower, "-", upper, ")"))

kable(lambda.table)
```

| trt                          | r.med | lower | upper | r                  |
|:-----------------------------|------:|------:|------:|:-------------------|
| feed                         |  0.12 | -0.29 |  0.46 | 0.12 (-0.29-0.46)  |
| feed-reducewolves            |  0.12 |  0.10 |  0.15 | 0.12 (0.1-0.15)    |
| reducemoose-reducewolves     |  0.11 |  0.06 |  0.17 | 0.11 (0.06-0.17)   |
| pen-reducewolves             |  0.10 |  0.05 |  0.14 | 0.1 (0.05-0.14)    |
| pen-reducemoose              |  0.06 | -0.09 |  0.19 | 0.06 (-0.09-0.19)  |
| reducewolves                 |  0.06 |  0.03 |  0.09 | 0.06 (0.03-0.09)   |
| reducewolves-sterilizewolves |  0.04 |  0.02 |  0.06 | 0.04 (0.02-0.06)   |
| pen-reducemoose-reducewolves | -0.02 | -0.20 |  0.18 | -0.02 (-0.2-0.18)  |
| Reference                    | -0.03 | -0.03 | -0.03 | -0.03 (-0.03–0.03) |
| transplant                   | -0.03 | -0.04 | -0.02 | -0.03 (-0.04–0.02) |
| reducemoose                  | -0.05 | -0.07 | -0.03 | -0.05 (-0.07–0.03) |

## Population growth Before-After

``` r
# Before-After ---------------------------------------------

## pull draws and organize into treatment and untreated (reference) timeframes for each herd
eff.draws <- demog.draws %>%
  ungroup() %>%
  filter(!trt %in% "Reference") %>%
  ungroup() %>%
  pivot_longer(r:R) %>%
  dplyr::select(.draw, trt, herd, yrs, name, eff = value) %>%
  ### add in reference
  left_join(
    demog.draws %>%
      ungroup() %>%
      filter(trt %in% "Reference") %>%
      group_by(.draw, herd) %>%
      summarise(across(r:R, ~ mean(.x, na.rm = TRUE))) %>%
      pivot_longer(r:R) %>%
      dplyr::select(.draw, herd, name, ref = value),
    by = c(".draw", "herd", "name")
  ) %>%
  ## calculate delta
  mutate(
    delta.r = eff - ref
  ) %>%
  ## add a label that includes sample sizes for plot
  left_join(trt.n, by = "trt") %>%
  left_join(label.lookup, by = "trt") %>%
  mutate(trt.label = paste0(new, "\n", n_herds, " subpops, ", n_yrs, " yrs"))



order <- eff.draws %>%
  filter(name == "r") %>%
  group_by(trt.label) %>%
  summarise(med = median(delta.r)) %>%
  arrange(-med)

eff.draws <- eff.draws %>%
  left_join(order, by = "trt.label") %>%
  mutate(
    trt.label = fct_reorder(trt.label, med),
    name = case_when(
      name == "r" ~ "Rate of increase (r)",
      name == "R" ~ "Recruitment (R)",
      name == "S" ~ "Survival (S)",
      TRUE ~ name
    )
  )


ggplot() +
  geom_density_ridges(
    data = eff.draws %>%
      filter(trt != "transplant") %>%
      group_by(name, trt.label, .draw) %>%
      summarise(delta = median(delta.r)) %>%
      ## trim for plot, removes <1% of values,
      filter(delta > -0.25 & delta < 0.35),
    aes(x = delta, y = trt.label, fill = name),
    scale = .9,
    rel_min_height = .01,
    size = 0.15,
    alpha = 0.9
  ) +
  geom_point(
    data = eff.draws %>%
      filter(trt != "transplant") %>%
      group_by(herd, name, trt.label) %>%
      summarise(delta = median(delta.r)),
    aes(y = trt.label, x = delta),
    shape = "|"
  ) +
  theme_ipsum() +
  facet_wrap(vars(name)) +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
     plot.title = element_text(size = 12),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Change in value", y = "Recovery action(s)", title = "Before-After Assessment of Effectivness") +
  scale_fill_manual(values = cols[c(1:3)])
```

![](README_files/figure-gfm/trt%20eff-%20BA-1.png)<!-- -->

``` r
#ggsave(here::here("plots", "ba_all.png"), width = 10, height = 7, bg = "white", dpi=600)
ggsave(here::here("plots", "final","Fig4.tiff"), width = 18, height = 13, bg = "white", units="cm", dpi=600)


trt_eff_ba_table <- eff.draws %>%
  filter(name == "Rate of increase (r)") %>%
  group_by(new, .draw) %>%
  summarise(delta = median(delta.r)) %>% ## get mean effect for each treatment across herds for each draw
  group_by(new) %>%
  summarise(
    delta.l = median(delta, na.rm = TRUE) %>% round(2), ## mean effect of each treatment
    lower = quantile(delta, 0.05, na.rm = TRUE) %>% round(2),
    upper = quantile(delta, 0.95, na.rm = TRUE) %>% round(2)
  ) %>%
  arrange(-delta.l) %>%
  mutate(delta.r = paste0(delta.l, " [", lower, "-", upper, "]")) %>%
  filter(new != "transplant") %>%
  dplyr::select(`Recovery action` = new, `Change in instantaneous growth rate (r)` = delta.r)


trt_eff_ba_table %>%
  gt() %>%
  gtsave(here::here("tables", "trt_eff_ba.rtf"))

trt_eff_ba_table %>%
  write_csv(here::here("tables", "trt_eff_ba.csv"))

kable(trt_eff_ba_table)
```

| Recovery action                  | Change in instantaneous growth rate (r) |
|:---------------------------------|:----------------------------------------|
| penning + wolf reduction         | 0.16 \[0.12-0.2\]                       |
| feeding                          | 0.15 \[-0.27-0.49\]                     |
| feeding + wolf reduction         | 0.14 \[0.09-0.2\]                       |
| moose & wolf reduction           | 0.11 \[0.04-0.18\]                      |
| wolf reduction                   | 0.08 \[0.02-0.13\]                      |
| penning + moose reduction        | 0.07 \[-0.17-0.3\]                      |
| wolf reduction + sterilization   | 0.07 \[0.01-0.13\]                      |
| penning + moose & wolf reduction | -0.01 \[-0.2-0.19\]                     |
| moose reduction                  | -0.04 \[-0.1-0.01\]                     |

``` r
eff.draws %>%
  filter(name == "Rate of increase (r)") %>%
  group_by(herd, name, trt) %>%
  summarise(delta = median(delta.r)) %>%
  ungroup() %>%
  arrange(delta) %>%
  count(delta > 0) %>%
  kable()
```

| delta \> 0 |   n |
|:-----------|----:|
| FALSE      |   5 |
| TRUE       |  26 |

## Summarize results by recovery action, including application intensity

``` r
eff.draws.app <- demog.draws %>%
  group_by(herd, yrs, trt) %>%
  mutate(totNMF.median = median(totNMF)) %>% # get average pop size so popsize threshold doesnt split low/standard in some years due to draws being above/below threshold
  ungroup() %>%
  mutate(
    application = case_when(intensity == "low" | totNMF.median < 30 ~ "low", TRUE ~ "standard")
  ) %>%
  filter(!trt %in% "Reference") %>%
  group_by(.draw, herd, trt, application) %>%
  summarise(across(r:R, ~ mean(.x, na.rm = TRUE))) %>% ### mean posterior per herd-treatment
  ungroup() %>%
  pivot_longer(r:R) %>%
  ungroup() %>%
  dplyr::select(.draw, application, trt, herd, name, eff = value) %>%
  left_join(
    demog.draws %>%
      ungroup() %>%
      filter(trt %in% "Reference") %>%
      group_by(.draw, herd) %>%
      summarise(across(r:R, ~ mean(.x, na.rm = TRUE))) %>% ### mean posterior per herd-treatment
      pivot_longer(r:R) %>%
      dplyr::select(.draw, herd, name, ref = value),
    by = c(".draw", "herd", "name")
  ) %>%
  ## calculate delta
  mutate(
    delta.r = eff - ref
  )

# eff.draws.app %>%
#   filter(
#     trt != "transplant",
#     name == "r"
#   ) %>%
#   filter(trt %in% "reducewolves") %>%
#   group_by(herd, trt, application) %>%
#   summarize(mean = mean(eff)) %>%
#   arrange(mean)

eff.draws.app %>%
  filter(
    trt != "transplant",
    name == "r"
  ) %>%
  filter(trt %in% "reducewolves") %>%
  group_by(trt, application, herd) %>%
  summarize(mean = mean(eff)) %>%
  group_by(trt, application) %>%
  summarize(
    n = n(),
    mean = median(mean)
  ) %>%
  arrange(mean)

demog.draws %>%
  filter(trt %in% "reducewolves") %>%
  group_by(herd, yrs, intensity, trt) %>%
  summarize(mean = mean(r)) %>%
  group_by(trt, intensity) %>%
  summarize(
    n = n(),
    mean = median(mean)
  )

eff.draws.app %>% write_csv(here::here("tables", "draws", "eff.draws.app.csv"))

## model individual effects by application intensity

eff.draws.app.model <- eff.draws.app %>%
  mutate(
    reducewolves = case_when(str_detect(trt, "reducewolves") ~ 1, TRUE ~ 0),
    sterilizewolves = case_when(str_detect(trt, "sterilizewolves") ~ 1, TRUE ~ 0),
    reducemoose = case_when(str_detect(trt, "reducemoose") ~ 1, TRUE ~ 0),
    pen = case_when(str_detect(trt, "pen") ~ 1, TRUE ~ 0),
    feed = case_when(str_detect(trt, "feed") ~ 1, TRUE ~ 0),
    transplant = case_when(str_detect(trt, "transplant") ~ 1, TRUE ~ 0)
  ) %>%
  filter(
    name == "r",
    application == "standard"
  )

ind.eff.app <- eff.draws.app.model %>%
  filter(name == "r") %>%
  group_by(.draw) %>%
  do(tidy(lm(delta.r ~ reducewolves + sterilizewolves + reducemoose + pen + feed, data = .)))

ind.eff.app <- ind.eff.app %>%
  filter(!term %in% c("(Intercept)")) %>%
  left_join(
    ind.eff.app %>%
      filter(term == "(Intercept)") %>%
      select(.draw, intercept = estimate),
    by = ".draw"
  ) %>%
  mutate(eff = intercept + estimate)

ind.eff.app %>%
  group_by(.draw) %>%
  group_by(term) %>%
  summarise(eff = median(eff))

ind.eff.app %>%
  group_by(term) %>%
  summarise(
    delta.l = median(eff, na.rm = TRUE) %>% round(2),
    lower = quantile(eff, 0.05, na.rm = TRUE) %>% round(2),
    upper = quantile(eff, 0.95, na.rm = TRUE) %>% round(2)
  ) %>%
  arrange(-delta.l) %>%
  mutate(delta.lambda = paste0(delta.l, " (", lower, "-", upper, ")")) %>%
  dplyr::select(Treatment = term, delta.lambda)

eff.draws.app.model %>% write_csv(here::here("tables", "draws", "eff.draws.app.model.csv"))

## do again but with all data included
eff.draws <- eff.draws %>%
  mutate(
    reducewolves = case_when(str_detect(trt, "reducewolves") ~ 1, TRUE ~ 0),
    sterilizewolves = case_when(str_detect(trt, "sterilizewolves") ~ 1, TRUE ~ 0),
    reducemoose = case_when(str_detect(trt, "reducemoose") ~ 1, TRUE ~ 0),
    pen = case_when(str_detect(trt, "pen") ~ 1, TRUE ~ 0),
    feed = case_when(str_detect(trt, "feed") ~ 1, TRUE ~ 0),
    transplant = case_when(str_detect(trt, "transplant") ~ 1, TRUE ~ 0)
  )

ind.eff <- eff.draws %>%
  filter(name == "Rate of increase (r)") %>%
  group_by(.draw) %>%
  do(tidy(lm(delta.r ~ reducewolves + sterilizewolves + reducemoose + pen + feed, data = .)))


ind.eff <- ind.eff %>%
  filter(!term %in% c("(Intercept)")) %>%
  left_join(
    ind.eff %>%
      filter(term == "(Intercept)") %>%
      select(.draw, intercept = estimate),
    by = ".draw"
  ) %>%
  mutate(eff = intercept + estimate)

ind.eff %>%
  group_by(term) %>%
  summarise(eff = median(estimate))

eff.draws %>% write_csv(here::here("tables", "draws", "eff.draws.csv"))
```

## Individual Treatment Effects

``` r
ind.eff.plot <- ggplot(ind.eff.app %>% left_join(label.lookup %>% rename(term = trt), by = "term"), aes(x = eff, y = fct_relevel(new, "wolf sterlization", "moose reduction", "feeding", "penning", "wolf reduction"), fill = new)) +
  geom_density_ridges(
    scale = .9,
    rel_min_height = .01,
    size = 0.25,
    alpha = 0.5
  ) +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Change in rate of increase", y = "Recovery action(s)", title = "a) Individual Treatment Effects", subtitle = "Partitioned using regression analysis, assuming effects are additive") +
  scale_fill_manual(values = cols[c(1:6)]) +
  xlim(-0.2, 0.25)
ind.eff.plot
```

![](README_files/figure-gfm/ind%20trt%20eff-1.png)<!-- -->

``` r
# ggsave(here::here("plots","ind_effects.png"), width=5, height=6, bg="white")


ind.eff.table <- ind.eff.app %>%
  group_by(term) %>%
  summarise(
    delta.l = median(eff, na.rm = TRUE) %>% round(2),
    lower = quantile(eff, 0.05, na.rm = TRUE) %>% round(2),
    upper = quantile(eff, 0.95, na.rm = TRUE) %>% round(2)
  ) %>%
  arrange(-delta.l) %>%
  mutate(delta.lambda = paste0(delta.l, " (", lower, "-", upper, ")")) %>%
  dplyr::select(Treatment = term, delta.lambda)

ind.eff.table %>%
  gt() %>%
  gtsave(here::here("tables", "ind_trt_eff_ba.rtf"))

kable(ind.eff.table)
```

| Treatment       | delta.lambda      |
|:----------------|:------------------|
| feed            | 0.11 (-0.13-0.31) |
| reducewolves    | 0.1 (0.03-0.17)   |
| pen             | 0.08 (0-0.17)     |
| reducemoose     | 0.03 (-0.03-0.08) |
| sterilizewolves | -0.01 (-0.1-0.07) |

## Simulate Conservation Intervention

``` r
n.sims <- 1000
start.pop <- 100

sim.ref <- demog.draws %>%
  filter(trt == "Reference" & totNMF < 150) %>%
  group_by(.draw, trt) %>%
  summarise(r = median(r)) %>% ## median lambda across herds
  ungroup() %>%
  sample_n(n.sims) %>%
  pull(r)

median(sim.ref)
```

    ## [1] -0.07485089

``` r
sim.trt <- ind.eff.app %>%
  group_by(term) %>%
  ungroup() %>%
  dplyr::select(.draw, term, eff) %>%
  pivot_wider(names_from = "term", values_from = "eff") %>%
  drop_na() %>%
  sample_n(n.sims)

median(sim.ref)
```

    ## [1] -0.07485089

``` r
year.end <- 9
sim.df <- list()
for (i in 1:n.sims) {
  ## effects
  wolf.sim.i <- exp(sim.ref[i] + sim.trt$reducewolves[i])
  feed.sim.i <- exp((sim.ref[i] + sim.trt$feed[i]))
  pen.sim.i <- exp((sim.ref[i] + sim.trt$pen[i]))
  moose.sim.i <- exp((sim.ref[i] + sim.trt$reducemoose[i]))
  steril.sim.i <- exp((sim.ref[i] + sim.trt$sterilizewolves[i]))
  wolffeed.sim.i <- exp((sim.ref[i] + sim.trt$feed[i] + sim.trt$reducewolves[i]))
  wolfpen.sim.i <- exp((sim.ref[i] + sim.trt$pen[i] + sim.trt$reducewolves[i]))

  ## sub out values that exceed the maximal sustained growth rate observed for caribou is r=0.26 (Heard and Oullet 1994)
  max.l <- exp(0.26)

  wolf.sim.i <- ifelse(wolf.sim.i > max.l, max.l, wolf.sim.i)
  feed.sim.i <- ifelse(feed.sim.i > max.l, max.l, feed.sim.i)
  pen.sim.i <- ifelse(pen.sim.i > max.l, max.l, pen.sim.i)
  moose.sim.i <- ifelse(moose.sim.i > max.l, max.l, moose.sim.i)
  steril.sim.i <- ifelse(steril.sim.i > max.l, max.l, steril.sim.i)
  wolffeed.sim.i <- ifelse(wolffeed.sim.i > max.l, max.l, wolffeed.sim.i)
  wolfpen.sim.i <- ifelse(wolfpen.sim.i > max.l, max.l, wolfpen.sim.i)
  ## project population change
  cons.sim.i <- tibble(yr = 1:(11 + year.end)) %>%
    dplyr::mutate(
      reference = c(start.pop * exp(sim.ref[i])^yr),
      reducewolves = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * wolf.sim.i^(yr - 10))
      ),
      feed = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * feed.sim.i^(yr - 10))
      ),
      pen = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * pen.sim.i^(yr - 10))
      ),
      reducemoose = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * moose.sim.i^(yr - 10))
      ),
      sterilizewolves = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * steril.sim.i^(yr - 10))
      ),
      `reducewolves+feed` = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * wolffeed.sim.i^(yr - 10))
      ),
      `reducewolves+pen` = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * wolfpen.sim.i^(yr - 10))
      ),
      sim = i
    )

  sim.df[[i]] <- cons.sim.i
}

sim.df <- do.call(rbind, sim.df)

nudge <- 0.15
sim.df.plot <- sim.df %>%
  pivot_longer(reference:`reducewolves+pen`) %>%
  group_by(name, sim) %>%
  mutate(
    inc = sum(last(value) > value[10]),
    ext = sum(last(value) <= 20),
    abund = last(value)
  ) %>%
  group_by(name, yr) %>%
  summarise(
    median = median(value),
    se = sd(value),
    upper = quantile(abund, 0.95),
    lower = quantile(abund, 0.05),
    inc = (mean(inc) * 100) %>% round(0),
    ext = (mean(ext) * 100) %>% round(0)
  ) %>%
  filter(
    !(yr < 10 & name != "reference"),
    yr <= (11 + year.end)
  ) %>%
  left_join(label.lookup %>% rename(name = trt), by = "name") %>%
  mutate(
    trt = paste0(new, " (", inc, ", ", ext, ", ", lower %>% round(0), "-", upper %>% round(0), ")"),
    yr = yr - 10
  ) %>%
  mutate(yr.shift = case_when(
    new %in% "sterilizewolves" & yr == year.end ~ yr + (2.5 * nudge),
    new %in% "reducemoose" & yr == year.end ~ yr + (3.5 * nudge),
    new %in% "reference" & yr == year.end ~ yr - nudge,
    new %in% "pen" & yr == year.end ~ yr,
    new %in% "feed" & yr == year.end ~ yr + nudge,
    new %in% "reducewolves" & yr == year.end ~ yr + (2 * nudge),
    new %in% "reducewolves+pen" & yr == year.end ~ yr + (3 * nudge),
    new %in% "reducewolves+feed" & yr == year.end ~ yr + (4 * nudge),
    TRUE ~ yr
  )) %>%
  mutate(trt = case_when(
    str_detect(trt, "\\+ wolf reduction") ~ str_wrap(trt, width = 26),
    TRUE ~ trt
  ))

# a <-sim.df.plot%>%
#   filter(yr == last(yr))%>%pull()

recov.cols <- c(cols[1:3], "black", cols[5:8])
recov.sims.plot <- ggplot() +
  annotate("rect",
    xmin = -10, xmax = year.end + 1, ymin = 0, ymax = 20,
    alpha = .1, fill = "black"
  ) +
  annotate("text", x = -4, y = 10, label = "Functionally extirpated") +
  geom_line(data = sim.df.plot %>% filter(name != "reference"), aes(x = yr, y = median, color = trt)) +
  geom_line(data = sim.df.plot %>% filter(name == "reference"), aes(x = yr, y = median, color = trt), size = 1.5, linetype = "dashed") +
  # geom_linerange(data=sim.df.plot%>%
  #                  group_by(name)%>%
  #                  filter(yr==max(yr)),aes(x=yr.shift,y=median,color=trt, ymin=median-se,ymax=median+se), alpha=0.7)+
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Years since intervention", y = "Abundance", title = "b) Simulated Options to Avert Caribou Extirpation", subtitle = "Labels = treatment (% samples increased, % samples extirpated, 90% end abundance interval)") +
  # geom_hline(yintercept = 10, linetype="dashed")+
  geom_text_repel(
    data = sim.df.plot %>%
      filter(
        yr == last(yr),
        name == "feed"
      ),
    aes(color = trt, label = trt, x = yr, y = median),
    size = 4,
    direction = "y",
    xlim = c(year.end + 2, 35),
    vjust = 2.5,
    segment.size = .7,
    segment.alpha = .5,
    segment.linetype = "dotted",
    box.padding = .8,
    seed = 999
  ) +
  geom_text_repel(
    data = sim.df.plot %>%
      filter(
        yr == last(yr),
        name != "feed"
      ),
    aes(color = trt, label = trt, x = yr, y = median),
    size = 4,
    direction = "y",
    xlim = c(year.end + 2, 35),
    hjust = 0,
    segment.size = .7,
    segment.alpha = .5,
    segment.linetype = "dotted",
    box.padding = .8,
    seed = 999
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(-10, year.end + 1),
    breaks = seq(-10, year.end + 1, by = 5)
  ) +
  theme(plot.margin = margin(1, 7, 1, 1.2, "cm")) +
  scale_color_manual(values = recov.cols) +
  geom_text(data = tibble(lab = "Status quo", x = -3, y = 90), aes(x = x, y = y, label = "Status quo"), hjust = "left", size = 5) +
  geom_curve(
    data = tibble(x = -2, y = 80, xend = -3, yend = 60),
    aes(x = x, y = y, yend = yend, xend = xend), inherit.aes = FALSE, curvature = -.3, arrow = arrow(length = unit(2, "mm"))
  )
recov.sims.plot
```

![](README_files/figure-gfm/cons%20sims-1.png)<!-- -->

``` r
## plot individual treatments together
recov.together <- ind.eff.plot + recov.sims.plot + plot_layout(widths = c(1.3, 1.6))
ggsave(plot = recov.together, here::here("plots/recov.together.png"), width = 13, height = 6, bg = "white", dpi=600)
```

## Map

``` r
## Prep Herd Bounds
herd.bounds <- st_read(here::here("data/Spatial/herds/u_bc_herds_2021_CL.shp")) %>%
  st_transform(3005) %>%
  dplyr::select(herd = HERD_NAME) %>%
  mutate(herd = case_when(
    herd %in% "Narraway" ~ "Bearhole Redwillow",
    herd %in% "Moberly" ~ "Klinse-Za",
    herd %in% "Scott" ~ "Scott West",
    herd %in% "Frisby Boulder" ~ "Frisby-Boulder",
    herd %in% "Purcell Central" ~ "Purcells Central",
    TRUE ~ herd
  )) %>%
  rbind(st_read(here::here("data/Spatial/herds/Caribou_Range.shp")) %>%
    st_transform(3005) %>%
    dplyr::select(herd = SUBUNIT) %>%
    mutate(herd = case_when(
      herd %in% "Narraway" ~ "Narraway AB",
      herd %in% "Redrock-Prairie Creek" ~ "Redrock/Prairie Creek",
      TRUE ~ herd
    ))) %>%
  rbind(st_read(here::here("data/Spatial/herds/rrpc/BC_RRPC.shp")) %>%
    st_transform(3005) %>%
    mutate(herd = "Redrock/Prairie Creek") %>%
    dplyr::select(herd)) %>%
  st_simplify(
    preserveTopology = FALSE,
    dTolerance = 1000
  ) %>%
  group_by(herd) %>%
  summarise() %>%
  ungroup() %>%
  filter(herd %in% c(herds, "Scott West"))
```

    ## Reading layer `u_bc_herds_2021_CL' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/Work/WSC/CaribouIPM_BCAB/data/Spatial/herds/u_bc_herds_2021_CL.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 57 features and 19 fields
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: -165343.7 ymin: 5422045 xmax: 1031821 ymax: 6709569
    ## Projected CRS: NAD83 / UTM zone 10N
    ## Reading layer `Caribou_Range' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/Work/WSC/CaribouIPM_BCAB/data/Spatial/herds/Caribou_Range.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 25 features and 9 fields
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 170844.4 ymin: 5689840 xmax: 819119.1 ymax: 6659319
    ## Projected CRS: NAD83 / Alberta 10-TM (Forest)
    ## Reading layer `BC_RRPC' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/Work/WSC/CaribouIPM_BCAB/data/Spatial/herds/rrpc/BC_RRPC.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 1 feature and 17 fields
    ## Geometry type: POLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 1345869 ymin: 959823.3 xmax: 1416126 ymax: 1036210
    ## Projected CRS: NAD83 / BC Albers

``` r
library(units)
disturb <- st_read(here::here("data/Spatial/disturbance/from_Emily/2022-01-19/anthrodisturbance_75to15incl0/anthrodisturbance_75to15incl0.shp")) %>%
  st_intersection(herd.bounds) %>%
  mutate(dist = st_area(.) %>% set_units("km^2") %>% as.numeric()) %>%
  tibble() %>%
  dplyr::select(herd, dist)
```

    ## Reading layer `anthrodisturbance_75to15incl0' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/Work/WSC/CaribouIPM_BCAB/data/Spatial/disturbance/from_Emily/2022-01-19/anthrodisturbance_75to15incl0/anthrodisturbance_75to15incl0.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 1 feature and 1 field
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 818361.6 ymin: 474487.5 xmax: 1806688 ymax: 1392210
    ## Projected CRS: NAD83 / BC Albers

``` r
herd.bounds <- herd.bounds %>%
  mutate(herd.area = st_area(.) %>% set_units("km^2") %>% as.numeric()) %>%
  left_join(disturb, by = "herd") %>%
  mutate(human = (dist / herd.area)) %>%
  dplyr::select(herd, human) %>%
  left_join(
    demog %>%
      filter(trt %in% c("Reference", "transplant")) %>%
      group_by(herd) %>%
      filter(yrs %in% (max(yrs) - 9):max(yrs)) %>% ## filter to 10 years before intervention started
      summarise(r = mean(log(lambda), na.rm = TRUE)), ## geo mean per herd
    by = "herd"
  ) %>%
  left_join(ext.yr %>% mutate(ext = 1) %>% ungroup() %>% dplyr::select(herd, ext) %>% rbind(tibble(herd = "Scott West", ext = 1)),
    by = "herd"
  ) %>%
  left_join(ecotype,
    by = "herd"
  ) %>%
  mutate(
    r = case_when(herd %in% "Scott West" ~ -0.05, TRUE ~ r),
    ECCC = case_when(herd %in% "Scott West" ~ "Central Group", TRUE ~ ECCC)
  ) %>% ## don't have Scott west. Extirpated
  mutate(
    human = round(human * 100, 0) %>% as.integer(),
    r.class = case_when(
      r <= -0.05 ~ "< -0.05",
      r > -0.05 & r <= -0.01 ~ "-0.05--0.01",
      r > -0.01 & r < 0.01 ~ "-0.01-0.01",
      r >= 0.01 ~ ">0.01"
    ),
    r.class2 = case_when(
      r <= -0.01 ~ "declining",
      r > -0.01 & r <= 0.01 ~ "stable",
      r > 0.01 ~ "increasing",
      is.na(r) & ext == 1 ~ "declining"
    )
  ) %>%
  st_as_sf() %>%
  st_transform(4326)

st_write(herd.bounds, here::here("data/Spatial/herds/ipm_herds.shp"), delete_dsn = TRUE)
```

    ## Deleting source `/Users/claytonlamb/Dropbox/Documents/University/Work/WSC/CaribouIPM_BCAB/data/Spatial/herds/ipm_herds.shp' using driver `ESRI Shapefile'
    ## Writing layer `ipm_herds' to data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/Work/WSC/CaribouIPM_BCAB/data/Spatial/herds/ipm_herds.shp' using driver `ESRI Shapefile'
    ## Writing 41 features with 9 fields and geometry type Polygon.

``` r
# mapview(herds,zcol="human")
# mapview(herds,zcol="lambda.class2")

## prep herd # labels
# herd.bounds %>%
#   cbind(st_centroid(.) %>% st_coordinates()) %>%
#   tibble() %>%
#   mutate(ECCC=fct_relevel(ECCC,"Northern Group","Central Group","Southern Group"))%>%
#   arrange(ECCC,-Y) %>%
#   mutate(number_label = 1:n()) %>%
#   dplyr::select (herd, human, number_label) %>%
#   write_csv(here::here("data", "clean", "labels.csv"))

herd.bounds <- herd.bounds %>%
  left_join(labels %>% dplyr::select(herd, number_label), by = "herd")

## Prep inset map
## Cities
cities <- st_read(here::here("data/Spatial/administrative/places.shp")) %>%
  filter(NAME %in% c("C", "VANCOUVER", "Prince George", "Fort St John")) %>%
  mutate(Name = str_to_title(NAME)) %>%
  dplyr::select(Name, geometry) %>%
  st_transform(4326) %>%
  rbind(tibble(Name = "Banff", Y = 51.18, X = -115.56) %>% st_as_sf(coords = c("X", "Y"), crs = 4326)) %>%
  st_transform(3005)
```

    ## Reading layer `places' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/Work/WSC/CaribouIPM_BCAB/data/Spatial/administrative/places.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 787 features and 3 fields
    ## Geometry type: POINT
    ## Dimension:     XY
    ## Bounding box:  xmin: -134.9996 ymin: 48.37586 xmax: -114.6224 ymax: 59.92801
    ## Geodetic CRS:  NAD27

``` r
## Pacific northwest for inset
pnw <- st_read(here::here("data/Spatial/administrative/North_America.shp")) %>%
  filter(FID_usa %in% c(1, 2, 8, 11) | FID_canada %in% c(8, 9)) %>%
  mutate(
    prov = c("WA", "MT", "ID", "OR", "AB", "BC", "BC"),
    country = c(rep("USA", times = 4), rep("CAN", times = 3))
  ) %>%
  st_make_valid() %>%
  group_by(prov, country) %>%
  summarise(id = mean(FID_canada)) %>%
  ungroup() %>%
  st_transform(3005)
```

    ## Reading layer `North_America' from data source 
    ##   `/Users/claytonlamb/Dropbox/Documents/University/Work/WSC/CaribouIPM_BCAB/data/Spatial/administrative/North_America.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 70 features and 2 fields
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: -178.2176 ymin: 18.92179 xmax: -52.62783 ymax: 83.11506
    ## Geodetic CRS:  WGS 84

``` r
library(ggspatial)
## inset
## custom crs to keep things straight
cust.crs <- "+proj=aea +lat_0=50 +lon_0=-114.9 +lat_1=49 +lat_2=50.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs"
inset <- ggplot() +
  geom_sf(data = pnw %>% st_transform(cust.crs), fill = "grey80", color = NA) +
  geom_sf(data = pnw %>% st_transform(cust.crs), size = 0.5, fill = NA, color = "grey30") +
  layer_spatial(sf::st_bbox(herd.bounds %>% st_buffer(100000) %>% st_transform(cust.crs)), fill = NA, linetype = "dashed", color = "grey99", linewidth = 1) +
  geom_sf_text(data = pnw %>% st_transform(cust.crs), aes(label = prov), inherit.aes = FALSE, size = 3) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = "transparent"),
    panel.border = element_rect(fill = NA, color = NA),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    legend.position = c(0.65, 0.075),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )




## get basemap
#
# register_google("Add your token here")
# bmap.big <- basemaps::basemap(
#   #ext = herd.bounds %>% st_buffer(200000) %>% group_by%>%summarise%>%st_transform(3857)%>%st_bbox,
#   ext=st_bbox(c(xmin = -15443769, xmax = -12714182, ymin = 6085017, ymax = 7962966), crs = st_crs(3857)),
#   map_res = 1, map_type = "terrain_bg", class="raster"
# ) %>% projectRaster(crs = cust.crs)
#
# writeRaster(bmap.big, here::here("data","Spatial","basemap.tif"))

bmap.big <- rast(here::here("data", "Spatial", "basemap.tif"))

bmap.big <- bmap.big %>% as.data.frame(xy = TRUE)



map <- ggplot() +
  geom_raster(
    data = bmap.big, aes(x = x, y = y),
    fill = rgb(
      red = bmap.big$basemap_1,
      green = bmap.big$basemap_2,
      blue = bmap.big$basemap_3,
      maxColorValue = 255
    ),
    show.legend = FALSE
  ) +
  theme_bw() +
  geom_sf(data = pnw %>%
    st_transform(cust.crs), size = 1, fill = NA, linetype = "dashed") +
  geom_sf(data = herd.bounds, aes(fill = r.class2), inherit.aes = FALSE, alpha = 0.7) +
  geom_sf(data = herd.bounds %>%
    filter(ext %in% 1), aes(color = "fnl extirpation"), fill = NA, inherit.aes = FALSE, alpha = 0.7, linewidth = 0.75) +
  geom_sf_text(data = st_centroid(herd.bounds), aes(label = number_label), inherit.aes = FALSE, size = 3, color = "white") +
  geom_sf_text(data = st_centroid(herd.bounds %>% filter(r.class2 == "stable")), aes(label = number_label), inherit.aes = FALSE, size = 2.5, color = "black") +
  geom_sf(data = cities %>% st_transform(cust.crs), inherit.aes = FALSE, size = 3, pch = 21, fill = "white", color = "black") +
  geom_sf_text(data = cities %>% filter(!Name %in% c("Fort St John", "Prince George")) %>% st_transform(cust.crs), aes(label = Name), inherit.aes = FALSE, size = 3, hjust = -0.1, vjust = -1, color = "white") +
  geom_sf_text(data = cities %>% filter(Name %in% c("Fort St John", "Prince George")) %>% st_transform(cust.crs), aes(label = Name), inherit.aes = FALSE, size = 3, hjust = 0.6, vjust = -1, color = "white") +
  theme_ipsum() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 15, color = "black"),
    legend.title = element_text(size = 15, color = "black"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.position = "top"
  ) +
  ggsn::scalebar(x.min = 10E4, x.max = 105E4, y.min = -19E4, y.max = 85E4, dist = 150, height = 0.03, dist_unit = "km", transform = FALSE, location = "bottomleft", st.color = "white", st.bottom = FALSE) +
  annotation_custom(ggplotGrob(inset), xmin = 65E4, xmax = 108E4, ymin = 52E4, ymax = 85E4) +
  scale_y_continuous(expand = c(0, 0), limits = c(-20E4, 85E4)) +
  scale_x_continuous(expand = c(0, 0), limits = c(5E4, 105E4)) +
  labs(fill = "Population growth\nwithout intervention", title = "a) Southern Mountain Caribou", color = "") +
  scale_fill_viridis_d() +
  scale_color_manual(values = c("fnl extirpation" = "red")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))


## plot together
map_combo <- map + abundance.all.plot + plot_layout(widths = c(1.7, 1.3))
map_combo
```

![](README_files/figure-gfm/map-1.png)<!-- -->

``` r
ggsave(plot = map_combo, here::here("plots/map.together.png"), width = 12, height = 8, bg = "white", dpi=600)
```
