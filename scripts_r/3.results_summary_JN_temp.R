## ----render, eval=FALSE,include=FALSE----------------------------------------------------------------------------------------------------------------------
## rmarkdown::render(here::here('CaribouIPM_BCAB.Rmd'),
##                   output_file =  "README.md")
##
## knitr::purl(input=here::here('CaribouIPM_BCAB.Rmd'),
##                   output =  here::here("scripts_r", '3.results_summary.r'))


## ----Load packages and data, results='hide', message=FALSE, warning=FALSE----------------------------------------------------------------------------------
library(ggmap)
library(raster)
library(RStoolbox)
library(ggsn)
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
library(ggeffects)
library(ggridges)
library(ggrepel)
library(tidylog)
library(gt)
library(patchwork)
library(sf)
library(mapview)
library(basemaps)
library(ggtext)
library(rlang)
library(tidyverse)

# Load data ---------------------------------------------------------------

## IPM Output
out <- readRDS(file = here::here("jags/output/BCAB_CaribouIPM_posteriors.rds"))

# Had to change to a local copy since the posteriors are not part of the GitHub
#  repo. I would generally suggest data be stored in Google Drive or similar.
#  That would provide a single source of truth (assuming no copies are made) and
#  version control the file while being accessible from anywhere.

# In the future using rjags returns a very simple data structure without
# redundancy and can really help reduce the memory load, but it does come at
# the cost of not having as much support from packages
# out <- readRDS(
#   file = "/Users/jamesnowak/Downloads/BCAB_CaribouIPM_posteriors.rds"
# )

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

#  Years of study
yrs <- seq(from = min(trt$year), to = max(trt$year), by = 1)
nyr <- length(yrs)
yr_idx <- seq(from = 1, to = nyr, by = 1)

# This seems like too many operations to get to the solution and a loss of
#  control of the output.
yr_df <- as.data.frame(cbind(yrs, yr_idx))


## ----Data housekeeping, results='hide', message=FALSE, warning=FALSE---------------------------------------------------------------------------------------
## set up colors for plotting
display.brewer.pal(8, "Accent")
cols <- RColorBrewer::brewer.pal(8, "Accent")

### Change Narraway BC to Bearhole Redwillow
# This was hard to read, maybe the styler Addin in RStudio if the tidyverse style
#  guide seems like a good fit
# Curious why %in% if only matching a single value
trt <- trt %>% mutate(
  herd = case_when(
    herd %in% "Narraway BC" ~ "Bearhole Redwillow",
    TRUE ~ herd
  )
)

hd <- hd %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
counts <- counts %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
afr <- afr %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
afs <- afs %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
ecotype <- ecotype %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))

# Pull out some values, summarize to demog data frame ---------------------

#### TREATMENT COMBOS
# Refactor to use summarize for clarity
treatment.combos <- trt |>
  filter(applied == 1) |>
  group_by(herd, year) |>
  summarize(
    trt = paste(
      gsub(" ", "", treatment), # Remove spaces
      collapse = "-" # Add - between treatments
    ),
    .groups = "drop"
  ) %>%
  select(herd, yrs = year, trt)

# treatment.combos <- trt %>%
#   filter(applied %in% 1) %>%
#   select(-intensity) %>%
#   pivot_wider(names_from = treatment, values_from = treatment) %>%
#   mutate(
#     trt = paste(
#       `reduce wolves`, `sterilize wolves`, `reduce moose`, pen, feed,
#       transplant,
#       sep = "-"
#     ) %>%
#     str_replace_all(" ", "") %>%
#     str_replace_all("-NA", "") %>%
#     str_replace_all("NA-", "")
#   ) %>%
#   select(herd, yrs = year, trt)

## pull posterior draws, add in herd, year, and treatment to each herd-year

# I did not find the parameter SR in the model file CaribouIPM_BCAB_jags_offset.txt .

# demog is grouped by i, just a note

demog <- out %>%
  spread_draws(
    c(totNMF, totN, totAdults, S, R_adj, lambda, R.ad, SR)[i, j],
    ndraws = 1000
  ) %>%
  median_qi(.width = 0.9) %>% # 0.9 is somewhat unconventional, just be explicit in writing, reminder that we can't do math on the median
  left_join(yr_df %>% rename(j = yr_idx), by = "j") %>%
  left_join(hd %>% select(herd, i = herd_num), by = "i") %>%
  left_join(treatment.combos, by = c("herd", "yrs")) %>%
  mutate(trt = replace_na(trt, "Reference")) %>%
  rename(R = R_adj, R.lower = R_adj.lower, R.upper = R_adj.upper)


## filter herds to a period where not fully extirpated

# No longer need this if using refactor below???
extirpated.yr <- demog %>%
  filter(round(totN, 0) <= 1) %>%
  group_by(herd) %>%
  filter(yrs == min(yrs)) %>%
  ungroup()


herds <- unique(demog$herd)

demog.trim <- tibble()
# A nice alternative, not sure if this is what you want, but it sets up the
# columns in the correct type
# demog.trim <- demog |> slice(0)

# I think all of this is because herds come back into exsitence or something?
# An alternative to below that sticks with dplyr syntax
# d_trim <- demog %>%
#   group_by(herd) |>
#   arrange(herd, yrs) |>
#   dplyr::mutate(
#     extrpt = cumsum(round(totN) < 2)
#   ) |>
#   dplyr::filter(extrpt == 0)

for (i in 1:length(herds)) {
  if (herds[i] %in% extirpated.yr$herd) {
    yr <- extirpated.yr %>%
      dplyr::filter(herd == !!herds[i])
    a <- demog %>%
      dplyr::filter(herd == !!herds[i] & yrs < yr$yrs)
  }

  if (!herds[i] %in% extirpated.yr$herd) {
    a <- demog %>%
      dplyr::filter(herd == !!herds[i])
  }

  demog.trim <- bind_rows(a, demog.trim)
}

demog <- demog.trim # Try not to reuse object names
rm(demog.trim)
rm(a)
rm(yr)

# write_csv(demog,"tables/demog.csv")

# Guess it doesn't matter, but the & keeps catching my eye. Using a comma
#  accomplishes the same thing. Maybe reserve the & for cases when () are
# required? This is very personal
demog.mod <- demog %>%
  filter(totAdults > 10 & totNMF > 20) ## modelling data that doesn't include functionally extirpated herds, demography gets unstable


## ----Plot herd abundance, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------

## Prep data and layout for plot

# Given that the names of the treatments were changed in treatment combos maybe
# it would be worthwhile to go back to trt and change the names there. Also, the
# multiple treatments per year combinations are lost in this trt.plot df. I
# think all of this is to help set the y-axis.
# e.g. trt |> filter(herd == "A La Peche", year == 1973)

trt.plot <- trt %>%
  filter(applied == 1) %>%
  left_join(
    demog %>%
      group_by(herd) %>%
      summarize(
        max = max(totNMF.upper)
      )
  ) %>%
  mutate(
    y = case_when(
      treatment %in% "reduce wolves" ~ max - (max * 0.1),
      treatment %in% "transplant" ~ max - (max * 0.15),
      treatment %in% "reduce moose" ~ max - (max * 0.2),
      treatment %in% "sterilize wolves" ~ max - (max * 0.25),
      treatment %in% "pen" ~ max - (max * 0.3),
      treatment %in% "feed" ~ max - (max * 0.35)
    )
  )

# Not sure why this is extracted from the package???
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

# I am running out of memory at this point, setting the ndraws to 1k
# Not sure if it is just me, but sims has a yrs and a year column and they do
# not align
sims_tmp <- out %>%
  gather_draws(c(pred_totNMF, totNMF)[i, j], ndraws = 1000)

# sims_tmp is grouped by i, j, variable

################################################################################
# Logic of approach below explained in auxillary script
#  /review/deriving_parameters.R

################################################################################
# Commented out old code that caused issues
# sim <- sims_tmp %>% # Changed in an attemp to save memory
#   mean_qi(.width = 0.9) %>% # Changed to mean, see comment on line 291
#   left_join(yr_df %>% rename(j = yr_idx)) %>%
#   left_join(hd %>% select(herd, i = herd_num)) %>%
#   left_join(treatment.combos) %>%
#   mutate(trt = case_when(is.na(trt) ~ "Reference", TRUE ~ trt))
#
# # Compare to previous work where we found the mean for i, j, variable
# sim |> filter(i == 1, j == 1, .variable == "totNMF")

# Something seems wrong here. I would expect 1 value like above, but we got 15.
# Also, I know I wrote it elsewhere, but the yrs and year column don't match and
# seems a likely reason for the error.

# Refactor of 278:283
sim <- sims_tmp %>%
  dplyr::mutate(
    herd = hd$herd[i],
    year = yr_df$yrs[j]
  ) |>
  dplyr::left_join(treatment.combos %>% rename(year = yrs), by = c("herd", "year")) |>
  dplyr::mutate(
    trt = tidyr::replace_na(trt, "Reference")
  )

# # Total caribou in the study???
# sims.plot <- sims %>%
#   # Seems easier to write this in the positive as opposed to the negative also
#   # probably makes the code more fragile this way e.g. add new column???
#   filter(!.variable %in% c("totAdults", "totCalvesMF", "pred_totCalvesMF")) %>%
#   group_by(yrs, .variable) %>%
#   summarise(across(.value:.upper, ~ sum(.x))) %>% ######### The .point column has the value median, but the median cannot be summarized like the mean, so this could be wrong !!!!!!!!!!!!! If it were the mean this would all be fine, but then why do 264:270? The summary object in out could be used to do the same thing, no?????
#   mutate(.variable = case_when(
#     .variable %in% "pred_totNMF" ~ "Status quo",
#     TRUE ~ "Observed"
#   )) %>%
#   mutate(.variable = fct_relevel(.variable, "Status quo", "Observed"))

# Refactor - get total population size
# This is an attempt to follow the logic in the deriving_parameters.R script
sims.plot <- sim |>
  dplyr::group_by(year, .variable, .draw) |>
  dplyr::summarize(
    mean = sum(.value),
    .groups = "drop"
  ) |>
  dplyr::group_by(year, .variable) |>
  dplyr::summarize(
    LCL = quantile(mean, 0.025),
    UCL = quantile(mean, 0.975),
    mean = mean(mean),
    .groups = "drop"
  )

### An alternative using the summary object
tst <- out$summary |>
  tibble::as_tibble(
    rownames = "param"
  ) |>
  dplyr::filter(
    grepl("totNMF", param)
  ) |>
  dplyr::mutate(
    param = gsub("\\[", ",", param),
    param = gsub("\\]", "", param)
  ) |>
  tidyr::separate(param, into = c("param", "herd", "year"), sep = ",") |>
  dplyr::mutate(
    herd_num = as.integer(herd), # for sorting
    herd = hd$herd[herd_num],
    year = yrs[as.integer(year)], # yrs to match sims.plot
    var = sd^2 # Variance can be summed
  )

# Population sum over all populations
tot_bou <- tst |>
  dplyr::group_by(param, year) |>
  summarize(
    dplyr::across(c(mean, var), sum),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    .variable = c("Observed", "Status quo")[(param == "pred_totNMF") + 1],
    LCL = mean - 1.96 * sqrt(var),
    UCL = mean + 1.96 * sqrt(var)
  )

# Compare - Works pretty well given short posteriors from subsetting gather
tot_bou |>
  filter(param == "pred_totNMF", year == 1973) |>
  select(param, year, LCL, mean, UCL)

sims.plot |>
  filter(.variable == "pred_totNMF", year == 1973) |>
  select(.variable, year, LCL, mean, UCL)

tot_bou |>
  filter(param == "pred_totNMF", year == 2021) |>
  select(param, year, LCL, mean, UCL)

sims.plot |>
  filter(.variable == "pred_totNMF", year == 2021) |>
  select(.variable, year, LCL, mean, UCL)


################################################################################
# The fast way to summarize would be to use the summary object.
# If for some reason we have to use the posteriors directly then we have to be
# very mindful of the dimensions and the functions being applied. The mean is
# the summary function to use for most every case, but the median can be the
# final summary if no other math needs to be done. We can't do math on the
# median.

################################################################################
ext.yr <- sim %>%
  dplyr::ungroup() |>
  dplyr::filter(.variable %in% c("totAdults", "totNMF")) %>%
  dplyr::select(.variable, .value, year, herd) %>%
  # pivot_wider(names_from=".variable", values_from=".value")%>%
  dplyr::group_by(herd, year, .variable) %>%
  dplyr::summarise(
    median = round(median(.value)),
    .groups = "drop"
  ) %>%
  dplyr::group_by(herd) %>%
  dplyr::mutate(
    # viable = case_when(
    #   yrs == 2021 & round(totAdults, 0) > 8 & round(totNMF,0) > 20 ~ 1,
    #   TRUE~0
    # ),
    nmf_viable = year == 2021 & (median > 8 & .variable == "totNMF"),
    adult_viable = year == 2021 & (median > 20 & .variable == "totAdults"),
    viable = max(nmf_viable, adult_viable)
  ) %>%
  # Doesn't viable == 0 cover the other two cases as well? Not sure I am
  # understanding this line
  dplyr::filter(viable == 0) %>% ## pull if 8 or fewer females present, or 20 or fewer total
  dplyr::filter(year == (min(year))) %>%
  dplyr::left_join(
    sims.plot %>%
      dplyr::filter(.variable %in% "Observed") %>%
      dplyr::select(year, mean),
    by = "year"
  )

## Herd Abundance
# ggplot() +
#   geom_ribbon(data = demog, aes(x = yrs, y = totNMF, ymin = totNMF.lower, ymax = totNMF.upper), alpha = 0.4, color = NA, fill = cols[3]) +
#   geom_line(data = demog, aes(x = yrs, y = totNMF, ymin = totNMF.lower, ymax = totNMF.upper), size = 1, color = cols[3]) +
#   geom_rug(
#     data = rbind(afr %>% select(herd, year, est) %>% mutate(type = "Recruit"), afs %>% select(herd, year, est) %>% mutate(type = "Surv"), counts %>% select(herd, year, est = Est_CL) %>% mutate(type = "Count")),
#     aes(x = year), sides = "t", length = unit(0.05, "npc"), alpha = 0.5
#   ) +
#   theme_ipsum() +
#   theme(legend.position = "none") +
#   ylab("") +
#   xlab("Year") +
#   facet_wrap(vars(herd), scales = "free_y", ncol = 8) +
#   labs(x = "Year", title = "Abundance") +
#   expand_limits(y = 0) +
#   scale_x_continuous(
#     limits = c(1974, 2026),
#     breaks = seq(1980, 2020, by = 20)
#   ) +
#   theme(
#     axis.title.x = element_text(size = 15),
#     axis.title.y = element_text(size = 15),
#     strip.text = element_textbox_highlight(
#       size = 14,
#       fill = "white", box.color = "white",
#       halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
#       padding = margin(4, 0, 4, 0), margin = margin(4, 1, 4, 1),
#       hi.labels = ext.yr$herd,
#       hi.fill = "firebrick", hi.box.col = "firebrick", hi.col = "white"
#     ),
#     legend.text = element_text(size = 13),
#     legend.title = element_text(size = 15)
#   ) +
#   geom_point(data = counts, aes(x = year, y = Est_CL), size = 0.5, alpha = 0.7) +
#   geom_linerange(data = counts %>% mutate(Est_CL.max = case_when(Est_CL.max > 5000 ~ 5000, TRUE ~ Est_CL.max)), aes(x = year, ymin = Est_CL.min, ymax = Est_CL.max), alpha = 0.5) +
#   geom_point(data = trt.plot, aes(x = year, y = y, group = treatment, color = treatment), size = 0.5) +
#   scale_color_manual(values = cols[-4]) +
#   geom_text(
#     data = trt.plot %>% distinct(herd, treatment, y) %>% mutate(t = str_remove(treatment, "reduce ") %>% str_sub(1, 1)), aes(label = t, colour = treatment, x = 2025, y = y),
#     direction = "y"
#   ) +
#   coord_cartesian(clip = "off")


# ggsave(here::here("plots", "abundance.png"), width = 15, height = 11, bg = "white")


## ----Plot total abundance, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------

#### Summarize bou pop in '73 vs 2021####
# sims.summary <- sim %>%
#   dplyr::group_by(herd) %>%
#   dplyr::filter(year == min(year) | year == max(year)) %>%
#   arrange(yrs)%>%
#   mutate(time=case_when(yrs==min(yrs)~"1973",
#                         TRUE~"2021"))%>%
#   group_by(time,.variable)%>%
#   summarise(abund=sum(.value),
#             lower=sum(.lower),
#             upper=sum(.upper))%>%
#   ungroup

# I think this is close to what was desired
sims.summary <- tst |>
  dplyr::group_by(herd) |>
  dplyr::filter(year %in% c(min(year), max(year))) |>
  dplyr::group_by(year, param) |>
  dplyr::summarise(
    abund = sum(mean),
    var = sum(var),
    lower = abund - 1.96 * sqrt(var),
    upper = abund + 1.96 * sqrt(var)
  ) |>
  dplyr::relocate(lower, .before = abund) |>
  dplyr::relocate(upper, .after = abund)


## what was the decline?
p.decline <- sims.summary %>%
  dplyr::filter(param == "totNMF") %>%
  dplyr::select(year, abund) %>%
  tidyr::pivot_wider(names_from = year, values_from = abund) %>%
  dplyr::summarise(
    dif = round(((`1973` - `2021`) / `1973`) * 100)
  )

## how many more caribou now than status quo w/ error
sims.draws <- out %>%
  gather_draws(
    pred_totNMF[i, j], totNMF[i, j], totCalvesMF[i, j], pred_totCalvesMF[i, j],
    ndraws = 1000
  )

n.recovery.all <- sims.draws %>%
  filter(j == 49, .variable %in% c("totNMF", "pred_totNMF")) %>%
  group_by(.draw, .variable) %>%
  summarise(across(.value, ~ sum(.x))) %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  mutate(dif = totNMF - pred_totNMF) %>%
  pull(dif)

quantile(n.recovery.all, c(0.05, 0.5, 0.95)) %>% round(0)

n.recovery <- median(n.recovery.all) %>% round(0)

## do again but just for calves
calves.recovered.all <- sims.draws %>%
  filter(.variable %in% c("totCalvesMF", "pred_totCalvesMF")) %>%
  group_by(.draw, .variable, j) %>%
  summarise(across(.value, ~ sum(.x))) %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  mutate(dif = totCalvesMF - pred_totCalvesMF) %>%
  group_by(.draw) %>% # I don't think I can wrap my head around this grouping. I really expected to take the mean here or similar, like with the n.recovery.all above
  summarise(dif = sum(dif)) %>%
  pull(dif)

quantile(calves.recovered.all, c(0.05, 0.5, 0.95)) %>% round(0)

calves.recovered <- median(calves.recovered.all) %>% round(0)

sum(out$mean$totCalvesMF - out$mean$pred_totCalvesMF)


## with error


median(n.recovery.all)
quantile(n.recovery.all, c(0.05, 0.5, 0.95))


#### Total Abundance####
# abundance.all.plot <- ggplot(data = sims.plot) +
#   # ggplot(aes(x=yrs,  y=.value, color=.variable,fill=.variable)) +
#   geom_ribbon(alpha = 0.3, aes(x = yrs, y = .value, ymin = .lower, ymax = .upper, color = .variable, fill = .variable)) +
#   geom_line(size = 1, aes(x = yrs, y = .value, ymin = .lower, ymax = .upper, color = .variable, fill = .variable)) +
#   geom_text(data = sims.plot %>% filter(yrs == 2021), aes(label = .variable, colour = .variable, x = Inf, y = .value), hjust = 0) +
#   geom_jitter(data = ext.yr, size = 1, aes(x = yrs, y = .value), alpha = 0.5) +
#   theme_ipsum() +
#   theme(legend.position = "none") +
#   ylab("") +
#   xlab("Year") +
#   labs(
#     x = "Year", y = "Abundance", title = "b) Population Trend",
#     subtitle = paste0(p.decline, "% decline since 1973, +", n.recovery, " caribou from recovery measures")
#   ) +
#   expand_limits(y = 0) +
#   scale_y_continuous(expand = c(0, 100)) +
#   scale_x_continuous(breaks = seq(1980, 2020, by = 20)) +
#   theme(
#     axis.title.x = element_text(size = 15),
#     axis.title.y = element_text(size = 15),
#     strip.text.x = element_text(size = 15),
#     strip.text.y = element_text(size = 15),
#     legend.text = element_text(size = 13),
#     legend.title = element_text(size = 15),
#     plot.margin = unit(c(1, 5, 1, 1), "lines")
#   ) +
#   geom_rug(data = trt %>% filter(applied %in% 1), aes(x = year), sides = "b", length = unit(0.05, "npc"), alpha = 0.05) +
#   scale_color_manual(values = cols[c(3, 1)]) +
#   annotate(geom = "text", x = 1974, y = 2200, label = "Recovery\nmeasures", hjust = "left") +
#   annotate(
#     geom = "curve", x = 1976, y = 1200, xend = 1985, yend = 300,
#     curvature = .3, arrow = arrow(length = unit(2, "mm"))
#   ) +
#   annotate(geom = "text", x = 1985, y = 8500, label = "Extirpation\nevent", hjust = "left") +
#   annotate(
#     geom = "curve", x = 1985, y = 8500, xend = (ext.yr %>% ungroup() %>% filter(yrs == min(yrs)) %>% pull(yrs)) + 1, yend = (ext.yr %>% ungroup() %>% filter(yrs == min(yrs)) %>% pull(.value)) + 100,
#     curvature = .3, arrow = arrow(length = unit(2, "mm"))
#   ) +
#   coord_cartesian(
#     clip = "off"
#   )
# abundance.all.plot
#
# #ggsave(plot = abundance.all.plot, here::here("plots", "abundance_all.png"), width = 6, height = 6, bg = "white")


## ----Plot total abundance ecotype, message=FALSE, warning=FALSE,eval=FALSE,include=FALSE-------------------------------------------------------------------
## ##for AB only:filter(!herd%in%c("A La Peche","Banff","Brazeau","Maligne","Narraway AB","Tonquin","Redrock/Prairie Creek"))%>%
##
## sims.plot.grp <- sims%>%
##   filter(.variable%in%c("pred_totNMF","totNMF"))%>%
##   left_join(ecotype)%>%
##   group_by(yrs,.variable,ECCC)%>%
##   summarise(across(.value:.upper,~sum(.x)))%>%
##   mutate(.variable=case_when(.variable%in%"pred_totNMF"~"Status quo*",
##                              TRUE~"Observed"))%>%
##   ungroup%>%
##   mutate(.variable=fct_relevel(.variable,"Status quo*","Observed"),
##          ECCC=fct_relevel(ECCC,"Southern Group","Central Group", "Northern Group"))
##
## ext.yr.grp <- sims%>%
##   filter(.variable=="totNMF")%>%
##   group_by(herd,yrs)%>%
##   summarise(abund=median(.value))%>%
##   filter(abund<20)%>%
##   filter(yrs==(min(yrs)))%>%
##   left_join(ecotype)%>%
##   left_join(sims.plot.grp%>%filter(.variable%in%"Observed")%>%select(yrs,.value,ECCC))%>%
##   select(herd,yrs,abund,ECCC,.value)%>%
##   ungroup%>%
##   mutate(ECCC=fct_relevel(ECCC,"Southern Group","Central Group", "Northern Group"))
##
##
## ggplot(data=sims.plot.grp)+
##   geom_ribbon(alpha=0.3, aes(x=yrs,  y=.value, ymin=.lower, ymax=.upper, color=.variable,fill=.variable))+
##   geom_path(size=1, aes(x=yrs,  y=.value, ymin=.lower, ymax=.upper, color=.variable,fill=.variable)) +
##   geom_text(data = sims.plot.grp%>%filter(yrs==2021), aes(label = .variable, colour = .variable, x = Inf, y = .value), hjust = 0) +  theme_ipsum()+
##   theme(legend.position = "none")+
##   ylab("")+
##   xlab("Year")+
##   labs(x="Year", y="Abundance",title="SMC Population Trend",
##        subtitle=paste0("",p.decline,"% decline since 1974, +", n.recovery," caribou from recovery measures"),
##        caption ="*Predicted population trend under a status quo scenario where no recovery measures were implemented")+
##   expand_limits(y=0)+
##   scale_y_continuous(expand = c(0, 100))+
##   scale_x_continuous(breaks = seq(1980, 2020, by = 20))+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=12),
##         strip.text.y = element_text(size=15),
##         plot.caption = element_text(color=cols[c(3)],size=12),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15),
##         plot.margin = unit(c(1,5,1,1), "lines"))+
##   geom_rug(data=trt%>%left_join(ecotype)%>%filter(applied%in%1), aes(x=year), sides = "b",length = unit(0.05, "npc"),alpha=0.05)+
##   scale_color_manual(values=cols[c(3,1)])+
##   geom_text(data=tibble(lab="Recovery\nmeasures", yrs = 1990, .value = 1300,ECCC="Northern Group"), aes(x=yrs,  y=.value, label="Recovery\nmeasures"),hjust = "left")+
##   geom_curve(data=tibble(x = 2002, y = 1000, xend = 2010, yend = 500,ECCC="Northern Group"),
##              aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE,curvature = -.3, arrow = arrow(length = unit(2, "mm")))+
##   coord_cartesian(
##     clip = "off"
##   )+
##   facet_wrap(vars(ECCC),scales="free_y")
##
##
## #ggsave(here::here("plots","abundance_bygroup.png"), width=13, height=6, bg="white")


## ----Plot vital rates, message=FALSE, warning=FALSE,eval=FALSE,include=FALSE-------------------------------------------------------------------------------
## ####Survival####
## ggplot() +
##   geom_ribbon(data=demog,aes(x=yrs, y=S, ymin=S.lower, ymax=S.upper, fill=herd),alpha=0.4)+
##   geom_point(data=demog,aes(x=yrs, y=S, ymin=S.lower, ymax=S.upper, fill=herd),alpha=0.7) +
##   theme_ipsum()+
##   theme(legend.position = "none")+
##   ylab("")+
##   xlab("Year")+
##   facet_wrap(vars(herd), scales="free_y")+
##   labs(x="Year",title="Survival")+
##   #expand_limits(y=0)+
##   scale_x_continuous(breaks = seq(1980, 2020, by = 20))+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15))+
##   geom_point(data=afs,aes(x=year, y=est, color=herd),alpha=0.7)+
##   geom_linerange(data=afs%>%mutate(ucl=case_when((est+sd)>1~1,
##                                                  TRUE~est+sd),
##                                    lcl=case_when((est-sd)<0~0,
##                                                  TRUE~est-sd)),aes(x=year, ymin=lcl,  ymax=ucl, color=herd),alpha=0.7)
## #ggsave(here::here("plots","survival.png"), width=14, height=11, bg="white")
##
## ####Recruitment####
## ggplot() +
##   geom_ribbon(data=demog,aes(x=yrs, y=R, ymin=R.lower, ymax=R.upper, fill=herd),alpha=0.4)+
##   geom_line(data=demog,aes(x=yrs, y=R, ymin=R.lower, ymax=R.upper, fill=herd)) +
##   geom_point(data=demog,aes(x=yrs, y=R, ymin=R.lower, ymax=R.upper, fill=herd)) +
##   theme_ipsum()+
##   theme(legend.position = "none")+
##   ylab("")+
##   xlab("Year")+
##   #expand_limits(y=0)+
##   scale_x_continuous(breaks = seq(1980, 2020, by = 20))+
##   facet_wrap(vars(herd), scales="free_y")+
##   labs(x="Year",title="Recruitment")+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15))+
##   geom_point(data=afr,aes(x=year, y=est, color=herd))+
##   geom_linerange(data=afr%>%mutate(ucl=case_when((est+sd)>1~1,
##                                                  TRUE~est+sd),
##                                    lcl=case_when((est-sd)<0~0,
##                                                  TRUE~est-sd)),aes(x=year, ymin=lcl,  ymax=ucl, color=herd),alpha=0.7)
## #ggsave(here::here("plots","recruitment.png"), width=14, height=11, bg="white")
##
##
## #write_csv(demog.trim, here::here("tables","demog.trim.csv"))
##
##
## ####Adult Recruitment####
## ggplot() +
## geom_ribbon(data=demog,aes(x=yrs, y=R.ad, ymin=R.ad.lower, ymax=R.ad.upper, fill=herd),alpha=0.4)+
##   geom_line(data=demog,aes(x=yrs, y=R.ad, ymin=R.ad.lower, ymax=R.ad.upper, fill=herd)) +
##   geom_point(data=demog,aes(x=yrs, y=R.ad, ymin=R.ad.lower, ymax=R.ad.upper, fill=herd)) +
##   theme_ipsum()+
##   theme(legend.position = "none")+
##   ylab("")+
##   xlab("Year")+
##   facet_wrap(vars(herd), scales="free_y")+
##   labs(x="Year",title="Recruitment (calves/adult)")+
##   expand_limits(y=0)+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15))
##
## ####Lambda####
## ggplot() +
##   geom_ribbon(data=demog,aes(x=yrs, y=lambda, ymin=lambda.lower, ymax=lambda.upper, fill=herd),alpha=0.4)+
##   geom_line(data=demog,aes(x=yrs, y=lambda, ymin=lambda.lower, ymax=lambda.upper, fill=herd)) +
##   geom_point(data=demog,aes(x=yrs, y=lambda, ymin=lambda.lower, ymax=lambda.upper, fill=herd)) +
##   theme_ipsum()+
##   theme(legend.position = "none")+
##   ylab("")+
##   xlab("Year")+
##   facet_wrap(vars(herd), scales="free_y")+
##   labs(x="Year",title="Population Growth")+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15))+
##   geom_hline(data=data.frame(herd=unique(demog$herd)),
##              aes(yintercept =1),linetype="dashed")
##


## ----trt eff- r, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------

## Gather draws
demog.draws <- out %>%
  gather_draws(lambda[i, j], S[i, j], R_adj[i, j], R.ad[i, j], totNMF[i, j], ndraws = 1000) %>%
  mutate(
    # .variable = case_when(
    #   .variable %in% "R_adj" ~ "R",
    #   TRUE ~ .variable
    # )
    .variable = replace(.variable, .variable == "R_adj", "R") # An alternative
  )

### pivot wider
demog.draws <- demog.draws %>%
  ungroup() %>%
  pivot_wider(names_from = .variable, values_from = .value)


## add in years and herd names
herd_df <- trt %>%
  distinct(herd_num, herd) %>%
  arrange(herd_num)

# AHHH, trying to troubleshoot but demog.draws relies on demog.draws...
demog.draws <- demog.draws %>%
  mutate(
    year = yr_df$yrs[match(j, yr_df$yr_idx)],
    herd = herd_df$herd[match(i, herd_df$herd_num)]
  ) %>%
  # left_join(
  #   trt %>%
  #     distinct(herd_num, herd) %>%
  #     select(i = herd_num, herd),
  #   #by = c("i") # This is wrong and is responsible for the yrs and year comments below
  #   by = "year"
  # ) %>%
  left_join(
    treatment.combos %>% rename(year = yrs),
    by = c("herd", "year")
  ) %>%
  left_join(
    trt %>%
      distinct(herd, year, intensity) %>%
      select(herd, yrs = year, intensity) %>%
      filter(intensity %in% "low"),
    by = c("herd", c("year" = "yrs")) # Year column name changed again
  ) %>%
  replace_na(
    list(
      trt = "Reference"
    )
  ) |>
  mutate(
    trt_active = trt != "Reference"
  )
# mutate(trt = case_when(is.na(trt) ~ "Reference", TRUE ~ trt))


# I don't get the next chunk (797-801), it may work, but it doesn't appear to do
# anything at first blush. If the desire is to remove lambda it would be really
# easy without the pivot_wider, but when wide you have to deal with other values
# too. Maybe set it to NA or something so the plotting code omits it?
# IF the first year is 1973 for all groups
demog.draws |>
  dplyr::mutate(
    lambda = replace(lambda, yrs == 1973, NA_real_)
  )

# Otherwise
demog.draws |>
  dplyr::group_by(herd) |>
  dplyr::mutate(
    lambda = replace(lambda, yrs == min(yrs), NA_real_)
  )

## remove first year lambda for each herd, as lambda==1
demog.draws <- demog.draws %>%
  group_by(herd, .draw) %>%
  arrange(year) %>% # Note that we have the yrs v year problem here as well
  slice(-1) %>%
  ungroup()
#
# # I offered an alternative above, I can't always get the code below to run
# ## filter herds to a period where not extirpated
# demog.draws.trim <- tibble()
# for (i in 1:length(herds)) {
#   if (herds[i] %in% ext.yr$herd) {
#     yr <- ext.yr %>%
#       dplyr::filter(herd == !!herds[i])
#     a <- demog.draws %>%
#       dplyr::filter(herd == !!herds[i] & yrs < yr$yrs)
#   }
#
#   if (!herds[i] %in% ext.yr$herd) {
#     a <- demog.draws %>%
#       dplyr::filter(herd == !!herds[i])
#   }
#
#   demog.draws.trim <- bind_rows(a, demog.draws.trim)
# }
#
# demog.draws <- demog.draws.trim
# rm(demog.draws.trim)

# Only two treatment types have transplant in the name, it would be more general
# to filter on trt, no? I can't make it work with yrs though, every year is
# present

### Remove years where transplant was being done (impossibly high lambda due to adding individuals, not population response)
# At this point demog.draws has been defined ~8 times...some of which could very
# well be my fault, consider consdensing after reading

demog.draws <- demog.draws %>%
  filter(
    !(herd %in% "Charlotte Alplands" & yrs %in% c(1984:1991)),
    !(herd %in% "Telkwa" & yrs %in% c(1997:1999)),
    !(herd %in% "South Selkirks" & yrs %in% c(1987:1989)),
    !(herd %in% "Purcell South" & yrs %in% c(2012:2013))
  )

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

# This retains each MCMC draw throughout, so not sure what is happening, refactor offered below might be closer to what is desired
# In my mind we want to transform or derive the new quantity at each draw and then summarize the resulting distribution over all the draws

# demog.draws.combotreat <- demog.draws %>%
#   # filter(!(intensity%in%"low"|totNMF<20))%>%
#   group_by(.draw, trt, herd) %>%
#   summarise(
#     lambda = exp(mean(log(lambda)))
#   ) %>% ## geo mean per herd-treatment-draw
#   group_by(.draw, trt) %>% # Seems like you would drop .draw at this point and calculate the mean or median to summarize the point estimate????
#   summarise(
#     lambda = exp(mean(log(lambda))) # Doing this twice seems odd, not sure what is happening here or the intent
#   ) %>% ## geo mean lambda per treatment-draw
#   mutate(
#     trt = case_when(
#       is.na(trt) ~ "Reference",
#       TRUE ~ trt
#     ),
#     group = case_when(
#       trt %in% "Reference" ~ "Reference",
#       TRUE ~ "Treatment"
#     )
#   )

geo_mean <- function(x) {
  prod(x)^(1 / length(x))
}

geo_alt <- demog.draws |>
  group_by(.draw, trt, herd) |> # Compute transformation at each draw
  summarise(
    lambda = geo_mean(lambda),
    .groups = "drop"
  ) |>
  group_by(trt, herd) |> # Summarise across draws
  summarise(
    LCL = quantile(lambda, 0.025),
    UCL = quantile(lambda, 0.975),
    lambda = mean(lambda),
    .groups = "drop"
  ) |>
  dplyr::relocate(UCL, .after = lambda)


demog.draws.combotreat.rug <- demog.draws %>%
  # filter(!(intensity%in%"low"|totNMF<20))%>%
  group_by(trt, herd) %>%
  summarise(
    lambda = exp(mean(log(lambda))) # This gives equivalent to chunk above, cool
  ) %>% ## geo mean per herd-treatment-draw
  mutate(
    trt = case_when(
      is.na(trt) ~ "Reference",
      TRUE ~ trt
    ),
    group = case_when(
      trt %in% "Reference" ~ "Reference",
      TRUE ~ "Treatment"
    )
  )

# Not sure the intent here, this seems like something that needed to be done
# because all the draws were still present in the original code???
# order <- demog.draws.combotreat%>%
#   group_by(trt)%>%
#   summarise(
#     med = median(lambda)
#   )


# demog.draws.combotreat %>%
#   left_join(order) %>%
# Here I show the failure using "fully" summarized data, but I think we want to
# retain the posteriors for the geom_density ridges, if we pass a single value
# to that function there is no density???
geo_alt |>
  filter(trt != "transplant") %>%
  ggplot(aes(x = log(lambda), y = fct_reorder(trt, lambda), fill = trt)) +
  geom_density_ridges( # scale = 1.5,
    # scale = 1.3,
    rel_min_height = .01,
    size = 0.25,
    alpha = 0.9
  ) +
  theme_ipsum()

# No need to create a new object, just start with demog.draws and move forward
demog.draws |>
  group_by(.draw, trt, herd) |> # Compute transformation at each draw
  summarise(
    lambda = geo_mean(lambda),
    .groups = "drop"
  ) |>
  mutate(
    group = "Treatment",
    group = replace(group, trt == "Reference", "Reference")
  ) |>
  filter(trt != "transplant") %>%
  ggplot(aes(x = log(lambda), y = fct_reorder(trt, lambda), fill = group)) +
  geom_density_ridges(
    rel_min_height = .01,
    size = 0.25,
    alpha = 0.9
  ) +
  theme_ipsum() +
  geom_point(
    data = demog.draws.combotreat.rug %>%
      # left_join(order) %>%
      filter(trt != "transplant"),
    aes(y = fct_reorder(trt, lambda), x = log(lambda)),
    shape = "|"
  ) +
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
  xlim(-0.25, 0.3) +
  labs(
    x = "Instantaneous rate of increase (r)",
    y = "Recovery measure(s)",
    fill = "",
    title = "Instantaneous Rate of Increase by Recovery Measure"
  ) +
  scale_fill_manual(values = c(cols[c(3, 1)]))


# ggsave(here::here("plots", "lambda_treatments.png"), width = 8, height = 7, bg = "white")


# Would be nice to show the years by herd that each prescription was applied
demog.draws |>
  group_by(.draw, trt, herd) |> # Compute transformation at each draw
  summarise(
    lambda = geo_mean(lambda),
    .groups = "drop"
  ) |>
  group_by(trt) |>
  summarise(
    LCL = round(quantile(lambda, 0.025), 2),
    UCL = round(quantile(lambda, 0.975), 2),
    lambda = round(median(lambda), 2)
  ) |>
  relocate(UCL, .after = lambda) |>
  arrange(desc(lambda)) |>
  mutate(
    lambda = paste0(lambda, " (", LCL, "-", UCL, ")")
  ) %>%
  select(-c(LCL, UCL)) |>
  gt()


## ----trt eff- BA, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------------
# Before-After ---------------------------------------------

# The following didn't work until I refactored some of the demog.draws code, not
# sure what you think about it

## pull draws and organize into treatment and untreated (reference) timeframes for each herd
eff.draws <- demog.draws %>%
  mutate(
    lambda = log(lambda)
  ) %>% ## little r
  ungroup() %>%
  filter(trt != "Reference") %>%
  group_by(.draw, herd, trt) %>%
  summarise(across(lambda:R, ~ mean(.x))) %>% # Don't see the exp portion of the geo mean ##little R
  pivot_longer(lambda:R) %>%
  ungroup() %>%
  select(.draw, trt, herd, name, eff = value) %>%
  left_join(
    demog.draws %>% # This is a big lift, seems like we should be able to leave Reference in the data and work with the data that way
      mutate(lambda = log(lambda)) %>% ## little r
      ungroup() %>%
      filter(trt %in% "Reference") %>%
      group_by(.draw, herd) %>%
      summarise(across(lambda:R, ~ mean(.x))) %>%
      pivot_longer(lambda:R) %>%
      select(.draw, herd, name, ref = value),
    by = c(".draw", "herd", "name")
  ) %>%
  mutate(
    delta.i = exp(eff - ref) - 1,
    delta.i.lr = eff - ref
  )

# This is at least the second object called order, not sure why this is needed
#  though
order <- eff.draws %>%
  filter(name == "lambda") %>%
  group_by(trt) %>%
  summarise(med = median(delta.i.lr)) %>%
  arrange(-med)

# Again, because of the reuse of the names I can't see the result and then
#  go back and explore the code. eff.draws depends on eff.draws, which was just
#  changed so now this code doesn't mean much and can't be run again. Consider
#  reducing the number of objects created and writing a few functions to reduce
#  the amount of code.
eff.draws <- eff.draws %>%
  left_join(order) %>%
  mutate(
    trt = fct_reorder(trt, med),
    name = case_when(
      name == "lambda" ~ "Rate of increase (r)",
      name == "R" ~ "Recruitment (R)",
      name == "S" ~ "Survival (S)",
      TRUE ~ name
    )
  )

# What is the scale of the x-axis here?
# I would think you could do this by making the units 1 SD, (x - mean)/sd kind
#  of approach, but if the units aren't the same it could be very misleading
ggplot() +
  geom_density_ridges(
    data = eff.draws %>%
      filter(trt != "transplant") %>%
      group_by(name, trt, .draw) %>%
      summarise(delta = median(delta.i.lr)),
    aes(x = delta, y = trt, fill = name),
    scale = .9,
    rel_min_height = .01,
    size = 0.25,
    alpha = 0.9
  ) +
  geom_point(
    data = eff.draws %>%
      filter(trt != "transplant") %>%
      group_by(herd, name, trt) %>%
      summarise(delta = median(delta.i.lr)),
    aes(y = trt, x = delta),
    shape = "|"
  ) +
  theme_ipsum() +
  facet_wrap(vars(name)) +
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
  labs(x = "Change in value", y = "Recovery measure(s)", title = "Before-After Assessment of Effectivness") +
  xlim(-0.2, 0.4) +
  scale_fill_manual(values = cols[c(1:3)])

# ggsave(here::here("plots", "ba_all.png"), width = 10, height = 7, bg = "white")


trt_eff_ba_table <- eff.draws %>%
  filter(name == "Rate of increase (r)") %>%
  group_by(trt) %>% # Don't need to group by draw, just summarise
  # summarise(delta = median(delta.i)) %>% # this is not needed and could be wrong if the grouping were not maintined by default the calculation in the next summarise would just be wrong
  summarise(
    lower = quantile(delta.i, 0.05, na.rm = TRUE) %>% round(2),
    delta.l = median(delta.i, na.rm = TRUE) %>% round(2),
    upper = quantile(delta.i, 0.95, na.rm = TRUE) %>% round(2)
  ) %>%
  arrange(-delta.l) %>%
  mutate(
    delta.lambda = paste0(delta.l, " (", lower, "-", upper, ")")
  ) %>%
  select(Treatment = trt, delta.lambda) %>%
  filter(Treatment != "transplant")

# Please separate code that saves from regular interactive evaluation, I don't
# want to have to hunt down too many files
trt_eff_ba_table %>%
  gt() # %>%
# gtsave(here::here("tables", "trt_eff_ba.rtf"))


eff.draws %>%
  filter(name == "Rate of increase (r)") %>%
  group_by(herd, name, trt) %>%
  summarise(delta = median(delta.i.lr)) %>%
  ungroup() %>%
  count(delta > 0)


## ----trt eff- BA w intensity, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------

# Would be nice if there were only one of these, not sure if it is possible

eff.draws.app <- demog.draws %>%
  group_by(herd, year, trt) %>%
  mutate( # This "should be" a summarize, not a mutate. We can't group by herd, year and treat, summarize and still have each draw available. Anyway, be sure not to use this quantity later except for the logical.
    totNMF.median = median(totNMF)
  ) %>% # get average pop size so popsize threshold doesnt split low/standard in some years due to draws being above/below threshold
  ungroup() %>%
  mutate(
    lambda = log(lambda), ## little r
    application = case_when(
      intensity == "low" | totNMF.median < 30 ~ "low",
      TRUE ~ "standard"
    )
  ) %>%
  ungroup() %>%
  # Since this is done so many times I added a column called trt_active, so you could do
  # filter(!trt_active) |>
  filter(!trt %in% "Reference") %>%
  group_by(.draw, herd, trt, application) %>%
  summarise(across(lambda:R, ~ mean(.x))) %>%
  pivot_longer(lambda:R) %>%
  ungroup() %>%
  select(.draw, application, trt, herd, name, eff = value) %>%
  left_join(
    demog.draws %>%
      mutate(lambda = log(lambda)) %>% ## little r
      ungroup() %>%
      filter(trt %in% "Reference") %>%
      group_by(.draw, herd) %>%
      summarise(across(lambda:R, ~ mean(.x))) %>%
      pivot_longer(lambda:R) %>%
      select(.draw, herd, name, ref = value),
    by = c(".draw", "herd", "name")
  ) %>%
  mutate(
    delta.i = exp(eff - ref) - 1,
    delta.i.lr = eff - ref
  )

# Given that the first one is not filtered, doesn't it include everything that
#  is in the second object? Seems like maybe the first should be
#  application != "standard"
eff.draws.app.bind <- eff.draws.app %>%
  mutate(group = "all") %>%
  dplyr::bind_rows( # The dplyr version is faster and will not coerce your data types
    eff.draws.app %>%
      filter(application == "standard") %>%
      mutate(group = "standard")
  )



order <- eff.draws.app.bind %>%
  filter(name == "lambda") %>%
  group_by(trt) %>%
  summarise(med = median(delta.i.lr)) %>%
  arrange(-med)

eff.draws.app.bind <- eff.draws.app.bind %>%
  left_join(order) %>% # I still don't get this approach. We are attaching fully summarized data to raw posteriors and just repeating it over and over?
  mutate(
    trt = fct_reorder(trt, med),
    name = case_when(
      name == "lambda" ~ "Rate of increase (r)",
      name == "R" ~ "Recruitment (R)",
      name == "S" ~ "Survival (S)",
      TRUE ~ name
    )
  )


ggplot() +
  geom_density_ridges(
    data = eff.draws.app.bind %>%
      filter(trt != "transplant") %>%
      group_by(name, trt, .draw, group) %>%
      summarise(delta = median(delta.i.lr)), # Huh? The median within a draw...I read this as you want to summarize over the herds or reduce some other dimension. I think I would use the mean here, not sure. Do you mean for the data to enter the function grouped?
    aes(x = delta, y = trt, fill = name, linetype = group),
    scale = .9,
    rel_min_height = .01,
    size = 0.25,
    alpha = 0.7
  ) +
  geom_point(
    data = eff.draws.app.bind %>%
      filter(trt != "transplant", group == "all", !(herd == "Frisby-Boulder" & application == "standard")) %>%
      group_by(name, herd, trt, application) %>% # No draw here, same data though, right?
      summarise(delta = median(delta.i.lr)), # Grouped data, fyi
    aes(y = trt, x = delta, color = application),
    shape = "|", size = 4
  ) +
  theme_ipsum() +
  facet_wrap(vars(name)) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Change in value", y = "Recovery measure(s)", title = "Before-After Assessment of Effectivness") +
  xlim(-0.2, 0.4) +
  scale_fill_manual(values = cols[c(1:3)]) +
  scale_color_manual(values = cols[c(8, 6)]) +
  guides(fill = "none")

# ggsave(here::here("plots", "ba_all_intensity.png"), width = 10, height = 7, bg = "white")


# ggplot()+
#   geom_density_ridges(data= eff.draws.app.bind%>%
#                         filter(trt=="reducemoose",name=="Rate of increase (r)", application=="low"),
#                       aes(x=delta.i.lr, y=trt,fill=herd,linetype=herd),
#                       scale = .9,
#                       rel_min_height = .01,
#                       size = 0.25,
#                       alpha=0.7)



eff.draws.app.bind %>%
  filter(trt != "transplant", name == "Rate of increase (r)") %>%
  # group_by(.draw, name, trt, group) %>%
  # summarise(delta = median(delta.i.lr, na.rm = TRUE)) %>%
  group_by(trt, group) %>%
  summarise(
    delta.l = median(delta.i.lr, na.rm = TRUE) %>% round(2),
    lower = quantile(delta.i.lr, 0.05, na.rm = TRUE) %>% round(2),
    upper = quantile(delta.i.lr, 0.95, na.rm = TRUE) %>% round(2)
  ) %>%
  arrange(-delta.l) %>%
  arrange(trt, group)


## ----trt eff- BACI w ecotype, message=FALSE, warning=FALSE,eval=FALSE,include=FALSE------------------------------------------------------------------------
## ####ECCC grouping
## demog.mod.baci <- demog.draws%>%
##   ungroup%>%
##   left_join(ecotype)%>%
##   select(.draw,herd,yrs,lambda,totNMF,trt,ECCC)
##
##
## trt.herds <- demog.mod.baci%>%
##   filter(!trt%in%c("Reference","transplant"),.draw==1)%>%
##   group_by(herd,trt,ECCC)%>%
##   summarise(trt.start=min(yrs),
##             trt.end=max(yrs),
##             dur=n(),
##             baci.start=trt.start-dur)%>%
##   ungroup
##
## baci <- tibble()
## for(i in 1:nrow(trt.herds)){
##
##   post.trt <- demog.mod.baci%>%
##     dplyr::filter(herd%in%!!trt.herds[i,"herd"][[1]],
##                   yrs%in%!!c(trt.herds[i,"trt.start"][[1]]:trt.herds[i,"trt.end"][[1]]))%>%
##     dplyr::mutate(CI=1)
##
##   pre.trt <- demog.mod.baci%>%
##     dplyr::filter(trt%in%c("Reference"),
##                   herd%in%!!trt.herds[i,"herd"][[1]],
##                   yrs%in%!!c((trt.herds[i,"trt.start"][[1]]-11):(trt.herds[i,"trt.start"][[1]]-1)))%>%
##     dplyr::mutate(CI=1)
##
##   control <-demog.mod.baci%>%
##     ungroup%>%
##     dplyr::filter(trt%in%c("Reference"),
##                   !herd%in%!!trt.herds[i,"herd"][[1]],
##                   ECCC%in%!!trt.herds[i,"ECCC"][[1]],
##                   yrs%in%!!c((trt.herds[i,"trt.start"][[1]]-11):trt.herds[i,"trt.end"][[1]]))%>%
##     dplyr::mutate(CI=0)
##
##   ##only keep herds with full data
## if(all(nrow(pre.trt)>0,nrow(post.trt)>0,nrow(control)>0)){
##   a <- rbind(pre.trt, post.trt,control)%>%
##     dplyr::mutate(herd.grp=trt.herds[i,"herd"][[1]],
##                   trt.grp=trt.herds[i,"trt"][[1]],
##                   BA=case_when(yrs<trt.herds[i,"trt.start"][[1]]~0,
##                                yrs>=trt.herds[i,"trt.start"][[1]]~1),
##                   BA.break=trt.herds[i,"trt.start"][[1]])
## }
##
##   baci <- rbind(baci,a)
## }
##
## baci <- baci%>%mutate(lambda=log(lambda))%>% ##little r
##   mutate(grp=paste(herd.grp,trt, sep="_"),
##          facet=paste(herd.grp,trt.grp, sep="_"))
##
## ggplot() +
##   geom_line(data=baci%>%
##               group_by(herd,yrs,CI,facet,trt.grp,herd.grp)%>%
##               summarise(lambda=median(lambda),
##                         totNMF=median(totNMF))%>%
##               filter(CI==1),aes(x=yrs, y=lambda,  group=herd), color="red",alpha=1, size=1) +
##   geom_line(data=baci%>%
##               group_by(herd,yrs,CI,facet,trt.grp,herd.grp)%>%
##               summarise(lambda=median(lambda),
##                         totNMF=median(totNMF)),aes(x=yrs, y=lambda, group=herd), alpha=0.3) +
##   theme_ipsum()+
##   ylab("")+
##   xlab("Year")+
##   labs(x="Year",y="Instantaneous rate of increase (r)", title="BACI set up", subtitle="control herds shown in grey")+
##   expand_limits(y=0)+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15))+
##   geom_vline(data=baci%>%
##                distinct(BA.break, facet,herd.grp,trt.grp),aes(xintercept = BA.break-0.5), linetype="dashed") +
##   facet_wrap(vars(trt.grp,herd.grp))
## #ggsave(here::here("plots", "appendix","baci_setup_ECCC.png"), width=13, height=9, bg="white")
##
## baci.draws <- baci%>%
##   group_by(.draw,trt.grp) %>%
##   do(tidy(lmer(lambda~BA+CI+BA*CI + (1|grp), data=.)))%>%
##   filter(term%in%"BA:CI")
##
## ggplot(data=baci.draws, aes(x=estimate, y=fct_reorder(trt.grp,estimate)))+
##   geom_density_ridges(scale = .9,
##                       rel_min_height = .01,
##                       size = 0.25,
##                       alpha=0.9,
##                       fill=cols[c(1)])+
##   theme_ipsum()+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15),
##         legend.position="none")+
##   geom_vline(xintercept = 0, linetype="dashed")+
##   labs(x="Change in rate of increase", y="Recovery measure(s)", title="Before-After-Control-Impact", subtitle="w/ spatio-temporal matching using ECCC ecotype")+
##   xlim(-0.3,0.6)
## #ggsave(here::here("plots", "appendix","baci_all_ECCC.png"), width=10, height=7, bg="white")
##
## baci.draws%>%
##   group_by(trt.grp)%>%
##   summarise(delta.l=median(estimate,na.rm=TRUE)%>%round(2),
##             lower=quantile(estimate,0.05,na.rm=TRUE)%>%round(2),
##             upper=quantile(estimate,0.95,na.rm=TRUE)%>%round(2))%>%
##   arrange(-delta.l)%>%
##   mutate(delta.lambda=paste0(delta.l," (",lower,"-",upper,")"))%>%
##   select(Treatment=trt.grp,`Delta r (BACI)`=delta.lambda)%>%
##   left_join(trt_eff_ba_table%>% select(Treatment,`Delta r (BA)`=delta.lambda))%>%
##   gt()%>%
##  gtsave(here::here("tables", "appendix","trt_eff_baci_ECCC.rtf"))
##
##
##
##
## ####Heard and Vagt grouping
## demog.mod.baci <- demog.draws%>%
##   ungroup%>%
##   left_join(ecotype)%>%
##   select(.draw,herd,yrs,lambda,totNMF,trt,Heard_Vagt1998)
##
##
## trt.herds <- demog.mod.baci%>%
##   filter(!trt%in%c("Reference","transplant"),.draw==1)%>%
##   group_by(herd,trt,Heard_Vagt1998)%>%
##   summarise(trt.start=min(yrs),
##             trt.end=max(yrs),
##             dur=n(),
##             baci.start=trt.start-dur)%>%
##   ungroup
##
## baci <- tibble()
## for(i in 1:nrow(trt.herds)){
##
##   post.trt <- demog.mod.baci%>%
##     dplyr::filter(herd%in%!!trt.herds[i,"herd"][[1]],
##                   yrs%in%!!c(trt.herds[i,"trt.start"][[1]]:trt.herds[i,"trt.end"][[1]]))%>%
##     dplyr::mutate(CI=1)
##
##   pre.trt <- demog.mod.baci%>%
##     dplyr::filter(trt%in%c("Reference"),
##                   herd%in%!!trt.herds[i,"herd"][[1]],
##                   yrs%in%!!c((trt.herds[i,"trt.start"][[1]]-11):(trt.herds[i,"trt.start"][[1]]-1)))%>%
##     dplyr::mutate(CI=1)
##
##   control <-demog.mod.baci%>%
##     ungroup%>%
##     dplyr::filter(trt%in%c("Reference"),
##                   !herd%in%!!trt.herds[i,"herd"][[1]],
##                   Heard_Vagt1998%in%!!trt.herds[i,"Heard_Vagt1998"][[1]],
##                   yrs%in%!!c((trt.herds[i,"trt.start"][[1]]-11):trt.herds[i,"trt.end"][[1]]))%>%
##     dplyr::mutate(CI=0)
##
##   ##only keep herds with full data
##   if(all(nrow(pre.trt)>0,nrow(post.trt)>0,nrow(control)>0)){
##     a <- rbind(pre.trt, post.trt,control)%>%
##       dplyr::mutate(herd.grp=trt.herds[i,"herd"][[1]],
##                     trt.grp=trt.herds[i,"trt"][[1]],
##                     BA=case_when(yrs<trt.herds[i,"trt.start"][[1]]~0,
##                                  yrs>=trt.herds[i,"trt.start"][[1]]~1),
##                     BA.break=trt.herds[i,"trt.start"][[1]])
##   }
##
##   baci <- rbind(baci,a)
## }
##
## baci <- baci%>%mutate(lambda=log(lambda))%>% ##little r
##   mutate(grp=paste(herd.grp,trt, sep="_"),
##          facet=paste(herd.grp,trt.grp, sep="_"))
##
## ggplot() +
##   geom_line(data=baci%>%
##               group_by(herd,yrs,CI,facet,trt.grp,herd.grp)%>%
##               summarise(lambda=median(lambda),
##                         totNMF=median(totNMF))%>%
##               filter(CI==1),aes(x=yrs, y=lambda,  group=herd), color="red",alpha=1, size=1) +
##   geom_line(data=baci%>%
##               group_by(herd,yrs,CI,facet,trt.grp,herd.grp)%>%
##               summarise(lambda=median(lambda),
##                         totNMF=median(totNMF)),aes(x=yrs, y=lambda, group=herd), alpha=0.3) +
##   theme_ipsum()+
##   ylab("")+
##   xlab("Year")+
##   labs(x="Year",y="Instantaneous rate of increase (r)", title="BACI set up", subtitle="control herds shown in grey")+
##   expand_limits(y=0)+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15))+
##   geom_vline(data=baci%>%
##                distinct(BA.break, facet,herd.grp,trt.grp),aes(xintercept = BA.break-0.5), linetype="dashed") +
##   facet_wrap(vars(trt.grp,herd.grp))
## #ggsave(here::here("plots", "appendix","baci_setup_Heard_Vagt1998.png"), width=13, height=9, bg="white")
##
## baci.draws <- baci%>%
##   group_by(.draw,trt.grp) %>%
##   do(tidy(lmer(lambda~BA+CI+BA*CI + (1|grp), data=.)))%>%
##   filter(term%in%"BA:CI")
##
## ggplot(data=baci.draws, aes(x=estimate, y=fct_reorder(trt.grp,estimate)))+
##   geom_density_ridges(scale = .9,
##                       rel_min_height = .01,
##                       size = 0.25,
##                       alpha=0.9,
##                       fill=cols[c(1)])+
##   theme_ipsum()+
##   theme(axis.title.x = element_text(size=15),
##         axis.title.y = element_text(size=15),
##         strip.text.x = element_text(size=15),
##         strip.text.y = element_text(size=15),
##         legend.text = element_text(size=13),
##         legend.title=element_text(size=15),
##         legend.position="none")+
##   geom_vline(xintercept = 0, linetype="dashed")+
##   labs(x="Change in rate of increase", y="Recovery measure(s)", title="Before-After-Control-Impact", subtitle="w/ spatio-temporal matching using Heard & Vagt 1998")+
##   xlim(-0.3,0.6)
## #ggsave(here::here("plots", "appendix","baci_all_Heard_Vagt1998.png"), width=10, height=7, bg="white")
##
## baci.draws%>%
##   group_by(trt.grp)%>%
##   summarise(delta.l=median(estimate,na.rm=TRUE)%>%round(2),
##             lower=quantile(estimate,0.05,na.rm=TRUE)%>%round(2),
##             upper=quantile(estimate,0.95,na.rm=TRUE)%>%round(2))%>%
##   arrange(-delta.l)%>%
##   mutate(delta.lambda=paste0(delta.l," (",lower,"-",upper,")"))%>%
##   select(Treatment=trt.grp,`Delta r (BACI)`=delta.lambda)%>%
##   left_join(trt_eff_ba_table%>% select(Treatment,`Delta r (BA)`=delta.lambda))%>%
##   gt()%>%
##   gtsave(here::here("tables","trt_eff_baci_Heard_Vagt1998.rtf"))


## ----trt eff- BACI w intensity, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------
## prep data with individual treatments 1/0
# This is starting to look like case_when is your favorite function=)
# Maybe something like:
# reducewolves = as.integer(grepl("reducewolves", trt))
eff.draws <- eff.draws %>%
  mutate(
    reducewolves = case_when(str_detect(trt, "reducewolves") ~ 1, TRUE ~ 0),
    sterilizewolves = case_when(str_detect(trt, "sterilizewolves") ~ 1, TRUE ~ 0),
    reducemoose = case_when(str_detect(trt, "reducemoose") ~ 1, TRUE ~ 0),
    pen = case_when(str_detect(trt, "pen") ~ 1, TRUE ~ 0),
    feed = case_when(str_detect(trt, "feed") ~ 1, TRUE ~ 0),
    transplant = case_when(str_detect(trt, "transplant") ~ 1, TRUE ~ 0)
  )

# I think you want to write this as a means parameterization
# ~ -1 + reducewolves + ... would get it done and give you the mean per category

# I am struggling a bit with the inference here, could be personal, but just not
# sure what is being gained over the original analysis
ind.eff <- eff.draws %>%
  filter(name == "Rate of increase (r)") %>%
  group_by(.draw) %>%
  do(tidy(lm(delta.i.lr ~ -1 + reducewolves + sterilizewolves + reducemoose + pen + feed, data = .))) # %>%
# filter(!term %in% "(Intercept)")

ind.eff <- eff.draws %>%
  filter(name == "Rate of increase (r)") %>%
  group_by(.draw) %>%
  do(tidy(lm(delta.i.lr ~ reducewolves + sterilizewolves + reducemoose + pen + feed, data = .))) %>%
  filter(!term %in% "(Intercept)")

ind.eff %>%
  group_by(term) %>%
  summarise(
    eff = median(estimate)
  )

order <- ind.eff %>%
  group_by(term) %>%
  summarize(med = median(estimate))

# Yikes, bimodal and all. This doesn't look consistent with earlier graphs, but
# I may just remembering it wrong
ind.eff.plot <- ggplot(
  ind.eff %>%
    left_join(order),
  aes(x = estimate, y = fct_reorder(term, med), fill = term)
) +
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
  labs(x = "Change in rate of increase", y = "Recovery measure(s)", title = "a) Individual Treatment Effects", subtitle = "Partitioned using regression analysis, assuming effects are additive") +
  xlim(-0.19, 0.22) +
  scale_fill_manual(values = cols[c(1:6)])
ind.eff.plot

# #ggsave(here::here("plots","ind_effects.png"), width=5, height=6, bg="white")


ind.eff %>%
  group_by(term) %>%
  summarise(
    delta.l = median(estimate, na.rm = TRUE) %>% round(2),
    lower = quantile(estimate, 0.05, na.rm = TRUE) %>% round(2),
    upper = quantile(estimate, 0.95, na.rm = TRUE) %>% round(2)
  ) %>%
  arrange(-delta.l) %>%
  mutate(delta.lambda = paste0(delta.l, " (", lower, "-", upper, ")")) %>%
  select(Treatment = term, delta.lambda) %>%
  gt()


## ----cons sims, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------
n.sims <- 900 # 1000
start.pop <- 100
# sim.ref   <- demog.draws.combotreat%>%ungroup%>%filter(trt=="Reference")%>%filter(lambda<quantile(lambda,0.95),lambda>quantile(lambda,0.05))%>%sample_n(n.sims)%>%pull(lambda)

# Why not just use random draws of lambda/treatment? Curious why the geometric mean transformation is desirable.
sim.ref <- demog.draws %>%
  # filter(trt=="Reference" & totNMF<100)%>%
  filter(trt == "Reference" & totNMF < 150) %>%
  group_by(.draw, trt, herd) %>%
  summarise(
    lambda = exp(mean(log(lambda)))
  ) %>% ## geo mean per herd-treatment-draw
  group_by(.draw, trt) %>%
  summarise(
    lambda = exp(mean(log(lambda))) # Not sure this needs to be done a second time, probably just call mean and summarize the thing? Anyway, I don't think doing this twice is correct, but could be wrong
  ) %>% ## geo mean lambda per treatment-draw
  ungroup() %>%
  filter(lambda < quantile(lambda, 0.95), lambda > quantile(lambda, 0.05)) %>% # I would ask how this is defensible. Just a devil's advocate type comment
  sample_n(min(c(dplyr::n(), n.sims))) %>%
  pull(lambda)

sim.trt <- ind.eff %>%
  mutate(estimate = exp(estimate) - 1) %>% ## back to lambda
  group_by(term) %>%
  filter(estimate < quantile(estimate, 0.95), estimate > quantile(estimate, 0.05)) %>%
  ungroup() %>%
  select(.draw, term, estimate) %>%
  pivot_wider(names_from = "term", values_from = "estimate") %>%
  drop_na() %>%
  sample_n(min(c(dplyr::n(), n.sims))) # This was easy to make too big

sim.trt$reducewolves %>% median()
sim.trt$feed %>% median()

year.end <- 14
sim.df <- list()
for (i in 1:n.sims) {
  cons.sim.i <- tibble(yr = 1:(11 + year.end)) %>%
    dplyr::mutate(
      reference = c(start.pop * (sim.ref[i])^yr),
      reducewolves = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * (sim.ref[i] + sim.trt$reducewolves[i])^(yr - 10))
      ),
      feed = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * (sim.ref[i] + sim.trt$feed[i])^(yr - 10))
      ),
      pen = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * (sim.ref[i] + sim.trt$pen[i])^(yr - 10))
      ),
      reducemoose = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * (sim.ref[i] + sim.trt$reducemoose[i])^(yr - 10))
      ),
      sterilizewolves = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * (sim.ref[i] + sim.trt$sterilizewolves[i])^(yr - 10))
      ),
      `reducewolves+feed` = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * (sim.ref[i] + sim.trt$feed[i] + sim.trt$reducewolves[i])^(yr - 10))
      ),
      `reducewolves+pen` = case_when(
        yr %in% 1:10 ~ reference,
        yr %in% 11:30 ~ c(reference[10] * (sim.ref[i] + sim.trt$pen[i] + sim.trt$reducewolves[i])^(yr - 10))
      ),
      sim = 1 # Did you mean to hardcode this value?
    )

  cons.sim.i$sim <- i # Oh, ok, this could have been done above
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
    median = median(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE),
    upper = quantile(abund, 0.95, na.rm = T), # Code didn't work for me without adding na.rm
    lower = quantile(abund, 0.05, na.rm = T),
    inc = (mean(inc, na.rm = T) * 100) %>% round(0),
    ext = (mean(ext, na.rm = T) * 100) %>% round(0)
  ) %>%
  filter(
    !(yr < 10 & name != "reference"),
    yr <= (11 + year.end)
  ) %>%
  mutate(
    trt = paste0(name, " (", inc, ", ", ext, ", ", lower %>% round(0), "-", upper %>% round(0), ")"),
    yr = yr - 10
  ) %>%
  mutate(yr.shift = case_when(
    name %in% "sterilizewolves" & yr == year.end ~ yr + (2.5 * nudge),
    name %in% "reducemoose" & yr == year.end ~ yr + (3.5 * nudge),
    name %in% "reference" & yr == year.end ~ yr - nudge,
    name %in% "pen" & yr == year.end ~ yr,
    name %in% "feed" & yr == year.end ~ yr + nudge,
    name %in% "reducewolves" & yr == year.end ~ yr + (2 * nudge),
    name %in% "reducewolves+pen" & yr == year.end ~ yr + (3 * nudge),
    name %in% "reducewolves+feed" & yr == year.end ~ yr + (4 * nudge),
    TRUE ~ yr
  ))

# a <-sim.df.plot%>%
#   filter(yr == last(yr))%>%pull()

recov.cols <- c(cols[1:3], "black", cols[5:8])
recov.sims.plot <- ggplot() +
  annotate("rect",
    xmin = -10, xmax = year.end + 1, ymin = 0, ymax = 20,
    alpha = .1, fill = "black"
  ) +
  annotate("text", x = -4, y = 5, label = "Functionally extirpated") +
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
      filter(yr == last(yr)),
    aes(color = trt, label = trt, x = yr, y = median),
    size = 4,
    direction = "y",
    xlim = c(year.end + 2, 35),
    hjust = 0,
    segment.size = .7,
    segment.alpha = .5,
    segment.linetype = "dotted",
    box.padding = .4,
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


# I am pretty skeptical that we will realize such fast growth if shoot wolves and feed
# It would be fun to chat about the shape of that line, the interpretation of a categorical effect in relation to time, etc
# This will be splashy for sure
recov.sims.plot



## plot individual treatments together
recov.together <- ind.eff.plot + recov.sims.plot + plot_layout(widths = c(1.3, 1.6))
# ggsave(plot = recov.together, here::here("plots/recov.together.png"), width = 13, height = 6, bg = "white")


## ----map, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------------
# I don't have spatial data in the repo, which is probably a good thing

## Prep Herd Bounds
herds <- st_read(here::here("data/Spatial/herds/u_bc_herds_2021_CL.shp")) %>%
  st_transform(3005) %>%
  select(herd = HERD_NAME) %>%
  mutate(herd = case_when(
    herd %in% "Narraway" ~ "Bearhole Redwillow",
    herd %in% "Moberly" ~ "Klinse-Za",
    herd %in% "Scott" ~ "Scott West",
    herd %in% "Frisby Boulder" ~ "Frisby-Boulder",
    herd %in% "Purcell Central" ~ "Purcells Central",
    TRUE ~ herd
  )) %>%
  rbind(st_read(here::here("data/Spatial/herds/Caribou_Range.shp")) %>% # Still very scared of rbind, but it only worked on matrices when I started this=)
    st_transform(3005) %>%
    select(herd = SUBUNIT) %>%
    mutate(herd = case_when(
      herd %in% "Narraway" ~ "Narraway AB",
      herd %in% "Jasper" ~ "Brazeau",
      herd %in% "Redrock-Prairie Creek" ~ "Redrock/Prairie Creek",
      TRUE ~ herd
    ))) %>%
  st_simplify(
    preserveTopology = FALSE,
    dTolerance = 1000
  )

dist.herd <- read_csv("data/Spatial/disturbance/from_Emily/2022-01-27/herds_propdist_20220126_forCL.csv") %>%
  mutate(herd = case_when(
    herd %in% c("Hart South", "Hart North") ~ "Hart Ranges",
    herd %in% "Narraway BC" ~ "Bearhole Redwillow",
    herd %in% "Central Selkirks" ~ "Nakusp",
    TRUE ~ herd
  )) %>%
  group_by(herd) %>%
  summarise(across(anthro_prop_dist:all_prop_dist, mean))

# This line reads like my comment above about is the same data in both
# I dist.herd and all of its data being combined with itself with the herd name
#  changed, is that the intent, two copies of the data with one called Nakusp
#  and one called Duncan?
dist.herd <- dist.herd %>%
  rbind(dist.herd %>% filter(herd %in% "Nakusp") %>% mutate(herd = "Duncan"))



herds <- dist.herd %>%
  dplyr::select(herd, human = anthro_prop_dist) %>%
  mutate(herd = case_when(
    herd %in% "Frisby Boulder" ~ "Frisby-Boulder",
    herd %in% "Purcell Central" ~ "Purcells Central",
    TRUE ~ herd
  )) %>%
  left_join(demog.draws %>%
    filter(trt %in% c("Reference")) %>%
    mutate(
      herd = case_when(
        herd %in% c("Hart South", "Hart North") ~ "Hart Ranges",
        TRUE ~ herd
      )
    ) %>%
    group_by(herd) %>%
    filter(yrs %in% (max(yrs) - 19):max(yrs)) %>%
    summarise(
      lambda = exp(mean(log(lambda)))
    )) %>% ## geo mean per herd
  left_join(
    ext.yr %>%
      mutate(ext = 1) %>%
      ungroup() %>%
      select(herd, ext) %>%
      rbind(tibble(herd = "Scott West", ext = 1))
  ) %>%
  mutate(lambda = case_when(herd %in% "Scott West" ~ 0.95, TRUE ~ lambda)) %>% ## don't have Scott west. Extirpated
  left_join(herds) %>% # This doesn't seem necessary, but really not sure
  drop_na(human) %>%
  mutate(
    human = round(human * 100, 0) %>% as.integer(),
    lambda.class = case_when(
      lambda <= 0.95 ~ "<0.95",
      lambda > 0.95 & lambda <= 0.99 ~ "0.95-0.99",
      lambda > 0.99 & lambda < 1.01 ~ "0.99-1.01",
      lambda >= 1.01 ~ ">1.01"
    ),
    lambda.class2 = case_when(
      lambda <= 0.98 ~ "declining",
      lambda > 0.98 & lambda <= 1.02 ~ "stable",
      lambda > 1.02 ~ "increasing",
      is.na(lambda) & ext == 1 ~ "declining"
    )
  ) %>%
  st_as_sf() %>%
  st_transform(4326)


# mapview(herds,zcol="human")
# mapview(herds,zcol="lambda.class2")

## Prep inset map
## Cities
cities <- st_read(here::here("data/Spatial/administrative/places.shp")) %>%
  filter(NAME %in% c("C", "VANCOUVER", "Prince George", "Fort St John")) %>%
  mutate(Name = str_to_title(NAME)) %>%
  select(Name, geometry) %>%
  st_transform(4326) %>%
  rbind(tibble(Name = "Banff", Y = 51.18, X = -115.56) %>% st_as_sf(coords = c("X", "Y"), crs = 4326)) %>%
  st_transform(3005)

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


library(ggspatial)
## inset
## custom crs to keep things straight
cust.crs <- "+proj=aea +lat_0=50 +lon_0=-114.9 +lat_1=49 +lat_2=50.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs"
inset <- ggplot() +
  geom_sf(data = pnw %>% st_transform(cust.crs), fill = "grey80", color = NA) +
  geom_sf(data = pnw %>% st_transform(cust.crs), size = 0.5, fill = NA, color = "grey30") +
  layer_spatial(sf::st_bbox(herds %>% st_buffer(100000) %>% st_transform(cust.crs)), fill = NA, linetype = "dashed", color = "grey99") +
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

register_google("AIzaSyCOwGx2D77XOqRgGhKmcb5F4Kt_S61tCLI")
# set_defaults(map_service = "osm", map_type = "terrain_bg")

bmap.big <- basemap_raster(
  ext = herds %>% st_buffer(200000) %>% st_transform(3857),
  map_res = 1, map_type = "terrain_bg"
) %>% projectRaster(crs = cust.crs)


map <- ggRGB(bmap.big, r = 1, g = 2, b = 3) +
  theme_bw() +
  geom_sf(data = pnw %>% st_transform(cust.crs), size = 1, fill = NA, linetype = "dashed") +
  geom_sf(data = herds, aes(fill = lambda.class2), inherit.aes = FALSE, alpha = 0.7) +
  geom_sf(data = herds %>% filter(ext %in% 1), aes(color = "Extirpated"), fill = NA, inherit.aes = FALSE, alpha = 0.7) +
  geom_sf_text(data = st_centroid(herds), aes(label = human), inherit.aes = FALSE, size = 2.5, color = "white") +
  geom_sf_text(data = st_centroid(herds %>% filter(lambda.class2 == "stable")), aes(label = human), inherit.aes = FALSE, size = 2.5, color = "black") +
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
  labs(fill = "Population growth\nw/o intervention", title = "a) Southern Mountain Caribou", color = "") +
  scale_fill_viridis_d() +
  scale_color_manual(values = c("Extirpated" = "red")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))


## plot together
map_combo <- map + abundance.all.plot + plot_layout(widths = c(1.7, 1.3))

# ggsave(plot = map_combo, here::here("plots/map.together.png"), width = 12, height = 8, bg = "white")




ggplot(herds %>% tibble(), aes(x = human, fill = lambda.class2)) +
  geom_histogram(alpha = 0.5)


## ----ecotype, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------------

demog.draws.combotreat.ecotype <- demog.draws %>%
  # filter(!(intensity%in%"low"|totNMF<20))%>%
  filter(year >= 2010) %>%
  # group_by(.draw, trt, herd) %>% # Seems unnecessary
  # summarise(lambda = exp(mean(log(lambda)))) %>% ## geo mean per herd-treatment-draw
  left_join(ecotype) %>%
  group_by(.draw, trt, ECCC) %>%
  summarise(
    lambda = exp(mean(log(lambda))) # Probably only need to do this once
  ) %>% ## geo mean lambda per treatment-draw
  mutate(trt = case_when(is.na(trt) ~ "Reference", TRUE ~ trt)) %>% # One mutate is fine, no reason to expand like this and this is another case where case_when probably isn't the right tool
  mutate(group = case_when(trt %in% "Reference" ~ "Reference", TRUE ~ "Treatment"))

order <- demog.draws.combotreat.ecotype %>%
  group_by(trt) %>%
  summarize(med = median(lambda))

demog.draws.combotreat.ecotype %>%
  filter(!trt %in% "transplant") %>%
  left_join(order) %>%
  ggplot(aes(x = lambda, y = fct_reorder(trt, med), fill = ECCC)) +
  geom_density_ridges( # scale = 1.5,
    # scale = 1.3,
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
    legend.title = element_text(size = 15)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  xlim(0.8, 1.3) +
  labs(
    x = "Population growth rate",
    y = "Recovery measure(s)",
    fill = "ECCC Ecotype",
    title = "Population Growth Rate by Recovery Measure"
  ) +
  scale_fill_manual(values = cols)


# ggsave(here::here("plots", "appendix", "lambda_treatments_ecotype.png"), width = 8, height = 7, bg = "white")

eff.draws.ecotype <- eff.draws %>%
  dplyr::select(.draw:delta.i.lr) %>%
  filter(name == "Rate of increase (r)") %>%
  left_join(ecotype) %>%
  pivot_longer(ECCC:Heard_Vagt1998, names_to = "ecotype")


ggplot() +
  geom_density_ridges(
    data = eff.draws.ecotype %>%
      filter(!trt %in% "transplant") %>%
      filter(ecotype == "ECCC") %>%
      group_by(name, trt, .draw, value) %>%
      summarise(delta = median(delta.i.lr)),
    aes(x = delta, y = trt, fill = value),
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
    legend.title = element_text(size = 15)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Delta r", y = "Recovery measure(s)", title = "Before-After Assessment of Effectivness", fill = "ECCC Ecotype") +
  xlim(-0.2, 0.4) +
  scale_fill_manual(values = cols)

# ggsave(here::here("plots", "appendix", "baci_treatments_ecotype.png"), width = 8, height = 7, bg = "white")

ggplot() +
  geom_density_ridges(
    data = eff.draws.ecotype %>%
      filter(!trt %in% "transplant") %>%
      filter(ecotype == "COSEWIC") %>%
      group_by(name, trt, .draw, value) %>%
      summarise(delta = median(delta.i.lr)),
    aes(x = delta, y = trt, fill = value),
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
    legend.title = element_text(size = 15)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "delta", y = "Recovery measure(s)", title = "COSEWIC Grouping: Before-After Assessment of Effectivness", fill = "COSEWIC Ecotype") +
  xlim(-0.2, 0.4)



ggplot() +
  geom_density_ridges(
    data = eff.draws.ecotype %>%
      filter(!trt %in% "transplant") %>%
      filter(ecotype == "Heard_Vagt1998") %>%
      group_by(name, trt, .draw, value) %>%
      summarise(delta = median(delta.i.lr)),
    aes(x = delta, y = trt, fill = value),
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
    legend.title = element_text(size = 15)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "delta", y = "Recovery measure(s)", title = "Before-After Assessment of Effectivness", fill = "Heard_Vagt1998 Ecotype") +
  xlim(-0.2, 0.4)
# ggsave(here::here("plots", "appendix", "baci_treatments_HeardVagt1998.png"), width = 8, height = 7, bg = "white")








# Again, this is an effects parameterization and so is relative to the intercept
# We could do the means parameterization again with the -1 leading
ind.eff.ecotype <- eff.draws %>%
  filter(name == "Rate of increase (r)" & !trt %in% "transplant") %>%
  left_join(ecotype) %>%
  group_by(.draw) %>%
  do(tidy(lm(delta.i ~ -1 + reducewolves + sterilizewolves + reducemoose + pen + feed + ECCC, data = .))) # %>%
# filter(!term %in% "(Intercept)")

# Not sure if the means is going to get it done for you this time. You may want to go with the effects and then calculate the various
#  contrats of interest. For example, the effect of reducing wolves in the Northern Group.

order <- ind.eff.ecotype %>%
  group_by(term) %>%
  summarize(med = median(estimate))

ggplot(ind.eff.ecotype %>% left_join(order), aes(x = estimate, y = fct_reorder(term, med), fill = term)) +
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
  labs(x = "Delta population growth", y = "Recovery measure", title = "Individual Treatment Effects w/ Ecotype", subtitle = "Partitioned using regression analysis, assuming effects are additive") +
  xlim(-0.15, 0.25) +
  scale_fill_manual(values = cols[c(1:8)])

# ggsave(here::here("plots", "appendix", "ind_eff_ecotype.png"), width = 8, height = 7, bg = "white")


ind.eff.ecotype %>%
  group_by(term) %>%
  summarise(
    eff = median(estimate) %>% round(3),
    lower = quantile(estimate, 0.05) %>% round(3),
    upper = quantile(estimate, 0.95) %>% round(3)
  ) %>%
  arrange(-eff) %>%
  mutate(delta.lambda = paste0(eff, " (", lower, "-", upper, ")")) %>%
  select(Treatment = term, delta.lambda) %>%
  gt() # %>%
# gtsave(here::here("tables", "appendix", "ind_eff_ecotype.rtf"))




ind.eff.ecotype.raw <- eff.draws %>%
  filter(name == "Rate of increase (r)" & !trt %in% "transplant") %>%
  left_join(ecotype) %>%
  dplyr::select(herd, trt, ECCC, COSEWIC, Heard_Vagt1998, reducewolves:feed, delta.i.lr) %>%
  dplyr::group_by(across(herd:feed)) %>%
  summarise(delta.i.lr = median(delta.i.lr)) %>%
  ungroup() %>%
  mutate(ECCC = str_split(ECCC, " ", simplify = TRUE)[, 1])


ggplot() +
  geom_boxplot(data = ind.eff.ecotype.raw, aes(x = ECCC, y = delta.i.lr), outlier.alpha = 0) +
  geom_jitter(data = ind.eff.ecotype.raw, aes(x = ECCC, y = delta.i.lr, color = trt), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "ECCC Recovery Ecotype", y = "Delta r")


ggplot() +
  geom_boxplot(data = ind.eff.ecotype.raw, aes(x = ECCC, y = delta.i.lr), outlier.alpha = 0) +
  geom_jitter(data = ind.eff.ecotype.raw, aes(x = ECCC, y = delta.i.lr, color = trt), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(vars(trt)) +
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
  labs(x = "ECCC Recovery Ecotype", y = "Delta r")
# ggsave(here::here("plots", "appendix", "trt_eff_boxplot_ecotype.png"), width = 8, height = 7, bg = "white")

ggplot() +
  geom_boxplot(data = ind.eff.ecotype.raw, aes(x = COSEWIC, y = delta.i.lr), outlier.alpha = 0) +
  geom_jitter(data = ind.eff.ecotype.raw, aes(x = COSEWIC, y = delta.i.lr, color = trt), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(vars(trt)) +
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
  labs(x = "COSEWIC Recovery Ecotype", y = "Delta r")

ggplot() +
  geom_boxplot(data = ind.eff.ecotype.raw, aes(x = Heard_Vagt1998, y = delta.i.lr), outlier.alpha = 0) +
  geom_jitter(data = ind.eff.ecotype.raw, aes(x = Heard_Vagt1998, y = delta.i.lr, color = trt), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(vars(trt)) +
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
  labs(x = "Heard_Vagt1998 Recovery Ecotype", y = "Delta r")


ind.eff.ecotype.raw2 <- eff.draws %>%
  filter(name == "Rate of increase (r)" & !trt %in% "transplant") %>%
  left_join(ecotype) %>%
  dplyr::select(herd, trt, ECCC, COSEWIC, Heard_Vagt1998, reducewolves:feed, delta.i.lr) %>%
  dplyr::group_by(across(herd:feed)) %>%
  summarise(delta.i.lr = median(delta.i.lr)) %>%
  pivot_longer(reducewolves:feed) %>%
  filter(value == 1) %>%
  mutate(ECCC = str_split(ECCC, " ", simplify = TRUE)[, 1])
ggplot() +
  geom_boxplot(data = ind.eff.ecotype.raw2, aes(x = ECCC, y = delta.i.lr), outlier.alpha = 0) +
  geom_jitter(data = ind.eff.ecotype.raw2, aes(x = ECCC, y = delta.i.lr, color = trt), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(vars(name)) +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "ECCC Recovery Ecotype", y = "Delta r") +
  guides(colour = guide_legend(nrow = 3))
# ggsave(here::here("plots", "appendix", "ind.trt_eff_boxplot_ecotype.png"), width = 8, height = 7, bg = "white")


## ----for pc, message=FALSE, warning=FALSE, eval=FALSE, include=FALSE---------------------------------------------------------------------------------------
## ##For PC contract..don't review
##
## ##Gather draws
## pc <- out %>%
##   gather_draws(S[i,j],R_adj[i,j],totAdultsMF[i,j])%>%
##   ungroup%>%
##   pivot_wider(names_from=.variable, values_from=.value)
##
##
## pc <- pc%>%
##   left_join(out%>%
##               gather_draws(N[i,j,a])%>%
##                ungroup%>%
##               pivot_wider(names_from=.variable, values_from=.value)%>%
##               select(i,j,.draw,a,N)%>%
##               pivot_wider(names_from=a, values_from=N,names_prefix = "N"), by=c("i","j",".draw"))%>%
##   left_join(yr_df%>%rename(j=yr_idx), by="j")%>%
##   left_join(trt%>%distinct(herd_num, herd)%>%select(i=herd_num,herd),by=c("i"))%>%
##   left_join(treatment.combos)%>%
##   left_join(trt%>%distinct(herd, year, intensity)%>%select(herd,yrs=year, intensity)%>%filter(intensity%in%"low"),by=c("herd","yrs"))%>%
##   mutate(trt=case_when(is.na(trt)~"Reference", TRUE~trt))
##
##
## ##remove first year lambda for each herd, as lambda==1
## pc  <- pc %>%
##   group_by(herd, .draw)%>%
##   arrange(yrs)%>%
##   slice(-1)%>%
##   ungroup
##
##
## ##filter herds to a period where not extirpated
## pc.trim <- tibble()
## herds<- unique(demog$herd)
## for(i in 1:length(herds)){
##
##   if(herds[i]%in%ext.yr$herd){
##     yr <-ext.yr%>%
##       dplyr::filter(herd==!!herds[i])
##     a <- pc %>%
##       dplyr::filter(herd==!!herds[i] & yrs<yr$yrs)}
##
##   if(!herds[i]%in%ext.yr$herd){
##     a <- pc %>%
##       dplyr::filter(herd==!!herds[i])}
##
##   pc.trim <- bind_rows(a,pc.trim)
## }
##
## pc  <- pc.trim
## rm(pc.trim)
##
## ##calculate N calves corrected to March
## pc <- pc%>%
##   mutate(N1_March=R_adj*totAdultsMF*0.5)
##
##
## ##export for Tal
## write_csv(pc%>%filter(herd%in%c("A La Peche", "Brazeau", "Maligne","Tonquin"),
##                       yrs>=2000)%>%select(.draw:intensity),"IPM_posteriors_AB.csv")






#### CL SIMPLE EXAMPLE

### JN APPROACH
# No need to create a new object, just start with demog.draws and move forward
demog.draws |>
  group_by(.draw, trt, herd) |> # Compute transformation at each draw
  summarise(
    lambda = geo_mean(lambda),
    .groups = "drop"
  ) |>
  mutate(
    group = "Treatment",
    group = replace(group, trt == "Reference", "Reference")
  ) |>
  filter(trt != "transplant") %>%
  ggplot(aes(x = log(lambda), y = fct_reorder(trt, lambda), fill = group)) +
  geom_density_ridges(
    rel_min_height = .01,
    size = 0.25,
    alpha = 0.9
  ) +
  theme_ipsum() +
  geom_point(
    data = demog.draws.combotreat.rug %>%
      # left_join(order) %>%
      filter(trt != "transplant"),
    aes(y = fct_reorder(trt, lambda), x = log(lambda)),
    shape = "|"
  ) +
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
  xlim(-0.25, 0.3) +
  labs(
    x = "Instantaneous rate of increase (r)",
    y = "Recovery measure(s)",
    fill = "",
    title = "Instantaneous Rate of Increase by Recovery Measure"
  ) +
  scale_fill_manual(values = c(cols[c(3, 1)]))



### CL APPROACH
# No need to create a new object, just start with demog.draws and move forward
demog.draws |>
  group_by(.draw, trt, herd) |> # Compute transformation at each draw
  summarise(
    lambda = geo_mean(lambda),
    .groups = "drop"
  ) |>
  group_by(.draw, trt) |> # Compute transformation at each draw
  summarise(
    lambda = median(lambda),
    .groups = "drop"
  ) |>
  mutate(
    group = "Treatment",
    group = replace(group, trt == "Reference", "Reference")
  ) |>
  filter(trt != "transplant") %>%
  ggplot(aes(x = log(lambda), y = fct_reorder(trt, lambda), fill = group)) +
  geom_density_ridges(
    rel_min_height = .01,
    size = 0.25,
    alpha = 0.9
  ) +
  theme_ipsum() +
  geom_point(
    data = demog.draws.combotreat.rug %>%
      # left_join(order) %>%
      filter(trt != "transplant"),
    aes(y = fct_reorder(trt, lambda), x = log(lambda)),
    shape = "|"
  ) +
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
  xlim(-0.25, 0.3) +
  labs(
    x = "Instantaneous rate of increase (r)",
    y = "Recovery measure(s)",
    fill = "",
    title = "Instantaneous Rate of Increase by Recovery Measure"
  ) +
  scale_fill_manual(values = c(cols[c(3, 1)]))



### JN approach visualized
demog.draws |>
  group_by(.draw, trt, herd) |> # Compute transformation at each draw
  summarise(
    lambda = geo_mean(lambda),
    .groups = "drop"
  ) |>
  mutate(
    group = "Treatment",
    group = replace(group, trt == "Reference", "Reference")
  ) |>
  filter(trt != "transplant") %>%
  ggplot(aes(x = log(lambda), y = fct_reorder(trt, lambda), color = herd, fill = herd)) +
  geom_density_ridges(
    rel_min_height = .01,
    size = 0.25,
    alpha = 0.1
  ) +
  theme_ipsum() +
  geom_point(
    data = demog.draws.combotreat.rug %>%
      # left_join(order) %>%
      filter(trt != "transplant"),
    aes(y = fct_reorder(trt, lambda), x = log(lambda)),
    shape = "|"
  ) +
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
  xlim(-0.25, 0.3) +
  labs(
    x = "Instantaneous rate of increase (r)",
    y = "Recovery measure(s)",
    fill = "",
    title = "Instantaneous Rate of Increase by Recovery Measure"
  )
