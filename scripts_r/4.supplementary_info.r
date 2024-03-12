## ----render, eval=FALSE,include=FALSE-----------------------------------------------------------------------------------------------
## rmarkdown::render(here::here("plots", "appendix", "CaribouIPM_BCAB_SM.Rmd"),
##   output_file = "README.md"
## )
## 
## knitr::purl(
##   input = here::here("plots", "appendix", "CaribouIPM_BCAB_SM.Rmd"),
##   output = here::here("scripts_r", "4.supplementary_info.r")
## )


## ----Load packages and data, results='hide', message=FALSE, warning=FALSE-----------------------------------------------------------
library(renv)
library(here)
library(RColorBrewer)
library(ggridges)
library(hrbrthemes)
library(broom.mixed)
library(lme4)
library(knitr)
library(ggrepel)
library(gt)
library(tidyverse)

# display.brewer.pal(8, "Accent")
cols <- RColorBrewer::brewer.pal(8, "Accent")
hd <- read.csv(here::here("data/clean/blueprint.csv"))
herds.keep <- hd %>%
  filter(herd != "Quintette Full") %>%
  distinct(herd) %>%
  pull()

demog.draws <- read_csv(here::here("tables/draws/demog.draws.csv")) %>%
  filter(herd %in% herds.keep)
demog <- read_csv(here::here("tables/demog.csv")) %>%
  filter(herd %in% herds.keep)
eff.draws <- read_csv(here::here("tables/draws/eff.draws.csv")) %>%
  filter(herd %in% herds.keep, name == "Rate of increase (r)")
eff.draws.app <- read_csv(here::here("tables/draws/eff.draws.app.csv")) %>%
  filter(herd %in% herds.keep)
ecotype <- read.csv(here::here("data/raw/treatment.csv")) %>%
  dplyr::select(herd = Herd, ECCC = ECCC_Recov_Grp, COSEWIC = COSEWIC_Grp, Heard_Vagt1998 = Heard.and.Vagt.1998.grouping) %>%
  distinct() %>%
  mutate(herd = case_when(
    herd %in% "Narraway BC" ~ "Bearhole Redwillow",
    TRUE ~ herd
  )) %>%
  filter(herd %in% herds.keep)
trt_eff_ba_table <- read_csv(here::here("tables", "trt_eff_ba.csv"))
sims <- read_csv(here::here("tables/draws/sims.draws.csv")) %>%
  filter(herd %in% herds.keep)
label.lookup <- read.csv(here::here("data/clean/label_lookup.csv"))

## for first plots
hn <- hd %>%
  dplyr::select(herd, herd_num) %>%
  filter(herd %in% herds.keep)
hn <- hd %>%
  dplyr::select(herd, herd_num) %>%
  filter(herd %in% herds.keep)
trt <- read.csv(here::here("data/clean/treatments.csv")) %>%
  arrange(herd) %>%
  left_join(hn, by = "herd") %>%
  filter(herd %in% herds.keep)

afs <- read.csv(here::here("data/clean/survival.csv")) %>%
  arrange(herd) %>%
  left_join(hn, by = "herd") %>%
  filter(herd %in% herds.keep)
afr <- read.csv(here::here("data/clean/recruitment.csv")) %>%
  arrange(herd) %>%
  left_join(hn, by = "herd") %>%
  filter(herd %in% herds.keep)
counts <- read.csv(here::here("data/clean/counts.csv")) %>%
  arrange(herd) %>%
  left_join(hn, by = "herd") %>%
  filter(herd %in% herds.keep)
ecotype <- read.csv(here::here("data/raw/treatment.csv")) %>%
  dplyr::select(herd = Herd, ECCC = ECCC_Recov_Grp, COSEWIC = COSEWIC_Grp, Heard_Vagt1998 = Heard.and.Vagt.1998.grouping) %>%
  distinct() %>%
  filter(herd %in% herds.keep)

## bearhole naming
counts <- counts %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
trt <- trt %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
afr <- afr %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
afs <- afs %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))
ecotype <- ecotype %>% mutate(herd = case_when(herd %in% "Narraway BC" ~ "Bearhole Redwillow", TRUE ~ herd))


## ----Firstyr, message=FALSE, warning=FALSE------------------------------------------------------------------------------------------
raw.demog <- rbind(
  afr %>% dplyr::select(herd, year, est) %>% mutate(type = "Recruit"),
  afs %>% dplyr::select(herd, year, est) %>% mutate(type = "Surv"),
  counts %>% dplyr::select(herd, year, est = Est_CL) %>% mutate(type = "Count")
) %>%
  left_join(ecotype, by = "herd")

### Identify first year of demographic data for each herd
first.yr <- raw.demog %>%
  dplyr::select(herd, ECCC, first.year = year) %>%
  distinct() %>%
  group_by(herd, ECCC) %>%
  filter(first.year == min(first.year)) %>%
  mutate(count = 1) %>%
  group_by(ECCC) %>%
  arrange(first.year) %>%
  mutate(cumsum = cumsum(count))

ggplot(first.yr, aes(x = first.year, y = cumsum)) +
  geom_step() +
  theme_ipsum() +
  labs(
    x = "Year", y = "Cumulative count",
    title = "First year of demographic data for the 41 SMC herds between 1973-2023",
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
  ) +
  facet_wrap(vars(ECCC), scales = "free_y")

first.yr %>%
  summarize(medianyr = median(first.year))

first.yr.summary <- first.yr %>%
  mutate(prop = cumsum / max(cumsum)) %>%
  filter(prop >= 0.5) %>%
  filter(first.year == min(first.year)) %>%
  ungroup()

kable(first.yr.summary)


## ----r2 plot, message=FALSE, warning=FALSE------------------------------------------------------------------------------------------
trt.plot <- trt %>%
  filter(applied == 1) %>%
  mutate(intensity = replace_na(intensity, "standard")) %>%
  group_by(herd, treatment, intensity) %>%
  count() %>%
  pivot_wider(names_from = "intensity", values_from = "n") %>%
  mutate(
    low = replace_na(low, 0),
    total = sum(low, standard, na.rm = TRUE),
    label = paste0(total, " (", low, ")")
  )

ggplot(trt.plot, aes(y = herd, x = treatment)) +
  geom_tile(aes(fill = total)) +
  geom_text(aes(label = label), color = "white") +
  scale_fill_viridis_c() +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  ) +
  labs(x = "Recovery action", y = "Subpopulation", fill = "Years applied")

ggsave(here::here("plots", "appendix", "treatment_matrix.png"), width = 10, height = 10, bg = "white")


## ----trajectory by ECCC ecotype, echo=TRUE, fig.height=6, fig.width=18, message=FALSE, warning=FALSE--------------------------------
sims.plot <- sims %>%
  filter(.variable %in% c("totNMF", "pred_totNMF")) %>%
  left_join(ecotype, by = "herd") %>%
  dplyr::group_by(yrs, ECCC, .variable, .draw) %>%
  dplyr::summarize(
    mean = sum(.value),
    .groups = "drop"
  ) %>%
  dplyr::group_by(ECCC, yrs, .variable) %>%
  dplyr::summarize(
    LCL = quantile(mean, 0.05) %>% round(0),
    UCL = quantile(mean, 0.95) %>% round(0),
    mean = mean(mean) %>% round(0),
    .groups = "drop"
  ) %>%
  left_join(first.yr.summary %>% dplyr::select(ECCC, first.year), by = "ECCC") %>%
  filter(yrs >= first.year) %>%
  mutate(ECCC = fct_relevel(ECCC, "Northern Group", "Central Group", "Southern Group"))

abundance.all.plot <- ggplot(data = sims.plot %>%
  mutate(.variable = case_when(
    .variable == "totNMF" ~ "Observed",
    TRUE ~ "Status quo"
  ))) +
  geom_ribbon(alpha = 0.3, aes(x = yrs, y = mean, ymin = LCL, ymax = UCL, fill = fct_relevel(.variable, "Status quo", "Observed")), color = NA) +
  geom_line(size = 1, aes(x = yrs, y = mean, ymin = LCL, ymax = UCL, color = fct_relevel(.variable, "Status quo", "Observed"), fill = fct_relevel(.variable, "Status quo", "Observed"))) +
  geom_text(data = sims.plot %>% filter(yrs == 2023) %>%
    mutate(.variable = case_when(
      .variable == "totNMF" ~ "Observed",
      TRUE ~ "Status quo"
    )), aes(label = fct_relevel(.variable, "Status quo", "Observed"), colour = .variable, x = Inf, y = mean), hjust = 0) +
  # geom_jitter(data = ext.yr, size = 1, aes(x = yrs, y = mean), alpha = 0.5) +
  theme_ipsum() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("Year") +
  labs(
    x = "Year", y = "Abundance", title = "Population Trend"
  ) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 100)) +
  scale_x_continuous(breaks = seq(1980, 2020, by = 20)) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    plot.margin = unit(c(1, 5, 1, 1), "lines")
  ) +
  geom_rug(data = trt %>% filter(applied %in% 1) %>% left_join(ecotype, by = "herd"), aes(x = year), sides = "b", length = unit(0.05, "npc"), alpha = 0.05) +
  scale_color_manual(values = cols[c(3, 1)]) +
  # annotate(geom = "text", x = 1974, y = 2600, label = "Recovery\nactions", hjust = "left") +
  # annotate(
  #   geom = "curve", x = 1978, y = 1800, xend = 1985, yend = 300,
  #   curvature = .4, arrow = arrow(length = unit(2, "mm"))
  # ) +
  #   annotate(geom = "text", x = 1975, y = 10000, label = "Demographic\ndata collected", hjust = "left") +
  # annotate(
  #   geom = "curve", x = 1988, y = 10000, xend = 1992, yend = 11300,
  #   curvature = .4, arrow = arrow(length = unit(2, "mm"))
  # ) +
  # annotate(geom = "text", x = 1990, y = 5500, label = "Subpopulation\nextirpation event", hjust = "left") +
  # annotate(
  #   geom = "curve", x =1989, y = 5800, xend = (ext.yr %>% ungroup() %>% filter(herd == "Central Rockies") %>% pull(yrs)), yend = (ext.yr %>% ungroup() %>% filter(herd == "Central Rockies") %>% pull(mean)) - 200,
  #   curvature = -.3, arrow = arrow(length = unit(2, "mm"))
  # ) +
  coord_cartesian(
    clip = "off"
  ) +
  geom_rug(
    data = rbind(
      afr %>% dplyr::select(year, herd, est) %>% mutate(type = "Recruit"),
      afs %>% dplyr::select(year, herd, est) %>% mutate(type = "Surv"),
      counts %>% dplyr::select(year, herd, est = Est_CL) %>% mutate(type = "Count")
    ) %>%
      left_join(ecotype, by = "herd") %>%
      ungroup(),
    aes(x = year), sides = "t", length = unit(0.05, "npc"), alpha = 0.02
  ) +
  facet_wrap(vars(fct_relevel(ECCC, "Northern Group", "Central Group", "Southern Group")), scales = "free_y") +
  theme(panel.spacing = unit(4, "lines"))
abundance.all.plot

ggsave(plot = abundance.all.plot, here::here("plots", "appendix", "abundance_byecotype.png"), width = 18, height = 6, bg = "white")


current.yr <- 2023
cosewic.abund <- sims.plot %>%
  group_by(ECCC) %>%
  dplyr::select(-first.year) %>%
  filter(.variable == "totNMF") %>%
  dplyr::select(-.variable) %>%
  filter(yrs >= current.yr - (9 * 3) - 1) %>% ## start at 3 generations back
  filter(yrs %in% c(current.yr, current.yr - 9, min(yrs))) ## use 3 generations back, or first value if not back that far

cosewic.abund %>%
  gt() %>%
  gtsave(here::here("tables", "appendix", "ecotype_abundance.rtf"))

cosewic.abund %>%
  dplyr::select(ECCC, yrs, mean) %>%
  pivot_wider(values_from = mean, names_from = yrs) %>%
  dplyr::select(ECCC, `1995`, `1996`, `2014`, `2023`) %>%
  mutate(
    decline.10yr = (((`2023` - `2014`) / `2014`) * 100) %>% round(0),
    decline.3generation = case_when(
      ECCC != "Central Group" ~ (((`2023` - `1995`) / `1995`) * 100) %>% round(0),
      ECCC == "Central Group" ~ (((`2023` - `1996`) / `1996`) * 100) %>% round(0)
    )
  ) %>%
  gt() %>%
  gtsave(here::here("tables", "appendix", "ecotype_change.rtf"))


## ----plotbyherd, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------
## prep data
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




## loop through all herds and plot
herds <- unique(demog$herd)
for (i in 1:length(herds)) {
  trt.dat.i <- trt.plot %>%
    filter(herd == !!herds[i]) %>%
    distinct(herd, treatment, y) %>%
    mutate(t = str_sub(treatment, 1, 1))
  if (nrow(trt.dat.i) == 0) {
    trt.dat.i <- add_row(trt.dat.i, herd = herds[i], treatment = NA)
  }

  ggplot() +
    geom_ribbon(data = demog %>% filter(herd == !!herds[i]), aes(x = yrs, y = totNMF, ymin = totNMF.lower, ymax = totNMF.upper), alpha = 0.4, color = NA, fill = cols[3]) +
    geom_line(data = demog %>% filter(herd == !!herds[i]), aes(x = yrs, y = totNMF, ymin = totNMF.lower, ymax = totNMF.upper), size = 1, color = cols[3]) +
    geom_rug(
      data = rbind(afr %>% dplyr::select(herd, year, est) %>% mutate(type = "Recruit") %>% filter(herd == !!herds[i]), afs %>% dplyr::select(herd, year, est) %>% mutate(type = "Surv") %>% filter(herd == !!herds[i]), counts %>% dplyr::select(herd, year, est = Est_CL) %>% mutate(type = "Count") %>% filter(herd == !!herds[i])),
      aes(x = year), sides = "t", length = unit(0.05, "npc"), alpha = 0.5
    ) +
    theme_ipsum() +
    theme(legend.position = "none") +
    ylab("") +
    xlab("Year") +
    labs(x = "Year", y = "Abundance", title = herds[i]) +
    expand_limits(y = 0) +
    scale_x_continuous(
      limits = c(1974, 2026),
      breaks = seq(1980, 2020, by = 20)
    ) +
    theme(
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 15)
    ) +
    geom_point(data = counts %>% filter(herd == !!herds[i]), aes(x = year, y = Est_CL), size = 0.5, alpha = 0.7) +
    geom_linerange(data = counts %>% filter(herd == !!herds[i]) %>% mutate(Est_CL.max = case_when(Est_CL.max > 5000 ~ 5000, herd == "Rainbows" & Est_CL.max > 500 ~ 500, herd == "South Selkirks" & Est_CL.max > 300 ~ 300, TRUE ~ Est_CL.max)), aes(x = year, ymin = Est_CL.min, ymax = Est_CL.max), alpha = 0.5) +
    geom_point(data = trt.plot %>% filter(herd == !!herds[i]), aes(x = year, y = y, group = treatment, color = treatment), size = 0.5) +
    scale_color_manual(values = cols[-4]) +
    geom_text_repel(
      data = trt.dat.i, aes(label = treatment, colour = treatment, x = 2023, y = y * 0.97),
      direction = "y",
      seed = 999,
      force = 0.5,
      nudge_x = 20,
      segment.size = .7,
      segment.alpha = .3,
      segment.linetype = "dotted",
      segment.curvature = -0.1,
      segment.ncp = 3,
      segment.angle = 20
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    geom_text(
      data = trt.plot %>% filter(herd == !!herds[i]) %>%
        group_by(herd, treatment, y) %>%
        summarize(year = mean(year)) %>%
        mutate(t = str_remove(treatment, "reduce ") %>% str_sub(1, 1)),
      aes(label = treatment, x = year, y = y),
      direction = "y",
      size = 3,
      vjust = -1
    )

  ggsave(here::here("plots", "by_herd", "with_treatments", paste0(herds[i] %>% str_replace_all("[:punct:]", " "), ".png")), width = 4, height = 4, bg = "white")


  ggplot() +
    geom_ribbon(data = demog %>% filter(herd == !!herds[i]), aes(x = yrs, y = totNMF, ymin = totNMF.lower, ymax = totNMF.upper), alpha = 0.4, color = NA, fill = cols[3]) +
    geom_line(data = demog %>% filter(herd == !!herds[i]), aes(x = yrs, y = totNMF, ymin = totNMF.lower, ymax = totNMF.upper), size = 1, color = cols[3]) +
    geom_rug(
      data = rbind(afr %>% dplyr::select(herd, year, est) %>% mutate(type = "Recruit") %>% filter(herd == !!herds[i]), afs %>% dplyr::select(herd, year, est) %>% mutate(type = "Surv") %>% filter(herd == !!herds[i]), counts %>% dplyr::select(herd, year, est = Est_CL) %>% mutate(type = "Count") %>% filter(herd == !!herds[i])),
      aes(x = year), sides = "t", length = unit(0.05, "npc"), alpha = 0.5
    ) +
    theme_ipsum() +
    theme(legend.position = "none") +
    ylab("") +
    xlab("Year") +
    labs(x = "Year", y = "Abundance", title = herds[i]) +
    expand_limits(y = 0) +
    scale_x_continuous(
      limits = c(1974, 2026),
      breaks = seq(1980, 2020, by = 20)
    ) +
    theme(
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 15)
    ) +
    geom_point(data = counts %>% filter(herd == !!herds[i]), aes(x = year, y = Est_CL), size = 0.5, alpha = 0.7) +
    geom_linerange(data = counts %>% filter(herd == !!herds[i]) %>% mutate(Est_CL.max = case_when(Est_CL.max > 5000 ~ 5000, herd == "Rainbows" & Est_CL.max > 500 ~ 500, herd == "South Selkirks" & Est_CL.max > 300 ~ 300, TRUE ~ Est_CL.max)), aes(x = year, ymin = Est_CL.min, ymax = Est_CL.max), alpha = 0.5) +
    coord_cartesian(
      clip = "off"
    )

  ggsave(here::here("plots", "by_herd", "without_treatments", paste0(herds[i] %>% str_replace_all("[:punct:]", " "), ".png")), width = 4, height = 4, bg = "white")



  ggplot() +
    geom_ribbon(data = demog %>% filter(herd == !!herds[i]), aes(x = yrs, y = totNMF, ymin = totNMF.lower, ymax = totNMF.upper), alpha = 0.4, color = NA, fill = cols[3]) +
    geom_line(data = demog %>% filter(herd == !!herds[i]), aes(x = yrs, y = totNMF, ymin = totNMF.lower, ymax = totNMF.upper), size = 1, color = cols[3]) +
    theme_ipsum() +
    theme(legend.position = "none") +
    ylab("") +
    xlab("Year") +
    labs(x = "Year", y = "Abundance", title = herds[i]) +
    expand_limits(y = 0) +
    scale_x_continuous(
      limits = c(1974, 2026),
      breaks = seq(1980, 2020, by = 20)
    ) +
    theme(
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 15)
    ) +
    geom_point(data = counts %>% filter(herd == !!herds[i]), aes(x = year, y = Est_CL), size = 0.5, alpha = 0.7) +
    geom_linerange(data = counts %>% filter(herd == !!herds[i]) %>% mutate(Est_CL.max = case_when(Est_CL.max > 5000 ~ 5000, herd == "Rainbows" & Est_CL.max > 500 ~ 500, herd == "South Selkirks" & Est_CL.max > 300 ~ 300, TRUE ~ Est_CL.max)), aes(x = year, ymin = Est_CL.min, ymax = Est_CL.max), alpha = 0.5) +
    coord_cartesian(
      clip = "off"
    )

  ggsave(here::here("plots", "by_herd", "without_rug or treatment", paste0(herds[i] %>% str_replace_all("[:punct:]", " "), ".png")), width = 4, height = 4, bg = "white")
}


## ----summarize by ECCC, echo=TRUE, message=FALSE, warning=FALSE---------------------------------------------------------------------
trt.yrs.app <- demog.draws %>%
  group_by(herd, yrs, trt, intensity) %>%
  mutate(totNMF.median = median(totNMF)) %>% # get average pop size so popsize threshold doesnt split low/standard in some years due to draws being above/below threshold
  ungroup() %>%
  mutate(
    application = case_when(intensity == "low" | totNMF.median < 30 ~ "low", TRUE ~ "standard")
  ) %>%
  ungroup() %>%
  filter(!trt %in% "Reference") %>%
  left_join(ecotype, by = "herd") %>%
  group_by(trt, application, ECCC) %>%
  summarize(
    n_herds = n_distinct(herd),
    n_herd_years = n_distinct(herd, yrs)
  )


eff.draws.app %>%
  filter(name == "r") %>%
  left_join(ecotype, by = "herd") %>%
  group_by(trt, application, ECCC, .draw) %>%
  summarise(eff = mean(eff)) %>%
  group_by(trt, application, ECCC) %>%
  summarise(
    eff.median = median(eff) %>% round(2),
    eff.lower = quantile(eff, 0.05) %>% round(2),
    eff.upper = quantile(eff, 0.95) %>% round(2)
  ) %>%
  left_join(trt.yrs.app, by = c("trt", "application", "ECCC")) %>%
  arrange(trt, application, ECCC) %>%
  ungroup() %>%
  gt() %>%
  gtsave(here::here("tables", "appendix", "ecotype_trt_eff_summary_ECCC.rtf"))


## ----summarize by HV, echo=TRUE, message=FALSE, warning=FALSE-----------------------------------------------------------------------
trt.yrs.app.hv <- demog.draws %>%
  group_by(herd, yrs, trt, intensity) %>%
  mutate(totNMF.median = median(totNMF)) %>% # get average pop size so popsize threshold doesnt split low/standard in some years due to draws being above/below threshold
  ungroup() %>%
  mutate(
    application = case_when(intensity == "low" | totNMF.median < 30 ~ "low", TRUE ~ "standard")
  ) %>%
  ungroup() %>%
  filter(!trt %in% "Reference") %>%
  left_join(ecotype, by = "herd") %>%
  group_by(trt, application, Heard_Vagt1998) %>%
  summarize(
    n_herds = n_distinct(herd),
    n_herd_years = n_distinct(herd, yrs)
  )


eff.draws.app %>%
  filter(name == "r") %>%
  left_join(ecotype, by = "herd") %>%
  group_by(trt, application, Heard_Vagt1998, .draw) %>%
  summarise(eff = mean(eff)) %>%
  group_by(trt, application, Heard_Vagt1998) %>%
  summarise(
    eff.median = median(eff) %>% round(2),
    eff.lower = quantile(eff, 0.05) %>% round(2),
    eff.upper = quantile(eff, 0.95) %>% round(2)
  ) %>%
  left_join(trt.yrs.app.hv, by = c("trt", "application", "Heard_Vagt1998")) %>%
  arrange(trt, application, Heard_Vagt1998) %>%
  ungroup() %>%
  gt() %>%
  gtsave(here::here("tables", "appendix", "ecotype_trt_eff_summary_Heard_Vagt1998.rtf"))


## ----trt eff- BACI w ECCC ecotype, echo=TRUE, fig.height=9, fig.width=13, message=FALSE, warning=FALSE------------------------------
## add ecotype from ECCC and from Heard Vagt 1998
demog.mod.baci <- demog.draws %>%
  ungroup() %>%
  left_join(ecotype) %>%
  dplyr::select(.draw, herd, yrs, r, totNMF, trt, ECCC, Heard_Vagt1998)

#### prep BACI with ECCC grouping
trt.herds.eccc <- demog.mod.baci %>%
  filter(!trt %in% c("Reference", "transplant"), .draw == min(.draw)) %>%
  group_by(herd, trt, ECCC) %>%
  summarise(
    trt.start = min(yrs),
    trt.end = max(yrs),
    dur = n(),
    baci.start = trt.start - dur
  ) %>%
  ungroup()

baci.eccc <- tibble()
for (i in 1:nrow(trt.herds.eccc)) {
  post.trt <- demog.mod.baci %>%
    dplyr::filter(
      herd %in% !!trt.herds.eccc[i, "herd"][[1]],
      yrs %in% !!c(trt.herds.eccc[i, "trt.start"][[1]]:trt.herds.eccc[i, "trt.end"][[1]])
    ) %>%
    dplyr::mutate(CI = 1)

  pre.trt <- demog.mod.baci %>%
    dplyr::filter(
      trt %in% c("Reference"),
      herd %in% !!trt.herds.eccc[i, "herd"][[1]],
      yrs %in% !!c((trt.herds.eccc[i, "trt.start"][[1]] - 11):(trt.herds.eccc[i, "trt.start"][[1]] - 1))
    ) %>%
    dplyr::mutate(CI = 1)

  control <- demog.mod.baci %>%
    ungroup() %>%
    dplyr::filter(
      trt %in% c("Reference"),
      !herd %in% !!trt.herds.eccc[i, "herd"][[1]],
      ECCC %in% !!trt.herds.eccc[i, "ECCC"][[1]],
      yrs %in% !!c((trt.herds.eccc[i, "trt.start"][[1]] - 11):trt.herds.eccc[i, "trt.end"][[1]])
    ) %>%
    dplyr::mutate(CI = 0)

  ## only keep herds with full data
  if (all(nrow(pre.trt) > 0, nrow(post.trt) > 0, nrow(control) > 0)) {
    a <- rbind(pre.trt, post.trt, control) %>%
      dplyr::mutate(
        herd.grp = trt.herds.eccc[i, "herd"][[1]],
        trt.grp = trt.herds.eccc[i, "trt"][[1]],
        BA = case_when(
          yrs < trt.herds.eccc[i, "trt.start"][[1]] ~ 0,
          yrs >= trt.herds.eccc[i, "trt.start"][[1]] ~ 1
        ),
        BA.break = trt.herds.eccc[i, "trt.start"][[1]]
      )
  }

  baci.eccc <- rbind(baci.eccc, a)
}

baci.eccc <- baci.eccc %>%
  mutate(
    grp = paste(herd.grp, trt, sep = "_"),
    facet = paste(herd.grp, trt.grp, sep = "_")
  )

## view data
ggplot() +
  geom_line(data = baci.eccc %>%
    group_by(herd, yrs, CI, facet, trt.grp, herd.grp) %>%
    summarise(
      r = median(r),
      totNMF = median(totNMF)
    ) %>%
    filter(CI == 1), aes(x = yrs, y = r, group = herd), color = "red", alpha = 1, size = 1) +
  geom_line(data = baci.eccc %>%
    group_by(herd, yrs, CI, facet, trt.grp, herd.grp) %>%
    summarise(
      r = median(r),
      totNMF = median(totNMF)
    ), aes(x = yrs, y = r, group = herd), alpha = 0.3) +
  theme_ipsum() +
  ylab("") +
  xlab("Year") +
  labs(x = "Year", y = "Instantaneous rate of increase (r)", title = "BACI set up using ECCC ecotype", subtitle = "control herds shown in grey") +
  expand_limits(y = 0) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  ) +
  geom_vline(data = baci.eccc %>%
    distinct(BA.break, facet, herd.grp, trt.grp), aes(xintercept = BA.break - 0.5), linetype = "dashed") +
  facet_wrap(vars(trt.grp, herd.grp))
ggsave(here::here("plots", "appendix", "baci_setup_ECCC.png"), width = 13, height = 9, bg = "white")


### plot baci delta r
baci.eccc %>%
  group_by(CI, BA, facet, trt.grp, herd.grp, ECCC, .draw) %>%
  summarise(
    r = mean(r)
  ) %>%
  mutate(baci = paste(BA, CI, sep = "_")) %>%
  ungroup() %>%
  dplyr::select(facet:baci) %>%
  pivot_wider(names_from = baci, values_from = r) %>%
  mutate(
    basechange = `1_0` - `0_0`,
    delta.r = `1_1` - `0_1`,
    delta.r.baci = delta.r - basechange
  ) %>%
  group_by(trt.grp, ECCC, .draw) %>%
  summarize(delta.r.baci = median(delta.r.baci)) %>%
  ggplot(aes(x = delta.r.baci, y = trt.grp, fill = ECCC)) +
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
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Rate of increase (r)",
    y = "Recovery measure(s)",
    fill = "ECCC Ecotype",
    title = "Population Growth Rate by Recovery Measure"
  ) +
  scale_fill_manual(values = cols)
ggsave(here::here("plots", "appendix", "baci_effmanual_ECCC.png"), width = 10, height = 7, bg = "white")



## RUN ECCC BACI
baci.eccc.draws <- baci.eccc %>%
  group_by(.draw, trt.grp) %>%
  do(tidy(lmer(r ~ BA + CI + BA * CI + (1 | grp), data = .))) %>%
  filter(term %in% "BA:CI")

ggplot(data = baci.eccc.draws, aes(x = estimate, y = fct_reorder(trt.grp, estimate))) +
  geom_density_ridges(
    scale = .9,
    rel_min_height = .01,
    size = 0.25,
    alpha = 0.9,
    fill = cols[c(1)]
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
  labs(x = "Change in rate of increase", y = "Recovery measure(s)", title = "Before-After-Control-Impact", subtitle = "w/ spatio-temporal matching using ECCC ecotype") +
  xlim(-0.3, 0.6)
ggsave(here::here("plots", "appendix", "baci_all_ECCC.png"), width = 10, height = 7, bg = "white")


## ----trt eff- BACI w HV1998 ecotype, echo=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------
#### Heard and Vagt grouping

trt.herds.hv <- demog.mod.baci %>%
  filter(!trt %in% c("Reference", "transplant"), .draw == min(.draw)) %>%
  group_by(herd, trt, Heard_Vagt1998) %>%
  summarise(
    trt.start = min(yrs),
    trt.end = max(yrs),
    dur = n(),
    baci.start = trt.start - dur
  ) %>%
  ungroup()

baci.hv <- tibble()
for (i in 1:nrow(trt.herds.hv)) {
  post.trt <- demog.mod.baci %>%
    dplyr::filter(
      herd %in% !!trt.herds.hv[i, "herd"][[1]],
      yrs %in% !!c(trt.herds.hv[i, "trt.start"][[1]]:trt.herds.hv[i, "trt.end"][[1]])
    ) %>%
    dplyr::mutate(CI = 1)

  pre.trt <- demog.mod.baci %>%
    dplyr::filter(
      trt %in% c("Reference"),
      herd %in% !!trt.herds.hv[i, "herd"][[1]],
      yrs %in% !!c((trt.herds.hv[i, "trt.start"][[1]] - 11):(trt.herds.hv[i, "trt.start"][[1]] - 1))
    ) %>%
    dplyr::mutate(CI = 1)

  control <- demog.mod.baci %>%
    ungroup() %>%
    dplyr::filter(
      trt %in% c("Reference"),
      !herd %in% !!trt.herds.hv[i, "herd"][[1]],
      Heard_Vagt1998 %in% !!trt.herds.hv[i, "Heard_Vagt1998"][[1]],
      yrs %in% !!c((trt.herds.hv[i, "trt.start"][[1]] - 11):trt.herds.hv[i, "trt.end"][[1]])
    ) %>%
    dplyr::mutate(CI = 0)

  ## only keep herds with full data
  if (all(nrow(pre.trt) > 0, nrow(post.trt) > 0, nrow(control) > 0)) {
    a <- rbind(pre.trt, post.trt, control) %>%
      dplyr::mutate(
        herd.grp = trt.herds.hv[i, "herd"][[1]],
        trt.grp = trt.herds.hv[i, "trt"][[1]],
        BA = case_when(
          yrs < trt.herds.hv[i, "trt.start"][[1]] ~ 0,
          yrs >= trt.herds.hv[i, "trt.start"][[1]] ~ 1
        ),
        BA.break = trt.herds.hv[i, "trt.start"][[1]]
      )
  }

  baci.hv <- rbind(baci.hv, a)
}

baci.hv <- baci.hv %>% mutate(
  grp = paste(herd.grp, trt, sep = "_"),
  facet = paste(herd.grp, trt.grp, sep = "_")
)

ggplot() +
  geom_line(data = baci.hv %>%
    group_by(herd, yrs, CI, facet, trt.grp, herd.grp) %>%
    summarise(
      r = median(r),
      totNMF = median(totNMF)
    ) %>%
    filter(CI == 1), aes(x = yrs, y = r, group = herd), color = "red", alpha = 1, size = 1) +
  geom_line(data = baci.hv %>%
    group_by(herd, yrs, CI, facet, trt.grp, herd.grp) %>%
    summarise(
      r = median(r),
      totNMF = median(totNMF)
    ), aes(x = yrs, y = r, group = herd), alpha = 0.3) +
  theme_ipsum() +
  ylab("") +
  xlab("Year") +
  labs(x = "Year", y = "Instantaneous rate of increase (r)", title = "BACI set up using Heard & Vagt 1998 ecotypes", subtitle = "control herds shown in grey") +
  expand_limits(y = 0) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  ) +
  geom_vline(data = baci.hv %>%
    distinct(BA.break, facet, herd.grp, trt.grp), aes(xintercept = BA.break - 0.5), linetype = "dashed") +
  facet_wrap(vars(trt.grp, herd.grp))
ggsave(here::here("plots", "appendix", "baci_setup_Heard_Vagt1998.png"), width = 13, height = 9, bg = "white")


### plot baci delta r
baci.hv %>%
  group_by(CI, BA, facet, trt.grp, herd.grp, Heard_Vagt1998, .draw) %>%
  summarise(
    r = mean(r)
  ) %>%
  mutate(baci = paste(BA, CI, sep = "_")) %>%
  ungroup() %>%
  dplyr::select(facet:baci) %>%
  pivot_wider(names_from = baci, values_from = r) %>%
  mutate(
    basechange = `1_0` - `0_0`,
    delta.r = `1_1` - `0_1`,
    delta.r.baci = delta.r - basechange
  ) %>%
  group_by(trt.grp, Heard_Vagt1998, .draw) %>%
  summarize(delta.r.baci = median(delta.r.baci)) %>%
  ggplot(aes(x = delta.r.baci, y = trt.grp, fill = Heard_Vagt1998)) +
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
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Rate of increase (r)",
    y = "Recovery measure(s)",
    fill = "Heard & Vagt 1998 Ecotype",
    title = "Population Growth Rate by Recovery Measure"
  ) +
  scale_fill_manual(values = cols)
ggsave(here::here("plots", "appendix", "baci_effmanual_Heard_Vagt1998.png"), width = 10, height = 7, bg = "white")



baci.hv.draws <- baci.hv %>%
  group_by(.draw, trt.grp) %>%
  do(tidy(lmer(r ~ BA + CI + BA * CI + (1 | grp), data = .))) %>%
  filter(term %in% "BA:CI")

ggplot(data = baci.hv.draws, aes(x = estimate, y = fct_reorder(trt.grp, estimate))) +
  geom_density_ridges(
    scale = .9,
    rel_min_height = .01,
    size = 0.25,
    alpha = 0.9,
    fill = cols[c(1)]
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
  labs(x = "Change in rate of increase", y = "Recovery measure(s)", title = "Before-After-Control-Impact", subtitle = "w/ spatio-temporal matching using Heard & Vagt 1998 ecotypes") +
  xlim(-0.3, 0.6)
ggsave(here::here("plots", "appendix", "baci_all_Heard_Vagt1998.png"), width = 10, height = 7, bg = "white")


## ----trt eff- BACI compare w BA, echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------------------
hv.summary <- baci.hv.draws %>%
  group_by(trt.grp) %>%
  summarise(
    delta.l = median(estimate, na.rm = TRUE) %>% round(2),
    lower = quantile(estimate, 0.05, na.rm = TRUE) %>% round(2),
    upper = quantile(estimate, 0.95, na.rm = TRUE) %>% round(2)
  ) %>%
  arrange(-delta.l) %>%
  mutate(delta.r = paste0(delta.l, " (", lower, "-", upper, ")")) %>%
  rename(trt = trt.grp) %>%
  left_join(label.lookup, by = "trt") %>%
  dplyr::select(Treatment = new, `Delta r (BACI): Heard Vagt (1998)` = delta.r)


eccc.summary <- baci.eccc.draws %>%
  group_by(trt.grp) %>%
  summarise(
    delta.l = median(estimate, na.rm = TRUE) %>% round(2),
    lower = quantile(estimate, 0.05, na.rm = TRUE) %>% round(2),
    upper = quantile(estimate, 0.95, na.rm = TRUE) %>% round(2)
  ) %>%
  arrange(-delta.l) %>%
  mutate(delta.r = paste0(delta.l, " (", lower, "-", upper, ")")) %>%
  rename(trt = trt.grp) %>%
  left_join(label.lookup, by = "trt") %>%
  dplyr::select(Treatment = new, `Delta r (BACI): ECCC` = delta.r)


trt_eff_ba_table %>%
  dplyr::select(Treatment = `Recovery action`, `Delta r (BA)` = `Change in instantaneous growth rate (r)`) %>%
  left_join(eccc.summary, by = "Treatment") %>%
  left_join(hv.summary, by = "Treatment") %>%
  rename(`Recovery action` = Treatment) %>%
  gt() %>%
  gtsave(here::here("tables", "appendix", "trt_eff_baci_compare.rtf"))


## ----r and delta r by ECCC ecotype, echo=TRUE, message=FALSE, warning=FALSE---------------------------------------------------------
demog.draws.combotreat.ecotype <- demog.draws %>%
  filter(yrs >= 2010) %>%
  group_by(.draw, trt, herd) %>%
  summarise(r = mean(r, na.rm = TRUE)) %>% ## mean per herd-treatment-draw
  left_join(ecotype) %>%
  group_by(.draw, trt, ECCC) %>%
  summarise(r = mean(r)) %>% ## mean r per treatment-draw
  mutate(trt = case_when(is.na(trt) ~ "Reference", TRUE ~ trt)) %>%
  mutate(group = case_when(trt %in% "Reference" ~ "Reference", TRUE ~ "Treatment"))

order <- demog.draws.combotreat.ecotype %>%
  group_by(trt) %>%
  summarize(med = median(r))

demog.draws.combotreat.ecotype %>%
  filter(!trt %in% "transplant") %>%
  left_join(order) %>%
  ggplot(aes(x = r, y = fct_reorder(trt, med), fill = ECCC)) +
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
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Rate of increase (r)",
    y = "Recovery measure(s)",
    fill = "ECCC Ecotype",
    title = "Population Growth Rate by Recovery Measure"
  ) +
  scale_fill_manual(values = cols)
ggsave(here::here("plots", "appendix", "r_treatments_ecotype.png"), width = 8, height = 7, bg = "white")

eff.draws.ecotype <- eff.draws %>%
  dplyr::select(.draw:delta.r) %>%
  filter(name == "Rate of increase (r)") %>%
  left_join(ecotype) %>%
  pivot_longer(ECCC:Heard_Vagt1998, names_to = "ecotype")

order <- eff.draws.ecotype %>%
  group_by(trt) %>%
  summarize(med = median(delta.r))

ggplot() +
  geom_density_ridges(
    data = eff.draws.ecotype %>%
      filter(!trt %in% "transplant") %>%
      filter(ecotype == "ECCC") %>%
      group_by(name, trt, .draw, value) %>%
      summarise(delta = median(delta.r)) %>%
      left_join(order),
    aes(x = delta, y = fct_reorder(trt, med), fill = value),
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
  scale_fill_manual(values = cols)
ggsave(here::here("plots", "appendix", "baci_treatments_ecotype.png"), width = 8, height = 7, bg = "white")


## ----individual by ECCC ecotype, echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------------------
ind.eff.ecotype <- eff.draws %>%
  filter(name == "Rate of increase (r)" & !trt %in% "transplant") %>%
  left_join(ecotype) %>%
  group_by(.draw) %>%
  do(tidy(lm(delta.r ~ reducewolves + sterilizewolves + reducemoose + pen + feed + ECCC, data = .))) %>%
  filter(!term %in% "(Intercept)")

# glmm gives analagous result
# ind.eff.ecotype <- eff.draws %>%
#   filter(name == "Rate of increase (r)" & !trt %in% "transplant") %>%
#   left_join(ecotype) %>%
#   group_by(.draw) %>%
#   do(tidy(lmer(delta.r ~ reducewolves + sterilizewolves + reducemoose + pen + feed + (1|ECCC), data = .)))%>%
#   filter(!term %in% "(Intercept)", effect=="fixed")


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
  scale_fill_manual(values = cols[c(1:7)])
ggsave(here::here("plots", "appendix", "ind_eff_ecotype_ECCC.png"), width = 8, height = 7, bg = "white")


ind.eff.ecotype %>%
  group_by(term) %>%
  summarise(
    eff = median(estimate) %>% round(3),
    lower = quantile(estimate, 0.05) %>% round(3),
    upper = quantile(estimate, 0.95) %>% round(3)
  ) %>%
  arrange(-eff) %>%
  mutate(delta.r = paste0(eff, " (", lower, "-", upper, ")")) %>%
  dplyr::select(Treatment = term, delta.r) %>%
  gt() %>%
  gtsave(here::here("tables", "appendix", "ind_eff_ecotype_ECCC.rtf"))




ind.eff.ecotype.raw <- eff.draws %>%
  filter(name == "Rate of increase (r)" & !trt %in% "transplant") %>%
  left_join(ecotype) %>%
  dplyr::select(herd, trt, ECCC, COSEWIC, Heard_Vagt1998, reducewolves:feed, delta.r) %>%
  dplyr::group_by(across(herd:feed)) %>%
  summarise(delta.r = median(delta.r)) %>%
  ungroup() %>%
  mutate(ECCC = str_split(ECCC, " ", simplify = TRUE)[, 1])


## ----individual by HV ecotype, echo=TRUE, fig.height=9, fig.width=13, message=FALSE, warning=FALSE----------------------------------
ind.eff.ecotype.hv <- eff.draws %>%
  filter(name == "Rate of increase (r)" & !trt %in% "transplant") %>%
  left_join(ecotype) %>%
  group_by(.draw) %>%
  do(tidy(lm(delta.r ~ reducewolves + sterilizewolves + reducemoose + pen + feed + Heard_Vagt1998, data = .))) %>%
  filter(!term %in% "(Intercept)")


order <- ind.eff.ecotype.hv %>%
  group_by(term) %>%
  summarize(med = median(estimate))

ggplot(ind.eff.ecotype.hv %>% left_join(order), aes(x = estimate, y = fct_reorder(term, med), fill = term)) +
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
  scale_fill_manual(values = cols[c(1:6)])
ggsave(here::here("plots", "appendix", "ind_eff_ecotype_Heard_Vagt1998.png"), width = 8, height = 7, bg = "white")


ind.eff.ecotype.hv %>%
  group_by(term) %>%
  summarise(
    eff = median(estimate) %>% round(3),
    lower = quantile(estimate, 0.05) %>% round(3),
    upper = quantile(estimate, 0.95) %>% round(3)
  ) %>%
  arrange(-eff) %>%
  mutate(delta.r = paste0(eff, " (", lower, "-", upper, ")")) %>%
  dplyr::select(Treatment = term, delta.r) %>%
  gt() %>%
  gtsave(here::here("tables", "appendix", "ind_eff_ecotype_Heard_Vagt1998.rtf"))


## ----raw results ECCC ecotype, echo=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------------
ggplot() +
  geom_boxplot(data = ind.eff.ecotype.raw, aes(x = ECCC, y = delta.r), outlier.alpha = 0) +
  geom_jitter(data = ind.eff.ecotype.raw, aes(x = ECCC, y = delta.r, color = trt), width = 0.2) +
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
  geom_boxplot(data = ind.eff.ecotype.raw, aes(x = ECCC, y = delta.r), outlier.alpha = 0) +
  geom_jitter(data = ind.eff.ecotype.raw, aes(x = ECCC, y = delta.r, color = trt), width = 0.2) +
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
ggsave(here::here("plots", "appendix", "trt_eff_boxplot_ecotype.png"), width = 10, height = 7, bg = "white")


ind.eff.ecotype.raw2 <- eff.draws %>%
  filter(name == "Rate of increase (r)" & !trt %in% "transplant") %>%
  left_join(ecotype) %>%
  dplyr::select(herd, trt, ECCC, COSEWIC, Heard_Vagt1998, reducewolves:feed, delta.r) %>%
  dplyr::group_by(across(herd:feed)) %>%
  summarise(delta.r = median(delta.r)) %>%
  pivot_longer(reducewolves:feed) %>%
  filter(value == 1) %>%
  mutate(ECCC = str_split(ECCC, " ", simplify = TRUE)[, 1])

ggplot() +
  geom_boxplot(data = ind.eff.ecotype.raw2, aes(x = ECCC, y = delta.r), outlier.alpha = 0) +
  geom_jitter(data = ind.eff.ecotype.raw2, aes(x = ECCC, y = delta.r, color = trt), width = 0.2) +
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
ggsave(here::here("plots", "appendix", "ind.trt_eff_boxplot_ecotype.png"), width = 8, height = 7, bg = "white")


## ----raw results HV1998 ecotype, echo=TRUE, message=FALSE, warning=FALSE------------------------------------------------------------
ggplot() +
  geom_boxplot(data = ind.eff.ecotype.raw, aes(x = Heard_Vagt1998, y = delta.r), outlier.alpha = 0) +
  geom_jitter(data = ind.eff.ecotype.raw, aes(x = Heard_Vagt1998, y = delta.r, color = trt), width = 0.2) +
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
  labs(x = "Heard & Vagt 1998 Recovery Ecotype", y = "Delta r")


ggplot() +
  geom_boxplot(data = ind.eff.ecotype.raw, aes(x = Heard_Vagt1998, y = delta.r), outlier.alpha = 0) +
  geom_jitter(data = ind.eff.ecotype.raw, aes(x = Heard_Vagt1998, y = delta.r, color = trt), width = 0.2) +
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
  labs(x = "Heard & Vagt 1998 Recovery Ecotype", y = "Delta r")
ggsave(here::here("plots", "appendix", "trt_eff_boxplot_ecotype_hv1998.png"), width = 10, height = 7, bg = "white")


ind.eff.ecotype.raw2 <- eff.draws %>%
  filter(name == "Rate of increase (r)" & !trt %in% "transplant") %>%
  left_join(ecotype) %>%
  dplyr::select(herd, trt, Heard_Vagt1998, reducewolves:feed, delta.r) %>%
  dplyr::group_by(across(herd:feed)) %>%
  summarise(delta.r = median(delta.r)) %>%
  pivot_longer(reducewolves:feed) %>%
  filter(value == 1)

ggplot() +
  geom_boxplot(data = ind.eff.ecotype.raw2, aes(x = Heard_Vagt1998, y = delta.r), outlier.alpha = 0) +
  geom_jitter(data = ind.eff.ecotype.raw2, aes(x = Heard_Vagt1998, y = delta.r, color = trt), width = 0.2) +
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
  labs(x = "Heard_Vagt1998 Recovery Ecotype", y = "Delta r") +
  guides(colour = guide_legend(nrow = 3))
ggsave(here::here("plots", "appendix", "ind.trt_eff_boxplot_ecotype_hv1998.png"), width = 8, height = 7, bg = "white")

