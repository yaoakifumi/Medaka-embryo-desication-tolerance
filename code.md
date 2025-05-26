# Medaka embryo desication tolerance
This analytic code was used to analyze and visualize data from egg dehydration experiment.

Raw data (two csv file) is deposited in Zenodo repository: https://zenodo.org/records/10784845

## Initial setting

```{R}
# library
library(tidyverse)
# version 2.0.0

library(MASS)
# version 7.3-60

library(coin)
# version 1.4-3 for wilcoxn signed rank test
```

## load data
```{R}
df_time <- read_csv("dehydration_final.csv")
df_plant <- read_csv("plant_final.csv")
```

## Analysis and Visualization for Figure 2A
```{R}
# Remove unused data point
df_time_sort <- df_time %>%
  filter(control_hatch == "OK") %>%
  dplyr::select(experiment_id, period, survival, hatch)

# make Figure 2a
df_time_sort %>%
  ggplot(aes(period, (hatch / 4))) +
  geom_smooth(method = "glm", method.args = list(family="binomial"(link="probit")), 
              se = FALSE, linewidth = 1.5, color = "darkcyan") +
  geom_jitter(width = 0.3, height = 0, alpha = 0.7) +
  xlab("Dehydration period (hour)") + 
  ylab("Hatching rate (%)") + 
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(-0.3, 24.3)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.15)) +
  theme_classic() +
  theme(text = element_text(size = 18, family = "Arial", colour = "black", face = "bold")) +
  theme(axis.ticks=element_line(colour = "black"),
        axis.text=element_text(colour = "black"))

ggsave("Fig1e_BiolLett.png", width = 4.8, height = 3.6)  
```

## calculate median lethal period using probit regression
```{R}
model_time <- glm((cbind(hatch, 4 - hatch)) ~ period, data = df_time_sort, family = binomial(link = "probit"))

summary(model_time)

dose.p(model_time, p = c(0.1, 0.5, 0.9))
```

## Analysis and Visualization for Figure 2B
```{R}
# Remove unused data point
df_plant_sort <- df_plant %>%
  filter(use == "OK")
  dplyr::select(experiment_id, plant, hatch, experiment)

# make Figure 2b
df_plant_sort %>%
  ggplot(aes(plant, hatch * 25)) + 
  geom_bar(stat = "summary", fun = "mean", fill = "darkcyan", width = 0.3) + 
  stat_summary(fun = "mean", fun.min = function(x)mean(x) - sd(x), 
               fun.max = function(x)mean(x)+sd(x), geom = "errorbar", width = 0.15) + 
  geom_jitter(width = 0.05, height = 0, alpha = 0.7) + 
  scale_x_discrete(limit=c("with", "without")) +
  xlab("With/without aquatic plants") +
  ylab("Hatching rate (%)") +
  theme_classic() +
  theme(text = element_text(size = 18, family = "Arial", colour = "black", face = "bold")) +
  theme(axis.ticks=element_line(colour = "black"),
        axis.text=element_text(colour = "black")) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 115))

ggsave("Fig1f_BiolLett.png", width = 4.8, height = 3.6)

```

## statistical analysis
```{R}
# make dataset with/without aquatic plants based on the experiment
with_plant <- c(4, 4, 4, 4, 4, 4, 3, 4)
without_plant <- c(0, 0, 0, 0, 0, 0, 0, 0)

# Wilcoxon singed-rank test 
# ref：https://zenn.dev/tmizuho/articles/e251a22d1a8e1e
coin::wilcoxsign_test(with_plant ~ without_plant, distribution = "exact", zero.method = "Wilcoxon")
```
