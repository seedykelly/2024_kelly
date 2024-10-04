rm(list=ls())
# if git is ahead by X commits do this: git reset --soft HEAD~8 (8=# of commits)

## ---- libraries ----
library(scales)
library(cowplot)
library(purrr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broman)
library(performance)
library(stringr)
library(bestNormalize)
library(brms)
library(rptR)
library(tidybayes)
library(MCMCglmm)
library(bayestestR)
library(modelr)
library(broom.mixed)
library(DHARMa)
library(bayesplot)
library(flextable)
library(modelsummary)
library(tinytable)
library(officer)
library(kableExtra)
library(knitr)
library(officedown)
library(smatr)
library(ggh4x)
library(scales)
## ---- end

## ---- setup ----
mobility_data<- read.csv(file="data/raw/giant_weta_behaviour.csv", header=TRUE, sep=",", dec=".")
morphological_data<- read.csv(file="data/raw/morphology_giant_weta.csv", header=TRUE, sep=",", dec=".")
condition_data<-read_csv(file="data/raw/giant_weta_composite_data.csv") #these are measures from general weta population, not individuals in the current study

# calculate condition (SMI)
mean_pronotum <- condition_data %>%
  group_by(sex) %>%
  summarise(mean_pro=mean(pronotum)) #population average pronotum length
#mean_pronotum[1,2]

allom_slope<-sma(mass~pronotum*sex, log="xy", method="SMA", data=condition_data)
summary(allom_slope)
coef(allom_slope)[2]

# # allometry of tibia
# allom_slope<-sma(Tibia~Pronotum*Sex, log="xy", method="SMA", data=morphology.Giant.weta)
# summary(allom_slope)

#calculate scaled mass index for each weta in this study and then add residual tibia length
morphology_long<-morphological_data %>%
  dplyr::select(ID, Pronotum, Tibia) %>%
  inner_join(mobility_data, by="ID") %>%
  mutate(smi = ifelse(sex == "f", mass*(12.6/Pronotum)^3.532828, mass*(10.9/Pronotum)^2.858071)) # SMI for each sex separately

# Subtracting 1 from observation number, so that intercept variance represents variation at the first trial (i.e. 0).
morphology_long$observation.n <- (morphology_long$observation - 1)

# Creating a mean-centered sex and time of day variable
morphology_long$sex.centred <- ifelse(morphology_long$sex == 'f', -0.5,
                                   ifelse(morphology_long$sex == 'm', 0.5, NA))

# Transform response variable
distanceBN<-bestNormalize(morphology_long$distance)
morphology_long$distance.tr <- distanceBN$x.t

#morphology_long$distance.t <- log(distance+1)

# Scale response variable
morphology_long$distance.z <- scale(as.numeric(morphology_long$distance.tr))
morphology_long$smi.z <- scale(as.numeric(morphology_long$smi))

# mean observations per sex
avg_obs <- mobility_data %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  group_by(sex, ID) %>%
  summarise(total_obs = n(), consecutive_obs=sum(!is.na(distance))) %>%
  summarise(mean_total=mean(total_obs), se.total=sd(total_obs)/sqrt(n()), n=n())

# Add mating success
morphology_long$ms <-ifelse(morphology_long$partner_id== "", 0,1)

mating_success_table <- morphology_long %>%
  group_by(ID, sex) %>%
  dplyr::summarise(mates=sum(ms), n=n(), ms=mates/n) #number of mates

# total number of observations by sex
obs_num <- mobility_data %>%
  group_by(sex) %>%
  dplyr::summarise(n=n())


## ---- end


#### PREDICTABILITY ####

# mobility_pred <- bf(distance.tr ~ sex.centred + observation.n + sex.centred:observation.n + (observation.n|ID), sigma ~ sex , family = gaussian)
mobility_pred <- bf(distance.tr ~ sex.centred + observation.n + sex.centred:observation.n + (observation.n|a|ID), sigma ~ sex + (1|a|ID), family = gaussian)

fit.model.brms.pred <- brm(mobility_pred, data = morphology_long, save_pars = save_pars(all = TRUE), 
                               warmup=6000, iter=10000, seed=12345, thin=4, chains=4, cores= 4, file = 'data/processed/fit.model.brms.pred')
summary(fit.model.brms.pred, prob = 0.95)

plot(fit.model.brms.pred, ask = F)
brms::pp_check(fit.model.brms.pred, resp = "distance.z", ndraws = 50)
conditional_effects(fit.model.brms.pred, prob = 0.95)

# get data in shape for plotting (see below)
mobility_female <-  exp(as_draws_df(fit.model.brms.pred)$"b_sigma_Intercept")^2
mobility_male <-  exp((as_draws_df(fit.model.brms.pred)$"b_sigma_Intercept") + (as_draws_df(fit.model.brms.pred)$"b_sigma_sexm"))^2

posterior_data_mobility = data.frame(mobility_female, mobility_male)

mobility.residual = posterior_data_mobility %>%
  dplyr::select(starts_with("mobility_")) %>%
  pivot_longer(cols = starts_with("mobility_"),
               names_to = 'Sex',
               names_prefix = "mobility_",
               values_to ="Estimate")


require(plyr)
mobility.mean <- ddply(mobility.residual, "Sex", summarise, grp.mean=mean(Estimate))
detach(package:plyr)

# males less predictable than females

## ---- data_analysis ----
fit.model.brms.pred = readRDS(file = "data/processed/fit.model.brms.pred.rds")

qmd.values <- fit.model.brms.pred %>%
  spread_draws(b_sigma_Intercept,b_sigma_sexm,sd_ID__sigma_Intercept, cor_ID__Intercept__sigma_Intercept) %>%
  median_qi(b_sigma_Intercept,b_sigma_sexm, sd_ID__sigma_Intercept, cor_ID__Intercept__sigma_Intercept) %>%
  select(-c(.width, .point,.interval)) %>%
  pivot_longer(cols=everything(), names_to = "parameter", values_to = "value")

get_variables(fit.model.brms.pred)

#### Behavioral predictability (individual rIIV's) ####

# Intercept Distance
bk.tr.dist <- exp(fixef(fit.model.brms.pred, pars = "sigma_Intercept")[1]) * sd(morphology_long$distance, na.rm = T)
# 7.66 m

#### Coefficient of variation in predictability (CVP) ####

log.norm.res.Dist <- exp(posterior_samples(fit.model.brms.pred)$"sd_ID__sigma_Intercept"^2)
CVP.long.Dist <- sqrt(log.norm.res.Dist - 1)
mean(CVP.long.Dist);HPDinterval(as.mcmc(CVP.long.Dist),0.95)

### Correlation between behavioural type and predictability ###

# COR.PERS.PRED <- 
#   as_draws_df(fit.model.brms.pred, 
#                     pars = c("cor_ID__Intercept__sigma_Intercept")) %>%
#   gather() %>%
#   separate(key,
#            c(NA,"Scale",NA,NA,NA,NA,NA,"Trait",NA),
#            sep = "([\\_\\__\\_\\__\\_\\_\\,])", fill = "right")

# Slope travel distance
cov.Trav <-
  posterior_samples(fit.model.brms.pred)[,11] *
  sqrt((posterior_samples(fit.model.brms.pred)[,7])^2) *
  sqrt((posterior_samples(fit.model.brms.pred)[,9])^2)
var.Trav <- (posterior_samples(fit.model.brms.pred)[,7])^2
slope_Trav <- cov.Trav / var.Trav
get_variables(fit.model.brms.pred)

mean_distance <- as_draws_df(fit.model.brms.pred, pars = "^r_ID")[1:47] %>%
  tidyr::gather(ID, value,
                "r_ID[W01,Intercept]" : 
                  "r_ID[W75,Intercept]") %>%
  dplyr::select(ID,value) %>%
  separate(ID,
           c(NA,NA,"ID",NA),
           sep = "([\\_\\[\\,])", fill = "right") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(Mean = mean(value),
         Up = mean(value) + 1.96 * sd(value),
         Lo = mean(value) - 1.96 * sd(value)) %>%
  dplyr::select (-value) %>%
  filter(!duplicated(ID))

names(mean_distance)<-c("ID","MeanDist","UpDist","LoDist")

sigma_distance <- as_draws_df(fit.model.brms.pred, pars = "^r_ID__sigma") %>%
  tidyr::gather(ID, value,
                "r_ID__sigma[W01,Intercept]" : 
                  "r_ID__sigma[W75,Intercept]")%>%
  dplyr::select(ID,value) %>%
  separate(ID,
           c(NA,"ID",NA),
           sep = "([\\[\\,])", fill = "right") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(Mean = mean(value),
         Up = mean(value) + 1.96 * sd(value),
         Lo = mean(value) - 1.96 * sd(value))%>%
  select (-value) %>%
  filter(!duplicated(ID))

names(sigma_distance)<-c("ID","predict_MeanDist","predict_UpDist","predict_LoDist")

### behavioural type vs mating success

mating_type_data <- right_join(mean_distance, mating_success_table, by ='ID') %>%
  filter(!is.na(MeanDist))

# x <-cor.test(mating_behaviour_data$Mean, mating_behaviour_data$ms)
# str(x)
# mating_behaviour_data %>%
#   group_by(sex) %>%
#   summarise(correlation = cor(Mean, ms))

mating_type_corr <- mating_type_data %>%
  group_by(sex) %>%
  dplyr::summarise(COR = stats::cor.test(MeanDist, ms)$estimate,
            pval = stats::cor.test(MeanDist, ms)$p.value,
            df = stats::cor.test(MeanDist, ms)$parameter
  ) %>%
  dplyr::mutate_at(vars(COR, pval, df), as.numeric)

### predictability vs mating success

mating_predictability_data <- right_join(sigma_distance, mating_success_table, by ='ID') %>%
  filter(!is.na(predict_MeanDist))

mating_predict_corr <- mating_predictability_data %>%
  group_by(sex) %>%
  dplyr::summarize(COR = stats::cor.test(predict_MeanDist, ms)$estimate,
            pval = stats::cor.test(predict_MeanDist, ms)$p.value,
            df = stats::cor.test(predict_MeanDist, ms)$parameter
  ) %>%
  dplyr::mutate_at(vars(COR, pval, df), as.numeric)

str(mating_predict_corr)

## ---- end ----

### FIGURES
sex.diff <- ggplot(mobility.residual, aes(x = Estimate, fill = Sex)) +
      geom_density(alpha=0.6)+ 
      geom_vline(data = mobility.mean, aes(xintercept=grp.mean, color = Sex),
             linetype="dashed", linewidth=1,show.legend=FALSE) +
      xlab('Residual within-individual variaiton in mobility\n (i.e. predictability)') +
      ylab('Density') +
      theme_bw() +
      theme_bw() +
      theme(legend.position = c(0.75,.75)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title.x = element_text(size=14)) +
      theme(axis.title.y = element_text(size=14)) +
      theme(axis.text.x = element_text(size=12)) +
      theme(axis.text.y = element_text(size=12))
ggsave(filename="figure_1.png", width=8, height=8, dpi=800,antialias="default")

correlation_data <- right_join(mean_distance, sigma_distance, by ='ID')

correlation_data_plot <- ggplot() +
  geom_segment(data = correlation_data[!duplicated(correlation_data$ID),],
               aes(x = LoDist,
                   xend = UpDist,
                   y = predict_MeanDist,
                   yend = predict_MeanDist),
               color = "#F8766D", alpha = 0.3) +

  geom_segment(data = correlation_data[!duplicated(correlation_data$ID),],
               aes(x = MeanDist,
                   xend = MeanDist,
                   y = predict_LoDist,
                   yend = predict_UpDist),
               color = "#F8766D", alpha = 0.3) +

  geom_point(data = correlation_data[!duplicated(correlation_data$ID),],
             aes(x = MeanDist, y = predict_MeanDist,
                 fill="#F8766D"),shape=22,
             size=3, color = "#F8766D", alpha = 0.9, stroke = 0) +

  geom_segment(aes(x = -1.5,
                   xend = 1.5,
                   y = 0 + mean(slope_Trav)*-1,
                   yend = 0 + mean(slope_Trav)*1),
               color="#F8766D",linewidth=1, alpha = 0.9) +

  ylab("Mobility predictability") +
  xlab("Mobility behavioral type") +
  scale_x_continuous(limits = c(-2.5,2.5), breaks = c(-2,0,2))+
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.75,0.75),
                     labels =c("-.5","0",".5"), expand = c(0,0)) +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  annotate("text",x = 1.3, y = 0.65, label = expression(paste(italic("r")*" = -0.10")), size = 4)
ggsave(filename="figure_2_supp.png", width=8, height=8, dpi=800,antialias="default")

# New facet label names for supp variable
sex.labs <- c("Females", "Males")
names(sex.labs) <- c("f", "m")

# distribution histograms
distributions <- ggplot(mobility_data, aes(x = distance)) +
  geom_histogram() + 
  facet_grid(~sex, labeller = labeller(sex = sex.labs)) +
  xlab('Distance travelled per night (m)') +
  ylab('Count') +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))
ggsave(filename="figure_1_supp.png", width=8, height=8, dpi=800,antialias="default")


##############
get_variables(fit.model.brms.pred)

# rename_explanator <- function(old_names) {
#   new_names <- gsub("b_|sd_ID__", "", old_names)
#   setNames(new_names, old_names)
# }
summary(fit.model.brms.pred)
coef(fit.model.brms.pred)
get_variables(fit.model.brms.pred)

