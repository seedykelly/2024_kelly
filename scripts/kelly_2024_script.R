#rm(list=ls())
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
library(forcats)
library(DHARMa)
library(bayesplot)
library(flextable)
library(modelsummary)
library(tinytable)
library(officer)
library(kableExtra)
library(knitr)
library(officedown)
library(gtsummary)
library(ggh4x)
## ---- end

check_brms <- function(model,             # brms model
                       integer = FALSE,   # integer response? (TRUE/FALSE)
                       plot = TRUE,       # make plot?
                       ...                # further arguments for DHARMa::plotResiduals 
) {
  
  mdata <- brms::standata(model)
  if (!"Y" %in% names(mdata))
    stop("Cannot extract the required information from this brms model")
  
  dharma.obj <- DHARMa::createDHARMa(
    simulatedResponse = t(brms::posterior_predict(model, ndraws = 1000)),
    observedResponse = mdata$Y, 
    fittedPredictedResponse = apply(
      t(brms::posterior_epred(model, ndraws = 1000, re.form = NA)),
      1,
      mean),
    integerResponse = integer)
  
  if (isTRUE(plot)) {
    plot(dharma.obj, ...)
  }
  
  invisible(dharma.obj)
  
}
## ---- data_analysis ----
behaviour_data<- read.csv(file="data/raw/bour_kelly_master_data.csv", header=TRUE, sep=",", dec=".")
beetle_data<- read.csv(file="data/raw/Beetle_data.csv", header=TRUE, sep=",", dec=".")

# rename column headings
behaviour_data_one <- behaviour_data %>%
  rename(distance=Distance.moved.Center.point.Total.cm, mobility = Mobility.Body.fill.Mean.., activity = Activity.Within.arena.Mean.., time_in_centre = centre, trial = Trial, arena=Arena, id=ID, sex=Sex, treatment=Treatment) %>%
  dplyr::select(-Subject.not.found, -Arena.settings) %>%
  filter(escaped==0) %>%
  dplyr::select(-escaped) %>%
  filter(id!="H6EOHO" & trial!=11)

behaviour_data_one[behaviour_data_one == "-"] <- NA

behaviour_data_one$trial = as.factor(behaviour_data_one$trial)
behaviour_data_one$sex = as.factor(behaviour_data_one$sex)
behaviour_data_one$arena = as.factor(behaviour_data_one$arena)
behaviour_data_one$treatment = as.factor(behaviour_data_one$treatment)
behaviour_data_one$id = as.factor(behaviour_data_one$id)
behaviour_data_one$latency_to_zone1 = as.numeric(behaviour_data_one$latency_to_zone1)
behaviour_data_one$latency_to_zone2 = as.numeric(behaviour_data_one$latency_to_zone2)
behaviour_data_one$latency_to_zone3 = as.numeric(behaviour_data_one$latency_to_zone3)
behaviour_data_one$latency_to_zone4 = as.numeric(behaviour_data_one$latency_to_zone4)

merged_data <-behaviour_data_one %>%
  arrange(id, trial) %>%
  group_by(id) %>%
  #mutate(assay = row_number()) %>%
  left_join(beetle_data %>% dplyr::select(id,age,mass_1,mass_2),
            by = "id") %>%
  ungroup() %>%
  mutate(treatment = recode_factor(treatment,
                                "High LPS" = "high",
                                "Low LPS" = "low",
                                "Saline" = "saline",
                                "Injured" = "injured",
                                "handled" = "handled"))

merged_data_2 <- merged_data %>%
  filter(!is.na(latency_to_zone1)) %>%
  filter(!is.na(latency_to_zone2)) %>%
  filter(!is.na(latency_to_zone3)) %>%
  filter(!is.na(latency_to_zone4)) %>%
  group_by(id, treatment) %>%
  summarise(count=n()) %>%
  filter(count==1) %>%
  print(n=1000) #36 ids with only one observation

merged_data_2 <- merged_data %>%
  group_by(treatment, id) %>%
  summarise(count=n()) %>%
  filter(count==1) %>%
  print(n=1000) #36 ids with only one observation

merged_data_3 <- merged_data %>%
  group_by(id) %>%
  filter(n()>1) %>%
  arrange(id, trial) %>%
  mutate(assay = row_number())

merged_data_3$age = as.factor(merged_data_3$age)

merged_data_four <- merged_data_3 %>%
  group_by(trial) %>%
  mutate(across(contains('to_zone'), replace_na, 600)) %>%
  mutate(explore_1 = pmax(latency_to_zone1, latency_to_zone2, latency_to_zone3, latency_to_zone4, na.rm = TRUE)) %>% 
  mutate(explore=600-explore_1) %>%
  dplyr::select(-latency_to_zone1, -latency_to_zone2, -latency_to_zone3, -latency_to_zone4) %>%
  mutate(treatment_pooled = factor(if_else(treatment == "saline" | treatment == "injured" | treatment == "handled" , "control", treatment))) %>%
  mutate(assay=factor(assay)) %>%
  mutate(assay_centred= if_else(assay == 1, 0, 1)) %>%
  mutate(assay_centred=factor(assay_centred))

# merged_data_five %>%
#   group_by(id) %>%
#   summarise(n=n()) %>%
#   filter(n==1) %>%
#   #count(treatment_two) %>%
#   print(n=1000) 

# sample sizes
# summary_1 <- merged_data_four %>%
#   group_by(id, sex,age,treatment) %>%
#   summarise(mean_activity=mean(distance), mean_explore=mean(explore), n=n()) %>%
#   group_by(sex,age,treatment) %>%
#   summarise(distance.x=mean(mean_activity),distance.sd= sd(mean_activity), explore.x=mean(mean_explore),explore.sd=sd(mean_activity), n=n()) %>%
#   arrange(desc(treatment)) %>%
#   pivot_wider(values_from = c(distance.x, distance.sd, explore.x, explore.sd, n), names_from = c("sex", "age")) %>%
#   relocate(treatment,distance.x_Female_Young,distance.sd_Female_Young,explore.x_Female_Young,explore.sd_Female_Young,n_Female_Young,
#            distance.x_Female_Old,distance.sd_Female_Old,explore.x_Female_Old,explore.sd_Female_Old,n_Female_Old,
#            distance.x_Male_Young,distance.sd_Male_Young,explore.x_Male_Young,explore.sd_Male_Young,n_Male_Young,
#            distance.x_Male_Old,distance.sd_Male_Old,explore.x_Male_Old,explore.sd_Male_Old,n_Male_Old) %>%
#   mutate(across(c(2:5,7:10, 12:15,17:20), ~ sprintf("%0.2f", .x))) %>%
#   unite(col='distance_female_young', c('distance.x_Female_Young', 'distance.sd_Female_Young'), sep=' ± ') %>%
#   unite(col='explore_female_young', c('explore.x_Female_Young', 'explore.sd_Female_Young'), sep=' ± ') %>%
#   unite(col='distance_female_old', c('distance.x_Female_Old', 'distance.sd_Female_Old'), sep=' ± ') %>%
#   unite(col='explore_female_old', c('explore.x_Female_Old', 'explore.sd_Female_Old'), sep=' ± ') %>%
#   unite(col='distance_male_young', c('distance.x_Male_Young', 'distance.sd_Male_Young'), sep=' ± ') %>%
#   unite(col='explore_male_young', c('explore.x_Male_Young', 'explore.sd_Male_Young'), sep=' ± ') %>%
#   unite(col='distance_male_old', c('distance.x_Male_Old', 'distance.sd_Male_Old'), sep=' ± ') %>%
#   unite(col='explore_male_old', c('explore.x_Male_Old', 'explore.sd_Male_Old'), sep=' ± ')

# sample sizes by assaytreatment# sample sizes by assay
# summary_1 <- merged_data_four %>%
#   group_by(id, age, treatment, sex, assay) %>%
#   summarise(mean_distance=mean(distance), mean_mobility=mean(mobility), mean_activity=mean(activity), mean_explore=mean(explore), mean_centre=mean(time_in_centre), n=n()) %>%
#   group_by(assay) %>%
#   summarise(distance=mean(mean_distance), mobility=mean(mean_mobility), activity=mean(mean_activity), explore=mean(mean_explore), boldness=mean(mean_centre),n=n())

# # distance data distributions
# ggplot(merged_data_four, aes(x = factor(treatment), y = log(distance+1), fill=factor(assay)) ) + 
#   ylab("Distance (cm)")+
#   xlab("Treatment") +
#   geom_boxplot() +
#   facet_wrap(vars(sex, age))
# 
# # exploration data distributions
# ggplot(merged_data_four, aes(x = factor(treatment), y = log(explore), fill=factor(assay)) ) + 
#   ylab("Time to visit all four zones (s)")+
#   xlab("Treatment") +
#   geom_boxplot() +
#   facet_wrap(vars(sex, age))

mass_data <- merged_data_four %>%
  group_by(id) %>%
  filter(row_number()==1)

mass <- aov(mass_1 ~ sex * age * treatment, data=mass_data)
summary.aov(mass) # no differences

mass_change <- glm(mass_2 ~ sex * age * treatment + mass_1, data=mass_data)
summary(mass_change)

## brms model results for rmarkdown. All analyses for these models are conducted below. 
fit_model.brms.activity.1 = readRDS(file = "data/processed/fit_model.brms.activity.1.rds")
fit_model.brms.activity.2 = readRDS(file = "data/processed/fit_model.brms.activity.2.rds")

fit_model.brms.explore.1 = readRDS(file = "data/processed/fit_model.brms.explore.1.rds")
fit_model.brms.explore.2 = readRDS(file = "data/processed/fit_model.brms.explore.2.rds")

activity_comp<-loo_compare(fit_model.brms.activity.1, fit_model.brms.activity.2, criterion = "loo") #  model 2 is best
explore_comp<-loo_compare(fit_model.brms.explore.1, fit_model.brms.explore.2,criterion = "loo") #  model 2 is best

load(file = "data/processed/contrasts.table.Va.RData")
load(file = "data/processed/contrasts.table.Vw.RData")
correlation.table <- readRDS(file = "data/processed/correlation.table.rda")
## ---- end

#### Bayesian analysis -- ACTIVITY ####

seed=123456
set.seed(seed)

distanceBN<-bestNormalize(merged_data_four$distance)
merged_data_four$distance.t <- distanceBN$x.t

time.lab <- c("Before transformation", "After transformation")
names(time.lab) <- c("before", "after")

plot.dist.distance <- merged_data_four %>%
  dplyr::select(distance,distance.t) %>%
  gather(treatment, value,distance:distance.t,factor_key=TRUE) %>%
  mutate(time=if_else(str_detect(treatment, 'distance.t'), 'after', 'before')) %>%
  mutate(time = factor(time, levels=c("before", "after"))) %>%
  ggplot(aes(x = value)) +
  geom_histogram(fill="lightblue") +
  facet_wrap(~ time, scales = "free",labeller = labeller(time=time.lab)) +
  ylab("Count") +
  xlab("Activity (Distance travelled cm)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

# model.brms=bf(distance.t ~ sex*treatment + (0+treatment||id) + (0+sex||id), sigma ~ 0+treatment, sigma ~ 0+sex) #Royaute-Dochterman
# model.brms=bf(distance.t ~ sex*treatment + (1|id))

model.brms.activity.1 <- bf(scale(distance.t)~ sex * treatment * age * assay + (0+sex:age||gr(id, by = treatment)), sigma ~ 0+age:sex:treatment, family = gaussian) 
fit_model.brms.activity.1 <- brm(model.brms.activity.1, data = merged_data_four, save_pars = save_pars(all = TRUE), warmup=500, iter=8000, seed=12345, thin=2, chains=4, cores= 4, file = 'data/processed/fit_model.brms.activity.1')
summary(fit_model.brms.activity.1)
fit_model.brms.activity.1 <- add_criterion(fit_model.brms.activity.1, "loo")
# # 
model.brms.activity.2 <- bf(scale(distance.t)~ sex+treatment+age + assay +(0+sex:age||gr(id, by = treatment)), sigma ~ 0+age:sex:treatment, family = gaussian) 
fit_model.brms.activity.2 <- brm(model.brms.activity.2, save_pars = save_pars(all = TRUE), data = merged_data_four, warmup=500, iter=8000, seed=12345, thin=2, chains=4, cores= 4, file = 'data/processed/fit_model.brms.activity.2')
summary(fit_model.brms.activity.2)
fit_model.brms.activity.2 <- add_criterion(fit_model.brms.activity.2, "loo")
# # 
fixef(fit_model.brms.activity.2)[8,1]

get_variables(fit_model.brms.activity.2)

# convergence and model check
# pp_check(fit_model.brms.activity.4)
# check_brms(fit_model.brms.activity.4)
# mcmc_trace(plot(fit_model.brms.activity.4))

# compare four models
activity_comp<-loo_compare(fit_model.brms.activity.1, fit_model.brms.activity.2, fit_model.brms.activity.3,fit_model.brms.activity.4,criterion = "loo") #  model 4 is best

# post <- 
#   fit_model.brms.activity.4 %>%
#   as_draws_df()
# 
# head(post)
# 
# post %>%
#   transmute(mu_High = b_Intercept,
#             mu_Low = b_Intercept + b_treatmentLowLPS,
#             mu_Saline = b_Intercept + b_treatmentSaline,
#             mu_Injured   = b_Intercept + b_treatmentInjured,
#             mu_handled = b_Intercept + b_treatmenthandled) %>%
#   gather() %>%
#   group_by(key) %>%
#   mean_hdi() %>% 
#   mutate_if(is.double, round, digits = 2)

#among
Va.activity.F.Y.handled <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexFemale:ageYoung:treatmenthandled"^2
Va.activity.F.Y.injured <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexFemale:ageYoung:treatmentinjured"^2
Va.activity.F.Y.saline <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexFemale:ageYoung:treatmentsaline"^2
Va.activity.F.Y.lo <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexFemale:ageYoung:treatmentlow"^2
Va.activity.F.Y.hi <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexFemale:ageYoung:treatmenthigh"^2

Va.activity.M.Y.handled <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexMale:ageYoung:treatmenthandled"^2
Va.activity.M.Y.injured <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexMale:ageYoung:treatmentinjured"^2
Va.activity.M.Y.saline <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexMale:ageYoung:treatmentsaline"^2
Va.activity.M.Y.lo <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexMale:ageYoung:treatmentlow"^2
Va.activity.M.Y.hi <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexMale:ageYoung:treatmenthigh"^2

Va.activity.F.O.handled <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexFemale:ageOld:treatmenthandled"^2
Va.activity.F.O.injured <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexFemale:ageOld:treatmentinjured"^2
Va.activity.F.O.saline <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexFemale:ageOld:treatmentsaline"^2
Va.activity.F.O.lo <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexFemale:ageOld:treatmentlow"^2
Va.activity.F.O.hi <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexFemale:ageOld:treatmenthigh"^2

Va.activity.M.O.handled <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexMale:ageOld:treatmenthandled"^2
Va.activity.M.O.injured <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexMale:ageOld:treatmentinjured"^2
Va.activity.M.O.saline <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexMale:ageOld:treatmentsaline"^2
Va.activity.M.O.lo <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexMale:ageOld:treatmentlow"^2
Va.activity.M.O.hi <- as_draws_df(fit_model.brms.activity.2)$"sd_id__sexMale:ageOld:treatmenthigh"^2

#within
Vw.activity.F.Y.handled <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageYoung:sexFemale:treatmenthandled")^2
Vw.activity.F.Y.injured <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageYoung:sexFemale:treatmentinjured")^2
Vw.activity.F.Y.saline <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageYoung:sexFemale:treatmentsaline")^2
Vw.activity.F.Y.lo <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageYoung:sexFemale:treatmentlow")^2
Vw.activity.F.Y.hi <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageYoung:sexFemale:treatmenthigh")^2

Vw.activity.M.Y.handled <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageYoung:sexMale:treatmenthandled")^2
Vw.activity.M.Y.injured <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageYoung:sexMale:treatmentinjured")^2
Vw.activity.M.Y.saline <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageYoung:sexMale:treatmentsaline")^2
Vw.activity.M.Y.lo <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageYoung:sexMale:treatmentlow")^2
Vw.activity.M.Y.hi <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageYoung:sexMale:treatmenthigh")^2

Vw.activity.F.O.handled <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageOld:sexFemale:treatmenthandled")^2
Vw.activity.F.O.injured <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageOld:sexFemale:treatmentinjured")^2
Vw.activity.F.O.saline <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageOld:sexFemale:treatmentsaline")^2
Vw.activity.F.O.lo <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageOld:sexFemale:treatmentlow")^2
Vw.activity.F.O.hi <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageOld:sexFemale:treatmenthigh")^2

Vw.activity.M.O.handled <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageOld:sexMale:treatmenthandled")^2
Vw.activity.M.O.injured <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageOld:sexMale:treatmentinjured")^2
Vw.activity.M.O.saline <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageOld:sexMale:treatmentsaline")^2
Vw.activity.M.O.lo <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageOld:sexMale:treatmentlow")^2
Vw.activity.M.O.hi <- exp(as_draws_df(fit_model.brms.activity.2)$"b_sigma_ageOld:sexMale:treatmenthigh")^2


post.data.all.activity = data.frame(Va.activity.F.Y.handled, Va.activity.F.Y.injured, Va.activity.F.Y.saline, Va.activity.F.Y.lo, Va.activity.F.Y.hi,
                                    Va.activity.M.Y.handled, Va.activity.M.Y.injured, Va.activity.M.Y.saline, Va.activity.M.Y.lo, Va.activity.M.Y.hi,
                                    Va.activity.F.O.handled, Va.activity.F.O.injured, Va.activity.F.O.saline, Va.activity.F.O.lo, Va.activity.F.O.hi, 
                                    Va.activity.M.O.handled, Va.activity.M.O.injured, Va.activity.M.O.saline, Va.activity.M.O.lo, Va.activity.M.O.hi, 
                                    Vw.activity.F.Y.handled, Vw.activity.F.Y.injured, Vw.activity.F.Y.saline, Vw.activity.F.Y.lo, Vw.activity.F.Y.hi, 
                                    Vw.activity.M.Y.handled, Vw.activity.M.Y.injured, Vw.activity.M.Y.saline, Vw.activity.M.Y.lo, Vw.activity.M.Y.hi, 
                                    Vw.activity.F.O.handled, Vw.activity.F.O.injured, Vw.activity.F.O.saline, Vw.activity.F.O.lo, Vw.activity.F.O.hi, 
                                    Vw.activity.M.Y.handled, Vw.activity.M.Y.injured, Vw.activity.M.Y.saline, Vw.activity.M.O.lo, Vw.activity.M.O.hi) 

all_activity.Va = post.data.all.activity %>%
  dplyr::select(starts_with("Va.activity.")) %>%
  pivot_longer(cols = starts_with("Va.activity."),
               names_to = 'treat',
               names_prefix = "Va.activity.",
               values_to ="Estimate")

Va_all_activity = all_activity.Va %>%
  dplyr::group_by(treat) %>%
  dplyr::summarise(Va_all_activity = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.89)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.89)[2], 3))%>%
  as.data.frame()

# Female-Young
post.data.all.activity$delta.va.activity.F.Y.m_i=
  with(post.data.all.activity, Va.activity.F.Y.handled-Va.activity.F.Y.injured)
post.data.all.activity$delta.va.activity.F.Y.m_s=
  with(post.data.all.activity, Va.activity.F.Y.handled-Va.activity.F.Y.saline)
post.data.all.activity$delta.va.activity.F.Y.m_l=
  with(post.data.all.activity, Va.activity.F.Y.handled-Va.activity.F.Y.lo)
post.data.all.activity$delta.va.activity.F.Y.m_h=
  with(post.data.all.activity, Va.activity.F.Y.handled-Va.activity.F.Y.hi)
post.data.all.activity$delta.va.activity.F.Y.i_s=
  with(post.data.all.activity, Va.activity.F.Y.injured-Va.activity.F.Y.saline)
post.data.all.activity$delta.va.activity.F.Y.i_l=
  with(post.data.all.activity, Va.activity.F.Y.injured-Va.activity.F.Y.lo)
post.data.all.activity$delta.va.activity.F.Y.i_h=
  with(post.data.all.activity, Va.activity.F.Y.injured-Va.activity.F.Y.hi)
post.data.all.activity$delta.va.activity.F.Y.s_l=
  with(post.data.all.activity, Va.activity.F.Y.saline-Va.activity.F.Y.lo)
post.data.all.activity$delta.va.activity.F.Y.s_h=
  with(post.data.all.activity, Va.activity.F.Y.saline-Va.activity.F.Y.hi)
post.data.all.activity$delta.va.activity.F.Y.l_h=
  with(post.data.all.activity, Va.activity.F.Y.lo-Va.activity.F.Y.hi)

# Male-Young
post.data.all.activity$delta.va.activity.M.Y.m_i=
  with(post.data.all.activity, Va.activity.M.Y.handled-Va.activity.M.Y.injured)
post.data.all.activity$delta.va.activity.M.Y.m_s=
  with(post.data.all.activity, Va.activity.M.Y.handled-Va.activity.M.Y.saline)
post.data.all.activity$delta.va.activity.M.Y.m_l=
  with(post.data.all.activity, Va.activity.M.Y.handled-Va.activity.M.Y.lo)
post.data.all.activity$delta.va.activity.M.Y.m_h=
  with(post.data.all.activity, Va.activity.M.Y.handled-Va.activity.M.Y.hi)
post.data.all.activity$delta.va.activity.M.Y.i_s=
  with(post.data.all.activity, Va.activity.M.Y.injured-Va.activity.M.Y.saline)
post.data.all.activity$delta.va.activity.M.Y.i_l=
  with(post.data.all.activity, Va.activity.M.Y.injured-Va.activity.M.Y.lo)
post.data.all.activity$delta.va.activity.M.Y.i_h=
  with(post.data.all.activity, Va.activity.M.Y.injured-Va.activity.M.Y.hi)
post.data.all.activity$delta.va.activity.M.Y.s_l=
  with(post.data.all.activity, Va.activity.M.Y.saline-Va.activity.M.Y.lo)
post.data.all.activity$delta.va.activity.M.Y.s_h=
  with(post.data.all.activity, Va.activity.M.Y.saline-Va.activity.M.Y.hi)
post.data.all.activity$delta.va.activity.M.Y.l_h=
  with(post.data.all.activity, Va.activity.M.Y.lo-Va.activity.M.Y.hi)

# Female-Old
post.data.all.activity$delta.va.activity.F.O.m_i=
  with(post.data.all.activity, Va.activity.F.O.handled-Va.activity.F.O.injured)
post.data.all.activity$delta.va.activity.F.O.m_s=
  with(post.data.all.activity, Va.activity.F.O.handled-Va.activity.F.O.saline)
post.data.all.activity$delta.va.activity.F.O.m_l=
  with(post.data.all.activity, Va.activity.F.O.handled-Va.activity.F.O.lo)
post.data.all.activity$delta.va.activity.F.O.m_h=
  with(post.data.all.activity, Va.activity.F.O.handled-Va.activity.F.O.hi)
post.data.all.activity$delta.va.activity.F.O.i_s=
  with(post.data.all.activity, Va.activity.F.O.injured-Va.activity.F.O.saline)
post.data.all.activity$delta.va.activity.F.O.i_l=
  with(post.data.all.activity, Va.activity.F.O.injured-Va.activity.F.O.lo)
post.data.all.activity$delta.va.activity.F.O.i_h=
  with(post.data.all.activity, Va.activity.F.O.injured-Va.activity.F.O.hi)
post.data.all.activity$delta.va.activity.F.O.s_l=
  with(post.data.all.activity, Va.activity.F.O.saline-Va.activity.F.O.lo)
post.data.all.activity$delta.va.activity.F.O.s_h=
  with(post.data.all.activity, Va.activity.F.O.saline-Va.activity.F.O.hi)
post.data.all.activity$delta.va.activity.F.O.l_h=
  with(post.data.all.activity, Va.activity.F.O.lo-Va.activity.F.O.hi)

# Male-Old
post.data.all.activity$delta.va.activity.M.O.m_i=
  with(post.data.all.activity, Va.activity.M.O.handled-Va.activity.M.O.injured)
post.data.all.activity$delta.va.activity.M.O.m_s=
  with(post.data.all.activity, Va.activity.M.O.handled-Va.activity.M.O.saline)
post.data.all.activity$delta.va.activity.M.O.m_l=
  with(post.data.all.activity, Va.activity.M.O.handled-Va.activity.M.O.lo)
post.data.all.activity$delta.va.activity.M.O.m_h=
  with(post.data.all.activity, Va.activity.M.O.handled-Va.activity.M.O.hi)
post.data.all.activity$delta.va.activity.M.O.i_s=
  with(post.data.all.activity, Va.activity.M.O.injured-Va.activity.M.O.saline)
post.data.all.activity$delta.va.activity.M.O.i_l=
  with(post.data.all.activity, Va.activity.M.O.injured-Va.activity.M.O.lo)
post.data.all.activity$delta.va.activity.M.O.i_h=
  with(post.data.all.activity, Va.activity.M.O.injured-Va.activity.M.O.hi)
post.data.all.activity$delta.va.activity.M.O.s_l=
  with(post.data.all.activity, Va.activity.M.O.saline-Va.activity.M.O.lo)
post.data.all.activity$delta.va.activity.M.O.s_h=
  with(post.data.all.activity, Va.activity.M.O.saline-Va.activity.M.O.hi)
post.data.all.activity$delta.va.activity.M.O.l_h=
  with(post.data.all.activity, Va.activity.M.O.lo-Va.activity.M.O.hi)

va.delta.all.activity = post.data.all.activity %>%
  dplyr::select(starts_with("delta.va.activity.")) %>%
  pivot_longer(cols = starts_with("delta.va.activity."),
               names_to = 'Contrast',
               names_prefix = "delta.va.activity.",
               values_to ="Estimate")

all.activity.va.delta = va.delta.all.activity %>%
  dplyr::group_by(Contrast)%>%
  dplyr::summarise(va.delta.all.activity = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.89)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.89)[2], 3)) %>%
  as.data.frame()

# within-individual
all_activity.Vw = post.data.all.activity %>%
  dplyr::select(starts_with("Vw.activity.")) %>%
  pivot_longer(cols = starts_with("Vw.activity."),
               names_to = 'treat',
               names_prefix = "Vw.activity.",
               values_to ="Estimate")

Vw_all_activity = all_activity.Vw %>%
  dplyr::group_by(treat)%>%
  dplyr::summarise(Vw_all_activity = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[2], 3))%>%
  as.data.frame()

# Female-Young
post.data.all.activity$delta.vw.activity.F.Y.m_i=
  with(post.data.all.activity, Vw.activity.F.Y.handled-Vw.activity.F.Y.injured)
post.data.all.activity$delta.vw.activity.F.Y.m_s=
  with(post.data.all.activity, Vw.activity.F.Y.handled-Vw.activity.F.Y.saline)
post.data.all.activity$delta.vw.activity.F.Y.m_l=
  with(post.data.all.activity, Vw.activity.F.Y.handled-Vw.activity.F.Y.lo)
post.data.all.activity$delta.vw.activity.F.Y.m_h=
  with(post.data.all.activity, Vw.activity.F.Y.handled-Vw.activity.F.Y.hi)
post.data.all.activity$delta.vw.activity.F.Y.i_s=
  with(post.data.all.activity, Vw.activity.F.Y.injured-Vw.activity.F.Y.saline)
post.data.all.activity$delta.vw.activity.F.Y.i_l=
  with(post.data.all.activity, Vw.activity.F.Y.injured-Vw.activity.F.Y.lo)
post.data.all.activity$delta.vw.activity.F.Y.i_h=
  with(post.data.all.activity, Vw.activity.F.Y.injured-Vw.activity.F.Y.hi)
post.data.all.activity$delta.vw.activity.F.Y.s_l=
  with(post.data.all.activity, Vw.activity.F.Y.saline-Vw.activity.F.Y.lo)
post.data.all.activity$delta.vw.activity.F.Y.s_h=
  with(post.data.all.activity, Vw.activity.F.Y.saline-Vw.activity.F.Y.hi)
post.data.all.activity$delta.vw.activity.F.Y.l_h=
  with(post.data.all.activity, Vw.activity.F.Y.lo-Vw.activity.F.Y.hi)

# Male-Young
post.data.all.activity$delta.vw.activity.M.Y.m_i=
  with(post.data.all.activity, Vw.activity.M.Y.handled-Vw.activity.M.Y.injured)
post.data.all.activity$delta.vw.activity.M.Y.m_s=
  with(post.data.all.activity, Vw.activity.M.Y.handled-Vw.activity.M.Y.saline)
post.data.all.activity$delta.vw.activity.M.Y.m_l=
  with(post.data.all.activity, Vw.activity.M.Y.handled-Vw.activity.M.Y.lo)
post.data.all.activity$delta.vw.activity.M.Y.m_h=
  with(post.data.all.activity, Vw.activity.M.Y.handled-Vw.activity.M.Y.hi)
post.data.all.activity$delta.vw.activity.M.Y.i_s=
  with(post.data.all.activity, Vw.activity.M.Y.injured-Vw.activity.M.Y.saline)
post.data.all.activity$delta.vw.activity.M.Y.i_l=
  with(post.data.all.activity, Vw.activity.M.Y.injured-Vw.activity.M.Y.lo)
post.data.all.activity$delta.vw.activity.M.Y.i_h=
  with(post.data.all.activity, Vw.activity.M.Y.injured-Vw.activity.M.Y.hi)
post.data.all.activity$delta.vw.activity.M.Y.s_l=
  with(post.data.all.activity, Vw.activity.M.Y.saline-Vw.activity.M.Y.lo)
post.data.all.activity$delta.vw.activity.M.Y.s_h=
  with(post.data.all.activity, Vw.activity.M.Y.saline-Vw.activity.M.Y.hi)
post.data.all.activity$delta.vw.activity.M.Y.l_h=
  with(post.data.all.activity, Vw.activity.M.Y.lo-Vw.activity.M.Y.hi)

# Female-Old
post.data.all.activity$delta.vw.activity.F.O.m_i=
  with(post.data.all.activity, Vw.activity.F.O.handled-Vw.activity.F.O.injured)
post.data.all.activity$delta.vw.activity.F.O.m_s=
  with(post.data.all.activity, Vw.activity.F.O.handled-Vw.activity.F.O.saline)
post.data.all.activity$delta.vw.activity.F.O.m_l=
  with(post.data.all.activity, Vw.activity.F.O.handled-Vw.activity.F.O.lo)
post.data.all.activity$delta.vw.activity.F.O.m_h=
  with(post.data.all.activity, Vw.activity.F.O.handled-Vw.activity.F.O.hi)
post.data.all.activity$delta.vw.activity.F.O.i_s=
  with(post.data.all.activity, Vw.activity.F.O.injured-Vw.activity.F.O.saline)
post.data.all.activity$delta.vw.activity.F.O.i_l=
  with(post.data.all.activity, Vw.activity.F.O.injured-Vw.activity.F.O.lo)
post.data.all.activity$delta.vw.activity.F.O.i_h=
  with(post.data.all.activity, Vw.activity.F.O.injured-Vw.activity.F.O.hi)
post.data.all.activity$delta.vw.activity.F.O.s_l=
  with(post.data.all.activity, Vw.activity.F.O.saline-Vw.activity.F.O.lo)
post.data.all.activity$delta.vw.activity.F.O.s_h=
  with(post.data.all.activity, Vw.activity.F.O.saline-Vw.activity.F.O.hi)
post.data.all.activity$delta.vw.activity.F.O.l_h=
  with(post.data.all.activity, Vw.activity.F.O.lo-Vw.activity.F.O.hi)

# Male-Old
post.data.all.activity$delta.vw.activity.M.O.m_i=
  with(post.data.all.activity, Vw.activity.M.O.handled-Vw.activity.M.O.injured)
post.data.all.activity$delta.vw.activity.M.O.m_s=
  with(post.data.all.activity, Vw.activity.M.O.handled-Vw.activity.M.O.saline)
post.data.all.activity$delta.vw.activity.M.O.m_l=
  with(post.data.all.activity, Vw.activity.M.O.handled-Vw.activity.M.O.lo)
post.data.all.activity$delta.vw.activity.M.O.m_h=
  with(post.data.all.activity, Vw.activity.M.O.handled-Vw.activity.M.O.hi)
post.data.all.activity$delta.vw.activity.M.O.i_s=
  with(post.data.all.activity, Vw.activity.M.O.injured-Vw.activity.M.O.saline)
post.data.all.activity$delta.vw.activity.M.O.i_l=
  with(post.data.all.activity, Vw.activity.M.O.injured-Vw.activity.M.O.lo)
post.data.all.activity$delta.vw.activity.M.O.i_h=
  with(post.data.all.activity, Vw.activity.M.O.injured-Vw.activity.M.O.hi)
post.data.all.activity$delta.vw.activity.M.O.s_l=
  with(post.data.all.activity, Vw.activity.M.O.saline-Vw.activity.M.O.lo)
post.data.all.activity$delta.vw.activity.M.O.s_h=
  with(post.data.all.activity, Vw.activity.M.O.saline-Vw.activity.M.O.hi)
post.data.all.activity$delta.vw.activity.M.O.l_h=
  with(post.data.all.activity, Vw.activity.M.O.lo-Vw.activity.M.O.hi)

vw.delta.all.activity = post.data.all.activity %>%
  dplyr::select(starts_with("delta.vw.activity.")) %>%
  pivot_longer(cols = starts_with("delta.vw.activity."),
               names_to = 'Contrast',
               names_prefix = "delta.vw.activity.",
               values_to ="Estimate")


all.activity.vw.delta = vw.delta.all.activity %>%
  dplyr::group_by(Contrast)%>%
  dplyr::summarise(vw.delta.all.activity = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[2], 3)) %>%
  as.data.frame()

# #tau
# tau.hi=Vi.hi/(Vi.hi+Vw.hi)
# tau.lo=Vi.lo/(Vi.lo+Vw.lo)
# 
# #delta
# var.brms$delta.Vi.F.Y.c_l=var.brms$Vi.F.Y.control-var.brms$Vi.F.Y.lo
# var.brms$delta.Vi.F.Y.c_h=var.brms$Vi.F.Y.control-var.brms$Vi.F.Y.hi
# var.brms$delta.Vi.F.Y.l_h=var.brms$Vi.F.Y.lo-var.brms$Vi.F.Y.hi
# 
# var.brms$delta.Vi.M.Y.c_l=var.brms$Vi.M.Y.control-var.brms$Vi.M.Y.lo
# var.brms$delta.Vi.M.Y.c_h=var.brms$Vi.M.Y.control-var.brms$Vi.M.Y.hi
# var.brms$delta.Vi.M.Y.l_h=var.brms$Vi.M.Y.lo-var.brms$Vi.M.Y.hi
# 
# var.brms$delta.Vi.F.O.c_l=var.brms$Vi.F.O.control-var.brms$Vi.F.O.lo
# var.brms$delta.Vi.F.O.c_h=var.brms$Vi.F.O.control-var.brms$Vi.F.O.hi
# var.brms$delta.Vi.F.O.l_h=var.brms$Vi.F.O.lo-var.brms$Vi.F.O.hi
# 
# var.brms$delta.Vi.M.O.c_l=var.brms$Vi.M.O.control-var.brms$Vi.M.O.lo
# var.brms$delta.Vi.M.O.c_h=var.brms$Vi.M.O.control-var.brms$Vi.M.O.hi
# var.brms$delta.Vi.M.O.l_h=var.brms$Vi.M.O.lo-var.brms$Vi.M.O.hi
# 
# var.brms$delta.Vw=var.brms$Vw.lo-var.brms$Vw.hi
# var.brms$delta.tau=var.brms$tau.lo-var.brms$tau.hi
# 
# t.2=var.brms %>% select(delta.Vi.F.Y.c_l,delta.Vi.F.Y.c_h,delta.Vi.F.Y.l_h) %>% stack() %>% 
#   group_by(ind) %>% summarise_if(is.numeric, describe_posterior)
# 
# t.2=t.2$values 
# t.2$level=c("Vi","Vw","tau")
# t.2=t.2 %>% select(level, Median, CI_low, CI_high, pd)
# colnames(t.2)=c(expression(Delta),"Median","lower_CI", "Upper_CI", "Pmcmc")
# 
# t.2 %>% knitr::kable(digits = 2)
# 
# 
# var.res <- exp(posterior_samples(fit_model.brms)$"sd_id__sexMale:ageYoung:treatmentLowLPS")^2 
# mean(var.brms$tau.hi);HPDinterval(as.mcmc(var.brms$tau.hi),0.95)
# 

#### Bayesian analysis -- EXPLORATION ####

exploreBN<-bestNormalize(merged_data_four$explore)
merged_data_four$explore.t <- exploreBN$x.t

plot.dist.explore <- merged_data_four %>%
  dplyr::select(explore,explore.t) %>%
  gather(treatment, value,explore:explore.t,factor_key=TRUE) %>%
  mutate(time=if_else(str_detect(treatment, 'explore.t'), 'after', 'before')) %>%
  mutate(time = factor(time, levels=c("before", "after"))) %>%
  ggplot(aes(x = value)) +
  geom_histogram(fill="orange") +
  facet_wrap(~ time, scales = "free",labeller = labeller(time=time.lab)) +
  ylab("Count") +
  xlab("Exploration (Time to visit all zones, s)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

model.brms.explore.1 <- bf(scale(explore.t)~ sex*treatment*age + assay +(0+sex:age||gr(id, by = treatment)), sigma ~ 0+age:sex:treatment, family = gaussian) 
fit_model.brms.explore.1 <- brm(model.brms.explore.1, data = merged_data_four, iter=4000, thin=2,seed=seed, cores= 4,file = 'data/processed/fit_model.brms.explore.1')
summary(fit_model.brms.explore.1)
fit_model.brms.explore.1 <- add_criterion(fit_model.brms.explore.1, "loo", moment_matching=TRUE)
# fit_model.brms.explore.1 = readRDS(file = "data/processed/fit_model.brms.explore.1.rds")

model.brms.explore.2 <- bf(scale(explore.t)~ sex+treatment+age + assay +(0+sex:age||gr(id, by = treatment)), sigma ~ 0+age:sex:treatment, family = gaussian) 
fit_model.brms.explore.2 <- brm(model.brms.explore.2, data = merged_data_four, thin=2,iter=4000, seed=seed, cores= 4, file = 'data/processed/fit_model.brms.explore.2')
summary(fit_model.brms.explore.2)
fit_model.brms.explore.2 <- add_criterion(fit_model.brms.explore.2, "loo", moment_matching=TRUE)

# fit_model.brms.explore.4 = readRDS(file = "data/processed/fit_model.brms.explore.4.rds")

#head(get_variables(fit_model.brms.explore.4),80)

# convergence and model check
# pp_check(fit_model.brms.explore.4)
# check_brms(fit_model.brms.explore.4)
# mcmc_trace(plot(fit_model.brms.explore.4))

# compare four models
explore_comp<-loo_compare(fit_model.brms.explore.1, fit_model.brms.explore.2,criterion = "loo") #  model 4 is best

# post <- 
#   fit_model.brms %>%
#   as_draws_df()
# 
# head(post)
# 
# post %>%
#   transmute(mu_High = b_Intercept,
#             mu_Low = b_Intercept + b_treatmentLowLPS,
#             mu_Saline = b_Intercept + b_treatmentSaline,
#             mu_Injured   = b_Intercept + b_treatmentInjured,
#             mu_handled = b_Intercept + b_treatmenthandled) %>%
#   gather() %>%
#   group_by(key) %>%
#   mean_hdi() %>% 
#   mutate_if(is.double, round, digits = 2)

#among
Va.explore.F.Y.handled <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexFemale:ageYoung:treatmenthandled"^2
Va.explore.F.Y.injured <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexFemale:ageYoung:treatmentinjured"^2
Va.explore.F.Y.saline <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexFemale:ageYoung:treatmentsaline"^2
Va.explore.F.Y.lo <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexFemale:ageYoung:treatmentlow"^2
Va.explore.F.Y.hi <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexFemale:ageYoung:treatmenthigh"^2

Va.explore.M.Y.handled <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexMale:ageYoung:treatmenthandled"^2
Va.explore.M.Y.injured <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexMale:ageYoung:treatmentinjured"^2
Va.explore.M.Y.saline <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexMale:ageYoung:treatmentsaline"^2
Va.explore.M.Y.lo <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexMale:ageYoung:treatmentlow"^2
Va.explore.M.Y.hi <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexMale:ageYoung:treatmenthigh"^2

Va.explore.F.O.handled <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexFemale:ageOld:treatmenthandled"^2
Va.explore.F.O.injured <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexFemale:ageOld:treatmentinjured"^2
Va.explore.F.O.saline <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexFemale:ageOld:treatmentsaline"^2
Va.explore.F.O.lo <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexFemale:ageOld:treatmentlow"^2
Va.explore.F.O.hi <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexFemale:ageOld:treatmenthigh"^2

Va.explore.M.O.handled <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexMale:ageOld:treatmenthandled"^2
Va.explore.M.O.injured <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexMale:ageOld:treatmentinjured"^2
Va.explore.M.O.saline <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexMale:ageOld:treatmentsaline"^2
Va.explore.M.O.lo <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexMale:ageOld:treatmentlow"^2
Va.explore.M.O.hi <- as_draws_df(fit_model.brms.explore.2)$"sd_id__sexMale:ageOld:treatmenthigh"^2

#within
Vw.explore.F.Y.handled <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageYoung:sexFemale:treatmenthandled")^2
Vw.explore.F.Y.injured <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageYoung:sexFemale:treatmentinjured")^2
Vw.explore.F.Y.saline <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageYoung:sexFemale:treatmentsaline")^2
Vw.explore.F.Y.lo <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageYoung:sexFemale:treatmentlow")^2
Vw.explore.F.Y.hi <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageYoung:sexFemale:treatmenthigh")^2

Vw.explore.M.Y.handled <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageYoung:sexMale:treatmenthandled")^2
Vw.explore.M.Y.injured <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageYoung:sexMale:treatmentinjured")^2
Vw.explore.M.Y.saline <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageYoung:sexMale:treatmentsaline")^2
Vw.explore.M.Y.lo <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageYoung:sexMale:treatmentlow")^2
Vw.explore.M.Y.hi <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageYoung:sexMale:treatmenthigh")^2

Vw.explore.F.O.handled <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageOld:sexFemale:treatmenthandled")^2
Vw.explore.F.O.injured <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageOld:sexFemale:treatmentinjured")^2
Vw.explore.F.O.saline <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageOld:sexFemale:treatmentsaline")^2
Vw.explore.F.O.lo <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageOld:sexFemale:treatmentlow")^2
Vw.explore.F.O.hi <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageOld:sexFemale:treatmenthigh")^2

Vw.explore.M.O.handled <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageOld:sexMale:treatmenthandled")^2
Vw.explore.M.O.injured <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageOld:sexMale:treatmentinjured")^2
Vw.explore.M.O.saline <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageOld:sexMale:treatmentsaline")^2
Vw.explore.M.O.lo <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageOld:sexMale:treatmentlow")^2
Vw.explore.M.O.hi <- exp(as_draws_df(fit_model.brms.explore.2)$"b_sigma_ageOld:sexMale:treatmenthigh")^2


post.data.all.explore = data.frame(Va.explore.F.Y.handled, Va.explore.F.Y.injured, Va.explore.F.Y.saline, Va.explore.F.Y.lo, Va.explore.F.Y.hi,
                                    Va.explore.M.Y.handled, Va.explore.M.Y.injured, Va.explore.M.Y.saline, Va.explore.M.Y.lo, Va.explore.M.Y.hi,
                                    Va.explore.F.O.handled, Va.explore.F.O.injured, Va.explore.F.O.saline, Va.explore.F.O.lo, Va.explore.F.O.hi, 
                                    Va.explore.M.O.handled, Va.explore.M.O.injured, Va.explore.M.O.saline, Va.explore.M.O.lo, Va.explore.M.O.hi, 
                                    Vw.explore.F.Y.handled, Vw.explore.F.Y.injured, Vw.explore.F.Y.saline, Vw.explore.F.Y.lo, Vw.explore.F.Y.hi, 
                                    Vw.explore.M.Y.handled, Vw.explore.M.Y.injured, Vw.explore.M.Y.saline, Vw.explore.M.Y.lo, Vw.explore.M.Y.hi, 
                                    Vw.explore.F.O.handled, Vw.explore.F.O.injured, Vw.explore.F.O.saline, Vw.explore.F.O.lo, Vw.explore.F.O.hi, 
                                    Vw.explore.M.Y.handled, Vw.explore.M.Y.injured, Vw.explore.M.Y.saline, Vw.explore.M.O.lo, Vw.explore.M.O.hi) 

all_explore.Va = post.data.all.explore %>%
  dplyr::select(starts_with("Va.explore.")) %>%
  pivot_longer(cols = starts_with("Va.explore."),
               names_to = 'treat',
               names_prefix = "Va.explore.",
               values_to ="Estimate")

Va_all_explore = all_explore.Va %>%
  dplyr::group_by(treat) %>%
  dplyr::summarise(Va_all_explore = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.89)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.89)[2], 3))%>%
  as.data.frame()

# Female-Young
post.data.all.explore$delta.va.explore.F.Y.m_i=
  with(post.data.all.explore, Va.explore.F.Y.handled-Va.explore.F.Y.injured)
post.data.all.explore$delta.va.explore.F.Y.m_s=
  with(post.data.all.explore, Va.explore.F.Y.handled-Va.explore.F.Y.saline)
post.data.all.explore$delta.va.explore.F.Y.m_l=
  with(post.data.all.explore, Va.explore.F.Y.handled-Va.explore.F.Y.lo)
post.data.all.explore$delta.va.explore.F.Y.m_h=
  with(post.data.all.explore, Va.explore.F.Y.handled-Va.explore.F.Y.hi)
post.data.all.explore$delta.va.explore.F.Y.i_s=
  with(post.data.all.explore, Va.explore.F.Y.injured-Va.explore.F.Y.saline)
post.data.all.explore$delta.va.explore.F.Y.i_l=
  with(post.data.all.explore, Va.explore.F.Y.injured-Va.explore.F.Y.lo)
post.data.all.explore$delta.va.explore.F.Y.i_h=
  with(post.data.all.explore, Va.explore.F.Y.injured-Va.explore.F.Y.hi)
post.data.all.explore$delta.va.explore.F.Y.s_l=
  with(post.data.all.explore, Va.explore.F.Y.saline-Va.explore.F.Y.lo)
post.data.all.explore$delta.va.explore.F.Y.s_h=
  with(post.data.all.explore, Va.explore.F.Y.saline-Va.explore.F.Y.hi)
post.data.all.explore$delta.va.explore.F.Y.l_h=
  with(post.data.all.explore, Va.explore.F.Y.lo-Va.explore.F.Y.hi)

# Male-Young
post.data.all.explore$delta.va.explore.M.Y.m_i=
  with(post.data.all.explore, Va.explore.M.Y.handled-Va.explore.M.Y.injured)
post.data.all.explore$delta.va.explore.M.Y.m_s=
  with(post.data.all.explore, Va.explore.M.Y.handled-Va.explore.M.Y.saline)
post.data.all.explore$delta.va.explore.M.Y.m_l=
  with(post.data.all.explore, Va.explore.M.Y.handled-Va.explore.M.Y.lo)
post.data.all.explore$delta.va.explore.M.Y.m_h=
  with(post.data.all.explore, Va.explore.M.Y.handled-Va.explore.M.Y.hi)
post.data.all.explore$delta.va.explore.M.Y.i_s=
  with(post.data.all.explore, Va.explore.M.Y.injured-Va.explore.M.Y.saline)
post.data.all.explore$delta.va.explore.M.Y.i_l=
  with(post.data.all.explore, Va.explore.M.Y.injured-Va.explore.M.Y.lo)
post.data.all.explore$delta.va.explore.M.Y.i_h=
  with(post.data.all.explore, Va.explore.M.Y.injured-Va.explore.M.Y.hi)
post.data.all.explore$delta.va.explore.M.Y.s_l=
  with(post.data.all.explore, Va.explore.M.Y.saline-Va.explore.M.Y.lo)
post.data.all.explore$delta.va.explore.M.Y.s_h=
  with(post.data.all.explore, Va.explore.M.Y.saline-Va.explore.M.Y.hi)
post.data.all.explore$delta.va.explore.M.Y.l_h=
  with(post.data.all.explore, Va.explore.M.Y.lo-Va.explore.M.Y.hi)

# Female-Old
post.data.all.explore$delta.va.explore.F.O.m_i=
  with(post.data.all.explore, Va.explore.F.O.handled-Va.explore.F.O.injured)
post.data.all.explore$delta.va.explore.F.O.m_s=
  with(post.data.all.explore, Va.explore.F.O.handled-Va.explore.F.O.saline)
post.data.all.explore$delta.va.explore.F.O.m_l=
  with(post.data.all.explore, Va.explore.F.O.handled-Va.explore.F.O.lo)
post.data.all.explore$delta.va.explore.F.O.m_h=
  with(post.data.all.explore, Va.explore.F.O.handled-Va.explore.F.O.hi)
post.data.all.explore$delta.va.explore.F.O.i_s=
  with(post.data.all.explore, Va.explore.F.O.injured-Va.explore.F.O.saline)
post.data.all.explore$delta.va.explore.F.O.i_l=
  with(post.data.all.explore, Va.explore.F.O.injured-Va.explore.F.O.lo)
post.data.all.explore$delta.va.explore.F.O.i_h=
  with(post.data.all.explore, Va.explore.F.O.injured-Va.explore.F.O.hi)
post.data.all.explore$delta.va.explore.F.O.s_l=
  with(post.data.all.explore, Va.explore.F.O.saline-Va.explore.F.O.lo)
post.data.all.explore$delta.va.explore.F.O.s_h=
  with(post.data.all.explore, Va.explore.F.O.saline-Va.explore.F.O.hi)
post.data.all.explore$delta.va.explore.F.O.l_h=
  with(post.data.all.explore, Va.explore.F.O.lo-Va.explore.F.O.hi)

# Male-Old
post.data.all.explore$delta.va.explore.M.O.m_i=
  with(post.data.all.explore, Va.explore.M.O.handled-Va.explore.M.O.injured)
post.data.all.explore$delta.va.explore.M.O.m_s=
  with(post.data.all.explore, Va.explore.M.O.handled-Va.explore.M.O.saline)
post.data.all.explore$delta.va.explore.M.O.m_l=
  with(post.data.all.explore, Va.explore.M.O.handled-Va.explore.M.O.lo)
post.data.all.explore$delta.va.explore.M.O.m_h=
  with(post.data.all.explore, Va.explore.M.O.handled-Va.explore.M.O.hi)
post.data.all.explore$delta.va.explore.M.O.i_s=
  with(post.data.all.explore, Va.explore.M.O.injured-Va.explore.M.O.saline)
post.data.all.explore$delta.va.explore.M.O.i_l=
  with(post.data.all.explore, Va.explore.M.O.injured-Va.explore.M.O.lo)
post.data.all.explore$delta.va.explore.M.O.i_h=
  with(post.data.all.explore, Va.explore.M.O.injured-Va.explore.M.O.hi)
post.data.all.explore$delta.va.explore.M.O.s_l=
  with(post.data.all.explore, Va.explore.M.O.saline-Va.explore.M.O.lo)
post.data.all.explore$delta.va.explore.M.O.s_h=
  with(post.data.all.explore, Va.explore.M.O.saline-Va.explore.M.O.hi)
post.data.all.explore$delta.va.explore.M.O.l_h=
  with(post.data.all.explore, Va.explore.M.O.lo-Va.explore.M.O.hi)

va.delta.all.explore = post.data.all.explore %>%
  dplyr::select(starts_with("delta.va.explore.")) %>%
  pivot_longer(cols = starts_with("delta.va.explore."),
               names_to = 'Contrast',
               names_prefix = "delta.va.explore.",
               values_to ="Estimate")

all.explore.va.delta = va.delta.all.explore %>%
  dplyr::group_by(Contrast)%>%
  dplyr::summarise(va.delta.all.explore = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.89)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.89)[2], 3)) %>%
  as.data.frame()

# within-individual
all_explore.Vw = post.data.all.explore %>%
  dplyr::select(starts_with("Vw.explore.")) %>%
  pivot_longer(cols = starts_with("Vw.explore."),
               names_to = 'treat',
               names_prefix = "Vw.explore.",
               values_to ="Estimate")

Vw_all_explore = all_explore.Vw %>%
  dplyr::group_by(treat)%>%
  dplyr::summarise(Vw_all_explore = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[2], 3))%>%
  as.data.frame()

# Female-Young
post.data.all.explore$delta.vw.explore.F.Y.m_i=
  with(post.data.all.explore, Vw.explore.F.Y.handled-Vw.explore.F.Y.injured)
post.data.all.explore$delta.vw.explore.F.Y.m_s=
  with(post.data.all.explore, Vw.explore.F.Y.handled-Vw.explore.F.Y.saline)
post.data.all.explore$delta.vw.explore.F.Y.m_l=
  with(post.data.all.explore, Vw.explore.F.Y.handled-Vw.explore.F.Y.lo)
post.data.all.explore$delta.vw.explore.F.Y.m_h=
  with(post.data.all.explore, Vw.explore.F.Y.handled-Vw.explore.F.Y.hi)
post.data.all.explore$delta.vw.explore.F.Y.i_s=
  with(post.data.all.explore, Vw.explore.F.Y.injured-Vw.explore.F.Y.saline)
post.data.all.explore$delta.vw.explore.F.Y.i_l=
  with(post.data.all.explore, Vw.explore.F.Y.injured-Vw.explore.F.Y.lo)
post.data.all.explore$delta.vw.explore.F.Y.i_h=
  with(post.data.all.explore, Vw.explore.F.Y.injured-Vw.explore.F.Y.hi)
post.data.all.explore$delta.vw.explore.F.Y.s_l=
  with(post.data.all.explore, Vw.explore.F.Y.saline-Vw.explore.F.Y.lo)
post.data.all.explore$delta.vw.explore.F.Y.s_h=
  with(post.data.all.explore, Vw.explore.F.Y.saline-Vw.explore.F.Y.hi)
post.data.all.explore$delta.vw.explore.F.Y.l_h=
  with(post.data.all.explore, Vw.explore.F.Y.lo-Vw.explore.F.Y.hi)

# Male-Young
post.data.all.explore$delta.vw.explore.M.Y.m_i=
  with(post.data.all.explore, Vw.explore.M.Y.handled-Vw.explore.M.Y.injured)
post.data.all.explore$delta.vw.explore.M.Y.m_s=
  with(post.data.all.explore, Vw.explore.M.Y.handled-Vw.explore.M.Y.saline)
post.data.all.explore$delta.vw.explore.M.Y.m_l=
  with(post.data.all.explore, Vw.explore.M.Y.handled-Vw.explore.M.Y.lo)
post.data.all.explore$delta.vw.explore.M.Y.m_h=
  with(post.data.all.explore, Vw.explore.M.Y.handled-Vw.explore.M.Y.hi)
post.data.all.explore$delta.vw.explore.M.Y.i_s=
  with(post.data.all.explore, Vw.explore.M.Y.injured-Vw.explore.M.Y.saline)
post.data.all.explore$delta.vw.explore.M.Y.i_l=
  with(post.data.all.explore, Vw.explore.M.Y.injured-Vw.explore.M.Y.lo)
post.data.all.explore$delta.vw.explore.M.Y.i_h=
  with(post.data.all.explore, Vw.explore.M.Y.injured-Vw.explore.M.Y.hi)
post.data.all.explore$delta.vw.explore.M.Y.s_l=
  with(post.data.all.explore, Vw.explore.M.Y.saline-Vw.explore.M.Y.lo)
post.data.all.explore$delta.vw.explore.M.Y.s_h=
  with(post.data.all.explore, Vw.explore.M.Y.saline-Vw.explore.M.Y.hi)
post.data.all.explore$delta.vw.explore.M.Y.l_h=
  with(post.data.all.explore, Vw.explore.M.Y.lo-Vw.explore.M.Y.hi)

# Female-Old
post.data.all.explore$delta.vw.explore.F.O.m_i=
  with(post.data.all.explore, Vw.explore.F.O.handled-Vw.explore.F.O.injured)
post.data.all.explore$delta.vw.explore.F.O.m_s=
  with(post.data.all.explore, Vw.explore.F.O.handled-Vw.explore.F.O.saline)
post.data.all.explore$delta.vw.explore.F.O.m_l=
  with(post.data.all.explore, Vw.explore.F.O.handled-Vw.explore.F.O.lo)
post.data.all.explore$delta.vw.explore.F.O.m_h=
  with(post.data.all.explore, Vw.explore.F.O.handled-Vw.explore.F.O.hi)
post.data.all.explore$delta.vw.explore.F.O.i_s=
  with(post.data.all.explore, Vw.explore.F.O.injured-Vw.explore.F.O.saline)
post.data.all.explore$delta.vw.explore.F.O.i_l=
  with(post.data.all.explore, Vw.explore.F.O.injured-Vw.explore.F.O.lo)
post.data.all.explore$delta.vw.explore.F.O.i_h=
  with(post.data.all.explore, Vw.explore.F.O.injured-Vw.explore.F.O.hi)
post.data.all.explore$delta.vw.explore.F.O.s_l=
  with(post.data.all.explore, Vw.explore.F.O.saline-Vw.explore.F.O.lo)
post.data.all.explore$delta.vw.explore.F.O.s_h=
  with(post.data.all.explore, Vw.explore.F.O.saline-Vw.explore.F.O.hi)
post.data.all.explore$delta.vw.explore.F.O.l_h=
  with(post.data.all.explore, Vw.explore.F.O.lo-Vw.explore.F.O.hi)

# Male-Old
post.data.all.explore$delta.vw.explore.M.O.m_i=
  with(post.data.all.explore, Vw.explore.M.O.handled-Vw.explore.M.O.injured)
post.data.all.explore$delta.vw.explore.M.O.m_s=
  with(post.data.all.explore, Vw.explore.M.O.handled-Vw.explore.M.O.saline)
post.data.all.explore$delta.vw.explore.M.O.m_l=
  with(post.data.all.explore, Vw.explore.M.O.handled-Vw.explore.M.O.lo)
post.data.all.explore$delta.vw.explore.M.O.m_h=
  with(post.data.all.explore, Vw.explore.M.O.handled-Vw.explore.M.O.hi)
post.data.all.explore$delta.vw.explore.M.O.i_s=
  with(post.data.all.explore, Vw.explore.M.O.injured-Vw.explore.M.O.saline)
post.data.all.explore$delta.vw.explore.M.O.i_l=
  with(post.data.all.explore, Vw.explore.M.O.injured-Vw.explore.M.O.lo)
post.data.all.explore$delta.vw.explore.M.O.i_h=
  with(post.data.all.explore, Vw.explore.M.O.injured-Vw.explore.M.O.hi)
post.data.all.explore$delta.vw.explore.M.O.s_l=
  with(post.data.all.explore, Vw.explore.M.O.saline-Vw.explore.M.O.lo)
post.data.all.explore$delta.vw.explore.M.O.s_h=
  with(post.data.all.explore, Vw.explore.M.O.saline-Vw.explore.M.O.hi)
post.data.all.explore$delta.vw.explore.M.O.l_h=
  with(post.data.all.explore, Vw.explore.M.O.lo-Vw.explore.M.O.hi)

vw.delta.all.explore = post.data.all.explore %>%
  dplyr::select(starts_with("delta.vw.explore.")) %>%
  pivot_longer(cols = starts_with("delta.vw.explore."),
               names_to = 'Contrast',
               names_prefix = "delta.vw.explore.",
               values_to ="Estimate")


all.explore.vw.delta = vw.delta.all.explore %>%
  dplyr::group_by(Contrast)%>%
  dplyr::summarise(vw.delta.all.explore = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[2], 3)) %>%
  as.data.frame()


## ---- end


# model comparison table
activity_comp.df<-data.frame(activity_comp)
explore_comp.df<-data.frame(explore_comp)
model.comp<-bind_rows(activity_comp.df, explore_comp.df) %>%
  select(elpd_diff, se_diff)

rownames(model.comp) <- c("(2) activity ~ sex + treatment + age + assay + (0 + sex:age||gr(id, by = treatment)), sigma ~ 0 + age:sex:treatment",
                          "(1) activity ~ sex * treatment * age + assay + (0 + sex:age||gr(id, by = treatment)), sigma ~ 0 + age:sex:treatment",
                          "(2) exploration ~ sex + treatment + age + assay + (0 + sex:age||gr(id, by = treatment)), sigma ~ 0 + age:sex:treatment",
                          "(1) exploration ~ sex * treatment * age + assay + (0 + sex:age||gr(id, by = treatment)), sigma ~ 0 + age:sex:treatment")

model.com2 <- model.comp %>%
  tibble::rownames_to_column(var = "model") %>%
  mutate_at(2:3, as.numeric) %>%
  dplyr::mutate_if(is.numeric, round, 1) %>%
  tidyr::unite("elpd difference ± SE", 2:3, sep=" ± ", remove=TRUE) %>% 
  huxtable::as_huxtable() %>%
  huxtable::insert_row(c("(a) Activity (total distance travelled, cm)", ""), after = 1) %>%
  huxtable::insert_row(c("(b) Exploration (latency to visit all quadrants, s)", ""), after = 4)



#### BIVARIATE CORRELATIONS ----

# factor-specific data frames

fem_y_man = merged_data_four[merged_data_four$sex=="Female" & merged_data_four$age=="Young" & merged_data_four$treatment=="handled",]
fem_y_inj = merged_data_four[merged_data_four$sex=="Female" & merged_data_four$age=="Young" & merged_data_four$treatment=="injured",]
fem_y_saline = merged_data_four[merged_data_four$sex=="Female" & merged_data_four$age=="Young" & merged_data_four$treatment=="saline",]
fem_y_low = merged_data_four[merged_data_four$sex=="Female" & merged_data_four$age=="Young" & merged_data_four$treatment=="low",]
fem_y_high = merged_data_four[merged_data_four$sex=="Female" & merged_data_four$age=="Young" & merged_data_four$treatment=="high",]
fem_o_man = merged_data_four[merged_data_four$sex=="Female" & merged_data_four$age=="Old" & merged_data_four$treatment=="handled",]
fem_o_inj = merged_data_four[merged_data_four$sex=="Female" & merged_data_four$age=="Old" & merged_data_four$treatment=="injured",]
fem_o_saline = merged_data_four[merged_data_four$sex=="Female" & merged_data_four$age=="Old" & merged_data_four$treatment=="saline",]
fem_o_low = merged_data_four[merged_data_four$sex=="Female" & merged_data_four$age=="Old" & merged_data_four$treatment=="low",]
fem_o_high = merged_data_four[merged_data_four$sex=="Female" & merged_data_four$age=="Old" & merged_data_four$treatment=="high",]
mal_y_man = merged_data_four[merged_data_four$sex=="Male" & merged_data_four$age=="Young" & merged_data_four$treatment=="handled",]
mal_y_inj = merged_data_four[merged_data_four$sex=="Male" & merged_data_four$age=="Young" & merged_data_four$treatment=="injured",]
mal_y_saline = merged_data_four[merged_data_four$sex=="Male" & merged_data_four$age=="Young" & merged_data_four$treatment=="saline",]
mal_y_low = merged_data_four[merged_data_four$sex=="Male" & merged_data_four$age=="Young" & merged_data_four$treatment=="low",]
mal_y_high = merged_data_four[merged_data_four$sex=="Male" & merged_data_four$age=="Young" & merged_data_four$treatment=="high",]
mal_o_man = merged_data_four[merged_data_four$sex=="Male" & merged_data_four$age=="Old" & merged_data_four$treatment=="handled",]
mal_o_inj = merged_data_four[merged_data_four$sex=="Male" & merged_data_four$age=="Old" & merged_data_four$treatment=="injured",]
mal_o_saline = merged_data_four[merged_data_four$sex=="Male" & merged_data_four$age=="Old" & merged_data_four$treatment=="saline",]
mal_o_low = merged_data_four[merged_data_four$sex=="Male" & merged_data_four$age=="Old" & merged_data_four$treatment=="low",]
mal_o_high = merged_data_four[merged_data_four$sex=="Male" & merged_data_four$age=="Old" & merged_data_four$treatment=="high",]

behav<-bf(mvbind(distance.t,explore.t)~1 + (1|p|id)) + set_rescor(TRUE)

#### female-young-handled ####
fem_y_man.brms <- brm(behav, data = fem_y_man,
            cores = 4,
            chains = 4,
            warmup = 1000,
            iter = 6000,
            thin = 2,
            seed = 12345,
            file = 'data/processed/fem_y_man.brms',
            control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(fem_y_man.brms)

draws.fem_y_man <- fem_y_man.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### female-young-injured ####
fem_y_inj.brms <- brm(behav, data = fem_y_inj,
                        cores = 4,
                        chains = 4,
                        warmup = 1000,
                        iter = 6000,
                        thin = 2,
                        seed = 12345,
                        file = 'data/processed/fem_y_inj.brms',
                        control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(fem_y_inj.brms)

draws.fem_y_inj <- fem_y_inj.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### female-young-saline ####
fem_y_sal.brms <- brm(behav, data = fem_y_saline,
                      cores = 4,
                      chains = 4,
                      warmup = 1000,
                      iter = 6000,
                      thin = 2,
                      seed = 12345,
                      file = 'data/processed/fem_y_sal.brms',
                      control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(fem_y_sal.brms)

draws.fem_y_sal <- fem_y_sal.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### female-young-low ####
fem_y_low.brms <- brm(behav, data = fem_y_low,
                       cores = 4,
                       chains = 4,
                       warmup = 1000,
                       iter = 6000,
                       thin = 2,
                       seed = 12345,
                       file = 'data/processed/fem_y_low.brms',
                       control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(fem_y_low.brms)

draws.fem_y_low <- fem_y_low.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable)  %>%
  mutate(category = rep(c('among','within'), 2, length.out = n())) 

#### female-young-high ####
fem_y_high.brms <- brm(behav, data = fem_y_high,
                      cores = 4,
                      chains = 4,
                      warmup = 1000,
                      iter = 6000,
                      thin = 2,
                      seed = 12345,
                      file = 'data/processed/fem_y_high.brms',
                      control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(fem_y_high.brms)

draws.fem_y_high <- fem_y_high.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n())) 

#### female-old-handled ####
fem_o_man.brms <- brm(behav, data = fem_o_man,
                       cores = 4,
                       chains = 4,
                       warmup = 1000,
                       iter = 6000,
                       thin = 2,
                       seed = 12345,
                       file = 'data/processed/fem_o_man.brms',
                       control = list(max_treedepth = 15, adapt_delta = 0.9999))

summary(fem_o_man.brms)

draws.fem_o_man <- fem_o_man.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable)  %>%
  mutate(category = rep(c('among','within'), 2, length.out = n())) 

#### female-old-injured ####
fem_o_inj.brms <- brm(behav, data = fem_o_inj,
                      cores = 4,
                      chains = 4,
                      warmup = 1000,
                      iter = 6000,
                      thin = 2,
                      seed = 12345,
                      file = 'data/processed/fem_o_inj.brms',
                      control = list(max_treedepth = 15, adapt_delta = 0.9999))

summary(fem_o_inj.brms)

draws.fem_o_inj <- fem_o_inj.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable)  %>%
  mutate(category = rep(c('among','within'), 2, length.out = n())) 

#### female-old-saline ####
fem_o_sal.brms <- brm(behav, data = fem_o_saline,
                      cores = 4,
                      chains = 4,
                      warmup = 1000,
                      iter = 6000,
                      thin = 2,
                      seed = 12345,
                      file = 'data/processed/fem_o_sal.brms',
                      control = list(max_treedepth = 15, adapt_delta = 0.9999))

summary(fem_o_sal.brms)

draws.fem_o_sal <- fem_o_sal.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable)  %>%
  mutate(category = rep(c('among','within'), 2, length.out = n())) 

#### female-old-low ####
fem_o_low.brms <- brm(behav, data = fem_o_low,
                      cores = 4,
                      chains = 4,
                      warmup = 1000,
                      iter = 6000,
                      thin = 2,
                      seed = 12345,
                      file = 'data/processed/fem_o_low.brms',
                      control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(fem_o_low.brms)

draws.fem_o_low <- fem_o_low.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### female-old-high ####
fem_o_high.brms <- brm(behav, data = fem_o_high,
                       cores = 4,
                       chains = 4,
                       warmup = 1000,
                       iter = 6000,
                       thin = 2,
                       seed = 12345,
                       file = 'data/processed/fem_o_high.brms',
                       control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(fem_o_high.brms)

draws.fem_o_high <- fem_o_high.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### male-young-handled ####
mal_y_man.brms <- brm(behav, data = mal_y_man,
                       cores = 4,
                       chains = 4,
                       warmup = 1000,
                       iter = 6000,
                       thin = 2,
                       seed = 12345,
                      file = 'data/processed/mal_y_man.brms',
                       control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(mal_y_man.brms)

draws.mal_y_man <- mal_y_man.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### male-young-injury ####
mal_y_inj.brms <- brm(behav, data = mal_y_inj,
                      cores = 4,
                      chains = 4,
                      warmup = 1000,
                      iter = 6000,
                      thin = 2,
                      seed = 12345,
                      file = 'data/processed/mal_y_inj.brms',
                      control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(mal_y_inj.brms)

draws.mal_y_inj <- mal_y_inj.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### male-young-saline ####
mal_y_sal.brms <- brm(behav, data = mal_y_saline,
                      cores = 4,
                      chains = 4,
                      warmup = 1000,
                      iter = 6000,
                      thin = 2,
                      seed = 12345,
                      file = 'data/processed/mal_y_sal.brms',
                      control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(mal_y_sal.brms)

draws.mal_y_sal <- mal_y_sal.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### male-young-low ####
mal_y_low.brms <- brm(behav, data = mal_y_low,
                      cores = 4,
                      chains = 4,
                      warmup = 1000,
                      iter = 6000,
                      thin = 2,
                      seed = 12345,
                      file = 'data/processed/mal_y_low.brms',
                      control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(mal_y_low.brms)

draws.mal_y_low <- mal_y_low.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### male-young-high ####
mal_y_high.brms <- brm(behav, data = mal_y_high,
                       cores = 4,
                       chains = 4,
                       warmup = 1000,
                       iter = 6000,
                       thin = 2,
                       seed = 12345,
                       file = 'data/processed/mal_y_high.brms',
                       control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(mal_y_high.brms)

draws.mal_y_high <- mal_y_high.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### male-old-handled ####
mal_o_man.brms <- brm(behav, data = mal_o_man,
                       cores = 4,
                       chains = 4,
                       warmup = 1000,
                       iter = 6000,
                       thin = 2,
                       seed = 12345,
                      file = 'data/processed/mal_o_man.brms',
                       control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(mal_o_man.brms)

draws.mal_o_man <- mal_o_man.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### male-old-injured ####
mal_o_inj.brms <- brm(behav, data = mal_o_inj,
                      cores = 4,
                      chains = 4,
                      warmup = 1000,
                      iter = 6000,
                      thin = 2,
                      seed = 12345,
                      file = 'data/processed/mal_o_inj.brms',
                      control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(mal_o_inj.brms)

draws.mal_o_inj <- mal_o_inj.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### male-old-saline ####
mal_o_sal.brms <- brm(behav, data = mal_o_saline,
                      cores = 4,
                      chains = 4,
                      warmup = 1000,
                      iter = 6000,
                      thin = 2,
                      seed = 12345,
                      file = 'data/processed/mal_o_sal.brms',
                      control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(mal_o_sal.brms)

draws.mal_o_sal <- mal_o_sal.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### male-old-low ####
mal_o_low.brms <- brm(behav, data = mal_o_low,
                      cores = 4,
                      chains = 4,
                      warmup = 1000,
                      iter = 6000,
                      thin = 2,
                      seed = 12345,
                      file = 'data/processed/mal_o_low.brms',
                      control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(mal_o_low.brms)

draws.mal_o_low <- mal_o_low.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

#### male-old-high ####
mal_o_high.brms <- brm(behav, data = mal_o_high,
                       cores = 4,
                       chains = 4,
                       warmup = 1000,
                       iter = 6000,
                       thin = 2,
                       seed = 12345,
                       file = 'data/processed/mal_o_high.brms',
                       control = list(max_treedepth = 15, adapt_delta = 0.9999))
summary(mal_o_high.brms)


draws.mal_o_high <- mal_o_high.brms %>%
  gather_draws(
    cor_id__distancet_Intercept__exploret_Intercept,
    rescor__distancet__exploret
  ) %>%
  median_qi %>%
  ungroup() %>%
  mutate(
    .variable = case_match(
      .variable,
      "cor_id__distancet_Intercept__exploret_Intercept" ~ "activity—exploration",
      "rescor__distancet__exploret" ~ "activity—exploration"
    )
  )  %>%
  dplyr::select(-.width, -.point, -.interval, -.variable) %>%
  mutate(category = rep(c('among','within'), 2, length.out = n()))

##### correlations #####

table_pooled_r<-data.frame(.value=c(cor.test(fem_y_man$distance.t, fem_y_man$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(fem_y_inj$distance.t, fem_y_inj$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(fem_y_saline$distance.t, fem_y_saline$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(fem_y_low$distance.t, fem_y_low$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(fem_y_high$distance.t, fem_y_high$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(fem_o_man$distance.t, fem_o_man$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(fem_o_inj$distance.t, fem_o_inj$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(fem_o_saline$distance.t, fem_o_saline$explore.t,  method = "pearson", use = "complete.obs")$estimate,                                   
                                    cor.test(fem_o_low$distance.t, fem_o_low$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(fem_o_high$distance.t, fem_o_high$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(mal_y_man$distance.t, mal_y_man$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(mal_y_inj$distance.t, mal_y_inj$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(mal_y_saline$distance.t, mal_y_saline$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(mal_y_low$distance.t, mal_y_low$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(mal_y_high$distance.t, mal_y_high$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(mal_o_man$distance.t, mal_o_man$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(mal_o_inj$distance.t, mal_o_inj$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(mal_o_saline$distance.t, mal_o_saline$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(mal_o_low$distance.t, mal_o_low$explore.t,  method = "pearson", use = "complete.obs")$estimate,
                                    cor.test(mal_o_high$distance.t, mal_o_high$explore.t,  method = "pearson", use = "complete.obs")$estimate),
                           .lower=c(cor.test(fem_y_man$distance.t, fem_y_man$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(fem_y_inj$distance.t, fem_y_inj$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(fem_y_saline$distance.t, fem_y_saline$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(fem_y_low$distance.t, fem_y_low$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(fem_y_high$distance.t, fem_y_high$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(fem_o_man$distance.t, fem_o_man$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(fem_o_inj$distance.t, fem_o_inj$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(fem_o_saline$distance.t, fem_o_saline$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],                                   
                                    cor.test(fem_o_low$distance.t, fem_o_low$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(fem_o_high$distance.t, fem_o_high$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(mal_y_man$distance.t, mal_y_man$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(mal_y_inj$distance.t, mal_y_inj$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(mal_y_saline$distance.t, mal_y_saline$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(mal_y_low$distance.t, mal_y_low$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(mal_y_high$distance.t, mal_y_high$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(mal_o_man$distance.t, mal_o_man$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(mal_o_inj$distance.t, mal_o_inj$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(mal_o_saline$distance.t, mal_o_saline$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(mal_o_low$distance.t, mal_o_low$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]],
                                    cor.test(mal_o_high$distance.t, mal_o_high$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[1]]),
                           .upper=c(cor.test(fem_y_man$distance.t, fem_y_man$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(fem_y_inj$distance.t, fem_y_inj$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(fem_y_saline$distance.t, fem_y_saline$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(fem_y_low$distance.t, fem_y_low$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(fem_y_high$distance.t, fem_y_high$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(fem_o_man$distance.t, fem_o_man$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(fem_o_inj$distance.t, fem_o_inj$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(fem_o_saline$distance.t, fem_o_saline$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],                                   
                                    cor.test(fem_o_low$distance.t, fem_o_low$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(fem_o_high$distance.t, fem_o_high$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(mal_y_man$distance.t, mal_y_man$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(mal_y_inj$distance.t, mal_y_inj$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(mal_y_saline$distance.t, mal_y_saline$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(mal_y_low$distance.t, mal_y_low$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(mal_y_high$distance.t, mal_y_high$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(mal_o_man$distance.t, mal_o_man$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(mal_o_inj$distance.t, mal_o_inj$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(mal_o_saline$distance.t, mal_o_saline$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(mal_o_low$distance.t, mal_o_low$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]],
                                    cor.test(mal_o_high$distance.t, mal_o_high$explore.t,  method = "pearson", use = "complete.obs")$conf.int[[2]])) %>%
                                    mutate(category = rep('raw')) %>%
                                    mutate(sex=rep(c("female", "male"), each = 10)) %>%
                                    mutate(age=rep(c("young", "old"), each = 5, times=2)) %>%
                                    mutate(treatment = rep(c("handled", "injured", "saline","low" ,"high"), each = 1, times=4)) 

partition.table <- rbind(draws.fem_y_man, draws.fem_y_inj, draws.fem_y_sal, draws.fem_y_low, draws.fem_y_high,
                         draws.fem_o_man, draws.fem_o_inj, draws.fem_o_sal, draws.fem_o_low, draws.fem_o_high,
                         draws.mal_y_man, draws.mal_y_inj, draws.mal_y_sal, draws.mal_y_low, draws.mal_y_high,
                         draws.mal_o_man, draws.mal_o_inj, draws.mal_o_sal, draws.mal_o_low, draws.mal_o_high) %>%
  mutate(sex=rep(c("female", "male"), each = 20)) %>%
  mutate(age=rep(c("young", "old"), each = 10, times=2)) %>%
  mutate(treatment = rep(c("handled", "injured", "saline", "low" ,"high"), each = 2, times=4)) 

correlation.table <- rbind(table_pooled_r,partition.table) %>%
  mutate(treatment = fct_relevel(treatment, 
                            "handled", "injured", "saline", "low", "high")) %>%
  mutate(age = factor(age, levels=c("young", "old")))

saveRDS(correlation.table, file = "data/processed/correlation.table.rda")
correlation.table <- readRDS(file = "data/processed/correlation.table.rda")

sex.lab <- c("Female", "Male")
names(sex.lab) <- c("female", "male")

age.lab <- c("Young", "Old")
names(age.lab) <- c("young", "old")

pooled.correlation.plot<- ggplot(correlation.table,aes(y = treatment, x = .value, xmin = .lower, xmax = .upper, shape=category)) +
  geom_vline(aes(xintercept=0),colour="grey", linetype='dotted') +
  geom_pointinterval(position=position_dodge(0.4), colour="slateblue2", fatten_point = 4) +
  xlab(expression(~italic(r)~" ± 95% CI")) +
  ylab("Treatment") +
  xlim(c(-1,1)) +
  facet_grid(sex ~ age,labeller = labeller(age=age.lab,sex = sex.lab)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_text(size = 12),
  axis.title.x = element_text(size=14,face="bold"),
  axis.title.y = element_text(size=14,face="bold"),
  axis.text.x = element_text(size=12, angle=45, vjust=0.5),
  axis.text.y = element_text(size=12), strip.background = element_rect(fill="white")) +
  scale_shape_manual(values=c(15,16,17),
                     name="Correlation",
                     breaks=c("raw", "among", "within"),
                     labels=c(expression(~italic(r)["p"]), expression(~italic(r)["ind"]), expression(~italic(r)["e"])))
ggsave(filename="figure_4.png", width=16.83, height=16.83, dpi=800,antialias="default")

#### Table ####

contrasts.va.table.activity <-all.activity.va.delta %>%
  rename(Va=va.delta.all.activity) %>%
  mutate(category=rep("activity",40)) %>%
  mutate(sex=rep(c("female","male"),each=20)) %>%
  mutate(age=rep(c("o", "y"), each = 10, times=2))

contrasts.va.table.explore <-all.explore.va.delta %>%
  rename(Va=va.delta.all.explore) %>%
  mutate(category=rep("explore",40)) %>%
  mutate(sex=rep(c("female","male"),each=20)) %>%
  mutate(age=rep(c("o", "y"), each = 10, times=2))

order<- c("Handled vs. Injured", "Handled vs. Saline", "Handled vs. Low", "Handled vs. High", "Injured vs. Saline", "Injured vs. Low", "Injured vs. High", "Saline vs. Low", "Saline vs. High")
contrasts.table.Va <-rbind(contrasts.va.table.activity, contrasts.va.table.explore) %>%
  mutate(Contrast = gsub("i_h", "Injured vs. High", Contrast)) %>%
  mutate(Contrast = gsub("i_s", "Injured vs. Saline", Contrast)) %>%
  mutate(Contrast = gsub("i_l", "Injured vs. Low", Contrast)) %>%
  mutate(Contrast = gsub("i_s", "Injured vs. Saline", Contrast)) %>%
  mutate(Contrast = gsub("l_h", "Low vs. High", Contrast)) %>%
  mutate(Contrast = gsub("m_h", "Handled vs. High", Contrast)) %>%
  mutate(Contrast = gsub("m_i", "Handled vs. Injured", Contrast)) %>%
  mutate(Contrast = gsub("m_l", "Handled vs. Low", Contrast)) %>%
  mutate(Contrast = gsub("m_s", "Handled vs. Saline", Contrast)) %>%
  mutate(Contrast = gsub("s_h", "Saline vs. High", Contrast)) %>%
  mutate(Contrast = gsub("s_l", "Saline vs. Low", Contrast)) %>%
  mutate(Contrast = str_remove(Contrast, "M.O.|F.O.|M.Y.|F.Y.")) %>%
  mutate(Contrast = fct_relevel(Contrast, order)) %>%
  arrange(category, sex,desc(age),Contrast) %>%
  pivot_wider(values_from = c(Va,lowerCI,upperCI), names_from = c("age","category")) %>%
  relocate(sex,Contrast,Va_y_activity,lowerCI_y_activity, upperCI_y_activity, Va_y_explore,lowerCI_y_explore,upperCI_y_explore,
           Va_o_activity,lowerCI_o_activity, upperCI_o_activity, Va_o_explore,lowerCI_o_explore,upperCI_o_explore) %>%
  mutate(across(where(is.numeric), ~ sprintf("%0.3f", .x))) %>%
  unite(col='Va_y_activity_CI', c('lowerCI_y_activity', 'upperCI_y_activity'), sep=', ') %>%
  unite(col='Va_y_explore_CI', c('lowerCI_y_explore', 'upperCI_y_explore'), sep=', ') %>%
  unite(col='Va_o_activity_CI', c('lowerCI_o_activity', 'upperCI_o_activity'), sep=', ') %>%
  unite(col='Va_o_explore_CI', c('lowerCI_o_explore', 'upperCI_o_explore'), sep=', ') %>%
  mutate(Va_y_activity_CI = paste0("(", Va_y_activity_CI, ")")) %>%
  mutate(Va_y_explore_CI = paste0("(", Va_y_explore_CI, ")")) %>%
  mutate(Va_o_activity_CI = paste0("(", Va_o_activity_CI, ")")) %>%
  mutate(Va_o_explore_CI = paste0("(", Va_o_explore_CI, ")")) %>%
  group_by(sex) %>%
  arrange(match(Contrast, c("Control vs. LPS Low", "Control vs. LPS High", "LPS Low vs. LPS High")),.by_group = TRUE) %>%
  mutate(sex= recode(sex, female="(a) Females", male="(b) Males")) %>%
  unite(col='Va_y_activity', c('Va_y_activity', 'Va_y_activity_CI'), sep=' ') %>%
  unite(col='Va_y_explore', c('Va_y_explore', 'Va_y_explore_CI'), sep=' ') %>%
  unite(col='Va_o_activity', c('Va_o_activity', 'Va_o_activity_CI'), sep=' ') %>%
  unite(col='Va_o_explore', c('Va_o_explore', 'Va_o_explore_CI'), sep=' ') %>%
  as_grouped_data(groups = c("sex"), columns=NULL) 
save(contrasts.table.Va, file = "data/processed/contrasts.table.Va.RData")

contrasts.vw.table.activity <-all.activity.vw.delta %>%
  rename(Vw=vw.delta.all.activity) %>%
  mutate(category=rep("activity",40)) %>%
  mutate(sex=rep(c("female","male"),each=20)) %>%
  mutate(age=rep(c("o", "y"), each = 10, times=2))

contrasts.vw.table.explore <-all.explore.vw.delta %>%
  rename(Vw=vw.delta.all.explore) %>%
  mutate(category=rep("explore",40)) %>%
  mutate(sex=rep(c("female","male"),each=20)) %>%
  mutate(age=rep(c("o", "y"), each = 10, times=2))

contrasts.table.Vw <-rbind(contrasts.vw.table.activity, contrasts.vw.table.explore) %>% 
  mutate(Contrast = gsub("i_h", "Injured vs. High", Contrast)) %>%
  mutate(Contrast = gsub("i_s", "Injured vs. Saline", Contrast)) %>%
  mutate(Contrast = gsub("i_l", "Injured vs. Low", Contrast)) %>%
  mutate(Contrast = gsub("i_s", "Injured vs. Saline", Contrast)) %>%
  mutate(Contrast = gsub("l_h", "Low vs. High", Contrast)) %>%
  mutate(Contrast = gsub("m_h", "Handled vs. High", Contrast)) %>%
  mutate(Contrast = gsub("m_i", "Handled vs. Injured", Contrast)) %>%
  mutate(Contrast = gsub("m_l", "Handled vs. Low", Contrast)) %>%
  mutate(Contrast = gsub("m_s", "Handled vs. Saline", Contrast)) %>%
  mutate(Contrast = gsub("s_h", "Saline vs. High", Contrast)) %>%
  mutate(Contrast = gsub("s_l", "Saline vs. Low", Contrast)) %>%
  mutate(Contrast = str_remove(Contrast, "M.O.|F.O.|M.Y.|F.Y.")) %>%
  mutate(Contrast = fct_relevel(Contrast, order)) %>%
  arrange(category, sex,desc(age),Contrast) %>%
  pivot_wider(values_from = c(Vw,lowerCI,upperCI), names_from = c("age","category")) %>%
  relocate(sex,Contrast,Vw_y_activity,lowerCI_y_activity, upperCI_y_activity, Vw_y_explore,lowerCI_y_explore,upperCI_y_explore,
           Vw_o_activity,lowerCI_o_activity, upperCI_o_activity, Vw_o_explore,lowerCI_o_explore,upperCI_o_explore) %>%
  mutate(across(where(is.numeric), ~ sprintf("%0.3f", .x))) %>%
  unite(col='Vw_y_activity_CI', c('lowerCI_y_activity', 'upperCI_y_activity'), sep=', ') %>%
  unite(col='Vw_y_explore_CI', c('lowerCI_y_explore', 'upperCI_y_explore'), sep=', ') %>%
  unite(col='Vw_o_activity_CI', c('lowerCI_o_activity', 'upperCI_o_activity'), sep=', ') %>%
  unite(col='Vw_o_explore_CI', c('lowerCI_o_explore', 'upperCI_o_explore'), sep=', ') %>%
  mutate(Vw_y_activity_CI = paste0("(", Vw_y_activity_CI, ")")) %>%
  mutate(Vw_y_explore_CI = paste0("(", Vw_y_explore_CI, ")")) %>%
  mutate(Vw_o_activity_CI = paste0("(", Vw_o_activity_CI, ")")) %>%
  mutate(Vw_o_explore_CI = paste0("(", Vw_o_explore_CI, ")")) %>%
  group_by(sex) %>%
  arrange(match(Contrast, c("Control vs. LPS Low", "Control vs. LPS High", "LPS Low vs. LPS High")),.by_group = TRUE) %>%
  mutate(sex= recode(sex, female="(a) Females", male="(b) Males")) %>%
  unite(col='Vw_y_activity', c('Vw_y_activity', 'Vw_y_activity_CI'), sep=' ') %>%
  unite(col='Vw_y_explore', c('Vw_y_explore', 'Vw_y_explore_CI'), sep=' ') %>%
  unite(col='Vw_o_activity', c('Vw_o_activity', 'Vw_o_activity_CI'), sep=' ') %>%
  unite(col='Vw_o_explore', c('Vw_o_explore', 'Vw_o_explore_CI'), sep=' ') %>%
  as_grouped_data(groups = c("sex"), columns=NULL) 
save(contrasts.table.Vw, file = "data/processed/contrasts.table.Vw.RData")

#### REPEATABILITy ####

# activity
# female - young
fem_y_man.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = fem_y_man, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_y_man.rpts.activity)

fem_y_inj.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = fem_y_inj, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_y_inj.rpts.activity)

fem_y_sal.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = fem_y_saline, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_y_sal.rpts.activity)

fem_y_low.rpts.activity <- rpt(distance.t ~ (1 | id), 
                              grname = "id", data = fem_y_low, 
                              datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_y_low.rpts.activity)

fem_y_high.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = fem_y_high, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_y_high.rpts.activity)

# female - old
fem_o_man.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = fem_o_man, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_o_man.rpts.activity)

fem_o_inj.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = fem_o_inj, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_o_inj.rpts.activity)

fem_o_sal.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = fem_o_saline, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_o_sal.rpts.activity)

fem_o_low.rpts.activity <- rpt(distance.t ~ (1 | id), 
                              grname = "id", data = fem_o_low, 
                              datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_o_low.rpts.activity)

fem_o_high.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = fem_o_high, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_o_high.rpts.activity)

# male-young
mal_y_man.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = mal_y_man, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_y_man.rpts.activity)

mal_y_inj.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = mal_y_inj, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_y_inj.rpts.activity)

mal_y_sal.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = mal_y_saline, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_y_sal.rpts.activity)

mal_y_low.rpts.activity <- rpt(distance.t ~ (1 | id), 
                              grname = "id", data = mal_y_low, 
                              datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_y_low.rpts.activity)

mal_y_high.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = mal_y_high, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_y_high.rpts.activity)

# male - old
mal_o_man.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = mal_o_man, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_o_man.rpts.activity)

mal_o_inj.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = mal_o_inj, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_o_inj.rpts.activity)

mal_o_sal.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = mal_o_saline, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_o_sal.rpts.activity)

mal_o_low.rpts.activity <- rpt(distance.t ~ (1 | id), 
                              grname = "id", data = mal_o_low, 
                              datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_o_low.rpts.activity)

mal_o_high.rpts.activity <- rpt(distance.t ~ (1 | id), 
                               grname = "id", data = mal_o_high, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_o_high.rpts.activity)

activity_R<-data.frame(.value=c(fem_y_man.rpts.activity$R[[1]], fem_y_inj.rpts.activity$R[[1]],fem_y_sal.rpts.activity$R[[1]],fem_y_low.rpts.activity$R[[1]], fem_y_high.rpts.activity$R[[1]],
                                fem_o_man.rpts.activity$R[[1]], fem_o_inj.rpts.activity$R[[1]],fem_o_sal.rpts.activity$R[[1]], fem_o_low.rpts.activity$R[[1]], fem_o_high.rpts.activity$R[[1]],
                                mal_y_man.rpts.activity$R[[1]], mal_y_inj.rpts.activity$R[[1]],mal_y_sal.rpts.activity$R[[1]],mal_y_low.rpts.activity$R[[1]], mal_y_high.rpts.activity$R[[1]],
                                mal_o_man.rpts.activity$R[[1]], mal_o_inj.rpts.activity$R[[1]],mal_o_sal.rpts.activity$R[[1]], mal_o_low.rpts.activity$R[[1]], mal_o_high.rpts.activity$R[[1]]),
                           .lower=c(fem_y_man.rpts.activity$CI_emp[[1]], fem_y_inj.rpts.activity$CI_emp[[1]],fem_y_sal.rpts.activity$CI_emp[[1]], fem_y_low.rpts.activity$CI_emp[[1]], fem_y_high.rpts.activity$CI_emp[[1]],
                                    fem_o_man.rpts.activity$CI_emp[[1]], fem_o_inj.rpts.activity$CI_emp[[1]],fem_o_sal.rpts.activity$CI_emp[[1]], fem_o_low.rpts.activity$CI_emp[[1]], fem_o_high.rpts.activity$CI_emp[[1]],
                                    mal_y_man.rpts.activity$CI_emp[[1]], mal_y_inj.rpts.activity$CI_emp[[1]],mal_y_sal.rpts.activity$CI_emp[[1]],mal_y_low.rpts.activity$CI_emp[[1]], mal_y_high.rpts.activity$CI_emp[[1]],
                                    mal_o_man.rpts.activity$CI_emp[[1]], mal_o_inj.rpts.activity$CI_emp[[1]], mal_o_sal.rpts.activity$CI_emp[[1]],mal_o_low.rpts.activity$CI_emp[[1]], mal_o_high.rpts.activity$CI_emp[[1]]),
                           .upper=c(fem_y_man.rpts.activity$CI_emp[[2]], fem_y_inj.rpts.activity$CI_emp[[2]],fem_y_sal.rpts.activity$CI_emp[[2]], fem_y_low.rpts.activity$CI_emp[[2]], fem_y_high.rpts.activity$CI_emp[[2]],
                                    fem_o_man.rpts.activity$CI_emp[[2]], fem_o_inj.rpts.activity$CI_emp[[2]],fem_o_sal.rpts.activity$CI_emp[[2]], fem_o_low.rpts.activity$CI_emp[[2]], fem_o_high.rpts.activity$CI_emp[[2]],
                                    mal_y_man.rpts.activity$CI_emp[[2]], mal_y_inj.rpts.activity$CI_emp[[2]],mal_y_sal.rpts.activity$CI_emp[[2]], mal_y_low.rpts.activity$CI_emp[[2]], mal_y_high.rpts.activity$CI_emp[[2]],
                                    mal_o_man.rpts.activity$CI_emp[[2]], mal_o_inj.rpts.activity$CI_emp[[2]], mal_o_sal.rpts.activity$CI_emp[[2]], mal_o_low.rpts.activity$CI_emp[[2]], mal_o_high.rpts.activity$CI_emp[[2]])) %>%
  mutate(sex=rep(c("female", "male"), each = 10)) %>%
  mutate(age=rep(c("young", "old"), each = 5, times=2)) %>%
  mutate(treatment = rep(c("handled","injured", "saline", "low" ,"high"), each = 1, times=4)) 

# exploration
# female - young
fem_y_man.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = fem_y_man, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_y_man.rpts.explore)

fem_y_inj.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = fem_y_inj, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_y_inj.rpts.explore)

fem_y_sal.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = fem_y_saline, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_y_sal.rpts.explore)

fem_y_low.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = fem_y_low, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_y_low.rpts.explore)

fem_y_high.rpts.explore <- rpt(explore.t ~ (1 | id), 
                                grname = "id", data = fem_y_high, 
                                datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_y_high.rpts.explore)

# female - old
fem_o_man.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = fem_o_man, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_o_man.rpts.explore)

fem_o_inj.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = fem_o_inj, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_o_inj.rpts.explore)

fem_o_sal.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = fem_o_saline, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_o_sal.rpts.explore)

fem_o_low.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = fem_o_low, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_o_low.rpts.explore)

fem_o_high.rpts.explore <- rpt(explore.t ~ (1 | id), 
                                grname = "id", data = fem_o_high, 
                                datatype = "Gaussian", nboot = 100, npermut = 100)
summary(fem_o_high.rpts.explore)

# male-young
mal_y_man.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = mal_y_man, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_y_man.rpts.explore)

mal_y_inj.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = mal_y_inj, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_y_inj.rpts.explore)

mal_y_sal.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = mal_y_saline, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_y_sal.rpts.explore)

mal_y_low.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = mal_y_low, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_y_low.rpts.explore)

mal_y_high.rpts.explore <- rpt(explore.t ~ (1 | id), 
                                grname = "id", data = mal_y_high, 
                                datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_y_high.rpts.explore)

# male - old
mal_o_man.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = mal_o_man, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_o_man.rpts.explore)

mal_o_inj.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = mal_o_inj, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_o_inj.rpts.explore)

mal_o_sal.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = mal_o_saline, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_o_sal.rpts.explore)

mal_o_low.rpts.explore <- rpt(explore.t ~ (1 | id), 
                               grname = "id", data = mal_o_low, 
                               datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_o_low.rpts.explore)

mal_o_high.rpts.explore <- rpt(explore.t ~ (1 | id), 
                                grname = "id", data = mal_o_high, 
                                datatype = "Gaussian", nboot = 100, npermut = 100)
summary(mal_o_high.rpts.explore)

explore_R<-data.frame(.value=c(fem_y_man.rpts.explore$R[[1]], fem_y_inj.rpts.explore$R[[1]],fem_y_sal.rpts.explore$R[[1]],fem_y_low.rpts.explore$R[[1]], fem_y_high.rpts.explore$R[[1]],
                                fem_o_man.rpts.explore$R[[1]], fem_o_inj.rpts.explore$R[[1]],fem_o_sal.rpts.explore$R[[1]], fem_o_low.rpts.explore$R[[1]], fem_o_high.rpts.explore$R[[1]],
                                mal_y_man.rpts.explore$R[[1]], mal_y_inj.rpts.explore$R[[1]],mal_y_sal.rpts.explore$R[[1]],mal_y_low.rpts.explore$R[[1]], mal_y_high.rpts.explore$R[[1]],
                                mal_o_man.rpts.explore$R[[1]], mal_o_inj.rpts.explore$R[[1]],mal_o_sal.rpts.explore$R[[1]], mal_o_low.rpts.explore$R[[1]], mal_o_high.rpts.explore$R[[1]]),
                       .lower=c(fem_y_man.rpts.explore$CI_emp[[1]], fem_y_inj.rpts.explore$CI_emp[[1]],fem_y_sal.rpts.explore$CI_emp[[1]], fem_y_low.rpts.explore$CI_emp[[1]], fem_y_high.rpts.explore$CI_emp[[1]],
                                fem_o_man.rpts.explore$CI_emp[[1]], fem_o_inj.rpts.explore$CI_emp[[1]],fem_o_sal.rpts.explore$CI_emp[[1]], fem_o_low.rpts.explore$CI_emp[[1]], fem_o_high.rpts.explore$CI_emp[[1]],
                                mal_y_man.rpts.explore$CI_emp[[1]], mal_y_inj.rpts.explore$CI_emp[[1]],mal_y_sal.rpts.explore$CI_emp[[1]],mal_y_low.rpts.explore$CI_emp[[1]], mal_y_high.rpts.explore$CI_emp[[1]],
                                mal_o_man.rpts.explore$CI_emp[[1]], mal_o_inj.rpts.explore$CI_emp[[1]], mal_o_sal.rpts.explore$CI_emp[[1]],mal_o_low.rpts.explore$CI_emp[[1]], mal_o_high.rpts.explore$CI_emp[[1]]),
                       .upper=c(fem_y_man.rpts.explore$CI_emp[[2]], fem_y_inj.rpts.explore$CI_emp[[2]],fem_y_sal.rpts.explore$CI_emp[[2]], fem_y_low.rpts.explore$CI_emp[[2]], fem_y_high.rpts.explore$CI_emp[[2]],
                                fem_o_man.rpts.explore$CI_emp[[2]], fem_o_inj.rpts.explore$CI_emp[[2]],fem_o_sal.rpts.explore$CI_emp[[2]], fem_o_low.rpts.explore$CI_emp[[2]], fem_o_high.rpts.explore$CI_emp[[2]],
                                mal_y_man.rpts.explore$CI_emp[[2]], mal_y_inj.rpts.explore$CI_emp[[2]],mal_y_sal.rpts.explore$CI_emp[[2]], mal_y_low.rpts.explore$CI_emp[[2]], mal_y_high.rpts.explore$CI_emp[[2]],
                                mal_o_man.rpts.explore$CI_emp[[2]], mal_o_inj.rpts.explore$CI_emp[[2]], mal_o_sal.rpts.explore$CI_emp[[2]], mal_o_low.rpts.explore$CI_emp[[2]], mal_o_high.rpts.explore$CI_emp[[2]])) %>%
  mutate(sex=rep(c("female", "male"), each = 10)) %>%
  mutate(age=rep(c("young", "old"), each = 5, times=2)) %>%
  mutate(treatment = rep(c("handled","injured", "saline", "low" ,"high"), each = 1, times=4)) 


repeatability_table <- rbind(activity_R, explore_R) %>%
  mutate(treatment = factor(treatment, levels=c("handled", "injured", "saline","low", "high"))) %>%
  mutate(behaviour=rep(c("Distance travelled", "Quadrant visitation latency"), each = 20)) 

#### FIGURES ####

##  Among individual

# activity
plot.activity.data = data.frame(Va.activity.F.Y.handled, Va.activity.F.Y.injured, Va.activity.F.Y.saline, Va.activity.F.Y.lo, Va.activity.F.Y.hi, 
                                Va.activity.M.Y.handled, Va.activity.M.Y.injured, Va.activity.M.Y.saline, Va.activity.M.Y.lo, Va.activity.M.Y.hi, 
                                Va.activity.F.O.handled, Va.activity.F.O.injured, Va.activity.F.O.saline, Va.activity.F.O.lo, Va.activity.F.O.hi,
                                Va.activity.M.O.handled, Va.activity.M.O.injured, Va.activity.M.O.saline, Va.activity.M.O.lo, Va.activity.M.O.hi,
                                Vw.activity.F.Y.handled, Vw.activity.F.Y.injured, Vw.activity.F.Y.saline, Vw.activity.F.Y.lo, Vw.activity.F.Y.hi, 
                                Vw.activity.M.Y.handled, Vw.activity.M.Y.injured, Vw.activity.M.Y.saline, Vw.activity.M.Y.lo, Vw.activity.M.Y.hi, 
                                Vw.activity.F.O.handled, Vw.activity.F.O.injured, Vw.activity.F.O.saline, Vw.activity.F.O.lo, Vw.activity.F.O.hi,
                                Vw.activity.M.O.handled, Vw.activity.M.O.injured, Vw.activity.M.O.saline, Vw.activity.M.O.lo, Vw.activity.M.O.hi) 

var.among.plot.activity <- plot.activity.data %>%
  select(Va.activity.F.Y.handled, Va.activity.F.Y.injured, Va.activity.F.Y.saline, Va.activity.F.Y.lo, Va.activity.F.Y.hi, 
         Va.activity.M.Y.handled, Va.activity.M.Y.injured, Va.activity.M.Y.saline, Va.activity.M.Y.lo, Va.activity.M.Y.hi,
         Va.activity.F.O.handled, Va.activity.F.O.injured, Va.activity.F.O.saline, Va.activity.F.O.lo, Va.activity.F.O.hi,
         Va.activity.M.O.handled, Va.activity.M.O.injured, Va.activity.M.O.saline, Va.activity.M.O.lo, Va.activity.M.O.hi) %>%
  gather(treatment, value,  factor_key=TRUE) %>%
  mutate(sex=if_else(str_detect(treatment, '.F.'), 'f', 'm')) %>%
  mutate(age=if_else(str_detect(treatment, '.O.'), 'old', 'young')) %>%
  mutate(treatment = str_remove(treatment, "Va.activity.F.Y.|Va.activity.F.O.|Va.activity.M.Y.|Va.activity.M.O.")) %>%
  mutate(treatment=recode(treatment,"handled"="handled"))

sex.lab <- c("Female", "Male")
names(sex.lab) <- c("f", "m")

age.lab <- c("Young", "Old")
names(age.lab) <- c("young", "old")

plot.Va.activity <- var.among.plot.activity %>%
  mutate(treatment = recode_factor(treatment,
                                   "hi" = "high",
                                   "lo" = "low")) %>%
  mutate(treatment = factor(treatment, levels=c("handled","injured", "saline", "low", "high"))) %>%
  mutate(age = factor(age, levels=c("young", "old"))) %>%
  ggplot(aes(x = treatment, y = value)) +
  #facet_grid(. ~ sex + age,labeller = labeller(age=age.lab,sex = sex.lab)) +
  facet_nested(.~sex+age,labeller = labeller(age=age.lab,sex = sex.lab)) +
  stat_dotsinterval(slab_fill="lightblue", slab_color="lightblue", point_interval = median_qi,
                    .width = 0.95, quantiles=100) +
  ylab("Distance travelled") +
  #xlab("Treatment") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12))
#ggsave(filename="figure_1.tiff", width=6.83, height=6.83, dpi=300,antialias="default")

## Exploration
plot.explore.data = data.frame(Va.explore.F.Y.handled, Va.explore.F.Y.injured, Va.explore.F.Y.saline, Va.explore.F.Y.lo, Va.explore.F.Y.hi, 
                                Va.explore.M.Y.handled, Va.explore.M.Y.injured, Va.explore.M.Y.saline, Va.explore.M.Y.lo, Va.explore.M.Y.hi, 
                                Va.explore.F.O.handled, Va.explore.F.O.injured, Va.explore.F.O.saline, Va.explore.F.O.lo, Va.explore.F.O.hi,
                                Va.explore.M.O.handled, Va.explore.M.O.injured, Va.explore.M.O.saline, Va.explore.M.O.lo, Va.explore.M.O.hi,
                                Vw.explore.F.Y.handled, Vw.explore.F.Y.injured, Vw.explore.F.Y.saline, Vw.explore.F.Y.lo, Vw.explore.F.Y.hi, 
                                Vw.explore.M.Y.handled, Vw.explore.M.Y.injured, Vw.explore.M.Y.saline, Vw.explore.M.Y.lo, Vw.explore.M.Y.hi, 
                                Vw.explore.F.O.handled, Vw.explore.F.O.injured, Vw.explore.F.O.saline, Vw.explore.F.O.lo, Vw.explore.F.O.hi,
                                Vw.explore.M.O.handled, Vw.explore.M.O.injured, Vw.explore.M.O.saline, Vw.explore.M.O.lo, Vw.explore.M.O.hi) 

var.among.plot.explore <- plot.explore.data %>%
  select(Va.explore.F.Y.handled, Va.explore.F.Y.injured, Va.explore.F.Y.saline, Va.explore.F.Y.lo, Va.explore.F.Y.hi, 
         Va.explore.M.Y.handled, Va.explore.M.Y.injured, Va.explore.M.Y.saline, Va.explore.M.Y.lo, Va.explore.M.Y.hi,
         Va.explore.F.O.handled, Va.explore.F.O.injured, Va.explore.F.O.saline, Va.explore.F.O.lo, Va.explore.F.O.hi,
         Va.explore.M.O.handled, Va.explore.M.O.injured, Va.explore.M.O.saline, Va.explore.M.O.lo, Va.explore.M.O.hi) %>%
  gather(treatment, value,  factor_key=TRUE) %>%
  mutate(sex=if_else(str_detect(treatment, '.F.'), 'f', 'm')) %>%
  mutate(age=if_else(str_detect(treatment, '.O.'), 'old', 'young')) %>%
  mutate(treatment = str_remove(treatment, "Va.explore.F.Y.|Va.explore.F.O.|Va.explore.M.Y.|Va.explore.M.O.")) %>%
  mutate(treatment=recode(treatment,"handled"="handled"))

sex.lab <- c("Female", "Male")
names(sex.lab) <- c("f", "m")

age.lab <- c("Young", "Old")
names(age.lab) <- c("young", "old")

plot.Va.explore <- var.among.plot.explore %>%
  mutate(treatment = recode_factor(treatment,
                                   "hi" = "high",
                                   "lo" = "low")) %>%
  mutate(treatment = factor(treatment, levels=c("handled","injured", "saline", "low", "high"))) %>%
  mutate(age = factor(age, levels=c("young", "old"))) %>%
  ggplot(aes(x = treatment, y = value)) +
  facet_grid(. ~ sex + age,labeller = labeller(age=age.lab,sex = sex.lab)) +
  stat_dotsinterval(slab_fill="orange", slab_color="orange", point_interval = median_qi,
                    .width = 0.95, quantiles=100) +
  ylab("Quadrant visitation latency") +
  xlab("Treatment") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text.x = element_text(size=12, angle=45, hjust=1),
        axis.text.y = element_text(size=12))
#ggsave(filename="figure_1.tiff", width=6.83, height=6.83, dpi=300,antialias="default")

figure_1<-plot_grid(plot.Va.activity, plot.Va.explore, ncol=1, nrow=2,
                    rel_heights = c(1,1, 1.2),hjust=-3.5,vjust=c(5.5,2,2),labels = c('(a)', '(b)'))
ggsave(filename="figure_2.tiff", width=10.83, height=10.83, dpi=300,antialias="default")    


### Within individual

# activity
Vwr.within.plot.long <- plot.activity.data %>%
  select(Vw.activity.F.Y.handled, Vw.activity.F.Y.injured, Vw.activity.F.Y.saline, Vw.activity.F.Y.lo, Vw.activity.F.Y.hi, 
         Vw.activity.M.Y.handled, Vw.activity.M.Y.injured, Vw.activity.M.Y.saline, Vw.activity.M.Y.lo, Vw.activity.M.Y.hi,
         Vw.activity.F.O.handled, Vw.activity.F.O.injured, Vw.activity.F.O.saline, Vw.activity.F.O.lo, Vw.activity.F.O.hi,
         Vw.activity.M.O.handled, Vw.activity.M.O.injured, Vw.activity.M.O.saline, Vw.activity.M.O.lo, Vw.activity.M.O.hi) %>%
  gather(treatment, value,  factor_key=TRUE) %>%
  mutate(sex=if_else(str_detect(treatment, '.F.'), 'f', 'm')) %>%
  mutate(age=if_else(str_detect(treatment, '.O.'), 'old', 'young')) %>%
  mutate(treatment = str_remove(treatment, "Vw.activity.F.Y.|Vw.activity.F.O.|Vw.activity.M.Y.|Vw.activity.M.O.")) %>%
  mutate(treatment=recode(treatment,"handled"="handled"))

sex.lab <- c("Female", "Male")
names(sex.lab) <- c("f", "m")

age.lab <- c("Young", "Old")
names(age.lab) <- c("young", "old")

plot.Vw.activity <- Vwr.within.plot.long %>%
  mutate(treatment = recode_factor(treatment,
                                   "hi" = "high",
                                   "lo" = "low")) %>%
  mutate(treatment = factor(treatment, levels=c("handled","injured", "saline", "low", "high"))) %>%
  mutate(age = factor(age, levels=c("young", "old"))) %>%
  ggplot(aes(x = treatment, y = value)) +
  facet_nested(.~sex+age,labeller = labeller(age=age.lab,sex = sex.lab)) +
  stat_dotsinterval(slab_fill="lightblue", slab_color="lightblue", point_interval = median_qi,
                    .width = 0.95, quantiles=100) +
  ylab("Distance travelled") +
  #xlab("Treatment") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12))
#ggsave(filename="figure_1.tiff", width=6.83, height=6.83, dpi=300,antialias="default")

## Exploration
Vwr.within.plot.long <- plot.explore.data %>%
  select(Vw.explore.F.Y.handled, Vw.explore.F.Y.injured, Vw.explore.F.Y.saline, Vw.explore.F.Y.lo, Vw.explore.F.Y.hi, 
         Vw.explore.M.Y.handled, Vw.explore.M.Y.injured, Vw.explore.M.Y.saline, Vw.explore.M.Y.lo, Vw.explore.M.Y.hi,
         Vw.explore.F.O.handled, Vw.explore.F.O.injured, Vw.explore.F.O.saline, Vw.explore.F.O.lo, Vw.explore.F.O.hi,
         Vw.explore.M.O.handled, Vw.explore.M.O.injured, Vw.explore.M.O.saline, Vw.explore.M.O.lo, Vw.explore.M.O.hi) %>%
  gather(treatment, value,  factor_key=TRUE) %>%
  mutate(sex=if_else(str_detect(treatment, '.F.'), 'f', 'm')) %>%
  mutate(age=if_else(str_detect(treatment, '.O.'), 'old', 'young')) %>%
  mutate(treatment = str_remove(treatment, "Vw.explore.F.Y.|Vw.explore.F.O.|Vw.explore.M.Y.|Vw.explore.M.O.")) %>%
  mutate(treatment=recode(treatment,"handled"="handled"))

sex.lab <- c("Female", "Male")
names(sex.lab) <- c("f", "m")

age.lab <- c("Young", "Old")
names(age.lab) <- c("young", "old")

plot.Vw.explore <- Vwr.within.plot.long %>%
  mutate(treatment = recode_factor(treatment,
                                   "hi" = "high",
                                   "lo" = "low")) %>%
  mutate(treatment = factor(treatment, levels=c("handled","injured", "saline", "low", "high"))) %>%
  mutate(age = factor(age, levels=c("young", "old"))) %>%
  ggplot(aes(x = treatment, y = value)) +
  facet_grid(. ~ sex + age,labeller = labeller(age=age.lab,sex = sex.lab)) +
  stat_dotsinterval(slab_fill="orange", slab_color="orange", point_interval = median_qi,
                    .width = 0.95, quantiles=100) +
  ylab("Quadrant visitation latency") +
  xlab("Treatment") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text.x = element_text(size=12, angle=45, hjust=1),
        axis.text.y = element_text(size=12))
#ggsave(filename="figure_1.tiff", width=6.83, height=6.83, dpi=300,antialias="default")

figure_2<-plot_grid(plot.Vw.activity, plot.Vw.explore, ncol=1, nrow=2,
                    rel_heights = c(1,1, 1.2),hjust=-3.5,vjust=c(5.5,2,2),labels = c('(a)', '(b)'))
ggsave(filename="figure_3.tiff", width=10.83, height=10.83, dpi=300,antialias="default")  

# REPEATABILITY
sex.lab <- c("Female", "Male")
names(sex.lab) <- c("female", "male")

repeat.plot<- ggplot(repeatability_table,aes(y = treatment, x = .value, xmin = .lower, xmax = .upper, shape=age)) +
  geom_vline(aes(xintercept=0),colour="grey", linetype='dotted') +
  geom_pointinterval(position=position_dodge(0.4), colour="slateblue2", fatten_point = 4) +
  facet_nested(sex~ behaviour,labeller = labeller(sex = sex.lab)) +
  theme(strip.text = element_text(
    size = 5, color = "dark green")) +
  xlab(expression(~italic(R)~" ± 95% CI")) +
  ylab("Treatment") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  strip.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        axis.title.x = element_text(size=18,face="bold"),
        axis.title.y = element_text(size=18,face="bold"),
        axis.text.x = element_text(size=12, angle=45, vjust=0.5),
        axis.text.y = element_text(size=12)) +
  scale_shape_manual(values=c(15,16,17),
                     name="Age",
                     breaks=c("young", "old"),
                     labels=c("young", "old"))
#ggsave(filename="figure_4.png", width=16.83, height=16.83, dpi=800,antialias="default")

 ##### ASSAY ANALYSIS

model.brms.activity.test <- bf(scale(distance.t)~ sex*treatment_pooled*age + assay + (0+assay:sex:age||gr(id, by = treatment_pooled)), sigma ~ 0+assay:age:sex:treatment_pooled, family = gaussian) 
fit_model.brms.activity.test <- brm(model.brms.activity.test, data = merged_data_four, iter=4000, seed=seed, thin=2,cores= 4)
summary(fit_model.brms.activity.test)

model.brms.activity.test <- bf(scale(distance.t)~ sex+treatment_pooled+age + assay + (0+assay||gr(id)), sigma ~ 0+assay:id, family = gaussian) 
fit_model.brms.activity.test <- brm(model.brms.activity.test, data = merged_data_four, iter=6000, seed=seed, thin=2,cores= 4)
summary(fit_model.brms.activity.test)


model.brms.activity.test2 <- bf(scale(distance.t)~ sex + treatment_pooled + age + assay_centred + (0+assay_centred||id), family = gaussian) 
fit_model.brms.activity.test2 <- brm(model.brms.activity.test2, data = merged_data_four, iter=6000, warmup=1000,seed=seed, thin=2,cores= 4)
summary(fit_model.brms.activity.test2)
fixef(fit_model.brms.activity.test2)

get_variables(fit_model.brms.activity.test2)
m <- fit_model.brms.activity.test2 %>%
  spread_draws(r_id[asssay,term]) %>%
  mean_qi(.width = 0.95)

BLUPS1 = fit_model.brms.activity.test2 %>% spread_draws(r_id__assay_centred0[ID,], r_id__assay_centred1[ID,]) %>% 
  mean_qi(.width = 0.95)

parnames(fit_model.brms.activity.test2)

ranef(fit_model.brms.activity.test2)
summary(fit_model.brms.activity.test)

assay2<-ranef(fit_model.brms.activity.test2)[[1]][, , "assay_centred1"]
assay1<-ranef(fit_model.brms.activity.test2)[[1]][, , "assay_centred0"] 

assay1x<- as.data.frame(assay1) %>%
  mutate(assay=rep(c("1"))) %>%
  summarise(mean=mean(Estimate))
assay1x <- data.frame(id = row.names(assay1x), assay1x)
assay2x<- as.data.frame(assay2) %>%
  mutate(assay=rep(c("2"))) %>%
  summarise(mean=mean(Estimate))
assay2x <- data.frame(id = row.names(assay2x), assay2x)


tablex <- bind_rows(assay1x,assay2x) %>%
  arrange(id)

ggplot(tablex, aes(x=assay, y=Estimate, group=id)) +
  geom_line(alpha=0.1)

ggplot(merged_data_four, aes(x=assay, y=distance, group=id)) +
  geom_line(alpha=0.1)

ggplot(m, aes(x=term, y=r_id, group=asssay)) +
  geom_line(alpha=0.1)

summy <-merged_data_four %>%
  group_by(assay) %>%
  summarise(mean=mean(distance))

Va.activity.all1 <- as_draws_df(fit_model.brms.activity.test2)$"sd_id__assay1"^2
Va.activity.all2 <- as_draws_df(fit_model.brms.activity.test2)$"sd_id__assay2"^2
post.data.all.assay = data.frame(Va.activity.all1,Va.activity.all2) 

all_activity.Va = post.data.all.assay %>%
  dplyr::select(starts_with("Va.activity.")) %>%
  pivot_longer(cols = starts_with("Va.activity."),
               names_to = 'treat',
               names_prefix = "Va.activity.",
               values_to ="Estimate")

Va_all_activity = all_activity.Va %>%
  dplyr::group_by(treat) %>%
  dplyr::summarise(Va_all_activity = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[2], 3))%>%
  as.data.frame()

Vw.activity.all1 <- exp(as_draws_df(fit_model.brms.activity.test2)$"sigma")^2
Vw.activity.all2 <- exp(as_draws_df(fit_model.brms.activity.test2)$"b_sigma_assay2")^2
post.data.avw.assay = data.frame(Vw.activity.all1) 
all_activity.Vw = post.data.avw.assay %>%
  dplyr::select(starts_with("Vw.activity.")) %>%
  pivot_longer(cols = starts_with("Vw.activity."),
               names_to = 'treat',
               names_prefix = "Vw.activity.",
               values_to ="Estimate")

Vw_all_activity = all_activity.Vw %>%
  dplyr::group_by(treat)%>%
  dplyr::summarise(Vw_all_activity = 
                     round(mean(Estimate), 3),
                   lowerCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[1], 3),
                   upperCI = 
                     round(rethinking::HPDI(Estimate, prob = 0.95)[2], 3))%>%
  as.data.frame()

model.brms.activity.test3 <- bf(scale(distance.t)~ sex + treatment_pooled + age + assay + (1+assay|id), family = gaussian)# gives correlation between intercepts and slope
fit_model.brms.activity.test3 <- brm(model.brms.activity.test3, data = merged_data_four, iter=6000, warmup=1000,seed=seed, thin=2,cores= 4)
summary(fit_model.brms.activity.test3)

