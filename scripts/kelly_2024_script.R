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

# Load Functions

# Load in the functions that will be needed to process output from the models to derive the various estimates needed for calculating the necessary parameters:
  
# {r funcs, include=TRUE, message=FALSE, warning=FALSE, class.source='klippy'}

# Population-level intercept   
func_Bp <- function(B0, B1){
  Bp <- (2*B0 + B1)/2 # assumes 50/50 ratio for B1 (in our example = sex)
  return(Bp)
}

# Total variance, equation 40. Function for summing dataframe of variance components 
func_sum_var <- function(vars){
  sum_var <- sapply(1:nrow(vars), function(x) sum(vars[x,], na.rm = T))
  return(sum_var)
}


# Variance for fixed effects
func_var_fixed <- function(B1, B2, data){
  # model matrix for fixed effects
  X <- model.matrix(~ 1 + sex.centred + sex.centred:observation.n , data = data) 
  X1 <- X[,2]
  X2 <- X[,3]
  var_fixed <- sapply(1:length(B1), function(x) var(X1*B1[x] + X2*B2[x]))
  return(var_fixed)
}		

# within-individual standard deviation
func_var_within <- function(Bpv, sigma_v0, var_fixed_v){
  
  # summed variance components
  var_within_exp <- func_sum_var(vars = data.frame(sigma_v0^2, var_fixed_v))
  # conversion back from log scale
  var_within <- func_ln_convert(mu_ln = Bpv, sigma_ln = sqrt(var_within_exp))$mu_raw
  
  return(var_within)
}

# Repeatability for the dispersion model for DHGLM    
func_Rp_var <- function(Bpv, sigma_v0, var_fixed_v, var_p){
  
  # sum of IDv0 and fixedV = total variance in residual variance on log-normal scale:
  var_residvar_exp <- func_sum_var(vars = data.frame(sigma_v0^2, var_fixed_v))
  # converting back to same scale as mean model: 
  sigma_residvar <- func_ln_convert(mu_ln = Bpv, sigma_ln = sqrt(var_residvar_exp))$sigma_raw
  
  # variance in phenotypic variance
  total_var_var <- 2*var_p^2 + 3*sigma_residvar^2
  
  # getting variance in individual component through the preservation of the proportionality (i.e. ratio method)
  ratio <- sigma_v0^2/(var_residvar_exp)
  var_ID <- sigma_residvar^2*ratio
  
  Rp_var <- var_ID/total_var_var
  return(Rp_var)
}

# Convert from ln to raw scale
func_ln_convert <- function(mu_ln, sigma_ln){
  
  mu_raw <- exp(mu_ln + sigma_ln^2/2)
  var_raw <- (exp(sigma_ln^2)-1)*exp(2*mu_ln + sigma_ln^2)
  
  x <- data.frame(mu_raw, sigma_raw = sqrt(var_raw))
  return(x)
}	 


# Function takes random slopes and intercepts and calculates the between individual correlation between them and returns this correlation for each row (i.e., posterior sampling iteration). It returns a posterior distribution of the correlation. 

cor_calc <- function(slopes, intercepts){
  cors <- c()
  for(i in 1:dim(intercepts)[1]){
    cors <- c(cors, cor(as.numeric(slopes[i,]),as.numeric(intercepts[i,])))
  }
  return(cors)
}





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

# Add mating success
morphology_long$ms <-ifelse(morphology_long$partner_id== "", 0,1)

mating_success_table <- mobility_data %>%
  mutate_all(funs(replace(., .=="", NA))) %>%
  group_by(ID, sex) %>%
  summarise(mates=sum(!is.na(partner_id)), n=n(), ms=mates/n) #number of mates
## ---- end


#### PREDICTABILITY ####

#mobility_pred <- bf(distance.tr ~ sex.centred + observation.n + sex.centred:observation.n + (observation.n|ID), sigma ~ sex , family = gaussian)
mobility_pred <- bf(distance.tr ~ sex.centred + observation.n + sex.centred:observation.n + (observation.n|a|ID), sigma ~ sex + (1|a|ID), family = gaussian)

fit.model.brms.pred <- brm(mobility_pred, data = morphology_long, save_pars = save_pars(all = TRUE), 
                               warmup=6000, iter=10000, seed=12345, thin=4, chains=4, cores= 4, file = 'data/processed/fit.model.brms.pred')
summary(fit.model.brms.pred, prob = 0.95)

plot(fit.model.brms.pred, ask = F)
brms::pp_check(fit.model.brms.pred, resp = "distance.z", ndraws = 50)
conditional_effects(fit.model.brms.pred, prob = 0.95)

#plot

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

ggplot(mobility.residual, aes(x = Estimate, fill = Sex)) +
  geom_density(alpha=0.6)+ 
  geom_vline(data = mobility.mean, aes(xintercept=grp.mean, color = Sex),
             linetype="dashed", linewidth=1) +
  xlab('Residual within-individual variaiton in mobility\n (i.e. predictability)') +
  ylab('Density') +
  theme_bw()

# males less predictable than females
# correlate rIIV with mating success = less predictable males have higher mating success

## ---- data_analysis ----
fit.model.brms.pred = readRDS(file = "data/processed/fit.model.brms.pred.rds")

qmd.values <- fit.model.brms.pred %>%
  spread_draws(b_sigma_Intercept,b_sigma_sexm,sd_ID__sigma_Intercept) %>%
  median_qi(b_sigma_Intercept,b_sigma_sexm, sd_ID__sigma_Intercept)

get_variables(fit.model.brms.pred)

#### Behavioral predictability (individual rIIV's) ####

# Intercept Distance
bk.tr.dist <- exp(fixef(fit.model.brms.pred, pars = "sigma_Intercept")[1]) * sd(morphology_long$distance, na.rm = T)
# 7.66 m

#### Coefficient of variation in predictability (CVP) ####

log.norm.res.Dist <- exp(posterior_samples(fit.model.brms.pred)$"sd_ID__sigma_Intercept"^2)
CVP.long.Dist <- sqrt(log.norm.res.Dist - 1)
mean(CVP.long.Dist);HPDinterval(as.mcmc(CVP.long.Dist),0.95)

## ---- end ----


COR.PERS.PRED <- 
  as_draws_df(fit.model.brms.pred, 
                    pars = c("cor_ID__Intercept__sigma_Intercept")) %>%
  gather() %>%
  separate(key,
           c(NA,"Scale",NA,NA,NA,NA,NA,"Trait",NA),
           sep = "([\\_\\__\\_\\__\\_\\_\\,])", fill = "right")

# Slope travel distance
cov.Trav <-
  posterior_samples(fit.model.brms.pred)[,11] *
  sqrt((posterior_samples(fit.model.brms.pred)[,7])^2) *
  sqrt((posterior_samples(fit.model.brms.pred)[,9])^2)
var.Trav <- (posterior_samples(fit.model.brms.pred)[,7])^2
slope_Trav <- cov.Trav / var.Trav

mean_distance <- as_draws_df(fit.model.brms.pred, pars = "^r_ID")[1:47] %>%
  tidyr::gather(ID, value,
                "r_ID[W01,Intercept]" : 
                  "r_ID[W75,Intercept]") %>%
  select(ID,value) %>%
  separate(ID,
           c(NA,NA,"ID",NA),
           sep = "([\\_\\[\\,])", fill = "right") %>%
  group_by(ID) %>%
  mutate(Mean = mean(value),
         Up = mean(value) + 1.96 * sd(value),
         Lo = mean(value) - 1.96 * sd(value))%>%
  select (-value)%>%
  filter(!duplicated(ID))

names(mean_distance)<-c("ID","MeanDist","UpDist","LoDist")

sigma_distance <- as_draws_df(fit.model.brms.pred, pars = "^r_ID__sigma") %>%
  tidyr::gather(ID, value,
                "r_ID__sigma[W01,Intercept]" : 
                  "r_ID__sigma[W75,Intercept]")%>%
  select(ID,value) %>%
  separate(ID,
           c(NA,"ID",NA),
           sep = "([\\[\\,])", fill = "right") %>%
  group_by(ID) %>%
  mutate(Mean = mean(value),
         Up = mean(value) + 1.96 * sd(value),
         Lo = mean(value) - 1.96 * sd(value))%>%
  select (-value) %>%
  filter(!duplicated(ID))

names(sigma_distance)<-c("ID","predict_MeanDist","predict_UpDist","predict_LoDist")

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
             size=1.4, color = "#F8766D", alpha = 0.9, stroke = 0) +
  
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
  #theme_classic(base_size = 10, base_family = "Georgia")+
  #mytheme +
  #myscale +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  annotate("text",x = 1.3, y = 0.65, label = expression(paste(italic("r")*" = -0.10")), size = 4)

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
  summarize(COR = stats::cor.test(MeanDist, ms)$estimate,
            pval = stats::cor.test(MeanDist, ms)$p.value,
            df = stats::cor.test(MeanDist, ms)$parameter
  ) %>%
  ungroup()

### predictability vs mating success

mating_predictability_data <- right_join(sigma_distance, mating_success_table, by ='ID') %>%
  filter(!is.na(predict_MeanDist))

mating_predict_corr <- mating_predictability_data %>%
  group_by(sex) %>%
  summarize(COR = stats::cor.test(predict_MeanDist, ms)$estimate,
            pval = stats::cor.test(predict_MeanDist, ms)$p.value,
            df = stats::cor.test(predict_MeanDist, ms)$parameter
  ) %>%
  ungroup()

exp(0.21)^2

## ---- end





install.packages("gitcreds")
library(gitcreds)
gitcreds_set()

##as_draws()#################################
# Mean and Dispersion Parameters
#################################

# (Step 1): Extract the population coefficients for mean model
B0_m <- brms::as_draws_matrix(fit.model.brms.pred, variable = "b_Intercept")	 
B1_m <- brms::as_draws_matrix(fit.model.brms.pred, variable = "b_sex.centred")
B2_m <- brms::as_draws_matrix(fit.model.brms.pred, variable = "b_sex.centred:observation.n")

mean(B0_m)

mean(B1_m)

mean(B2_m)

# (Step 2): Extract the population coefficients for dispersion model	 
B0_sd_exp <- brms::as_draws_matrix(fit.model.brms.pred, variable = "b_sigma_Intercept")	
B1_sd_exp <- brms::as_draws_matrix(fit.model.brms.pred, variable = "b_sigma_sexm")
#B2_sd_exp <- brms::as_draws_matrix(fit.model.brms.pred, variable = "b_sigma_age.Z")

mean(2*B0_sd_exp)

mean(2*B1_sd_exp)

#mean(2*B2_sd_exp)

# (Step 3): Extract the random ID intercept variance components for random intercept. Note that we model ln(sd) so we need to convert to variance scale by multipling by 2.
sigma_ID_m0 <- brms::as_draws_matrix(fit.model.brms.pred, variable = "sd_ID__Intercept")

mean(sigma_ID_m0)

sigma_ID_sd0 <- brms::as_draws_matrix(fit.model.brms.pred, variable = "sd_ID__observation.n")
sigma_ID_v0 <- 2*sigma_ID_sd0

mean(sigma_ID_v0)

# (Step 4): **New step as we now have a sd for the random slopes.**
sigma_ID_m2 <- brms::posterior_samples(brms_personality_plasticity_dhglm_aggression, pars = "sd_id__age.Z")[,1]

mean(sigma_ID_m2)

# (Step 5): Calculate the fixed effect variance for the mean model (var_fixed_m) and the dispersion model (B_pv). Note here that the calculations are the same because we assume intercept at an average age.
B_pm <- func_Bp(B0 = B0_m, B1 = B1_m)
B_pv <- 2*func_Bp(B0 = B0_sd_exp, B1 = B1_sd_exp)

mean(B_pm)

mean(B_pv)

# (Step 6): Calculate the population intercept for the mean model (var_fixed_m) and dispersion model (var_fixed_v). Note that a variance on the ln(sd) scale can be converted to a ln(var) by multipling by 4. See supplemental text.
var_fixed_m <- func_var_fixed(B1 = B1_m, B2 = B2_m, data = morphology_long)
var_fixed_v <- 4*func_var_fixed(B1 = B1_sd_exp, data = morphology_long)

mean(var_fixed_m)
mean(stanproc_personality_plasticity_dhglm_aggression$var1_fixed_m)

mean(var_fixed_v)
mean(stanproc_personality_plasticity_dhglm_aggression$var1_fixed_v)



###################################
# Mean Model - Repeatability and CV
###################################

# (Step 1): Calculate average within-individual / residual variance
sigma_w <- func_var_within(Bpv = B_pv, 
                           sigma_v0 = sigma_ID_v0, 
                           var_fixed_v = var_fixed_v)

mean(sigma_w)
mean(stanproc_personality_dhglm_aggression$sigma1_w)

# (Step 2): Calculate total phenotypic variance 
var_p <- func_sum_var(vars = data.frame(sigma_ID_m0^2, var_fixed_m, sigma_w^2))

mean(var_p)
mean(stanproc_personality_dhglm_aggression$var1_p)

# (Step 3): Calculate repeatability
Rp_mu <- sigma_ID_m0^2/var_p

mean(Rp_mu)
mean(stanproc_personality_dhglm_aggression$Rp1_mu)

# (Step 4): Convert back from z-scale to calculate CVm
mu_unscaled <- mean(dat$stim.dist_0_5cm_dur_Aggression)
sigma_unscaled <-   sd(dat$stim.dist_0_5cm_dur_Aggression)

Mu <- B_pm*sigma_unscaled + mu_unscaled
Sigma <- sqrt(var_p)*sigma_unscaled
indiv_sd <- sigma_ID_m0*sigma_unscaled

# (Step 5): Calculate CVm for mean model
CV_mu <- indiv_sd/Mu

mean(CV_mu)
mean(stanproc_personality_dhglm_aggression$CV1_mu)

###################################
# Dispersion Model - Repeatability and CV
###################################

# (Step 1): Calculate repeatability for the dispersion model
Rp_var <- func_Rp_var(Bpv = B_pv, 
                      sigma_v0 = sigma_ID_v0, 
                      var_fixed_v = var_fixed_v, 
                      var_p = var_p)

mean(Rp_var)
mean(stanproc_personality_dhglm_aggression$Rp1_var)

# (Step 2): Sum of IDv0 and fixedV = total variance in residual variance on log-normal scale
var_residvar_exp <- func_sum_var(vars = data.frame(sigma_ID_v0^2, var_fixed_v))


# (Step 3): Converting back to same scale as mean model  
sigma_residvar <- func_ln_convert(mu_ln = B_pv, 
                                  sigma_ln = sqrt(var_residvar_exp))$sigma_raw

# (Step 4): Getting variance in individual component through the preservation of the proportionality (i.e. ratio method)
ratio <- sigma_ID_v0^2/(var_residvar_exp)
var_ID <- sigma_residvar^2*ratio

# (Step 5): Calculate coefficient of variation 
CV_var <- sqrt(var_ID)/sigma_w^2   

mean(CV_var)
mean(stanproc_personality_dhglm_aggression$CV1_var)



