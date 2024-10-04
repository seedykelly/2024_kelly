
# R scripts to reproduce results for
# "Biologging reveals individual variation in behavioral predictability in the wild."
# by Anne G. Hertel, Raphael Royaute, Andreas Zedrosser, Thomas Mueller

# enable parallel processing in R 4.0 and RStudio
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() == "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}

library(data.table);library(brms);library(ggplot2);library(dplyr)
library(tidyr);library(coda);library(RVAideMemoire);library(tidybayes)
library(parallel);library(stringr);library(bayestestR)

ncores = parallel::detectCores()

PERS <- readRDS("data/raw/DATA_R1.rds")
PERS$Bearyear <- factor(PERS$Bearyear)

length(levels(droplevels(PERS$Bearyear[PERS$subDiurn == 1])))
length(levels(droplevels(PERS$Bearyear[PERS$subDist == 1])))
length(levels(droplevels(PERS$Bearyear[PERS$subDist == 1 & PERS$subDiurn == 1])))

# scale response variables
s.DailyDistance <- scale(PERS$DailyDistance)
s.DailyDiurnality <- scale(PERS$DailyDiurnality)

table(PERS$Status)

# 1. MODEL: Multivariate DHGLM diurnality & distance

# 2. Model fit, importance and interpretation of fixed effects in mean model

# 3. Repeatability, coefficient of variation

# 4. Coefficient of variation in predictability (CVP) 

# 5. Correlation mean x rIIV

# 6. Behavioral syndrome (Correlation mean x mean)

# 7. Predictability syndrome (Correlation rIIV x rIIV)

# 8. Non-sensible correlations

# 9. Behavioral types (individual means)

# 10. Behavioral predictability (individual rIIV's)

#----------------------------------------------------------------------

#### 1. MODEL: Multivariate DHGLM diurnality & distance ####

#----------------------------------------------------------------------

DiurnDHLGM = bf(scale(DailyDiurnality) | subset(subDiurn)  ~ 
                  day.scaled * Status + Young:age.scaled +
                  (day.scaled|a|Bearyear) + (day.scaled|b|BearID),
                sigma ~ Status + age.scaled + (1|a|Bearyear) + (1|b|BearID))

DistDHGLM = bf(scale(DailyDistance) | subset(subDist) ~
                 day.scaled * Status + age.scaled +
                 (day.scaled|a|Bearyear) + (day.scaled|b|BearID),
               sigma ~ Status + age.scaled + (1|a|Bearyear) + (1|b|BearID)) 

get_prior(DistDHGLM,data = PERS)

priors <- c(
  set_prior("normal(0, 1)", class = "b", resp = c("sxcaleDailyDistance","scaleDailyDiurnality")),
  set_prior("normal(0, 1)", class = "Intercept", resp = c("scaleDailyDistance","scaleDailyDiurnality")),
  set_prior("normal(0, 1)", class = "sd", resp = c("scaleDailyDistance","scaleDailyDiurnality")),
  set_prior("normal(0, 1)", class = "b", dpar = "sigma",resp = c("scaleDailyDistance","scaleDailyDiurnality")),
  set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma",resp = c("scaleDailyDistance","scaleDailyDiurnality")),
  set_prior("normal(0, 1)", class = "sd", dpar = "sigma", resp = c("scaleDailyDistance","scaleDailyDiurnality")),
  set_prior("lkj(2)", class = "cor"))

CompDHGLM <- brm(DiurnDHLGM + DistDHGLM + set_rescor(FALSE),
                 data = PERS,
                 iter  = 10000, warmup = 6000, thin = 4,
                 chains = 4, cores = 4, seed = 12345)

#saveRDS(CompDHGLM,"CompositeModel.rds")
CompDHGLM <- readRDS("CompositeModel.rds")

#stancode(CompDHGLM)

#----------------------------------------------------------------------

#### 2. Model fit, importance and interpretation of fixed effects in mean model ####

#----------------------------------------------------------------------

# Check mixing of chains
plot(CompDHGLM)

# Model validation - Diurnality
pp_check(CompDHGLM, resp = "scaleDailyDiurnality")+ 
  theme_bw(base_size = 20)
pp_check(CompDHGLM, 
         resp = "scaleDailyDiurnality",
         nsamples = 1e3, 
         type = "stat_2d") + 
  theme_bw(base_size = 20)

# Model validation - Distance
pp_check(CompDHGLM, resp = "scaleDailyDistance")
pp_check(CompDHGLM, 
         resp = "scaleDailyDistance",
         nsamples = 1e3, 
         type = "stat_2d") + 
  theme_bw(base_size = 20)

# Effects plot - mean model
head(get_variables(CompDHGLM))
DIUR.EFF <- posterior_samples(CompDHGLM, pars = "^b_scaleDailyDiurnality")%>% 
  gather() %>% 
  mutate(predictor = str_remove(key, "b_scaleDailyDiurnality_"),
         Trait = "Diurnality")%>%
  filter(predictor != "Intercept")%>% 
  droplevels()

DIST.EFF <- posterior_samples(CompDHGLM, pars = "^b_scaleDailyDistance")%>% 
  gather() %>% 
  mutate(predictor = str_remove(key, "b_scaleDailyDistance_"),
         Trait = "Daily Distance")%>%
  filter(predictor != "Intercept")%>% 
  droplevels()

EFF <- rbind(DIUR.EFF,DIST.EFF)
EFF$predictor <- factor(EFF$predictor)

l <- levels(EFF$predictor)

EFF$PRED <- NA

EFF[EFF$predictor == l[1],"PRED"] <- "Age"
EFF[EFF$predictor == l[2],"PRED"] <- "Ordinal day"
EFF[EFF$predictor == l[3],"PRED"] <- "Ordinal day : Cubs"
EFF[EFF$predictor == l[4],"PRED"] <- "Ordinal day : Yearling offspring"
EFF[EFF$predictor == l[5],"PRED"] <- "Cubs"
EFF[EFF$predictor == l[6],"PRED"] <- "Yearling offspring"
EFF[EFF$predictor == l[7],"PRED"] <- "Age : Adult(>4 yrs)"
EFF[EFF$predictor == l[8],"PRED"] <- "Age : Young(3-4 yrs)"

EFF$PRED <- factor(EFF$PRED, levels = c("Age : Young(3-4 yrs)" ,
                                                       "Age : Adult(>4 yrs)" ,
                                                       "Age",
                                                       "Ordinal day : Yearling offspring",
                                                       "Ordinal day : Cubs",
                                                       "Yearling offspring",
                                                       "Cubs",
                                                       "Ordinal day"))

ggplot(EFF, aes(x = value, y = PRED)) +
  geom_halfeyeh()+
  geom_vline(xintercept = 0)+
  labs(x = "Posterior distribution",y=" ")+
  facet_wrap( ~ Trait, ncol=2)+ 
  theme_classic(base_size = 10)

# Evaluate fixed effect backtransformed to original scale (for 
# ecological interpretation)
# Diurnality
    # population level mean diurnality for solitary females
    mean(fixef(CompDHGLM, pars = c("scaleDailyDiurnality_Intercept"), 
          resp = "scaleDailyDiurnality", summary = F)) * attr(s.DailyDiurnality, 'scaled:scale') + 
      attr(s.DailyDiurnality, 'scaled:center')
    
    # population level mean diurnality for females with cubs
    (mean(fixef(CompDHGLM, pars = c("scaleDailyDiurnality_Intercept"), 
               resp = "scaleDailyDiurnality", summary = F)) +
     mean(fixef(CompDHGLM, pars = c("scaleDailyDiurnality_StatusWithCOY"), 
                   resp = "scaleDailyDiurnality", summary = F))) * attr(s.DailyDiurnality, 'scaled:scale') + 
      attr(s.DailyDiurnality, 'scaled:center')
    
    # population level mean diurnality for fem with offspring (>= 1 yrs old)
    (mean(fixef(CompDHGLM, pars = c("scaleDailyDiurnality_Intercept"), 
                resp = "scaleDailyDiurnality", summary = F)) +
        mean(fixef(CompDHGLM, pars = c("scaleDailyDiurnality_StatusWithOFFSPRING"), 
                   resp = "scaleDailyDiurnality", summary = F))) * attr(s.DailyDiurnality, 'scaled:scale') + 
      attr(s.DailyDiurnality, 'scaled:center')

# Distance
    # population level mean daily distance for solitary females
    mean(fixef(CompDHGLM, pars = c("scaleDailyDistance_Intercept"), 
               resp = "scaleDailyDistance", summary = F)) * attr(s.DailyDistance, 'scaled:scale') + 
      attr(s.DailyDistance, 'scaled:center')
    
    # population level mean daily distance for females with cubs
    (mean(fixef(CompDHGLM, pars = c("scaleDailyDistance_Intercept"), 
                resp = "scaleDailyDistance", summary = F)) +
        mean(fixef(CompDHGLM, pars = c("scaleDailyDistance_StatusWithCOY"), 
                   resp = "scaleDailyDistance", summary = F))) * attr(s.DailyDistance, 'scaled:scale') + 
      attr(s.DailyDistance, 'scaled:center')
    
    # population level mean daily distance for fem with offspring (>= 1 yrs old)
    (mean(fixef(CompDHGLM, pars = c("scaleDailyDistance_Intercept"), 
                resp = "scaleDailyDistance", summary = F)) +
        mean(fixef(CompDHGLM, pars = c("scaleDailyDistance_StatusWithOFFSPRING"), 
                   resp = "scaleDailyDistance", summary = F))) * attr(s.DailyDistance, 'scaled:scale') + 
      attr(s.DailyDistance, 'scaled:center')

#----------------------------------------------------------------------

#### 3. Repeatability, coefficient of variation ####

#----------------------------------------------------------------------

## shortterm repeatability ("Bearyear")
## longterm repeatability  ("BearID")
## Intercept repeatability  ("BearID")
## shortterm CVP 
## longterm CVP 

    names(posterior_samples(CompDHGLM))[1:35]
    
# Variance components DIURNALITY 
var.BearID.Diurn <- posterior_samples(CompDHGLM)$"sd_BearID__scaleDailyDiurnality_Intercept"^2
var.Bearyear.Diurn <- posterior_samples(CompDHGLM)$"sd_Bearyear__scaleDailyDiurnality_Intercept"^2
var.res.Diurn <- exp(posterior_samples(CompDHGLM)$"b_sigma_scaleDailyDiurnality_Intercept")^2

# Variance components DISTANCE 
var.BearID.Dist <- posterior_samples(CompDHGLM)$"sd_BearID__scaleDailyDistance_Intercept"^2
var.Bearyear.Dist <- posterior_samples(CompDHGLM)$"sd_Bearyear__scaleDailyDistance_Intercept"^2
var.res.Dist <- exp(posterior_samples(CompDHGLM)$"b_sigma_scaleDailyDistance_Intercept")^2


### shortterm repeatability ("Bearyear")
RDiurn.short <- (var.Bearyear.Diurn + var.BearID.Diurn) / 
  (var.Bearyear.Diurn + var.BearID.Diurn + var.res.Diurn)
mean(RDiurn.short);HPDinterval(as.mcmc(RDiurn.short),0.95)
# 0.57 [0.49, 0.65]

RDist.short <- (var.Bearyear.Dist + var.BearID.Dist) / 
  (var.Bearyear.Dist + var.BearID.Dist + var.res.Dist)
mean(RDist.short);HPDinterval(as.mcmc(RDist.short),0.95)
# 0.21 [0.16, 0.26]


### longterm repeatability  ("BearID") 
RDiurn.long <- var.BearID.Diurn / (var.Bearyear.Diurn + var.BearID.Diurn + var.res.Diurn)
mean(RDiurn.long);HPDinterval(as.mcmc(RDiurn.long),0.95)
# 0.41 [0.31, 0.53]

RDist.long <- var.BearID.Dist / (var.Bearyear.Dist + var.BearID.Dist + var.res.Dist)
mean(RDist.long);HPDinterval(as.mcmc(RDist.long),0.95)
# 0.07 [0.02, 0.12]


### Intercept repeatability  ("BearID")
RDiurn.int <- var.BearID.Diurn / (var.Bearyear.Diurn + var.BearID.Diurn)
mean(RDiurn.int);HPDinterval(as.mcmc(RDiurn.int),0.95)
# 0.72 [0.60, 0.83]

RDist.int <- var.BearID.Dist / (var.Bearyear.Dist + var.BearID.Dist)
mean(RDist.int);HPDinterval(as.mcmc(RDist.int),0.95)
# 0.34 [0.15, 0.55]


### Slope repeatability  ("BearID")
var.BearID.slope.Diurn <- posterior_samples(CompDHGLM)$"sd_BearID__scaleDailyDiurnality_day.scaled"^2
var.Bearyear.slope.Diurn <- posterior_samples(CompDHGLM)$"sd_Bearyear__scaleDailyDiurnality_day.scaled"^2

var.BearID.slope.Dist <- posterior_samples(CompDHGLM)$"sd_BearID__scaleDailyDistance_day.scaled"^2
var.Bearyear.slope.Dist <- posterior_samples(CompDHGLM)$"sd_Bearyear__scaleDailyDistance_day.scaled"^2

### Slope repeatability  ("BearID")
RDiurn.slope <- var.BearID.slope.Diurn / 
  (var.Bearyear.slope.Diurn + var.BearID.slope.Diurn)
mean(RDiurn.slope);HPDinterval(as.mcmc(RDiurn.slope),0.95)
# 0.2 [0, 0.41]

RDist.slope <- var.BearID.slope.Dist / 
  (var.Bearyear.slope.Dist + var.BearID.slope.Dist)
mean(RDist.slope);HPDinterval(as.mcmc(RDist.slope),0.95)
# 0.19 [0, 0.41]

#----------------------------------------------------------------------

#### 4. Coefficient of variation in predictability (CVP) ####

#----------------------------------------------------------------------

### shortterm rIIV & CVP ("Bearyear)
log.norm.res.Diurn <- exp(posterior_samples(CompDHGLM)$"sd_Bearyear__sigma_scaleDailyDiurnality_Intercept"^2)
CVP.short.Diurn <- sqrt(log.norm.res.Diurn - 1)
mean(CVP.short.Diurn);HPDinterval(as.mcmc(CVP.short.Diurn),0.95)
# 0.22 [0.19, 0.26]

log.norm.res.Dist <- exp(posterior_samples(CompDHGLM)$"sd_Bearyear__sigma_scaleDailyDistance_Intercept"^2)
CVP.short.Dist <- sqrt(log.norm.res.Dist - 1)
mean(CVP.short.Dist);HPDinterval(as.mcmc(CVP.short.Dist),0.95)
# 0.24 [0.2, 0.28]

### longterm rIIV & CVP ("BearID")
log.norm.res.Diurn <- exp(posterior_samples(CompDHGLM)$"sd_BearID__sigma_scaleDailyDiurnality_Intercept"^2)
CVP.long.Diurn <- sqrt(log.norm.res.Diurn - 1)
mean(CVP.long.Diurn);HPDinterval(as.mcmc(CVP.long.Diurn),0.95)
# 0.16 [0.10, 0.21]

log.norm.res.Dist <- exp(posterior_samples(CompDHGLM)$"sd_BearID__sigma_scaleDailyDistance_Intercept"^2)
CVP.long.Dist <- sqrt(log.norm.res.Dist - 1)
mean(CVP.long.Dist);HPDinterval(as.mcmc(CVP.long.Dist),0.95)
# 0.12 [0.06, 0.19]

#----------------------------------------------------------------------

#### 5. Correlation personality x predictability ####

#----------------------------------------------------------------------

COR.PERS.PRED <- 
  posterior_samples(CompDHGLM, 
                    pars = c("cor_BearID__scaleDailyDiurnality_Intercept__sigma_scaleDailyDiurnality_Intercept",
               "cor_Bearyear__scaleDailyDiurnality_Intercept__sigma_scaleDailyDiurnality_Intercept",
               "cor_BearID__scaleDailyDistance_Intercept__sigma_scaleDailyDistance_Intercept",
               "cor_Bearyear__scaleDailyDistance_Intercept__sigma_scaleDailyDistance_Intercept")) %>%
  gather() %>%
  separate(key,
           c(NA,"Scale",NA,NA,NA,NA,NA,"Trait",NA),
           sep = "([\\_\\__\\_\\__\\_\\_\\,])", fill = "right")

COR.PERS.PRED %>%
    ggplot(aes(x = value, y = Scale)) +
    geom_halfeyeh() +
    geom_vline(xintercept =  0) +
    facet_wrap( ~ Trait, ncol=2)
  
# credible intervals (precision of correlation)
# pd (accurracy that effect exists)

mean(COR.PERS.PRED[COR.PERS.PRED$Scale == "BearID" & 
                     COR.PERS.PRED$Trait == "scaleDailyDiurnality","value"])
HPDinterval(as.mcmc(COR.PERS.PRED[COR.PERS.PRED$Scale == "BearID" & 
                                    COR.PERS.PRED$Trait == "scaleDailyDiurnality","value"]))
pd(COR.PERS.PRED[COR.PERS.PRED$Scale == "BearID" & 
                     COR.PERS.PRED$Trait == "scaleDailyDiurnality","value"])
# 0.7 [0.44, 0.93], 100%

mean(COR.PERS.PRED[COR.PERS.PRED$Scale == "Bearyear" & 
                     COR.PERS.PRED$Trait == "scaleDailyDiurnality","value"])
HPDinterval(as.mcmc(COR.PERS.PRED[COR.PERS.PRED$Scale == "Bearyear" & 
                                    COR.PERS.PRED$Trait == "scaleDailyDiurnality","value"]))
pd(COR.PERS.PRED[COR.PERS.PRED$Scale == "Bearyear" & 
                     COR.PERS.PRED$Trait == "scaleDailyDiurnality","value"])
# 0.27 [0.08, 0.47]; 99.6%

mean(COR.PERS.PRED[COR.PERS.PRED$Scale == "BearID" & 
                     COR.PERS.PRED$Trait == "scaleDailyDistance","value"])
HPDinterval(as.mcmc(COR.PERS.PRED[COR.PERS.PRED$Scale == "BearID" & 
                                    COR.PERS.PRED$Trait == "scaleDailyDistance","value"]))
pd(COR.PERS.PRED[COR.PERS.PRED$Scale == "BearID" & 
                     COR.PERS.PRED$Trait == "scaleDailyDistance","value"])
# 0.45 [0.01, 0.84]; pd = 95.8%

mean(COR.PERS.PRED[COR.PERS.PRED$Scale == "Bearyear" & 
                     COR.PERS.PRED$Trait == "scaleDailyDistance","value"])
HPDinterval(as.mcmc(COR.PERS.PRED[COR.PERS.PRED$Scale == "Bearyear" & 
                                    COR.PERS.PRED$Trait == "scaleDailyDistance","value"]))
pd(COR.PERS.PRED[COR.PERS.PRED$Scale == "Bearyear" & 
                     COR.PERS.PRED$Trait == "scaleDailyDistance","value"])
# 0.68 [0.55, 0.81]; 100%

#----------------------------------------------------------------------

#### 6. Behavioral syndrome (Correlation BT Movement x BT Diurnality) ####

#----------------------------------------------------------------------

COR.BEHSYN <- 
  posterior_samples(CompDHGLM, 
                    pars = c("cor_BearID__scaleDailyDiurnality_Intercept__scaleDailyDistance_Intercept",
                             "cor_Bearyear__scaleDailyDiurnality_Intercept__scaleDailyDistance_Intercept")) %>%
  gather() %>%
  separate(key,
           c(NA,"Scale",NA,NA,NA,NA,NA,NA),
           sep = "([\\_\\__\\_\\__\\_\\,])", fill = "right")

COR.BEHSYN %>%
  ggplot(aes(x = value, y = Scale)) +
  geom_halfeyeh() +
  geom_vline(xintercept =  0) +
  ggtitle("Behavioral syndrome diurnality - movement")

# long-term *
mean(COR.BEHSYN[COR.BEHSYN$Scale=="BearID","value"])
HPDinterval(as.mcmc(COR.BEHSYN[COR.BEHSYN$Scale=="BearID","value"]))
pd(COR.BEHSYN[COR.BEHSYN$Scale=="BearID","value"])
#0.41 [0.08, 0.76]; 98.75%

# short-term
mean(COR.BEHSYN[COR.BEHSYN$Scale=="Bearyear","value"])
HPDinterval(as.mcmc(COR.BEHSYN[COR.BEHSYN$Scale=="Bearyear","value"]))
pd(COR.BEHSYN[COR.BEHSYN$Scale=="Bearyear","value"])
#-0.11 [-0.31, 0.09]; 86.4%

#----------------------------------------------------------------------

#### 7. Predictability syndrome (Correlation rIIV x rIIV) ####

#----------------------------------------------------------------------

summary(CompDHGLM)
COR.PREDSYN <- 
  posterior_samples(CompDHGLM, 
                    pars = c("cor_BearID__sigma_scaleDailyDiurnality_Intercept__sigma_scaleDailyDistance_Intercept",
                             "cor_Bearyear__sigma_scaleDailyDiurnality_Intercept__sigma_scaleDailyDistance_Intercept")) %>%
  gather() %>%
  separate(key,c(NA,"Scale",NA,NA,NA,NA,NA,NA,NA,NA),
           sep = "([\\_\\__\\_\\_\\__\\_\\_\\,])", fill = "right")

COR.PREDSYN %>%
  ggplot(aes(x = value, y = Scale)) +
  geom_halfeyeh() +
  geom_vline(xintercept =  0) +
  ggtitle("Predictability syndrome diurnality - movement")

# long-term *
mean(COR.PREDSYN[COR.PREDSYN$Scale=="BearID","value"])
HPDinterval(as.mcmc(COR.PREDSYN[COR.PREDSYN$Scale=="BearID","value"]))
pd(COR.PREDSYN[COR.PREDSYN$Scale=="BearID","value"])
# 0.51 [0.09, 0.87]; pd = 98%

# short-term
mean(COR.PREDSYN[COR.PREDSYN$Scale=="Bearyear","value"])
HPDinterval(as.mcmc(COR.PREDSYN[COR.PREDSYN$Scale=="Bearyear","value"]))
pd(COR.PREDSYN[COR.PREDSYN$Scale=="Bearyear","value"])
# 0.07 [-0.13, 0.28]; pd = 0.62

#----------------------------------------------------------------------

#### 8. Nonsensible correlations ####

#----------------------------------------------------------------------

nons <- posterior_samples(CompDHGLM, pars = "^cor_")%>% 
  gather() %>%
  separate(key,
           c("Scale","Trait1","Trait2"),
           sep = "__", fill = "left") %>%
  separate(Scale,
           c(NA,"Scale"),
           sep = "_") %>%
  mutate(Trait1 = str_remove(Trait1, "scaleDaily")) %>%
  mutate(Trait2 = str_remove(Trait2, "scaleDaily")) %>%
  mutate(COR = paste(Trait1, Trait2, sep = "-"))

levels(factor(nons$COR))

nons$sensible <- "0"
sensible <- c("Distance_Intercept-Distance_day.scaled","Distance_Intercept-sigma_Distance_Intercept",
  "Diurnality_Intercept-Distance_Intercept","Diurnality_Intercept-Diurnality_day.scaled",
  "Diurnality_Intercept-sigma_Diurnality_Intercept","sigma_Diurnality_Intercept-sigma_Distance_Intercept")
nons[nons$COR %in% sensible, "sensible"] <- "1"

ggplot(nons, aes(x = value, y = COR, fill = sensible)) +
  geom_halfeyeh()+
  geom_vline(xintercept = 0)+
  labs(x = "Posterior distribution",y=" ")+
  facet_wrap( ~ Scale, ncol=2)+ 
  theme_classic(base_size = 10)+
  scale_fill_manual(name = "",
                    labels = c("Non-sensible","Sensible"),
                    values = viridis::inferno(6, alpha = .6)[c(2,4)])

nons <- nons %>% 
  group_by(COR, Scale) %>%
  summarize(mean = mean(value),
            HDPup = HPDinterval(as.mcmc(value))[1],
            HDPlo = HPDinterval(as.mcmc(value))[2])

#----------------------------------------------------------------------

#### 9. Behavioral types ####

#----------------------------------------------------------------------

s.DailyDistance <- scale(PERS$DailyDistance)
s.DailyDiurnality <- scale(PERS$DailyDiurnality)

#-----------------------------------------------------------------------------------------------
# Travel distance BearID (long-term)
#-----------------------------------------------------------------------------------------------
  BT_Dist_LT <- posterior_samples(CompDHGLM, pars = "^r_BearID__scaleDailyDistance")[1:62] %>%
  tidyr::gather(BearID, value,
                "r_BearID__scaleDailyDistance[Bear1,Intercept]" : 
                  "r_BearID__scaleDailyDistance[Bear62,Intercept]")%>%
  select(BearID,value)%>%
  separate(BearID,
           c(NA,NA,NA,NA,"BearID",NA),
           sep = "([\\_\\__\\[\\,])", fill = "right") %>%
  group_by(BearID) %>%
  mutate(Mean = mean(value),
         Up = mean(value) + 1.96 * sd(value),
         Lo = mean(value) - 1.96 * sd(value))%>%
  select (-value)%>%
  filter(!duplicated(BearID))

names(BT_Dist_LT)<-c("BearID","LT_MeanDist","LT_UpDist","LT_LoDist")
range(BT_Dist_LT$LT_MeanDist)
# -0.36 - 0.36


# what does this mean in terms of actual distances?
# Add population intercept and backtransform
LT_MeanDist_back <- BT_Dist_LT %>%
  mutate(LT_MeanDist  = LT_MeanDist + fixef(CompDHGLM, pars = "scaleDailyDistance_Intercept")[1]) %>%  
  mutate(LT_MeanDist = LT_MeanDist * attr(s.DailyDistance, 'scaled:scale') + 
  attr(s.DailyDistance, 'scaled:center')) 
range(LT_MeanDist_back$LT_MeanDist)

# 7415 - 9642

#-----------------------------------------------------------------------------------------------
# Travel distance BearYear (short-term)
#-----------------------------------------------------------------------------------------------

BT_Dist_AN <- posterior_samples(CompDHGLM, pars = "^r_Bearyear__scaleDailyDistance")[1:187] %>%
  tidyr::gather(BearID, value,
                "r_Bearyear__scaleDailyDistance[Bear1_2006,Intercept]" : 
                  "r_Bearyear__scaleDailyDistance[Bear9_2009,Intercept]")%>%
  select(BearID,value)%>%
  separate(BearID,
           c(NA,NA,NA,NA,"BearID","year",NA),
           sep = "([\\_\\__\\[\\_\\,])", fill = "right") %>%
  left_join(select(BT_Dist_LT, BearID, LT_MeanDist),
            by="BearID") %>%
  mutate(value = value + LT_MeanDist) %>%
  group_by(BearID,year) %>%
  mutate(Mean = mean(value),
         Up = mean(value) + 1.96 * sd(value),
         Lo = mean(value) - 1.96 * sd(value),
         Bearyear = paste(BearID,year,sep="_")) %>%
  select (-c(value,LT_MeanDist)) %>%
  filter(!duplicated(Bearyear)) %>%
  left_join(select(PERS[!duplicated(PERS$Bearyear),], Bearyear, Status),
            by="Bearyear")

names(BT_Dist_AN)<-c("BearID","year","AN_MeanDist","AN_UpDist","AN_LoDist",
                     "Bearyear","Status")
range(BT_Dist_AN$AN_MeanDist)
# -1.1 - 1.01

# what does this mean in terms of actual distances?
# Add population intercept and reproductive status coefficient and then backtransform
AN_MeanDist_back <- BT_Dist_AN %>%
mutate(AN_MeanDist  = ifelse(Status == "SolitaryF",
                       (AN_MeanDist + fixef(CompDHGLM, pars = "scaleDailyDistance_Intercept")[1]),
                       ifelse(Status == "WithCOY",
                              (AN_MeanDist +
                                 fixef(CompDHGLM, pars = "scaleDailyDistance_Intercept")[1] +
                                 fixef(CompDHGLM, pars = "scaleDailyDistance_StatusWithCOY")[1]),
                              (AN_MeanDist +
                                 fixef(CompDHGLM, pars = "scaleDailyDistance_Intercept")[1] +
                                 fixef(CompDHGLM, pars = "scaleDailyDistance_StatusWithOFFSPRING")[1])))) %>%
mutate(AN_MeanDist = (AN_MeanDist * attr(s.DailyDistance, 'scaled:scale')) + 
         attr(s.DailyDistance, 'scaled:center'))
range(AN_MeanDist_back$AN_MeanDist)
# 3878 - 11659 m

#-----------------------------------------------------------------------------------------------
# Diurnality BearID (long-term)
#-----------------------------------------------------------------------------------------------
BT_Diurn_LT <- posterior_samples(CompDHGLM, pars = "^r_BearID__scaleDailyDiurnality")[1:62] %>%
  tidyr::gather(BearID, value,
                "r_BearID__scaleDailyDiurnality[Bear1,Intercept]" : 
                "r_BearID__scaleDailyDiurnality[Bear62,Intercept]")%>%
  select(BearID,value)%>%
  separate(BearID,
           c(NA,NA,NA,NA,"BearID",NA),
           sep = "([\\_\\__\\[\\,])", fill = "right") %>%
  group_by(BearID) %>%
  mutate(Mean = mean(value),
         Up = mean(value) + 1.96 * sd(value),
         Lo = mean(value) - 1.96 * sd(value))%>%
  select (-value)%>%
  filter(!duplicated(BearID))
names(BT_Diurn_LT)<-c("BearID","LT_MeanDiurn","LT_UpDiurn","LT_LoDiurn")
range(BT_Diurn_LT$LT_MeanDiurn)

# what does this mean in terms of actual diurnality?
# Add population intercept and then backtransform
LT_MeanDiurn_back <- BT_Diurn_LT %>%
  mutate(LT_MeanDiurn  = LT_MeanDiurn + fixef(CompDHGLM, pars = "scaleDailyDiurnality_Intercept")[1]) %>%  
  mutate(LT_MeanDiurn = (LT_MeanDiurn * attr(s.DailyDiurnality, 'scaled:scale')) + 
           attr(s.DailyDiurnality, 'scaled:center')) 
range(LT_MeanDiurn_back$LT_MeanDiurn)
# - 0.45 - 0.32

#-----------------------------------------------------------------------------------------------
# Diurnality BearYear (short-term)
#-----------------------------------------------------------------------------------------------
BT_Diurn_AN <- posterior_samples(CompDHGLM, pars = "^r_Bearyear__scaleDailyDiurnality")[1:187] %>%
  tidyr::gather(BearID, value,
                "r_Bearyear__scaleDailyDiurnality[Bear1_2006,Intercept]" : 
                  "r_Bearyear__scaleDailyDiurnality[Bear9_2009,Intercept]")%>%
  select(BearID,value)%>%
  separate(BearID,
           c(NA,NA,NA,NA,"BearID","year",NA),
           sep = "([\\_\\__\\[\\_\\,])", fill = "right") %>%
  left_join(select(BT_Diurn_LT, BearID, LT_MeanDiurn),
            by="BearID") %>%
  mutate(value = value + LT_MeanDiurn) %>%
  group_by(BearID,year) %>%
  mutate(Mean = mean(value),
         Up = mean(value) + 1.96 * sd(value),
         Lo = mean(value) - 1.96 * sd(value),
         Bearyear = paste(BearID,year,sep="_")) %>%
  select (-c(value,LT_MeanDiurn))%>%
  filter(!duplicated(Bearyear))%>%
  left_join(select(PERS[!duplicated(PERS$Bearyear),], Bearyear, Status),
            by="Bearyear")
names(BT_Diurn_AN)<-c("BearID","year","AN_MeanDiurn","AN_UpDiurn","AN_LoDiurn","Bearyear","Status")
range(BT_Diurn_AN$AN_MeanDiurn)

# what does this mean in terms of actual distances?
# Add population intercept and reproductive status coefficient and then backtransform
AN_MeanDiurn_back <- BT_Diurn_AN %>%
  mutate(AN_MeanDiurn  = ifelse(Status == "SolitaryF",
                               (AN_MeanDiurn + fixef(CompDHGLM, pars = "scaleDailyDiurnality_Intercept")[1]),
                               ifelse(Status == "WithCOY",
                                      (AN_MeanDiurn +
                                         fixef(CompDHGLM, pars = "scaleDailyDiurnality_Intercept")[1] +
                                         fixef(CompDHGLM, pars = "scaleDailyDiurnality_StatusWithCOY")[1]),
                                      (AN_MeanDiurn +
                                         fixef(CompDHGLM, pars = "scaleDailyDiurnality_Intercept")[1] +
                                         fixef(CompDHGLM, pars = "scaleDailyDiurnality_StatusWithOFFSPRING")[1])))) %>%
  mutate(AN_MeanDiurn = (AN_MeanDiurn * attr(s.DailyDiurnality, 'scaled:scale')) + 
           attr(s.DailyDiurnality, 'scaled:center')) 
range(AN_MeanDiurn_back$AN_MeanDiurn)
# -0.49 - 0.6

#-----------------------------------------------------------------------------------------------

#### 10. Behavioral predictability (individual rIIV's) ####

#-----------------------------------------------------------------------------------------------

# Intercept Distance
exp(fixef(CompDHGLM, pars = "sigma_scaleDailyDistance_Intercept")[1]) * sd(PERS$DailyDistance, na.rm = T)
# 0 = 2500 m

# Intercept Diurnality
exp(fixef(CompDHGLM, pars = "sigma_scaleDailyDiurnality_Intercept")[1]) * sd(PERS$DailyDiurnality, na.rm = T)
# 0.19

#-----------------------------------------------------------------------------------------------
# Travel distance BearID (long-term)
#-----------------------------------------------------------------------------------------------

BS_Dist_LT <- posterior_samples(CompDHGLM, pars = "^r_BearID__sigma_scaleDailyDistance") %>%
  tidyr::gather(BearID, value,
                "r_BearID__sigma_scaleDailyDistance[Bear1,Intercept]" : "r_BearID__sigma_scaleDailyDistance[Bear49,Intercept]")%>%
  select(BearID,value)%>%
  separate(BearID,
           c(NA,"BearID",NA),
           sep = "([\\[\\,])", fill = "right") %>%
  group_by(BearID) %>%
  mutate(Mean = mean(value),
         Up = mean(value) + 1.96 * sd(value),
         Lo = mean(value) - 1.96 * sd(value))%>%
  select (-value)%>%
  filter(!duplicated(BearID))
names(BS_Dist_LT)<-c("BearID","BS_LT_MeanDist","BS_LT_UpDist","BS_LT_LoDist")
range(BS_Dist_LT$BS_LT_MeanDist)

# what does this mean in terms of actual diurnality?
# Add population intercept and then backtransform
BS_Dist_LT_back <- BS_Dist_LT %>%
  mutate(BS_LT_MeanDist  = BS_LT_MeanDist + fixef(CompDHGLM, pars = "sigma_scaleDailyDistance_Intercept")[1]) %>%  
  mutate(BS_LT_MeanDist = exp(BS_LT_MeanDist)) %>% 
  mutate(BS_LT_MeanDist = BS_LT_MeanDist * sd(PERS$DailyDistance)) 
range(BS_Dist_LT_back$BS_LT_MeanDist)
# 2152 - 3090 m

#-----------------------------------------------------------------------------------------------
# Travel distance BearYear (short-term)
#-----------------------------------------------------------------------------------------------
BS_Dist_AN <- posterior_samples(CompDHGLM, pars = "^r_Bearyear__sigma_scaleDailyDistance") %>%
  tidyr::gather(BearID, value,
                "r_Bearyear__sigma_scaleDailyDistance[Bear1_2006,Intercept]" :
                  "r_Bearyear__sigma_scaleDailyDistance[Bear9_2009,Intercept]")%>%
  select(BearID,value)%>%
  separate(BearID,
           c(NA,NA,NA,NA,NA,"BearID","year",NA),
           sep = "([\\[\\_\\,])", fill = "right") %>%
  left_join(select(BS_Dist_LT, BearID, BS_LT_MeanDist),
            by="BearID") %>%
  mutate(value = value + BS_LT_MeanDist) %>%
  group_by(BearID,year) %>%
  mutate(Mean = mean(value),
         Up = mean(value) + 1.96 * sd(value),
         Lo = mean(value) - 1.96 * sd(value),
         Bearyear = paste(BearID,year,sep="_"))%>%
  select (-c(value,BS_LT_MeanDist))%>%
  filter(!duplicated(Bearyear))%>% 
  left_join(select(PERS[!duplicated(PERS$Bearyear),], Bearyear, Status),
            by="Bearyear")
names(BS_Dist_AN)<-c("BearID","year","BS_AN_MeanDist","BS_AN_UpDist","BS_AN_LoDist","Bearyear", "Status")
range(BS_Dist_AN$BS_AN_MeanDist, na.rm=T)

# what does this mean in terms of actual distances?
# Add population intercept and reproductive status coefficient and then backtransform
BS_Dist_AN_back <- BS_Dist_AN %>%
  mutate(BS_AN_MeanDist  = ifelse(Status == "SolitaryF",
                               (BS_AN_MeanDist + fixef(CompDHGLM, pars = "sigma_scaleDailyDistance_Intercept")[1]),
                               ifelse(Status == "WithCOY",
                                      (BS_AN_MeanDist +
                                         fixef(CompDHGLM, pars = "sigma_scaleDailyDistance_Intercept")[1] +
                                         fixef(CompDHGLM, pars = "sigma_scaleDailyDistance_StatusWithCOY")[1]),
                                      (BS_AN_MeanDist +
                                         fixef(CompDHGLM, pars = "sigma_scaleDailyDistance_Intercept")[1] +
                                         fixef(CompDHGLM, pars = "sigma_scaleDailyDistance_StatusWithOFFSPRING")[1])))) %>%
  mutate(BS_AN_MeanDist = exp(BS_AN_MeanDist)) %>% 
  mutate(BS_AN_MeanDist = BS_AN_MeanDist * sd(PERS$DailyDistance, na.rm = T)) 

range(BS_Dist_AN_back$BS_AN_MeanDist, na.rm=T)

# 1039 - 5431 meters

#-----------------------------------------------------------------------------------------------
#  diurnality BearID (long-term)
#-----------------------------------------------------------------------------------------------
BS_Diurn_LT <- posterior_samples(CompDHGLM, pars = "^r_BearID__sigma_scaleDailyDiurnality") %>%
  tidyr::gather(BearID, value,
                "r_BearID__sigma_scaleDailyDiurnality[Bear1,Intercept]" : 
                "r_BearID__sigma_scaleDailyDiurnality[Bear62,Intercept]")%>%
  select(BearID,value)%>%
  separate(BearID,
           c(NA,"BearID",NA),
           sep = "([\\[\\,])", fill = "right") %>%
  group_by(BearID) %>%
  mutate(Mean = mean(value),
         Up = mean(value) + 1.96 * sd(value),
         Lo = mean(value) - 1.96 * sd(value))%>%
  select (-value)%>%
  filter(!duplicated(BearID))

names(BS_Diurn_LT)<-c("BearID","BS_LT_MeanDiurn","BS_LT_UpDiurn","BS_LT_LoDiurn")
range(BS_Diurn_LT$BS_LT_MeanDiurn)

# what does this mean in terms of actual diurnality?
# Add population intercept and then backtransform
BS_Diurn_LT_back <- BS_Diurn_LT %>%
  mutate(BS_LT_MeanDiurn  = BS_LT_MeanDiurn + fixef(CompDHGLM, pars = "sigma_scaleDailyDiurnality_Intercept")[1]) %>%  
  mutate(BS_LT_MeanDiurn = exp(BS_LT_MeanDiurn)) %>% 
  mutate(BS_LT_MeanDiurn = BS_LT_MeanDiurn * sd(PERS$DailyDiurnality, na.rm = T)) 
range(BS_Diurn_LT_back$BS_LT_MeanDiurn, na.rm=T)
# 0.15 - 0.25 Diurnality Index units

#-----------------------------------------------------------------------------------------------
#  diurnality BearYear (short-term)
#-----------------------------------------------------------------------------------------------

BS_Diurn_AN <- posterior_samples(CompDHGLM, pars = "^r_Bearyear__sigma_scaleDailyDiurnality") %>%
  tidyr::gather(BearID, value,
                "r_Bearyear__sigma_scaleDailyDiurnality[Bear1_2006,Intercept]" :
                  "r_Bearyear__sigma_scaleDailyDiurnality[Bear9_2009,Intercept]")%>%
  select(BearID,value)%>%
  separate(BearID,
           c(NA,NA,NA,NA,NA,"BearID","year",NA),
           sep = "([\\[\\_\\,])", fill = "right") %>%
  left_join(select(BS_Diurn_LT, BearID, BS_LT_MeanDiurn),
            by="BearID") %>%
  mutate(value = value + BS_LT_MeanDiurn) %>%
  group_by(BearID,year) %>%
  mutate(Mean = mean(value),
         Up = mean(value) + 1.96 * sd(value),
         Lo = mean(value) - 1.96 * sd(value),
         Bearyear = paste(BearID,year,sep="_"))%>%
  select (-c(value,BS_LT_MeanDiurn))%>%
  filter(!duplicated(Bearyear)) %>% 
  left_join(select(PERS[!duplicated(PERS$Bearyear),], Bearyear, Status),
            by="Bearyear")

names(BS_Diurn_AN) <- c("BearID","year","BS_AN_MeanDiurn","BS_AN_UpDiurn",
                        "BS_AN_LoDiurn","Bearyear", "Status")
range(BS_Diurn_AN$BS_AN_MeanDiurn)


# what does this mean in terms of actual distances?
# Add population intercept and reproductive status coefficient and then backtransform
BS_Diurn_AN_back <- BS_Diurn_AN %>%
  mutate(BS_AN_MeanDiurn  = ifelse(Status == "SolitaryF",
                                  (BS_AN_MeanDiurn + fixef(CompDHGLM, pars = "sigma_scaleDailyDiurnality_Intercept")[1]),
                                  ifelse(Status == "WithCOY",
                                         (BS_AN_MeanDiurn +
                                            fixef(CompDHGLM, pars = "sigma_scaleDailyDiurnality_Intercept")[1] +
                                            fixef(CompDHGLM, pars = "sigma_scaleDailyDiurnality_StatusWithCOY")[1]),
                                         (BS_AN_MeanDiurn +
                                            fixef(CompDHGLM, pars = "sigma_scaleDailyDiurnality_Intercept")[1] +
                                            fixef(CompDHGLM, pars = "sigma_scaleDailyDiurnality_StatusWithOFFSPRING")[1])))) %>%
  mutate(BS_AN_MeanDiurn = exp(BS_AN_MeanDiurn)) %>% 
  mutate(BS_AN_MeanDiurn = BS_AN_MeanDiurn * sd(PERS$DailyDiurnality, na.rm=T)) 

range(BS_Diurn_AN_back$BS_AN_MeanDiurn)
# 0.11 - 0.39


# Identify a more predictable individuals

head(LT_MeanDist_back[order(LT_MeanDist_back$LT_MeanDist),c("BearID","LT_MeanDist")],10)
head(BS_Dist_LT_back[order(BS_Dist_LT_back$BS_LT_MeanDist),c("BearID","BS_LT_MeanDist")],20)

head(LT_MeanDiurn_back[order(LT_MeanDiurn_back$LT_MeanDiurn),c("BearID","LT_MeanDiurn")],10)
head(BS_Diurn_LT_back[order(BS_Diurn_LT_back$BS_LT_MeanDiurn),c("BearID","BS_LT_MeanDiurn")],10)

# More predictable
# Bear 33
# 7415 +- 2259
# -0.272 +- 0.16

# Identify a less predictable individuals

head(LT_MeanDist_back[order(-LT_MeanDist_back$LT_MeanDist),c("BearID","LT_MeanDist")],10)
head(BS_Dist_LT_back[order(-BS_Dist_LT_back$BS_LT_MeanDist),c("BearID","BS_LT_MeanDist")],10)

head(LT_MeanDiurn_back[order(-LT_MeanDiurn_back$LT_MeanDiurn),c("BearID","LT_MeanDiurn")],10)
head(BS_Diurn_LT_back[order(-BS_Diurn_LT_back$BS_LT_MeanDiurn),c("BearID","BS_LT_MeanDiurn")],10)

# More unpredictable
# Bear 37
# 9640 +- 3090
# 0.256 +- 0.238
