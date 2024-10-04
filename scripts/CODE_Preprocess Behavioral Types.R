
#-----------------------------------------------------------------------------------------------
### BEHAVIORAL TYPE
#-----------------------------------------------------------------------------------------------
PERS <- readRDS("DATA_R1.rds")

CompDHGLM <- readRDS("CompositeModel.rds")

### BT - BearID - Diurnality
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

### BT - BearID - Distance

BT_Dist_LT <-   BT_Dist_LT <- posterior_samples(CompDHGLM, pars = "^r_BearID__scaleDailyDistance")[1:62] %>%
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

### BT - Bearyear - Diurnality
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


# BT - Bearyear - Distance

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
  filter(!duplicated(Bearyear)) 

names(BT_Dist_AN)<-c("BearID","year","AN_MeanDist","AN_UpDist","AN_LoDist",
                     "Bearyear")

#-----------------------------------------------------------------------------------------------
### BEHAVIORAL SPECIALIZATION
#-----------------------------------------------------------------------------------------------

# BS - BearID - Distance

BS_Dist_LT <- posterior_samples(CompDHGLM, pars = "^r_BearID__sigma_scaleDailyDistance") %>%
  tidyr::gather(BearID, value,
                "r_BearID__sigma_scaleDailyDistance[Bear1,Intercept]" : 
                  "r_BearID__sigma_scaleDailyDistance[Bear62,Intercept]")%>%
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

# BS - BearID - Diurnality
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

# BS - Bearyear - Diurnality

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
  filter(!duplicated(Bearyear))
names(BS_Diurn_AN) <- c("BearID","year","BS_AN_MeanDiurn","BS_AN_UpDiurn",
                        "BS_AN_LoDiurn","Bearyear")

# BS - Bearyear - Distance
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
  filter(!duplicated(Bearyear))
names(BS_Dist_AN)<-c("BearID","year","BS_AN_MeanDist","BS_AN_UpDist","BS_AN_LoDist","Bearyear")


nrow(BT_Dist_AN)
nrow(BT_Diurn_AN)
nrow(BT_Dist_LT)
nrow(BT_Diurn_LT)

AN_BT<-merge(BT_Diurn_LT,BT_Dist_LT,by="BearID",all.x=T)
AN_BT<-merge(AN_BT,BS_Diurn_LT,by="BearID",all.x=T)
AN_BT<-merge(AN_BT,BS_Dist_LT,by="BearID",all.x=T)
AN_BT<-merge(AN_BT,BT_Diurn_AN[,c(1,3:7)],by="BearID",all.x=T)
AN_BT<-merge(AN_BT,BT_Dist_AN[,c(1,3:6)],by=c("BearID","Bearyear"),all.x=T)
AN_BT<-merge(AN_BT,BS_Diurn_AN[,c(1,3:6)],by=c("BearID","Bearyear"),all.x=T)
AN_BT<-merge(AN_BT,BS_Dist_AN[,c(1,3:6)],by=c("BearID","Bearyear"),all.x=T)

head(AN_BT)
head(PERS)

#AN_BT<-merge(AN_BT,PERS[!duplicated(PERS$Bearyear),c(2,8)],by=c("Bearyear"),all.x=T)

saveRDS(AN_BT,"AN_BT.rds")
