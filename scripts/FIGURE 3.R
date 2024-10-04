library(pracma)
library(mvnormtest)
library(coda)
library("postMCMCglmm")
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(data.table)
library(extrafont)
library(ggpubr)
library(brms)

Personality.N<-readRDS("DATA_R1.rds")

mytheme<-theme(legend.justification = c(-0.1, 1.1), 
               legend.position = c(0, 1), 
               legend.title = element_text(size=8),
               legend.text = element_text(size=8),
               legend.background = element_rect(fill ="transparent"),
               legend.direction ="vertical",
               legend.box = "vertical",
               legend.margin = margin(0,0,-0.3,0, unit = "cm"),
               legend.box.just = c("top"), 
               axis.title.x = element_text(colour = "black",size = 10),
               axis.title.y = element_text(colour = "black",size = 10),
               axis.text.x = element_text(colour = "black",size = 10),
               axis.text.y = element_text(colour = "black",size = 10),
               panel.border = element_rect(fill = NA),
               plot.tag.position = c(0.015, 0.99), 
               plot.tag = element_text(face="bold", size = 10))

myscale <- scale_colour_manual(name = "",
                    values = c("#08306B"),
                    labels = c("Short-term scale"))
  

AN_BT <- readRDS("AN_BT.rds")

CompDHGLM <- readRDS("CompositeModel.rds")

#-----------------------------------------------------------------------------------------------

### CORRELATION OF Behavioral Type & Predictability

#-----------------------------------------------------------------------------------------------
get_variables(CompDHGLM)
# Slope Movement - BearID
cov.Trav <-
  posterior_samples(CompDHGLM)[,49] *
  sqrt((posterior_samples(CompDHGLM)[,27])^2) *
  sqrt((posterior_samples(CompDHGLM)[,29])^2)
var.Trav <- (posterior_samples(CompDHGLM)[,27])^2
LT_slope_Trav <- cov.Trav / var.Trav
mean(LT_slope_Trav)

# Slope Movement - Bearyear
cov.Trav_AN <-
  posterior_samples(CompDHGLM)[,64] *
  sqrt((posterior_samples(CompDHGLM)[,33])^2) *
  sqrt((posterior_samples(CompDHGLM)[,35])^2)
var.Trav_AN <- (posterior_samples(CompDHGLM)[,33])^2
AN_slope_Trav <- cov.Trav_AN / var.Trav_AN


# Panel a)

BT_SP_Trav_a <- ggplot() + 
  geom_segment(data = AN_BT[!duplicated(AN_BT$BearID),],
               aes(x = LT_LoDist,
                   xend = LT_UpDist,
                   y = BS_LT_MeanDist,
                   yend = BS_LT_MeanDist), 
               color = "#F8766D", alpha = 0.2) +
  
  geom_segment(data = AN_BT[!duplicated(AN_BT$BearID),], 
               aes(x = LT_MeanDist,
                   xend = LT_MeanDist, 
                   y = BS_LT_LoDist,
                   yend = BS_LT_UpDist), 
               color = "#F8766D", alpha = 0.2) +
 
   geom_point(data = AN_BT[!duplicated(AN_BT$BearID),],
             aes(x = LT_MeanDist, y = BS_LT_MeanDist,
                 fill="#F8766D"),shape=22,
             size=1.4, color = "#F8766D", alpha = 0.8, stroke = 0) +
  
  geom_segment(aes(x = -0.8,
                   xend = 0.8,
                   y = 0 + mean(LT_slope_Trav)*-1, 
                   yend = 0 + mean(LT_slope_Trav)*1),
  color="#F8766D",size=1, alpha = 0.8)+
  
  ylab("Predictability movement") + 
  xlab("Behavioral type movement")+ 
  scale_x_continuous(limits = c(-2.1,2.4), breaks = c(-2,0,2))+
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1,1.1),
                     labels =c("-.5","0",".5"), expand = c(0,0))+
  theme_classic(base_size = 10, base_family = "Georgia")+
  mytheme +
  myscale +
  scale_fill_manual("",values = "#F8766D",labels = "Long-term scale")+
  guides(
    fill=guide_legend(order = 2,
                      keywidth = 0.08,
                      keyheight = 0.1,
                      default.unit = "inch"),
    color=guide_legend(order = 1,
                       keywidth = 0.08,
                       keyheight = 0.1,
                       default.unit = "inch"))+
  annotate("text",x = 1.3, y = -0.5, 
           label = expression(paste(italic("r")" = 0.44")),
           size = 2.5)+
  labs(tag = "(a)")


# Panel b)

BT_SP_Trav_b <- ggplot() + 
  geom_segment(data = AN_BT,
               aes(x = AN_LoDist,
                   xend = AN_UpDist,
                   y = BS_AN_MeanDist,
                   yend = BS_AN_MeanDist), color = "#08306B", alpha = 0.2) +
  
  geom_segment(data=AN_BT, aes(x = AN_MeanDist,xend = AN_MeanDist,
                               y = BS_AN_LoDist,yend = BS_AN_UpDist), color = "#08306B", alpha = 0.2) +
  
  geom_point(data=AN_BT,aes(x = AN_MeanDist, y = BS_AN_MeanDist,color = "#08306B"), 
             shape=16,size=1.2, alpha = 0.8) +

  geom_segment(aes(x = -1.5, 
                   xend = 1.3, 
                   y = 0 + mean(AN_slope_Trav)*-1.5, 
                   yend = 0 + mean(AN_slope_Trav)*1.3),
               color = "#08306B",size = 1, alpha = 0.8)+
  
  ylab("Predictability movement") + 
  xlab("Behavioral type movement")+ 
  scale_x_continuous(limits = c(-2.1,2.4), breaks = c(-2,0,2))+
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1,1.1),
                     labels = c("-.5","0",".5"), expand = c(0,0))+
  theme_classic(base_size = 10, base_family = "Georgia")+
  mytheme +
  myscale+
  scale_fill_manual("",values = "#F8766D",labels = "Long-term scale")+
  guides(
    fill=guide_legend(order = 2,
                      keywidth = 0.08,
                      keyheight = 0.1,
                      default.unit = "inch"),
    color=guide_legend(order = 1,
                       keywidth = 0.08,
                       keyheight = 0.1,
                       default.unit = "inch"))+
  annotate("text",x = 1.3, y = -0.5, 
           label = expression(paste(italic("r")[short-term]*" = 0.68")),
           size = 2.5)+
  labs(tag = "(b)")

BT_SP_Trav_b


#-----------------------------------------------------------------------------------------------
# Diurnality
#-----------------------------------------------------------------------------------------------

# Slope Diurnality - BearID
cov.DIUR <-
  posterior_samples(CompDHGLM)[,37] *
  sqrt((posterior_samples(CompDHGLM)[,24])^2) *
  sqrt((posterior_samples(CompDHGLM)[,26])^2)
VAR.DIUR <- (posterior_samples(CompDHGLM)[,24])^2
LT_slope_DIURN <- cov.DIUR / VAR.DIUR

# Slope Diurnality - BearYear
cov.DIURN_BearYear <-
  posterior_samples(CompDHGLM)[,52] *
  sqrt((posterior_samples(CompDHGLM)[,30])^2) *
  sqrt((posterior_samples(CompDHGLM)[,32])^2)
var.diurn.BearYear <- (posterior_samples(CompDHGLM)[,30])^2
AN_slope_DIURN <- cov.DIURN_BearYear / var.diurn.BearYear


# Panel c)

BT_SP_Diur_c <- ggplot() + 
  geom_segment(data = AN_BT[!duplicated(AN_BT$BearID),],
               aes(x = LT_LoDiurn,
                   xend = LT_UpDiurn,
                   y = BS_LT_MeanDiurn,
                   yend = BS_LT_MeanDiurn), 
               color = "#F8766D", alpha = 0.2) +
  
  geom_segment(data = AN_BT[!duplicated(AN_BT$BearID),], 
               aes(x = LT_MeanDiurn, 
                   xend = LT_MeanDiurn, 
                   y = BS_LT_LoDiurn,
                   yend = BS_LT_UpDiurn), 
               color = "#F8766D", alpha = 0.2) +
  
  geom_point(data = AN_BT[!duplicated(AN_BT$BearID),],
             aes(x = LT_MeanDiurn, 
                 y = BS_LT_MeanDiurn,
                 fill = "#F8766D"),
             shape = 22, size = 1.4, color = "#F8766D",alpha = 0.8, stroke = 0) +
  
  geom_segment(aes(x = -1.5,
                   xend = 1.5,
                   y = 0 + mean(LT_slope_DIURN)*-1.5, 
                   yend = 0 + mean(LT_slope_DIURN)*1.5), 
               color = "#F8766D", size = 1, alpha = 0.8)+
  
  ylab("Predictability diurnality") + 
  xlab("Behavioral type diurnality")+ 
  scale_x_continuous(limits = c(-2.1,2.4), breaks = c(-2,0,2))+
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1,1.1),
                     labels = c("-.5","0",".5"), expand = c(0,0))+
  theme_classic(base_size = 10, base_family='Georgia')+
  mytheme +
  myscale +
  scale_fill_manual("", values = "#F8766D",labels = "Long-term scale")+
  guides(
    fill=guide_legend(order = 2,
                      keywidth = 0.15,
                      keyheight = 0.15,
                      default.unit = "inch"),
    color=guide_legend(order = 1,
                       keywidth = 0.15,
                       keyheight = 0.15,
                       default.unit = "inch")) +
  labs(tag = "(c)")+
  annotate("text",x = 1.3, y = -0.5, 
           label = expression(paste(italic("r")[long-term]*" = 0.69")),
           size = 2.5)


# Panel d)

BT_SP_Diur_d <- ggplot() + 
  geom_segment(data = AN_BT,
               aes(x = AN_LoDiurn,
                   xend = AN_UpDiurn,
                   y = BS_AN_MeanDiurn,
                   yend = BS_AN_MeanDiurn), 
               color = "#08306B", alpha = 0.2) +
  
  geom_segment(data = AN_BT, 
               aes(x = AN_MeanDiurn,
                   xend = AN_MeanDiurn,
                   y = BS_AN_LoDiurn,
                   yend = BS_AN_UpDiurn), 
               color = "#08306B", alpha = 0.2) +
  
  geom_point(data=AN_BT,
             aes(x = AN_MeanDiurn, 
                 y = BS_AN_MeanDiurn,
                 color = "#08306B"),
             shape=16,size=1, alpha = 0.8) +
 
   geom_segment(aes(x = -1.8,
                   xend = 1.8,
                   y = 0 + mean(AN_slope_DIURN)*-1.8,
                   yend = 0 + mean(AN_slope_DIURN)*1.8),
               color="#08306B",size=1, alpha = 0.8)+
  
  ylab("Predictability diurnality") + 
  xlab("Behavioral type diurnality")+ 
  scale_x_continuous(limits = c(-2.1,2.4), breaks = c(-2,0,2))+
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1,1.1),
                     labels = c("-.5","0",".5"), expand = c(0,0))+
  theme_classic(base_size = 10, base_family='Georgia')+
  mytheme+
  myscale+
  scale_fill_manual("",values = "#F8766D",labels = "Long-term scale")+
  guides(
    fill=guide_legend(order = 2,
                      keywidth = 0.15,
                      keyheight = 0.15,
                      default.unit = "inch"),
    color=guide_legend(order = 1,
                       keywidth = 0.15,
                       keyheight = 0.15,
                       default.unit = "inch")) +
  labs(tag = "(d)")+
  annotate("text",x = 1.3, y = -0.5, 
           label = expression(paste(italic("r")[short-term]*" = 0.27")),
           size = 2.5)


png("FIGURE3.png", width = 174, 
     height = 160, units = 'mm', 
     res = 300, pointsize = 8)

grid.arrange(arrangeGrob(BT_SP_Trav_a,BT_SP_Trav_b,
                         BT_SP_Diur_c, BT_SP_Diur_d,
                         nrow=2,ncol=2))
dev.off()


