rm(list = ls())

library(tidyverse)
library(wesanderson)
library(reshape2)
library(car)
library(scico)
library(png)
library(cowplot)



#setwd("/Users/mandskog/Documents/Assemble")
setwd("C:/Users/admin/ownCloud/Assemble/")
full  <-read.csv("20210805_Results.csv", dec=".",sep = ",",stringsAsFactors = FALSE)
PAR.algae  <-read.csv("PAR_cal_chambers_Algae", dec=".",sep = ",",stringsAsFactors = FALSE)
PAR.seagrass  <-read.csv("PAR_cal_chambers_Seagrass.csv", dec=".",sep = ",",stringsAsFactors = FALSE)

# O2, chamber and plant variables --------------------------------------------------

names(full)
shoot<- full %>%
  group_by(Chamber,Organism,Treatment,Shoot.ID) %>%
  summarise_if(is.numeric, mean)

#O2 is with a single volume (average) for all seagrass bags
#O2vol is with the original volume measured

names(shoot)
params=shoot[c("Chamber","Organism","Treatment","Volume","Plant.biomass","Shoot.Dry.weight","Root.Dry.weight",
               "Initial.O2.original","End.O2",
               "Initial.DOC.original","End.WC.DOC",
               "Initial.DON.original","End.DON","Initial.c.original","c.end",
               "Initial.CDOM.original","CDOM.end", "Fucoidan")]
Plant.biomass <- params %>% 
  dplyr::group_by(Chamber,Organism, Treatment) %>% 
  dplyr::summarise(Plant.biomass = sum(Plant.biomass, na.rm = T),
            New.biomass = sum(Shoot.Dry.weight, na.rm = T),
            Old.biomass = sum(Root.Dry.weight, na.rm = T))

chamber<- params %>%
  group_by(Chamber,Organism, Treatment) %>%
  summarise_all(mean) %>% 
  select(-c(Plant.biomass, Shoot.Dry.weight,Root.Dry.weight))

chamber <- merge(chamber, Plant.biomass, by=c("Chamber","Organism", "Treatment"))


# Initial vs End figures --------------------------------------------------

chamber.1.6 = chamber[!is.na(chamber$Fucoidan),]

#Extract all initial columns
initial.strings <- c("Initial","initial")
end.strings <- c("End", 'end')

initial = cbind(chamber.1.6[1:3],chamber.1.6[,grepl(paste0(initial.strings,collapse="|"),colnames(chamber.1.6))])
end = cbind(chamber.1.6[1:3],chamber.1.6[,grepl(paste0(end.strings,collapse="|"),colnames(chamber.1.6))])
end.light <- subset(end, end$Treatment == "light")
end.dark <- subset(end, end$Treatment == "dark")

initial$Treatment <- "initial"
names(initial)

new.names=c("Chamber","Organism","Treatment","O2","DOC","TDN","c","CDOM")
colnames(initial) <- new.names
colnames(end.light) <- new.names
colnames(end.dark) <- new.names

init.v.end = rbind(initial, end.light, end.dark)
init.v.end$Treatment <- factor(init.v.end$Treatment , levels=c("initial", "light", "dark"))
# init.v.end$Treatment = as.factor(init.v.end$Treatment)

##Figures
fonts <- list(
  sans = "Helvetica",
  mono = "Consolas",
  `Times New Roman` = "DejaVu Serif"
)
font_family<-"sans"
font_size<-10

time_cat_colors<-c(scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.7),
                   scico(2, palette = "vikO",begin=0.25,end=0.75, alpha=0.7))

ylabs = data.frame(matrix(nrow = 5))
ylabs[1,1] = expression('O'[2]~(mu*'mol'~L^-1))
ylabs[2,1] = expression('C'~('mg'~L^-1))
ylabs[3,1] = expression('N'~('mg'~L^-1))
ylabs[4,1] = expression('Fluorescence')
ylabs[5,1] = expression('Fluorescence')

#Algae figures
Fucus = data.frame(subset(init.v.end, init.v.end$Organism == "Algae.unfiltered"))
Fucus$Treatment<-factor(Fucus$Treatment, 
                        levels=c("initial", "dark", "light"))

#O2
Fuc_O2 <- ggplot(data = Fucus)+
  geom_boxplot(aes(x = Treatment, y = O2, fill=Treatment), position = position_dodge(width = 1),
           color="#262626")+
  geom_point(aes(x=Treatment, y=O2, group=Treatment, fill=Treatment),
             pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
  # annotate("text", x=0.5, y=max(Fucus$O2), label = "italic(Fucus~vesiculosus)", parse=T, vjust = 1, hjust = 0,size=3)+
  scale_fill_manual(name="", labels=c("initial", "dark", "light"),
                    values = time_cat_colors)+
  scale_y_continuous(name=expression('O'[2]~(mu*'mol'~L^-1)))+
  scale_x_discrete(name="")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))+
  labs(title = expression(O[2]))


#DOC
Fuc_DOC <- ggplot(data = Fucus)+
  geom_boxplot(aes(x = Treatment, y = DOC, fill=Treatment), position = position_dodge(width = 1),
           color="#262626")+
  geom_point(aes(x=Treatment, y=DOC, group=Treatment, fill=Treatment),
             pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
  # annotate("text", x=0.5, y=max(Fucus$DOC), label = "italic(Fucus~vesiculosus)", parse=T, vjust = 1, hjust = 0,size=3)+
  scale_fill_manual(name="", labels=c("initial", "dark", "light"),
                    values = time_cat_colors)+
  scale_y_continuous(name=expression('C'~('mg'~L^-1)))+
  scale_x_discrete(name="")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))+
  labs(title = expression(DOC))

#TDN
Fuc_TDN <- ggplot(data = Fucus)+
  geom_boxplot(aes(x = Treatment, y = TDN, fill=Treatment), position = position_dodge(width = 1),
           color="#262626")+
  geom_point(aes(x=Treatment, y=TDN, group=Treatment, fill=Treatment),
             pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
  # annotate("text", x=0.5, y=max(Fucus$TDN), label = "italic(Fucus~vesiculosus)", parse=T, vjust = 1, hjust = 0,size=3)+
  scale_fill_manual(name="", labels=c("initial", "dark", "light"),
                    values = time_cat_colors)+
  scale_y_continuous(name=expression('N'~('mg'~L^-1)))+
  scale_x_discrete(name="")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))+
  labs(title = expression(TDN))

#c
Fuc_c <- ggplot(data = Fucus)+
  geom_boxplot(aes(x = Treatment, y = c, fill=Treatment), position = position_dodge(width = 1),
           color="#262626")+
  geom_point(aes(x=Treatment, y=c, group=Treatment, fill=Treatment),
             pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
  # annotate("text", x=0.5, y=max(Fucus$c), label = "italic(Fucus~vesiculosus)", parse=T, vjust = 1, hjust = 0,size=3)+
  scale_fill_manual(name="", labels=c("initial", "dark", "light"),
                    values = time_cat_colors)+
  scale_y_continuous(name="Fluorescence")+
  scale_x_discrete(name="")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))+
  labs(title = "Humic-like DOM")

#CDOM
Fuc_CDOM <- ggplot(data = Fucus)+
  geom_boxplot(aes(x = Treatment, y = CDOM, fill=Treatment), position = position_dodge(width = 1),
           color="#262626")+
  geom_point(aes(x=Treatment, y=CDOM, group=Treatment, fill=Treatment),
             pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
  # annotate("text", x=0.5, y=max(Fucus$CDOM), label = "italic(Fucus~vesiculosus)", parse=T, vjust = 1, hjust = 0,size=3)+
  scale_fill_manual(name="", labels=c("initial", "dark", "light"),
                    values = time_cat_colors)+
  scale_y_continuous(name="Absorbance")+
  scale_x_discrete(name="")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))+
  labs(title = "CDOM 254 nm")

plot_grid(Fuc_O2, Fuc_DOC,Fuc_TDN,Fuc_c,Fuc_CDOM, nrow=1)

##Figures Zostera
Zostera = data.frame(subset(init.v.end, init.v.end$Organism == "Seagrass"))
Zostera$Treatment<-factor(Zostera$Treatment, 
                        levels=c("initial", "dark", "light"))

#O2
Zos_O2 <- ggplot(data = Zostera)+
  geom_boxplot(aes(x = Treatment, y = O2, fill=Treatment), position = position_dodge(width = 1),color="#262626")+
  geom_point(aes(x=Treatment, y=O2, group=Treatment, fill=Treatment),
             pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
  # annotate("text", x=0.5, y=max(Zostera$O2), label = "italic(Zostera~marina)", parse=T, vjust = 1, hjust = 0,size=3)+
  scale_fill_manual(name="", labels=c("initial", "dark", "light"),
                    values = time_cat_colors)+
  scale_y_continuous(name=expression('O'[2]~(mu*'mol'~L^-1)))+
  scale_x_discrete(name="")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))+
  labs(title = expression(O[2]))


#DOC
Zos_DOC <- ggplot(data = Zostera)+
  geom_boxplot(aes(x = Treatment, y = DOC, fill=Treatment), position = position_dodge(width = 1),color="#262626")+
  geom_point(aes(x=Treatment, y=DOC, group=Treatment, fill=Treatment),
             pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
  # annotate("text", x=0.5, y=max(Zostera$DOC), label = "italic(Zostera~marina)", parse=T, vjust = 1, hjust = 0,size=3)+
  scale_fill_manual(name="", labels=c("initial", "dark", "light"),
                    values = time_cat_colors)+
  scale_y_continuous(name=expression('C'~('mg'~L^-1)))+
  scale_x_discrete(name="")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))+
  labs(title = "DOC")

#TDN
Zos_TDN <- ggplot(data = Zostera)+
  geom_boxplot(aes(x = Treatment, y = TDN, fill=Treatment), position = position_dodge(width = 1),color="#262626")+
  geom_point(aes(x=Treatment, y=TDN, group=Treatment, fill=Treatment),
             pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
  # annotate("text", x=0.5, y=max(Zostera$TDN), label = "italic(Zostera~marina)", parse=T, vjust = 1, hjust = 0,size=3)+
  scale_fill_manual(name="", labels=c("initial", "dark", "light"),
                    values = time_cat_colors)+
  scale_y_continuous(name=expression('N'~('mg'~L^-1)))+
  scale_x_discrete(name="")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))+
  labs(title = "TDN")

#c
Zos_c <- ggplot(data = Zostera)+
  geom_boxplot(aes(x = Treatment, y = c, fill=Treatment), position = position_dodge(width = 1),color="#262626")+
  geom_point(aes(x=Treatment, y=c, group=Treatment, fill=Treatment),
             pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
  # annotate("text", x=0.5, y=max(Zostera$c), label = "italic(Zostera~marina)", parse=T, vjust = 1, hjust = 0,size=3)+
  scale_fill_manual(name="", labels=c("initial", "dark", "light"),
                    values = time_cat_colors)+
  scale_y_continuous(name="Fluorescence")+
  scale_x_discrete(name="")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))+
  labs(title = "Humic-like DOM")

#CDOM
Zos_CDOM <- ggplot(data = Zostera)+
  geom_boxplot(aes(x = Treatment, y = CDOM, fill=Treatment), position = position_dodge(width = 1),color="#262626")+
  geom_point(aes(x=Treatment, y=CDOM, group=Treatment, fill=Treatment),
             pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
  # annotate("text", x=0.5, y=max(Zostera$CDOM), label = "italic(Zostera~marina)", parse=T, vjust = 1, hjust = 0,size=3)+
  scale_fill_manual(name="", labels=c("initial", "dark", "light"),
                    values = time_cat_colors)+
  scale_y_continuous(name="Fluorescence")+
  scale_x_discrete(name="")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))+
  labs(title = "CDOM 254 nm")

plot_grid(Zos_O2, Zos_DOC,Zos_TDN,Zos_c,Zos_CDOM, nrow=1)


# save fig on incubations #####

tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210806_Fig_incubations.tiff",
     height = 16, width=30, units = "cm", res=300)
plot_grid(Fuc_O2, Fuc_DOC,Fuc_TDN,Fuc_c,Fuc_CDOM, Zos_O2, Zos_DOC,Zos_TDN,Zos_c,Zos_CDOM, nrow=2)
dev.off()


#########END##########

blue_env<-wes_palette("Zissou1")[1]
yellow_light<-wes_palette("Zissou1")[3]
red_dark<-wes_palette("Zissou1")[5]

names(PAR.algae)
PAR.algae$Raw_time = as.factor(PAR.algae$Raw_time)
PAR.algae$Raw_time = as.numeric(PAR.algae$Raw_time)

#For every timepoint calculate a percentage light in bag from ambient light. 

plot(PAR.algae$Raw_time, PAR.algae$Lightintensity_lux*0.03916453, type = "l", col=blue_env, ylim=c(0,1100),
     las =1,cex.main=1,cex.axis=1, ylab = "Light intensity (PAR)") 
lines(PAR.algae$Raw_time, PAR.algae$Lightintensity_light_chamber*0.03916453 , type = "l", col=yellow_light) 
lines(PAR.algae$Raw_time, PAR.algae$Lightintensity_dark_chamber*0.03916453 , type = "l", col=red_dark) 

# par(mfrow = c(2,5))
# par(mar = c(3,6,4,1))

for (i in 4:8){
  y = Fucus[,i]
  x = Fucus$Treatment
  cols = unique(Fucus$colour)
  
  if (names(Fucus[i]) == "c"){
    title = "Humic-like
DOM"
  } else if (names(Fucus[i]) == "CDOM"){
    title = "CDOM
254 nm"
    # } else if (names(Fucus[i]) == "O2"){
    #   title = expression('O'[2])
    #Didn't get O2 title to work with subscript, because it un-bolded it for some reason.
  } else {
    title = names(Fucus[i])
  }
  fucus = data.frame(x, y)
  if (names(Fucus[i])=="DOC"){
    boxplot(fucus$y ~ fucus$x, main = title,ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1)
  }
  else if (names(Fucus[i])=="TDN"){
    boxplot(fucus$y ~ fucus$x, main = title,ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1)
  }
  else if (names(Fucus[i])=="c"){
    boxplot(fucus$y ~ fucus$x, main = title,ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1)
  }
  else if (names(Fucus[i])=="CDOM"){
    boxplot(fucus$y ~ fucus$x, main = title,ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1)
  }
  else{
    boxplot(fucus$y ~ fucus$x, main = title,ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1)
  }
  title(ylab = ylabs[i-3,1], mgp = c(3.8, 1, 0), cex = 1)
  
  mylevels <- levels(fucus$x)
  levelProportions <- summary(fucus$x)/nrow(fucus)
  
  for(j in 1:length(mylevels)){
    thislevel <- mylevels[j]
    thisvalues <- fucus[fucus$x==thislevel, "y"]
    print(thisvalues)
    
    # take the x-axis indices and add a jitter, proportional to the N in each level
    myjitter <- jitter(rep(j, length(thisvalues)), amount=levelProportions[j]/2)
    points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9), cex = 1.5)
    
  }
  
}


#Seagrass figures
Zostera = data.frame(subset(init.v.end, init.v.end$Organism == "Seagrass"))

for (i in 4:8){
  y = Zostera[,i]
  x = Zostera$Treatment
  cols = unique(Zostera$colour)
  
  if (names(Zostera[i]) == "c"){
    title = "Humic-like
DOM"
  } else if (names(Zostera[i]) == "CDOM"){
    title = "CDOM
254 nm"
    # } else if (names(Zostera[i]) == "O2"){
    #   title = expression('O'[2])
    #Didn't get O2 title to work with subscript, because it un-bolded it for some reason.
  } else {
    title = names(Zostera[i])
  }
  zostera = data.frame(x, y)
  if (names(Zostera[i])=="DOC"){
    boxplot(zostera$y ~ zostera$x, main = "",ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1)
  }
  else if (names(Zostera[i])=="TDN"){
    boxplot(zostera$y ~ zostera$x, main = "",ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1)
  }
  else if (names(Zostera[i])=="c"){
    boxplot(zostera$y ~ zostera$x, main = "",ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1)
  }
  else if (names(Zostera[i])=="CDOM"){
    boxplot(zostera$y ~ zostera$x, main = "",ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1)
  }
  else{
    boxplot(zostera$y ~ zostera$x, main = "",ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1)
  }
  title(ylab = ylabs[i-3,1], mgp = c(3.8, 1, 0), cex = 1)
  
  mylevels <- levels(zostera$x)
  levelProportions <- summary(zostera$x)/nrow(zostera)
  
  for(j in 1:length(mylevels)){
    thislevel <- mylevels[j]
    thisvalues <- zostera[zostera$x==thislevel, "y"]
    print(thisvalues)
    
    # take the x-axis indices and add a jitter, proportional to the N in each level
    myjitter <- jitter(rep(j, length(thisvalues)), amount=levelProportions[j]/2)
    points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9), cex = 1.5)
    
  }
  
}


# Production figures -----------------------------------------------------

#Multiply each column to get to the exact mg/ug in the bags
volume_multiplied = chamber
for(i in 5:18) {
  # volume_multiplied[,i] <- chamber[, 4] * chamber[, i]
  
  volume_multiplied[,i] <- 1 * chamber[, i]
}
# write.csv(volume_multiplied, "volume_multiplied.csv")


#Then calculate how much the algae produced of the parameter in the bag
names(volume_multiplied)
production <- volume_multiplied %>%
  dplyr::group_by(Chamber,Organism, Treatment, Volume, Plant.biomass,New.biomass,Old.biomass) %>%  
  dplyr::mutate(O2 =  End.O2-Initial.O2.original,
         DOC =  End.WC.DOC-Initial.DOC.original,
         TDN =  End.DON-Initial.DON.original,
         c =  c.end-Initial.c.original,
         CDOM =  CDOM.end-Initial.CDOM.original) %>% 
  dplyr::select(-c("Initial.O2.original","End.O2",
            "Initial.DOC.original","End.WC.DOC",
            "Initial.DON.original","End.DON","Initial.c.original","c.end",
            "Initial.CDOM.original","CDOM.end"))

# production.per.gram = production
# for(i in 8:length(production)) {
#   production.per.gram[,i] <- (production[, i] / production[, 5])/5
# }


####Figure 2####
Fucus.production = data.frame(subset(production, production$Organism == "Algae.unfiltered"))
Fucus.production$Treatment = as.factor(Fucus.production$Treatment)
Fucus.production <- Fucus.production %>%
  select_if(~ !any(is.na(.)))

colours = scico(10, palette = 'vikO', alpha=0.7, begin=0.15, end=0.85)
colour <- function(x) {
  if(x == "light") return(colours[2])
  if (x == "dark") return(colours[7])
}
Fucus.production$colour = NA
Fucus.production$colour <- sapply(Fucus.production$Treatment,colour)

ylabs = data.frame(matrix(nrow = 5))

ylabs[1,1] = expression('O'[2]~(mu*'mol'~L^-1~h^-1))
ylabs[2,1] = expression('C'~('mg'~L^-1~h^-1))
ylabs[3,1] = expression('N'~('mg'~L^-1~h^-1))
ylabs[4,1] = expression('Fluorescence'~L^-1~(h^-1))
ylabs[5,1] = expression('Fluorescence'~L^-1~(h^-1))


par(mfrow = c(2,5))
par(mar = c(2,6,4,1))

for (i in 8:12){
  x = Fucus.production[,i]
  y = Fucus.production$Treatment
  cols = unique(Fucus.production$colour)
  
  if (names(Fucus.production[i]) == "c"){
    title = "Humic-like
DOM"
  } else if (names(Fucus.production[i]) == "CDOM"){
    title = "CDOM
254 nm"
    # } else if (names(Fucus.production[i]) == "O2"){
    #   title = expression('O'[2])
    #Didn't get O2 title to work with subscript, because it un-bolded it for some reason.
  } else {
    title = names(Fucus.production[i])
  }
  fucus = data.frame(x, y)
  if (names(Fucus.production[i])=="DOC"){
    boxplot(fucus$x ~ fucus$y, main = title,ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1, ylim=c(0, max(x)*1.1))
  }
  else if (names(Fucus.production[i])=="TDN"){
    boxplot(fucus$x ~ fucus$y, main = title,ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1, ylim=c(0, max(x)*1.1))
  }
  else if (names(Fucus.production[i])=="c"){
    boxplot(fucus$x ~ fucus$y, main = title,ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1, ylim=c(0, max(x)*1.1))
  }
  else if (names(Fucus.production[i])=="CDOM"){
    boxplot(fucus$x ~ fucus$y, main = title,ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1, ylim=c(0, max(x)*1.1))
  }
  else{
    boxplot(fucus$x ~ fucus$y, main = title,ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1, ylim=c(min(x)*1.1, max(x)*1.1))
  }
  abline(h = 0, lty = 2, lwd = 2)
  title(ylab = ylabs[i-7,1], mgp = c(3.8, 1, 0), cex = 1)
  
  mylevels <- levels(fucus$y)
  levelProportions <- summary(fucus$y)/nrow(fucus)
  
  for(j in 1:length(mylevels)){
    thislevel <- mylevels[j]
    thisvalues <- fucus[fucus$y==thislevel, "x"]
    
    # take the x-axis indices and add a jitter, proportional to the N in each level
    myjitter <- jitter(rep(j, length(thisvalues)), amount=levelProportions[j]/2)
    points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9), cex = 2)
    
  }
  
}



#3 x 15

#Seagrass#
Zostera.production = data.frame(subset(production, production$Organism == "Seagrass"))
Zostera.production = na.omit(Zostera.production)
Zostera.production$Treatment = as.factor(Zostera.production$Treatment)
Zostera.production <- Zostera.production %>%
  select_if(~ !any(is.na(.)))

colours = scico(10, palette = 'davos')
colour <- function(x) {
  if(x == "light") return(colours[2]) 
  if (x == "dark") return(colours[7]) 
}
Zostera.production$colour = NA
Zostera.production$colour <- sapply(Zostera.production$Treatment,colour)

for (i in 8:12){
  x = Zostera.production[,i]
  y = Zostera.production$Treatment
  cols = unique(Zostera.production$colour)
  
  if (names(Zostera.production[i]) == "c"){
    title = "Humic-like
DOM"
  } else if (names(Zostera.production[i]) == "CDOM"){
    title = "CDOM
254 nm"
    # } else if (names(Zostera.production[i]) == "O2"){
    #   title = expression('O'[2])
    #Didn't get O2 title to work with subscript, because it un-bolded it for some reason.
  } else {
    title = names(Zostera.production[i])
  }
  Zostera = data.frame(x, y)
  if (names(Zostera.production[i])=="DOC"){
    boxplot(Zostera$x ~ Zostera$y, main = "",ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1, ylim=c(min(x)*1.1, max(x)*1.1))
  }
  else if (names(Zostera.production[i])=="TDN"){
    boxplot(Zostera$x ~ Zostera$y, main = "",ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1, ylim=c(0, max(x)*1.1))
  }
  else if (names(Zostera.production[i])=="c"){
    boxplot(Zostera$x ~ Zostera$y, main = "",ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1, ylim=c(min(x)*1.1, max(x)*1.1))
  }
  else if (names(Zostera.production[i])=="CDOM"){
    boxplot(Zostera$x ~ Zostera$y, main = "",ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1, ylim=c(min(x)*1.1, 0))
  }
  else{
    boxplot(Zostera$x ~ Zostera$y, main = "",ylab = "", xlab = NULL, col = cols, las =1,
            cex.main=1,cex.axis=1, ylim=c(min(x)*1.1, 0))
  }
  abline(h = 0, lty = 2, lwd = 2)
  title(ylab = ylabs[i-7,1], mgp = c(3.8, 1, 0), cex = 1)
  
  mylevels <- levels(Zostera$y)
  levelProportions <- summary(Zostera$y)/nrow(Zostera)
  
  for(j in 1:length(mylevels)){
    thislevel <- mylevels[j]
    thisvalues <- Zostera[Zostera$y==thislevel, "x"]
    
    # take the x-axis indices and add a jitter, proportional to the N in each level
    myjitter <- jitter(rep(j, length(thisvalues)), amount=levelProportions[j]/2)
    points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9), cex = 2)
    
  }
  
}
# 3 x 15


