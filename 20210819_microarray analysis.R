# microarray analysis
# printed in June 2021
# ELISA analysis below

rm(list=ls())
require(ggplot2)
require(tidyr)
require(reshape2)
require(wesanderson)
require(dplyr)
require(tidyverse)
require(scico)
require(Rmisc)
require(readxl)
require(ggpubr)
require(rstatix)
require(cowplot)
require(FSA)
require(lme4)
require(lmerTest)
require(emmeans)


fonts <- list(
  sans = "Helvetica",
  mono = "Consolas",
  `Times New Roman` = "DejaVu Serif"
)
font_family<-"sans"
font_size<-10

time_cat_colors<-c(scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.7),
                   scico(2, palette = "vikO",begin=0.25,end=0.75, alpha=0.7))
time_cat_4colors<-c("grey",scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.7),
                   scico(2, palette = "vikO",begin=0.25,end=0.75, alpha=0.7))
time_cat_5colors<-c("white", "grey",scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.7),
                    scico(2, palette = "vikO",begin=0.25,end=0.75, alpha=0.7))

path_data<-"C:/Users/admin/ownCloud/Assemble/Data/microarray/20210712_microarraydata_forR.xlsx"

sheets<-excel_sheets(path_data)

raw<-read_excel(path=path_data, sheet = sheets[1])
names<-read_excel(path=path_data, sheet=sheets[2])
head(raw)

sample_names<-character()

for (line in seq(1,length(names$`Sample number`)*2-1,2)){
  sample_names[line]<-names$`Sample name`[(line+1)/2]
  sample_names[line+1]<-names$`Sample name`[(line+1)/2]
}

raw$Sample_name<-NA
raw$Sample_name[1:116]<-sample_names


meta_temp=data.frame(plant=character(),
                     timepoint=character(),
                     treatment=character(), 
                     time_treat=character())

plant_temp=character()
timepoint_temp=character()
treatment_temp=character()
time_treat_temp=character()

for (i in 1:length(raw$Sample)){
  if (!is.na(raw$Sample_name[i])){
    if (grepl(raw$Sample_name[i], pattern="dialysis_blank")){
      plant_temp<-"none"
      timepoint_temp<-"none"
      treatment_temp<-"dialysis_blank"
      time_treat_temp<-"dialysis_blank"
    }
  if (grepl(raw$Sample_name[i], pattern = "-")){
    name_parts_temp<-unlist(str_split(raw$Sample_name[i], "-"))
    if (name_parts_temp[1]=="A"){
      plant_temp<-"Fucus"  
    }
    
    if (name_parts_temp[1]=="S"){
      plant_temp<-"Zostera"
    }
    
    if (grepl(name_parts_temp[2], pattern="I$")){
      timepoint_temp<-"initial"
      time_treat_temp<-"initial"
    }
    
    if (grepl(name_parts_temp[2], pattern="E$")){
      timepoint_temp<-"end"
      name_subparts_temp<-unlist(strsplit(name_parts_temp[2], split=""))
      if(name_subparts_temp[3]<=3){
        treatment_temp="light"
        time_treat_temp<-"end_light"
      }
      if(name_subparts_temp[3]>=4){
        treatment_temp="dark"
        time_treat_temp<-"end_dark"
      }
      
    }
  }
  }
  else {
    if (grepl(raw$Sample[i], pattern = "^S_")){
      standard_parts_temp<-unlist(str_split(raw$Sample[i], "_"))
      plant_temp<-standard_parts_temp[2]
      timepoint_temp<-"none"
      treatment_temp<-"none"
      time_treat_temp<-"none"
    }
    else{
    plant_temp<-"none"
    timepoint_temp<-"none"
    treatment_temp<-"none"
    time_treat_temp<-"none"
    }
  }
  meta_temp<-rbind(meta_temp, c(as.character(plant_temp), as.character(timepoint_temp), 
                                as.character(treatment_temp), as.character(time_treat_temp)))
  colnames(meta_temp)<-c("plant", "timepoint", "treatment", "time_treat_temp")
}

raw$plant<-meta_temp$plant
raw$treatment<-meta_temp$treatment
raw$timepoint<-meta_temp$timepoint
raw$time_treat<-meta_temp$time_treat_temp


melted_raw<-melt(raw, measure.vars = colnames(raw)[3:9])

melted_raw$treatment[melted_raw$Sample_name=="seagrass_blank"]<-"none"
melted_raw$treatment[melted_raw$Sample_name=="algae_blank"]<-"none"
melted_raw$time_treat[melted_raw$Sample_name=="seagrass_blank"]<-"extraction_blank"
melted_raw$time_treat[melted_raw$Sample_name=="algae_blank"]<-"extraction_blank"
melted_raw$timepoint[melted_raw$Sample_name=="seagrass_blank"]<-"none"
melted_raw$timepoint[melted_raw$Sample_name=="algae_blank"]<-"none"
melted_raw$time_treat[which(melted_raw$Sample_name%in%c("BlankHB1", "BlankHB2", "BlankHB3",
  "BlankHB4", "BlankHB5", "BlankHB6"))]<-"extraction_blank"
melted_raw$plant[which(melted_raw$Sample_name%in%c("BlankHB1", "BlankHB2", "BlankHB3",
  "BlankHB4", "BlankHB5", "BlankHB6"))]<-"extraction_blank"

# exploratory plotting ----------------------------------------------------

unique(melted_raw$time_treat)

melted_raw$time_treat<-factor(melted_raw$time_treat, levels=c("none",
                                                              "dialysis_blank",
                                                              "extraction_blank",
                                                              "initial", 
                                                              "end_dark",
                                                              "end_light"
                                                              ),
                              ordered=T)

ggplot(data=melted_raw[!is.na(melted_raw$`Concentration, µg/ml`),])+
  geom_point(aes(x=as.numeric(`Concentration, µg/ml`), y=value, color=plant))+
  facet_wrap(facets = "variable")

ggplot(data=melted_raw[melted_raw$`Concentration, µg/ml`=="NA" & melted_raw$time_treat!="dialysis_blank",])+
  geom_boxplot(aes(x=variable, y=value, color=time_treat), outlier.shape = NA)+
  geom_point(aes(x=variable, y=value, color=time_treat), position = position_dodge2(0.8))+
  facet_wrap(facets = "plant", nrow = 3)+
  theme_light()


tiff(filename="C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210712_Fig3c_microarray.tiff",
     height=12, width=20, units = "cm", res=300)
ggplot(data=melted_raw[melted_raw$`Concentration, µg/ml`=="NA" & melted_raw$plant!="Zostera",])+
  geom_boxplot(aes(x=variable, y=value, fill=time_treat), pch=21, outlier.shape = NA)+
  geom_point(aes(x=variable, y=value, fill=time_treat),
             pch=21, position = position_dodge2(0.8), size=1.5)+
  scale_fill_manual(name="Treatment", labels=c("blank\ndialysis", "blank\nextraction", "initial", "dark", "light"),
                    values = time_cat_5colors)+
  scale_x_discrete(name="", labels=c("BAM1 #1\n(Fucoidan)", 
                                     "BAM1 #2\n(Fucoidan)",
                                     "BAM1 #25\n(Fucoidan)",
                                     "BAM1 #26\n(Fucoidan)",
                                     "BAM2\n(Fucoidan)",
                                     "BAM7\n(Alginate)",
                                     "Anti-rat\n(Control 2ary antibody)"))+
  scale_y_continuous(name="Fluorescence")+
#  facet_wrap(facets = "plant", nrow = 3)+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.y = element_rect(fill = "#ECECEC"))
dev.off()


tiff(filename="C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210712_Fig3c_microarray_Zostera.tiff",
     height=12, width=20, units = "cm", res=300)
ggplot(data=melted_raw[melted_raw$`Concentration, µg/ml`=="NA" & melted_raw$plant!="Zostera",])+
  geom_boxplot(aes(x=variable, y=value, fill=time_treat), pch=21, outlier.shape = NA)+
  geom_point(aes(x=variable, y=value, fill=time_treat),
             pch=21, position = position_dodge2(0.8), size=1.5)+
  scale_fill_manual(name="Treatment", labels=c("blank", "initial", "dark", "light"),
                    values = time_cat_5colors)+
  scale_x_discrete(name="", labels=c("BAM1 #1\n(Fucoidan)", 
                                     "BAM1 #2\n(Fucoidan)",
                                     "BAM1 #25\n(Fucoidan)",
                                     "BAM1 #26\n(Fucoidan)",
                                     "BAM2\n(Fucoidan)",
                                     "BAM7\n(Alginate)",
                                     "Anti-rat\n(Control 2ary antibody)"))+
  scale_y_continuous(name="Fluorescence")+
  #  facet_wrap(facets = "plant", nrow = 3)+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.y = element_rect(fill = "#ECECEC"))
dev.off()



# correlate with monosaccharides ------------------------------------------

alldata<-merge(triplicate_data_corr_means, melted_raw, by="Sample_name", all=T)

ggplot(data=alldata[alldata$Plant=="Zostera" & alldata$variable=="BAM7 (Alginate)" & !is.na(alldata$time_treat),])+
  geom_point(aes(x=mean_conc_corr, y=value, color=Monomer), size=2)+
  scale_y_continuous(name="Antibody normalized antibody signal")+
  scale_x_continuous(name="Monosaccharide concentration (microg per L)")+
  facet_wrap("time_treat", nrow = 1)+
  theme_bw()



# calculate means of duplicate prints -------------------------------------

melted_means<-melted_raw %>%
  dplyr::group_by(Sample_name, variable, plant, treatment, timepoint, time_treat)%>%
  dplyr::summarise(mean_duplicates=mean(value))

# correlate means with monosaccharides
allmeans<-merge(triplicate_data_corr_means, melted_means, by="Sample_name", all = T)

allmeans$time_treat<-factor(allmeans$time_treat, levels = c("none", "dialysis_blank",
                                                            "extraction_blank", "initial",
                                                            "end_dark", "end_light"),
                            ordered = T)

# group low abundant sugars under "other" 

max_vals<-allmeans%>%
  dplyr::group_by(Monomer)%>%
  dplyr::summarise(max_int=max(mean_env))
allmeans$sugars_plot<-as.character(0)

max_vals<-max_vals[!is.na(max_vals$Monomer),]
allmeans<-allmeans[!is.na(allmeans$Monomer),]

for (i in 1:length(allmeans$Monomer)){
  if (max_vals$max_int[max_vals$Monomer==allmeans$Monomer[i]]<170){
    allmeans$sugars_plot[i]<-"other"
  }
  else{
    allmeans$sugars_plot[i]<-as.character(allmeans$Monomer[i])
  }
}
allmeans$sugars_plot<-factor(allmeans$sugars_plot,
                         levels=c("Fucose", "Galactose", "Glucose", "Mannose", "Rhamnose", "Xylose", "other"),
                         ordered=T)


formula<-y~x
stats_text_size<-3

Fuc_BAM1<-filter(allmeans, Plant=="Fucus", variable=="BAM1 #1 (Fucoidan)") %>% 
  ggplot()+
  geom_point(aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot), size=2, pch=21,
             alpha=0.7)+
  scale_x_continuous(name="BAM1 normalized antibody signal")+
  scale_fill_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_color_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_y_continuous(name=expression(paste("Monosaccharide concentration (", mu,"g",  L^-1, ")")), 
                     limits = c(-300,2400))+
  annotate("text", x=5.5, y=-200, label = "italic(Fucus~vesiculosus)", parse=T)+
  geom_smooth(method="lm", aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, color=sugars_plot), 
              formula=formula, alpha=0.3)+
  stat_fit_glance(method="lm", method.args = list(formula=formula), geom="text", 
                  aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, 
                      label = sprintf('R²=%.3f p=%.2g',stat(r.squared), stat(p.value))), 
                  color=scico(7, alpha=1, end=0.8, palette = "vikO"),
                  label.x=3, label.y=seq(2400, 1000, -200), size=stats_text_size)+
#  ggtitle("Fucus, BAM1")+
  #  facet_wrap("time_treat", ncol = 3)+
  theme_bw()+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))

Fuc_BAM2<-filter(allmeans, Plant=="Fucus", variable=="BAM2 (Fucoidan)") %>% 
  ggplot()+
  geom_point(aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot), size=2, pch=21,
             alpha=0.7)+
  scale_x_continuous(name="BAM2 normalized antibody signal")+
  scale_y_continuous(name="", limits = c(-300,2400))+
  scale_fill_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_color_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  annotate("text", x=4, y=-200, label = "italic(Fucus~vesiculosus)", parse=T)+
  geom_smooth(method="lm", aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, color=sugars_plot), 
              formula=formula, alpha=0.3)+
  stat_fit_glance(method="lm", method.args = list(formula=formula), geom="text", 
                  aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, 
                      label = sprintf('R²=%.3f p=%.2g',stat(r.squared), stat(p.value))), 
                  color=scico(7, alpha=1, end=0.8, palette = "vikO"),
                  label.x=2, label.y=seq(2400, 1000, -200), size=stats_text_size)+
#  ggtitle("Fucus, BAM2")+
  #  facet_wrap("time_treat", ncol = 3)+
  theme_bw()+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))

Zos_BAM7<-filter(allmeans, Plant=="Zostera", variable=="BAM7 (Alginate)") %>% 
  ggplot()+
  geom_point(aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot), size=2, pch=21, 
             alpha=0.7)+
  scale_x_continuous(name="BAM7 normalized antibody signal")+
  scale_y_continuous(name="", limits = c(-150,1200))+
  scale_fill_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_color_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  annotate("text", x=15, y=-90, label = "italic(Zostera~marina)", parse=T)+
  geom_smooth(method="lm", aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, color=sugars_plot), 
              formula=formula, alpha=0.3, 
              show.legend=T)+
  stat_fit_glance(method="lm", method.args = list(formula=formula), geom="text", 
                  aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, 
                      label = sprintf('R²=%.3f p=%.2g',stat(r.squared), stat(p.value))),
                  label.x=8, label.y=seq(1200, 500, -100), color=scico(7, alpha=1, end=0.8, palette = "vikO"), 
                  size=stats_text_size, family=font_family
                  )+
#  ggtitle("Zostera, BAM7")+
  #  facet_wrap("time_treat", ncol = 3)+
  theme_bw()+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = "right",
        strip.background.y = element_rect(fill = "#ECECEC"))

Zos_BAM7_sq<-filter(allmeans, Plant=="Zostera", variable=="BAM7 (Alginate)") %>% 
  ggplot()+
  geom_point(aes(y=mean_conc_corr, x=mean_duplicates^2, fill=sugars_plot), size=2, pch=21, 
             alpha=0.7)+
  scale_x_continuous(name="BAM7 normalized antibody signal")+
  scale_y_continuous(name="", limits = c(-150,1200))+
  scale_fill_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_color_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  annotate("text", x=600, y=-90, label = "italic(Zostera~marina)", parse=T)+
  geom_smooth(method="lm", aes(y=mean_conc_corr, x=mean_duplicates^2, fill=sugars_plot, color=sugars_plot), 
              formula=formula, alpha=0.3, 
              show.legend=T)+
  stat_fit_glance(method="lm", method.args = list(formula=formula), geom="text", 
                  aes(y=mean_conc_corr, x=mean_duplicates^2, fill=sugars_plot, 
                      label = sprintf('R²=%.3f p=%.2g',stat(r.squared), stat(p.value))),
                  label.x=300, label.y=seq(1200, 500, -100), color=scico(7, alpha=1, end=0.8, palette = "vikO"), 
                  size=stats_text_size, family=font_family
  )+
  #  ggtitle("Zostera, BAM7")+
  #  facet_wrap("time_treat", ncol = 3)+
  theme_bw()+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = "right",
        strip.background.y = element_rect(fill = "#ECECEC"))

tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210713_Fig_correlation_fluo-monosac_sqZM.tiff",
     height=12, width=25, units = "cm", res=300)
plot_grid(Fuc_BAM1, Fuc_BAM2, Zos_BAM7_sq, align = "v", axis = "top", nrow=1, rel_widths = c(1,1,1.5))
dev.off()



# Plot BAM1 2 and 7 for both algae and sg ----------------------------------------------

Fuc_BAM1<-filter(allmeans, Plant=="Fucus", variable=="BAM1 #1 (Fucoidan)") %>% 
  ggplot()+
  geom_point(aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot), size=2, pch=21,
             alpha=0.7)+
  scale_x_continuous(name="BAM1 normalized antibody signal")+
  scale_fill_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_color_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_y_continuous(name=expression(paste("Monosaccharide concentration (", mu,"g",  L^-1, ")")), 
                     limits = c(-300,2400))+
  annotate("text", x=5.5, y=-200, label = "italic(Fucus~vesiculosus)", parse=T)+
  geom_smooth(method="lm", aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, color=sugars_plot), 
              formula=formula, alpha=0.3)+
  stat_fit_glance(method="lm", method.args = list(formula=formula), geom="text", 
                  aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, 
                      label = sprintf('R²=%.3f p=%.2g',stat(r.squared), stat(p.value))), 
                  color=scico(7, alpha=1, end=0.8, palette = "vikO"),
                  label.x=3, label.y=seq(2400, 1000, -200), size=stats_text_size)+
  #  ggtitle("Fucus, BAM1")+
  #  facet_wrap("time_treat", ncol = 3)+
  theme_bw()+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))

Fuc_BAM2<-filter(allmeans, Plant=="Fucus", variable=="BAM2 (Fucoidan)") %>% 
  ggplot()+
  geom_point(aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot), size=2, pch=21,
             alpha=0.7)+
  scale_x_continuous(name="BAM2 normalized antibody signal")+
  scale_y_continuous(name="", limits = c(-300,2400))+
  scale_fill_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_color_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  annotate("text", x=4, y=-200, label = "italic(Fucus~vesiculosus)", parse=T)+
  geom_smooth(method="lm", aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, color=sugars_plot), 
              formula=formula, alpha=0.3)+
  stat_fit_glance(method="lm", method.args = list(formula=formula), geom="text", 
                  aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, 
                      label = sprintf('R²=%.3f p=%.2g',stat(r.squared), stat(p.value))), 
                  color=scico(7, alpha=1, end=0.8, palette = "vikO"),
                  label.x=2, label.y=seq(2400, 1000, -200), size=stats_text_size)+
  #  ggtitle("Fucus, BAM2")+
  #  facet_wrap("time_treat", ncol = 3)+
  theme_bw()+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))

Fuc_BAM7<-filter(allmeans, Plant=="Fucus", variable=="BAM7 (Alginate)") %>% 
  ggplot()+
  geom_point(aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot), size=2, pch=21, 
             alpha=0.7)+
  scale_x_continuous(name="BAM7 normalized antibody signal")+
  scale_y_continuous(name="", limits = c(-300,2400))+
  scale_fill_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_color_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  annotate("text", x=5, y=-200, label = "italic(Fucus~vesiculosus)", parse=T)+
  geom_smooth(method="lm", aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, color=sugars_plot), 
              formula=formula, alpha=0.3, 
              show.legend=T)+
  stat_fit_glance(method="lm", method.args = list(formula=formula), geom="text", 
                  aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, 
                      label = sprintf('R²=%.3f p=%.2g',stat(r.squared), stat(p.value))),
                  label.x=2, label.y=seq(2400, 1000, -200), color=scico(7, alpha=1, end=0.8, palette = "vikO"), 
                  size=stats_text_size, family=font_family
  )+
  #  ggtitle("Zostera, BAM7")+
  #  facet_wrap("time_treat", ncol = 3)+
  theme_bw()+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = "right",
        strip.background.y = element_rect(fill = "#ECECEC"))


Zos_BAM1<-filter(allmeans, Plant=="Zostera", variable=="BAM1 #1 (Fucoidan)") %>% 
  ggplot()+
  geom_point(aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot), size=2, pch=21,
             alpha=0.7)+
  scale_x_continuous(name="BAM1 normalized antibody signal")+
  scale_fill_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_color_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_y_continuous(name=expression(paste("Monosaccharide concentration (", mu,"g",  L^-1, ")")), 
                     limits = c(-150,1200))+
  annotate("text", x=0.8, y=-90, label = "italic(Zostera~marina)", parse=T)+
  geom_smooth(method="lm", aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, color=sugars_plot), 
              formula=formula, alpha=0.3)+
  stat_fit_glance(method="lm", method.args = list(formula=formula), geom="text", 
                  aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, 
                      label = sprintf('R²=%.3f p=%.2g',stat(r.squared), stat(p.value))), 
                  color=scico(7, alpha=1, end=0.8, palette = "vikO"),
                  label.x=0.45, label.y=seq(1200, 500, -100), size=stats_text_size)+
  #  ggtitle("Fucus, BAM1")+
  #  facet_wrap("time_treat", ncol = 3)+
  theme_bw()+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))

Zos_BAM2<-filter(allmeans, Plant=="Zostera", variable=="BAM2 (Fucoidan)") %>% 
  ggplot()+
  geom_point(aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot), size=2, pch=21,
             alpha=0.7)+
  scale_x_continuous(name="BAM2 normalized antibody signal")+
  scale_y_continuous(name="", limits = c(-150,1200))+
  scale_fill_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_color_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  annotate("text", x=1.2, y=-90, label = "italic(Zostera~marina)", parse=T)+
  geom_smooth(method="lm", aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, color=sugars_plot), 
              formula=formula, alpha=0.3)+
  stat_fit_glance(method="lm", method.args = list(formula=formula), geom="text", 
                  aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, 
                      label = sprintf('R²=%.3f p=%.2g',stat(r.squared), stat(p.value))), 
                  color=scico(7, alpha=1, end=0.8, palette = "vikO"),
                  label.x=0.6, label.y=seq(1200, 500, -100), size=stats_text_size)+
  #  ggtitle("Fucus, BAM2")+
  #  facet_wrap("time_treat", ncol = 3)+
  theme_bw()+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.y = element_rect(fill = "#ECECEC"))

Zos_BAM7<-filter(allmeans, Plant=="Zostera", variable=="BAM7 (Alginate)") %>% 
  ggplot()+
  geom_point(aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot), size=2, pch=21, 
             alpha=0.7)+
  scale_x_continuous(name="BAM7 normalized antibody signal")+
  scale_y_continuous(name="", limits = c(-150,1200))+
  scale_fill_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_color_manual(name="Monosaccharide", values=scico(7, alpha=0.7, end=0.8, palette = "vikO"))+
  annotate("text", x=16, y=-90, label = "italic(Zostera~marina)", parse=T)+
  geom_smooth(method="lm", aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, color=sugars_plot), 
              formula=formula, alpha=0.3, 
              show.legend=T)+
  stat_fit_glance(method="lm", method.args = list(formula=formula), geom="text", 
                  aes(y=mean_conc_corr, x=mean_duplicates, fill=sugars_plot, 
                      label = sprintf('R²=%.3f p=%.2g',stat(r.squared), stat(p.value))),
                  label.x=10, label.y=seq(1200, 500, -100), color=scico(7, alpha=1, end=0.8, palette = "vikO"), 
                  size=stats_text_size, family=font_family
  )+
  #  ggtitle("Zostera, BAM7")+
  #  facet_wrap("time_treat", ncol = 3)+
  theme_bw()+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = "right",
        strip.background.y = element_rect(fill = "#ECECEC"))


tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210713_Fig_correlation_fluo-monosac_6plots.tiff",
     height=14, width=24, units = "cm", res=300)
plot_grid(Fuc_BAM1, Fuc_BAM2, Fuc_BAM7, Zos_BAM1, Zos_BAM2, Zos_BAM7, 
          align = "v", axis = "top", nrow=2, rel_widths = c(1,1,1.4))
dev.off()




Fuc_all<-ggplot(data=melted_means[melted_means$time_treat!="none" & melted_means$plant!="Zostera",])+
  geom_bar(aes(x=variable, y=mean_duplicates, fill=time_treat), stat="summary",  
               position = position_dodge(width = 1), color="#262626")+
  geom_point(aes(x=variable, y=mean_duplicates, group=time_treat, fill=time_treat),
             pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
  annotate("text", x=4, y=11.8, label = "italic(Fucus~vesiculosus)", parse=T)+
  scale_fill_manual(name="Treatment", labels=c("dialysis\nblank", "extraction\nblank", "initial", "dark", "light"),
                    values = time_cat_5colors)+
  scale_x_discrete(name="", labels=c("BAM1 #1\n(Fucoidan)", 
                                     "BAM1 #2\n(Fucoidan)",
                                     "BAM1 #25\n(Fucoidan)",
                                     "BAM1 #26\n(Fucoidan)",
                                     "BAM2\n(Fucoidan)",
                                     "BAM7\n(Alginate)",
                                     "Anti-rat\n(Control 2ary antibody)"))+
  scale_y_continuous(name="Normalized antibody signal")+
  #  facet_wrap(facets = "plant", nrow = 3)+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.y = element_rect(fill = "#ECECEC"))

Zos_all<-ggplot(data=melted_means[melted_means$time_treat!="none" & melted_means$plant!="Fucus",])+
  geom_bar(aes(x=variable, y=mean_duplicates, fill=time_treat), stat="summary", color="#262626", 
               position = position_dodge(width = 1))+
  geom_point(aes(x=variable, y=mean_duplicates, group=time_treat, fill=time_treat),
             pch=21, position = position_jitterdodge(jitter.width = .5, dodge.width = 1), size=1.5)+
  annotate("text", x=4, y=29, label = "italic(Zostera~marina)", parse=T)+
  scale_fill_manual(name="Treatment", labels=c("dialysis\nblank", "extraction\nblank", "initial", "dark", "light"),
                    values = time_cat_5colors)+
  scale_x_discrete(name="", labels=c("BAM1 #1\n(Fucoidan)", 
                                     "BAM1 #2\n(Fucoidan)",
                                     "BAM1 #25\n(Fucoidan)",
                                     "BAM1 #26\n(Fucoidan)",
                                     "BAM2\n(Fucoidan)",
                                     "BAM7\n(Alginate)",
                                     "Anti-rat\n(Control 2ary antibody)"))+
  scale_y_continuous(name="Normalized antibody signal")+
  
  #  facet_wrap(facets = "plant", nrow = 3)+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.y = element_rect(fill = "#ECECEC"))

tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210714_Fig_fluo-all.tiff",
     height=17.5, width=24, units = "cm", res=300)
plot_grid(Fuc_all, Zos_all, nrow=2)
dev.off()



# ELISA -------------------------------------------------------------------

elisa_raw<-read.table("C:/Users/admin/ownCloud/Assemble/Data/microarray/20210802_ELISA_Finland.txt",
                      h=T)
head(elisa_raw)

#blank_mean_bam1<-elisa_raw %>%
#  filter(sample_type=="buffer_blank") %>%
#  summarise(mean(BAM1))

#blank_mean_bam2<-elisa_raw %>%
#  filter(sample_type=="buffer_blank") %>%
#  summarise(mean(BAM2))

#elisa_raw$BAM1_net<-elisa_raw$BAM1-blank_mean_bam1[1,1]
#elisa_raw$BAM2_net<-elisa_raw$BAM2-blank_mean_bam2[1,1]


# Fig Elisa on raw samples

elisa_plot1_data<-elisa_raw[which(elisa_raw$sample_type %in% c("raw_sample", "blank_sample", "buffer_blank")),]
elisa_plot1_data<-melt(elisa_plot1_data, measure.vars = c("BAM1", "BAM2"))
elisa_plot1_means<-elisa_plot1_data%>%
  dplyr::group_by(Sample, sample_nr, sample_type, sample_group, variable)%>%
  dplyr::summarise(intensity=mean(value))

elisa_plot1_means$sample_group<-factor(elisa_plot1_means$sample_group,
                                       levels = c("buffer_blank", "extraction_blank", "initial",
                                                  "dark_end", "light_end"),
                                       ordered = F)
elisa_model_data<-elisa_plot1_means[which(elisa_plot1_means$sample_group %in% 
                                            c("extraction_blank", "initial","dark_end","light_end")),]
noise_BAM1<-elisa_plot1_means%>%
  filter(sample_group=="buffer_blank")%>%
  filter(variable=="BAM1")%>%
  dplyr::summarise(mean_intensity=mean(intensity))

elisa_model_BAM1<-data.frame(elisa_model_data[elisa_model_data$variable=="BAM1",])
elisa_model_BAM1$net_intensity<-elisa_model_BAM1$intensity-mean(noise_BAM1$mean_intensity)
model1_BAMs<-lm(log(intensity)~sample_group, data = elisa_model_BAM1)
anova(model1_BAMs)
plot(model1_BAMs)
levels(elisa_model_BAM1$sample_group)
kruskal.test(x=elisa_model_BAM1$net_intensity, g=elisa_model_BAM1$sample_group)
dunnTest(elisa_model_BAM1$net_intensity~as.factor(elisa_model_BAM1$sample_group), method = 'bh')

noise_BAM2<-elisa_plot1_means%>%
  filter(sample_group=="buffer_blank")%>%
  filter(variable=="BAM2")%>%
  dplyr::summarise(mean_intensity=mean(intensity))

elisa_model_BAM2<-elisa_model_data[elisa_model_data$variable=="BAM2",]
elisa_model_BAM2$net_intensity<-elisa_model_BAM2$intensity-mean(noise_BAM2$mean_intensity)
model1_BAMs<-lm(log(intensity)~sample_group, data = elisa_model_BAM2)
anova(model1_BAMs)
plot(model1_BAMs)
kruskal.test(x=elisa_model_BAM2$net_intensity, g=elisa_model_BAM2$sample_group)
dunnTest(x=elisa_model_BAM2$net_intensity, g=elisa_model_BAM2$sample_group)


# mixed effects models ------

elisa_model_data<-elisa_plot1_data[which(elisa_plot1_data$sample_group %in% 
                                            c("extraction_blank", "initial","dark_end","light_end")),]

# BAM1

noise_BAM1<-elisa_plot1_means%>%
  filter(sample_group=="buffer_blank")%>%
  filter(variable=="BAM1")%>%
  dplyr::summarise(mean_intensity=mean(intensity))

elisa_model_BAM1<-data.frame(elisa_model_data[elisa_model_data$variable=="BAM1",])
elisa_model_BAM1$net_intensity<-elisa_model_BAM1$value-mean(noise_BAM1$mean_intensity)



elisa_model_BAM1$group_normalized<-NA
for (sg in elisa_model_BAM1$sample_group){
  sg_mean<-mean(elisa_model_BAM1$value[elisa_model_BAM1$sample_group==sg])
  for (i in 1:length(elisa_model_BAM1$Sample)){
    if (elisa_model_BAM1$sample_group[1]==sg){
      elisa_model_BAM1$group_normalized=elisa_model_BAM1$net_intensity/sg_mean
    }
  }
}

fmBAM1<-lmer(log(value)~sample_group+(1|sample_nr), data=elisa_model_BAM1)

plot(fmBAM1)
qqnorm(resid(fmBAM1))
anova(fmBAM1)
emmeans(fmBAM1, list(pairwise ~ sample_group), adjust = "tukey")

# BAM2

noise_BAM2<-elisa_plot1_means%>%
  filter(sample_group=="buffer_blank")%>%
  filter(variable=="BAM2")%>%
  dplyr::summarise(mean_intensity=mean(intensity))

elisa_model_BAM2<-data.frame(elisa_model_data[elisa_model_data$variable=="BAM2",])
elisa_model_BAM2$net_intensity<-elisa_model_BAM2$value-mean(noise_BAM2$mean_intensity)



elisa_model_BAM2$group_normalized<-NA
for (sg in elisa_model_BAM2$sample_group){
  sg_mean<-mean(elisa_model_BAM2$value[elisa_model_BAM2$sample_group==sg])
  for (i in 1:length(elisa_model_BAM2$Sample)){
    if (elisa_model_BAM2$sample_group[1]==sg){
      elisa_model_BAM2$group_normalized=elisa_model_BAM2$net_intensity/sg_mean
    }
  }
}

fmBAM2<-lmer(log(value)~sample_group+(1|sample_nr), data=elisa_model_BAM2)

plot(fmBAM2)
qqnorm(resid(fmBAM2))
anova(fmBAM2)
emmeans(fmBAM2, list(pairwise ~ sample_group), adjust = "tukey")


# WORKS! mixed effects on elisa - finished --------------------------------
# assumptions okay



Fig_ELISA_rawsamples<-ggplot(data=elisa_plot1_means)+
  geom_bar(aes(x=sample_group, y=intensity, fill=sample_group), stat="summary", color="#262626")+
  geom_point(aes(x=sample_group, y=intensity, fill=sample_group, group=sample_group),
             pch=21, position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size=2)+
  scale_fill_manual(name="", labels=c("buffer\nblank", "procedure\nblank", "seawater", "dark\ninc.", "light\ninc."),
                    values=time_cat_5colors)+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Antibody signal")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "#262626"),
        axis.text.x = element_blank(),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.ticks.x = element_blank(),
        strip.background.x = element_rect(fill = "#ECECEC"))+
  facet_wrap(~variable, ncol=2)
  
ann_text <- data.frame(x= c(1,2,3,4,5), y = c(0.15, 0.15, 0.35,0.55,0.55),lab = c("a", "a", "b", "c", "c"),
                         variable = factor("BAM1", levels = c("BAM1","BAM2")))
Fig_ELISA_rawsamples<-Fig_ELISA_rawsamples + geom_text(data = ann_text, aes(x=x, y=y, label = lab), 
                                                      family=font_family, size=4.5, colour="#262626")
ann_text <- data.frame(x= c(1,2,3,4,5), y = c(0.15, 0.15, 0.15, 0.45, 0.45),lab = c("a", "ab", "b", "c", "c"),
                       variable = factor("BAM2", levels = c("BAM1","BAM2")))
Fig_ELISA_rawsamples<-Fig_ELISA_rawsamples + geom_text(data = ann_text, aes(x=x, y=y, label = lab), 
                                 family=font_family, size=4.55, colour="#262626")

tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210819_Fig_ELISA.tiff", 
     height=8, width = 12, units = "cm", res=300)
Fig_ELISA_rawsamples
dev.off()

pdf("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210819_Fig_ELISA.pdf", 
     height=3, width = 5)
Fig_ELISA_rawsamples
dev.off()


# Correlate ELISA BAM1 and monosaccharides raw samples --------------------

raw_samples<-read_xls("C:/Users/admin/ownCloud/Assemble/Data/raw samples AH/20210802_monosaccharides_raw_samples_AH_results_for_R.xls")

raw_samples$Group_for_plot<-factor(raw_samples$Group_for_plot,
                                   levels=c("dialysis_blank", "procedure_blank",
                                            "initial", "dark_end", "light_end"),
                                   ordered = T)

# for merging datasets
elisa_4melt<-elisa_plot1_means[which(elisa_plot1_means$sample_group %in% 
                                      c("extraction_blank",
                                        "initial", "dark_end", "light_end")),]
raw_4melt<-raw_samples[which(raw_samples$Group_for_plot %in%
                               c("initial", "light_end", "dark_end",
                                 "procedure_blank")),]
raw_4melt<-raw_4melt[raw_4melt$Sample_group!="ANX_flowthrough",]

raw_4melt$others<-NA
for (i in 1:length(raw_4melt$nr)){
  raw_4melt$others[i]<-sum(raw_4melt[i,c(8,9,10,11,14,16,17,18,20,21)])
}
rawmelted_4melt<-melt(raw_4melt, id.vars = colnames(raw_4melt)[1:6], measure.vars = c("Fucose", "Galactose", "Xylose", "Glucuronic_acid",
                                                  "Glucose", "others"))

elisa_and_raw<-merge(rawmelted_4melt, elisa_4melt, by = "Sample")
elisasum_and_raw<-elisa_and_raw %>%
  dplyr::group_by(Sample, variable.x)%>%
  dplyr::summarise(sum_int=sum(intensity), mono_conc=mean(value))
elisasum_and_raw$variable.x<-factor(elisasum_and_raw$variable.x, 
                                    levels=c("Fucose", "Glucuronic_acid", "Xylose", 
                                             "Galactose", "Glucose", "others"),
                                    labels=c("Fucose", "Glucuronic acid", "Xylose", 
                                             "Galactose", "Glucose", "others summed"))

ELISA_BAM1_corr<- 
  ggplot(data=elisasum_and_raw)+
  
  scale_x_continuous(name="BAM1 and BAM2 summed antibody signal")+
  scale_fill_manual(name="", values=scico(6, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_color_manual(name="", values=scico(6, alpha=0.7, end=0.8, palette = "vikO"))+
  scale_y_continuous(name=expression(sqrt(paste("Monosaccharide concentration (", mu,"g",  L^-1, ")"))))+
  #annotate("text", x=5.5, y=-200, label = "italic(Fucus~vesiculosus)", parse=T)+
  geom_smooth(method="lm", aes(y=sqrt(mono_conc), x=sum_int, fill=variable.x, color=variable.x), 
              formula=formula, alpha=0.3)+
  geom_point(aes(y=sqrt(mono_conc), x=sum_int, fill=variable.x), size=2, pch=21,
             alpha=0.7)+
  #stat_fit_glance(method="lm", method.args = list(formula=formula), geom="text", 
  #                aes(y=value, x=intensity, fill=variable.x, 
  #                    label = sprintf('R²=%.3f p=%.2g',stat(r.squared), stat(p.value), stat(slope))), 
  #                color=scico(6, alpha=1, end=0.8, palette = "vikO"),
  #                label.x=0.1, label.y=seq(16000, 10000, -1000), size=stats_text_size)+
  #  ggtitle("Fucus, BAM1")+
  #  facet_wrap("time_treat", ncol = 3)+
  theme_bw()+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        panel.grid = element_blank(),
        legend.position = c(0.2, 0.75),
        strip.background.y = element_rect(fill = "#ECECEC"), 
        legend.background = element_rect(fill=NA),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "#262626"))

tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210818_bam_mono_correlation.tiff", 
     height=10, width = 12, units = "cm", res=300)
ELISA_BAM1_corr
dev.off()

pdf("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210819_bam_mono_correlation.pdf", 
    height=4, width = 4.8)
ELISA_BAM1_corr
dev.off()





# Fig Elisa on procedure samples

elisa_plot2_data<-elisa_raw[which(elisa_raw$sample_group %in% c("abc_elution_light_end",  "abc_elution_intial",    
                                                                "NaCl_elution_light_end", "NaCl_elution_initial", 
                                                                "flow_through_light_end", "flow_through_initial")),]
elisa_plot2_data$sample_nr<-unlist(strsplit(elisa_plot2_data$sample_nr,split="_"))[seq(3,54,3)]
elisa_plot2_data_add<-elisa_raw[which(elisa_raw$sample_nr %in% c("1E", "1I")),]
elisa_plot2_data<-rbind(elisa_plot2_data, elisa_plot2_data_add)

elisa_plot2_melted<-melt(elisa_plot2_data, measure.vars = c("BAM1", "BAM2"))

elisa_plot2_means<-elisa_plot2_melted %>%
  group_by(Sample, sample_nr, sample_type, sample_group, variable)%>%
  summarise(intensity=mean(value))

ggplot(elisa_plot2_means)+
  geom_bar(aes(x=sample_group, y=intensity, fill=sample_nr, group=variable), stat="summary", color="#262626", 
           position = "stack")+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Normalized antibody signal")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.ticks.x = element_blank(),
        strip.background.x = element_rect(fill = "#262626"))
  facet_wrap(~variable, nrow = 2)
