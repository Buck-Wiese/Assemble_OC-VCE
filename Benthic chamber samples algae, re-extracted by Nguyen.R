##### Benthic chamber samples algae, re-extracted by Nguyen, dialyzed, acid hydrolyzed, HPAEC-PAD


#packages ####
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

fonts <- list(
  sans = "Helvetica",
  mono = "Consolas",
  `Times New Roman` = "DejaVu Serif"
)
font_family<-"sans"
font_size<-10

# Iv: NGUYEN#S EXTRACTIONS --------------------------------------


# 1 -  data import ####

NN_data<-read_xls("C:/Users/admin/ownCloud/Assemble_OC-VCE/AEX_data/20210603_AEX_Nguyentest/20210603_AEX_Nguyentest_forR.xls")
#rownames(triplicate_data)<-triplicate_data$nr


# replace n.a. in sugar measurements with 0
for (i in 1:length(NN_data$Injection_nr)){
  for (j in 4:33){
    if (NN_data[i,j]=="n.a."){
      NN_data[i,j]<-"0"
    }
  }
}

# change sugar measurements to numeric
for (i in 4:33){
  NN_data[, i]<-as.numeric(unlist(NN_data[,i]))
}


# Turn sample names into metadata -----------------------------------------
meta_temp=data.frame(plant=character(),
                     timepoint=character(),
                     treatment=character(), 
                     time_treat=character(),
                     elution=character(),
                     stringsAsFactors = F)

plant_temp=character()
timepoint_temp=character()
treatment_temp=character()
time_treat_temp=character()
elution_temp=character()

for (i in 1:length(NN_data$Sample)){
  if (grepl(NN_data$Sample[i], pattern = "_")){
    name_parts_temp<-unlist(str_split(NN_data$Sample[i], "_"))
# select by name structure    
    if (length(name_parts_temp)==3){
      elution_temp<-name_parts_temp[3]
      if(name_parts_temp[1]=="Std"){
        plant_temp<-"Standard"
        timepoint_temp<-"none"
        treatment_temp<-"none"
        time_treat_temp<-"Standard"
      }
      if (name_parts_temp[1]=="Blank"){
      plant_temp<-"Blank"
      timepoint_temp<-"none"
      treatment_temp<-"none"
      time_treat_temp<-"Blank"
    }
      if (name_parts_temp[1]=="A"){
      plant_temp<-"Fucus"
        if (grepl(name_parts_temp[2], pattern="I$")){
        timepoint_temp<-"initial"
        time_treat_temp<-"initial"
        treatment_temp<-"none"
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
  }
  else {
    plant_temp<-"none"
    timepoint_temp<-"none"
    treatment_temp<-"none"
    time_treat_temp<-"none"
    elution_temp<-"none"
    
  }
  meta_temp<-rbind(meta_temp, c(as.character(plant_temp), as.character(timepoint_temp), 
                                as.character(treatment_temp), as.character(time_treat_temp),
                                as.character(elution_temp)))
  colnames(meta_temp)<-c("plant", "timepoint", "treatment", "time_treat","elution")
}

NN_data$Plant<-meta_temp$plant
NN_data$Timepoint<-meta_temp$timepoint  
NN_data$Treatment<-meta_temp$treatment 
NN_data$Time_treat<-meta_temp$time_treat
NN_data$Elution<-meta_temp$elution


# plotting corrected means by plant ---------------------------------------
time_cat_colors<-c(scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.8),
                   scico(1, palette = "davos",begin=0.1,end=0.2, alpha=0.8),
                   scico(1, palette = "lajolla",begin=0.2,end=0.3, alpha=0.8),
                   "black")

NN_sub<-NN_data[,c(1:18, 34:38)]

NN_melt<-melt(NN_sub, id.vars = c("Injection_nr", "Injection_name", "Sample",
                                  "Plant", "Timepoint", "Treatment", "Time_treat", "Elution"),
              measure.vars = c("Fucose_ngpermL", "Rhamnose_ngpermL", "Galactosamine_ngpermL",
                               "Arabinose_ngpermL", "Glucosamine_ngpermL", "Galactose_ngpermL",
                               "Glucose_ngpermL", "Mannose_ngpermL", "Xylose_ngpermL",
                               "Gluconicacid_ngpermL", "Muramicacid_ngpermL", "Galacturonicacid_ngpermL",
                               "Glucuronicacid_ngpermL", "Mannuronicacid_ngpermL", "Iduronicacid_ngpermL"))
colnames(NN_melt)[9]<-"Monomer"

#NN_melt$Monomer<-factor(NN_melt$Monomer,
#                                           levels=c("Fucose", "Rhamnose","Galactosamine","Arabinose","Glucosamine",
#                                                    "Galactose","Glucose","Mannose","Xylose",
#                                                    "Gluconic_acid","Muramic_acid","Galacturonic_acid",
#                                                    "Glucuronic_acid","Mannuronic_acid","Iduronic_acid"))
NN_fucus<-NN_melt[which(NN_melt$Plant %in% c("Fucus", "Blank")),]
NN_fucus$Time_treat<-factor(NN_fucus$Time_treat, levels=c("Blank", "initial", "end_dark", "end_light"))

ggplot(data=NN_fucus)+
  geom_bar(aes(x=Time_treat, y=value, fill=Elution),
           position = "stack", stat = "summary", fun="mean", colour="black")+
  geom_errorbar(aes(x=Time_treat, y=value),
                position = "identity", stat = "summary", fun="sd", colour="black")+
  scale_x_discrete(name="", breaks=c("Blank", "initial", "end_dark", "end_light"),
                   labels=c("Blanks", "Initials", "Dark\nincubations", "Light\nincubations"))+
  scale_y_continuous(name=expression(paste(Concentration~(microg~mL^-1))))+
#  scale_fill_manual(name=NULL,  breaks=c("initial", "end_dark", "end_light"),
#                    labels=c("before incubation",
#                             "incubation in dark",
#                             "incubation in light"),
#                    values = time_cat_colors)+
  facet_grid(~ Monomer)+
  theme_bw()
#NN_plot<-NN_fucus[which(NN_fucus$Monomer %in% c("Fucose_ngpermL", "Galactosamine_ngpermL","Galactose_ngpermL", "Xylose_ngpermL",
#                                                "Mannose_ngpermL")),]

NN_plot<-NN_fucus

NN_plot$monomers_plot<-as.character("character")

for (i in 1:length(NN_plot$Monomer)){
  NN_plot$monomers_plot[i]<-as.character(NN_plot$Monomer[i])
}
  
#NN_plot$monomers_plot<-factor(NN_plot$monomers_plot, levels=c("Fucose_ngpermL", "Galactosamine_ngpermL","Galactose_ngpermL", "Xylose_ngpermL",
#                                                              "Mannose_ngpermL", "Glucuronicacid_ngpermL"))
NN_plot$monomers_plot<-factor(NN_plot$monomers_plot, levels=c("Fucose_ngpermL", "Rhamnose_ngpermL", "Galactosamine_ngpermL",
                                                              "Arabinose_ngpermL", "Glucosamine_ngpermL", "Galactose_ngpermL",
                                                              "Glucose_ngpermL", "Mannose_ngpermL", "Xylose_ngpermL",
                                                              "Gluconicacid_ngpermL", "Muramicacid_ngpermL", "Galacturonicacid_ngpermL",
                                                              "Glucuronicacid_ngpermL", "Mannuronicacid_ngpermL", "Iduronicacid_ngpermL"),
                              labels=c("Fucose", "Rhamnose", "Galactosamine", "Arabinose", "Glucosamine",
                                       "Galactose", "Glucose", "Mannose", "Xylose", "Gluconic acid", "Muramic acid",
                                       "Galacturonic acid", "Glucuronic acid", "Mannuronic acid", "Iduronic acid"))

fucose_lm<-lm(value~Time_treat+Elution+monomers_plot, NN_plot)
plot(fucose_lm)

elutions_colors<-scico(2, alpha=0.5, palette = "vikO", begin=0.3, end=0.7)

stat.test <- NN_plot %>%
  group_by(monomers_plot, Elution) %>%
  t_test(value ~ Time_treat) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

Fig3a<-ggbarplot(
  NN_plot,
  x = "Time_treat", y = "value", add = "mean_se",
  fill = "Elution", palette = "Paired", 
  xlab="", ylab="Intensity", nrow=1, facet.by = "monomers_plot", color = "#262626")+
  scale_x_discrete(labels=c("blank", "\ninitial", "dark", "\nlight"))+
  scale_y_continuous(name=expression(paste(Concentration~(microg~L^-1))))+
  scale_fill_manual(name="Consecutive elutions:", breaks=c("abc", "NaCl"), 
                    labels=c("2 M ammonium bicarbonate", "5 M sodium chloride"),
                    values = elutions_colors)+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family), 
        strip.background.x = element_rect(fill = "#262626"))
  

tiff(filename = "Figures/20210617_Fig3a_allMonosaccharides.tiff", width = 50, height = 20, units = "cm", res=200)
Fig3a  
dev.off()


means_for_writing_NN<-NN_plot %>%
  dplyr::group_by(Time_treat, Elution, monomers_plot)%>%
  dplyr::summarise(mean_reps=mean(value))

NN_plot2<-separate(NN_plot, "Sample", into=c("plant", "sample_name", "elution"), sep = "_")

NN_plot2$elution_vol<-25
NN_plot2$elution_vol[NN_plot2$Elution=="abc"]<-10
NN_plot2$extr_dial_corr<-NN_plot2$value*NN_plot2$elution_vol/100/10
NN_plot2$Elution<-factor(NN_plot2$Elution, levels=c("NaCl", "abc"), ordered = T)


#NN_plot<-NN_fucus[which(NN_fucus$Monomer %in% c("Fucose_ngpermL", "Galactosamine_ngpermL","Galactose_ngpermL", "Xylose_ngpermL",
#                                                "Mannose_ngpermL")),]
tiff(filename = "Figures/20210701_Fig3a_Monosaccharides.tiff", width = 18, height = 13, units = "cm", res=200)
ggbarplot(
  NN_plot2[which(NN_fucus$Monomer %in% c("Fucose_ngpermL", "Galactosamine_ngpermL","Galactose_ngpermL", "Xylose_ngpermL",
                                                                                         "Mannose_ngpermL")),],
  x = "Time_treat", y = "extr_dial_corr", add = "mean_se",
  fill = "Elution", palette = "Paired", 
  xlab="", ylab="Intensity", nrow=1, facet.by = "monomers_plot", color = "#262626")+
  scale_x_discrete(labels=c("blank", "\ninitial", "dark", "\nlight"))+
  scale_y_continuous(name=expression(paste(Concentration~(microg~L^-1))))+
  scale_fill_manual(name="Consecutive elutions:", breaks=c("abc", "NaCl"), 
                    labels=c("2 M ammonium bicarbonate", "5 M sodium chloride"),
                    values = elutions_colors)+
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family), 
        strip.background.x = element_rect(fill = "#262626"))
dev.off()

NN_summed<-NN_plot2 %>%
  dplyr::group_by(sample_name, monomers_plot, Time_treat)%>%
  dplyr::summarise(bothelutions=sum(extr_dial_corr), n=n())

NN_summed$ext_dial_corr<-NN_summed$bothelutions/10/10
NN_summed$ext_dial_corr[NN_summed$]

NN_summary<-NN_summed %>%
  dplyr::group_by(Time_treat, monomers_plot) %>%
  dplyr::summarise(means=mean(ext_dial_corr), SD=sd(ext_dial_corr), n=n())

