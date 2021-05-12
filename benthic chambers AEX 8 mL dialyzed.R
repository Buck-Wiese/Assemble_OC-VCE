### first analysis of extracts from incubations
#8 mL dialyzed from anionic extraction

#packages ####
require(ggplot2)
require(tidyr)
require(reshape2)
require(wesanderson)
require(dplyr)
require(tidyverse)
require(scico)
require(Rmisc)

# II. benthic chambers ####

# 1 -  data import ####
#hydrolysate from 4 different column types on HPAEC-PAD Nov 2020

chamber_data<-read.table("C:/Users/admin/ownCloud/Assemble/Data/AEX/20210204_Monosaccharide_results_AEX_incubations_forR_according_to_tube_labels.txt", 
                         h=T, stringsAsFactors = F, colClasses = character())

temp_raw<-sapply(chamber_data[,4:18], as.numeric)
chamber_data$Amount_ugperL_sum<-rowSums(temp_raw)

temp_colname<-NULL

for (i in (1:length(colnames(chamber_data)))){
  if (grepl("Amount_ugperL_",colnames(chamber_data)[i])){
    temp_colname<-unlist(strsplit(colnames(chamber_data)[i],"_"))
    colnames(chamber_data)[i]<-temp_colname[3]
  }
}

#2 - exploration ####

#2.1 blanks and standards ####
ggplot(data=chamber_data[which(chamber_data$Sample%in% c("Blank","Standard")), ])+
  geom_point(aes(x=Sample,y=sum))


blksnstds_chambers<-chamber_data[which(chamber_data$Sample%in% c("Blank","Standard")), ]

blksnstds_melt_chambers<-melt(blksnstds_chambers,id.vars = c("Number","Standard_concentration","Sample","sum"))

ggplot(data=blksnstds_melt_chambers)+
  geom_point(aes(x=variable,y=as.numeric(value),fill=Sample), 
             pch=21, position=position_dodge(0.3), size=2.5)+
  theme_bw()

ggplot()+
  geom_point(data=blksnstds_melt_chambers[blksnstds_melt_chambers$Sample=="Standard",] , 
             aes(x=Standard_concentration,y=as.numeric(value),fill=variable),pch=21, size=3)+
  geom_point(data=blksnstds_melt_chambers[blksnstds_melt_chambers$Sample=="Blank",],
             aes(x=1, y=as.numeric(value), fill=variable),pch=21, size=3)+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()


#2.2 dialysed samples ####

dialysed_chambers<-chamber_data[chamber_data$Sample!="Blank",]
dialysed_chambers<-dialysed_chambers[dialysed_chambers$Sample!="Standard",]
dialysed_chambers<-dialysed_chambers[,c(1,3:19)]


#2.2.1 controls ####
control_colors<-scico(4,palette = "bilbao",end=0.8,alpha=0.8)

dialysed_controls<-dialysed_chambers[which(dialysed_chambers$Sample %in% c("Dialysis_blank","BlankHB1",
                                                                           "BlankHB2", "BlankHB3", "BlankHB4",
                                                                           "BlankHB5","BlankHB6",
                                                                           "LowSTD1", "LowSTD2","LowSTD3",
                                                                           "HighSTD1", "HighSTD2", "HighSTD3",
                                                                           "maybe_extraction_blank")),]
dialysed_controls$Type<-c("Dialysis_blank","Extraction_blank",
                          rep("Dialysis_blank",5),rep("Extraction_blank",2),"Dialysis_blank",
                          rep("Extraction_blank",3),rep("Low_extraction_standard",3),
                          rep("High_extraction_standard",3), rep("Dialysis_blank",4))

dialysed_controls_melted<-melt(dialysed_controls, id.vars = c("Number", "Sample", "sum", "Type"))
dialysed_controls_melted$Type<-factor(dialysed_controls_melted$Type, 
                                      levels = c("Dialysis_blank","Extraction_blank",
                                                 "Low_extraction_standard","High_extraction_standard"))


SE_dialysed_controls_melted<-summarySE(dialysed_controls_melted, measurevar = "value",
                                       groupvars = c("Type", "variable"))

# plotting controls -------------------------------------------------------

# boxplot controls

ggplot(data=dialysed_controls_melted)+
  geom_boxplot(aes(x=variable, y=value, fill=Type))+
  scale_x_discrete(name="Monosaccharide")+
  scale_y_continuous(name=expression(paste(Concentration~(?g~L^-1))))+
  scale_fill_manual(name=NULL, labels=c("Dialysis blanks",
                                        "Extraction blanks",
                                        "Low-concentrated standards",
                                        "High-concentrated standards"),
                    values=control_colors)+
  theme_bw()

# barplot controls

tiff("C:/Users/admin/ownCloud/Assemble/Data/AEX/Monosaccharides/20210204_monosaccharides_barplot_controls.tiff",
     width=36, height=15, units = "cm", res=300)
ggplot(data=dialysed_controls_melted)+
  geom_bar(aes(x=variable, y=value, fill=Type),
           position=position_dodge(), stat= "summary",summary.fun.y="mean" , colour="black")+
  geom_point(aes(y=value, x=variable, fill=Type), 
             position = position_dodge(.9), pch=21, size=2)+
  scale_x_discrete(name="Monosaccharide", 
                   labels=c("fucose", "rhamnose","galactosamine","arabinose","glucosamine",
                            "galactose","glucose","mannose\n(bad separation)","xylose",
                            "gluconic\nacid","muramic\nacid","galacturonic\nacid",
                            "glucuronic\nacid","mannuronic\nacid","iduronic\nacid"))+
  scale_y_continuous(name=expression(paste(Concentration~(?g~L^-1))))+
  scale_fill_manual(name=NULL, labels=c("Dialysis blanks",
                                        "Extraction blanks",
                                        "Low-concentrated standards",
                                        "High-concentrated standards"),
                    values=control_colors)+
  theme_bw()

dev.off()

#2.2.2 metadata for samples ####

dialysed_samples<-dialysed_chambers %>%
  filter(!Sample %in% c("Dialysis_blank","BlankHB1",
                        "BlankHB2", "BlankHB3", "BlankHB4",
                        "BlankHB5","BlankHB6",
                        "LowSTD1", "LowSTD2","LowSTD3",
                        "HighSTD1", "HighSTD2", "HighSTD3",
                        "erroneous_run", "maybe_extraction_blank",
                        "maybe_extraction_blank_smalltube"))

temp_waterbody<-NULL
temp_plant<-NULL
temp_condition<-NULL
temp_timing<-NULL
temp_replicate<-NULL
temp_enriched13C<-NULL
ALGAE_temp<-NULL
EXP_temp<-NULL
EXP_temp2<-NULL
REP_COND_temp<-NULL

raw_metadata<-data.frame(plant=character(),
                         waterbody=character(),
                         condition=character(),
                         timing=character(),
                         replicate=character(),
                         enriched13C=character())

for (i in c(1:length(dialysed_samples$Sample))){
  if (grepl("A-",dialysed_samples$Sample[i])){
    ALGAE_temp<-unlist(strsplit(dialysed_samples$Sample[i],"-"))
    temp_plant<-"ALGAE"
    EXP_temp<-gsub("[^a-zA-Z]", "-", ALGAE_temp[2])
    EXP_temp2<-unlist(strsplit(EXP_temp, "-"))
    temp_waterbody<-EXP_temp2[1]
    temp_timing<-EXP_temp2[2]
    REP_COND_temp<-str_extract(ALGAE_temp[2],"[0-9]+")
    if (REP_COND_temp <= 3){
      temp_condition<-"light"
      temp_replicate<-REP_COND_temp
    }
    if (REP_COND_temp >=4){
      temp_condition<-"dark"
      temp_replicate<-as.numeric(REP_COND_temp)-3
    }
    temp_enriched13C<-"NO"
  }
  if (grepl("S-",dialysed_samples$Sample[i])){
    SEAGRASS_temp<-unlist(strsplit(dialysed_samples$Sample[i],"-"))
    temp_plant<-"SEAGRASS"
    EXP_temp<-gsub("[^a-zA-Z]", "-", SEAGRASS_temp[2])
    EXP_temp2<-unlist(strsplit(EXP_temp, "-"))
    temp_waterbody<-EXP_temp2[1]
    temp_timing<-EXP_temp2[length(EXP_temp2)]
    REP_COND_temp<-str_extract(SEAGRASS_temp[2],"[0-9]+")
    if (as.numeric(REP_COND_temp) <= 3){
      temp_condition<-"light"
      temp_replicate<-REP_COND_temp
      temp_enriched13C<-"NO"
    }
    if (as.numeric(REP_COND_temp) >=4 && as.numeric(REP_COND_temp) <=6){
      temp_condition<-"dark"
      temp_replicate<-as.numeric(REP_COND_temp)-3
      temp_enriched13C<-"NO"
    }
    if (as.numeric(REP_COND_temp) >=7 && as.numeric(REP_COND_temp) <=9){
      temp_condition<-"light"
      temp_replicate<-as.numeric(REP_COND_temp)-6
      temp_enriched13C<-"YES"
    }
    if (as.numeric(REP_COND_temp) >=10){
      temp_condition<-"dark"
      temp_replicate<-as.numeric(REP_COND_temp)-9
      temp_enriched13C<-"YES"
    }
  }
  raw_metadata<-rbind(raw_metadata, c(temp_plant, 
                                      temp_waterbody, 
                                      temp_timing,
                                      temp_condition,
                                      temp_replicate,
                                      temp_enriched13C))
  colnames(raw_metadata)<-c("Plant", "Waterbody", "Timing", "Condition", 
                            "Replicate", "Enrichment13C")
  temp_plant<-NULL
  temp_waterbody<-NULL
  temp_timting<-NULL
  temp_condition<-NULL
  temp_replicate<-NULL
  temp_enriched13C<-NULL
}

dialysed_samples$Plant<-raw_metadata$Plant
dialysed_samples$Waterbody<-raw_metadata$Waterbody
dialysed_samples$Timing<-raw_metadata$Timing
dialysed_samples$Condition<-raw_metadata$Condition
dialysed_samples$Replicate<-raw_metadata$Replicate
dialysed_samples$Enrichment13C<-raw_metadata$Enrichment13C

# 2.2.3 samples ####
time_cat<-NULL

for (i in (1:length(dialysed_samples$Number))){
  if (dialysed_samples$Timing[i]=="I"){
    time_cat[i]<-"egal"
  }
  if (dialysed_samples$Timing[i]=="E" & dialysed_samples$Condition[i]=="dark"){
    time_cat[i]<-"end_dark"
  }
  if (dialysed_samples$Timing[i]=="E" & dialysed_samples$Condition[i]=="light"){
    time_cat[i]<-"end_light"
  }
}

dialysed_samples$time_cat<-time_cat

time_cat_colors<-c(scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.8),
                   scico(1, palette = "davos",begin=0.1,end=0.2, alpha=0.8),
                   scico(1, palette = "lajolla",begin=0.2,end=0.3, alpha=0.8))

dialysed_samples_melted<-melt(dialysed_samples, id.vars = c("Number","Sample", "Plant", 
                                                            "Waterbody", "Timing", "Condition", 
                                                            "Replicate", "Enrichment13C", "sum", "time_cat"))

SE_dialysed_samples_melted<-summarySE(dialysed_samples_melted, measurevar = "value",
                                      groupvars = c("Plant","time_cat","variable"))



# plotting samples --------------------------------------------------------


# boxplot algae

ggplot(data=dialysed_samples_melted[dialysed_samples_melted$Plant=="ALGAE",])+
  geom_boxplot(aes(x=variable,y=value,fill=time_cat))+
  scale_x_discrete(name="Monosaccharide")+
  scale_y_continuous(name=expression(paste(Concentration~(?g~L^-1))))+
  scale_fill_manual(name="Algae\nincubation", labels=c("background",
                                                       "in dark",
                                                       "in light"),
                    values = time_cat_colors)+
  theme_bw()


# barplot algae

ggplot(data=SE_dialysed_samples_melted[SE_dialysed_samples_melted$Plant=="ALGAE",],
       aes(x=variable, y=value, fill=time_cat))+
  geom_bar(position = position_dodge(), stat = "identity", colour="black")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.2, position=position_dodge(.9))+
  scale_x_discrete(name="Monosaccharide")+
  scale_y_continuous(name=expression(paste(Concentration~(?g~L^-1))))+
  scale_fill_manual(name="Algae\nincubation", labels=c("background",
                                                       "in dark",
                                                       "in light"),
                    values = time_cat_colors)+
  theme_bw()



# boxplot seagrass

ggplot(data=dialysed_samples_melted[dialysed_samples_melted$Plant=="SEAGRASS",])+
  scale_x_discrete(name="Monosaccharide")+
  scale_y_continuous(name=expression(paste(Concentration~(?g~L^-1))))+
  geom_boxplot(aes(x=variable,y=value,fill=time_cat))+
  ggtitle("Seagrass incubations\nZostera marina, 5 h")+
  scale_fill_manual(name=NULL, labels=c("before incubation",
                                        "incubation in dark",
                                        "incubation in light"),
                    values = time_cat_colors)+
  theme_bw()

#barplot seagrass

ggplot(data=SE_dialysed_samples_melted[SE_dialysed_samples_melted$Plant=="SEAGRASS",],
       aes(x=variable, y=value, fill=time_cat))+
  geom_bar(position = position_dodge(), stat = "identity", colour="black")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.2, position=position_dodge(.9))+
  scale_x_discrete(name="Monosaccharide")+
  scale_y_continuous(name=expression(paste(Concentration~(?g~L^-1))))+
  scale_fill_manual(name=NULL,  labels=c("before incubation",
                                         "incubation in dark",
                                         "incubation in light"),
                    values = time_cat_colors)+
  theme_bw()


#barplot both

tiff("C:/Users/admin/ownCloud/Assemble/Data/AEX/Monosaccharides/20210209_monosaccharides_barplot_samples_labelling.tiff",
     width = 40, height=20, units = "cm", res=300)

ggplot(data=dialysed_samples_melted)+
  geom_bar(aes(x=variable, y=value, fill=time_cat),
           position = position_dodge(), stat = "summary", colour="black")+
  geom_point(aes(y=value, x=variable, fill=time_cat), 
             position = position_dodge(.9), pch=21, size=2)+
  scale_x_discrete(name="Monosaccharide", 
                   labels=c("fucose", "rhamnose","galactosamine","arabinose","glucosamine",
                            "galactose","glucose","mannose\n(bad separation)","xylose",
                            "gluconic\nacid","muramic\nacid","galacturonic\nacid",
                            "glucuronic\nacid","mannuronic\nacid","iduronic\nacid"))+
  scale_y_continuous(name=expression(paste(Concentration~(?g~L^-1))))+
  scale_fill_manual(name=NULL,  labels=c("before incubation",
                                         "incubation in dark",
                                         "incubation in light"),
                    values = time_cat_colors)+
  facet_grid(rows="Plant")+
  theme_bw()
dev.off()


######DOC#####
DOC  <-read.csv("C:/Users/admin/ownCloud/Assemble/DOC.csv", dec=".",sep = ",",stringsAsFactors = FALSE)
DOC$Organism = as.factor(DOC$Organism)
DOC$Types = as.factor(DOC$Types)
DOC$Timepoint = as.factor(DOC$Timepoint)
DOC$Treatment = as.factor(DOC$Treatment)
DOC$Timepoint  = factor(DOC$Timepoint, levels=c("I", "E"))



# correlation DOC fucose --------------------------------------------------
subset_fucose <-dialysed_samples_melted %>%
  filter(variable=="fucose")

name<-data.frame(name=character())
for (i in (1:length(DOC$Organism))){
  if (DOC$Types[i]=="PW"){
    next
  }
  if (DOC$Organism[i]=="Seagrass"){
    Organism<-"S"
  }
  if (DOC$Organism[i]=="Algae"){
    Organism<-"A"
  }
  Type<-"WC"
  Timepoint<-DOC$Timepoint[i]
  Chamber<-str_extract(DOC$Chamber[i],"[0-9]+")
  temp_name<-paste(Organism,"-",Type,Chamber, Timepoint,sep="")
  name<-rbind(name,temp_name)
  colnames(name)<-"Sample"
}

new_DOC<-data.frame(name,DOC$DOC[DOC$Types=="WC"],DOC$DON[DOC$Types=="WC"])

corr_table<-merge(new_DOC,subset_fucose[subset_fucose$Enrichment13C=="NO",], by = "Sample")

tiff("C:/Users/admin/ownCloud/Assemble/Data/AEX/Monosaccharides/DOCvsFucose_labelling.tiff",
     height=12, width=15, units="cm", res=300)
ggplot(corr_table[corr_table$Plant=="ALGAE",], aes(x = DOC.DOC.DOC.Types.....WC.., y = value, col = time_cat,
                                                   label = Sample))+
  geom_point(size = 3)+
  scale_y_continuous(name=expression(paste(Fucose~concentration~(?g~L^-1))))+
  scale_x_continuous(name=expression(paste(DOC~concentration~(mg~C~L^-1))),limits = c(3,7))+
  scale_color_manual(name=NULL,  labels=c("before incubation",
                                          "incubation in dark",
                                          "incubation in light"),
                     values = time_cat_colors)+
  geom_text(aes(label=Sample),size= 3, hjust=1, vjust=0)+
  theme_bw()
dev.off()

