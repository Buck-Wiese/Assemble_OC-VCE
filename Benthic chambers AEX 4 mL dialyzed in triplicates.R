##### AEX extract from Benthic chambers reanalyzed in triplicates


#packages ####
require(ggplot2)
require(tidyr)
require(reshape2)
require(wesanderson)
require(dplyr)
require(tidyverse)
require(scico)
require(Rmisc)


# III. Re: Benthic chambers (triplicates) --------------------------------------

# 1 -  data import ####

triplicate_data<-read.table("C:/Users/admin/ownCloud/Assemble/Data/AEX/20210428_PA10_3_Finland_AEX_with samplenames_forR.txt", 
                            h=T, stringsAsFactors = F, colClasses = character())
#rownames(triplicate_data)<-triplicate_data$nr

#temp_raw<-sapply(chamber_data[,4:18], as.numeric)
#chamber_data$Amount_ugperL_sum<-rowSums(temp_raw)

# replace n.a. in sugar measurements with 0
for (i in 1:length(triplicate_data$nr)){
  for (j in 4:18){
    if (triplicate_data[i,j]=="n.a."){
      triplicate_data[i,j]<-0
    }
  }
}

# change sugar measurements to numeric
for (i in 4:18){
  triplicate_data[, i]<-as.numeric(triplicate_data[,i])
}



# standards - visualizing signal decrease over time -----------------------


quant_stds<-triplicate_data[complete.cases(triplicate_data$Standard_concentration),]


plot(x=as.numeric(quant_stds[,'Standard_concentration']),y=as.numeric(quant_stds[,'Fucose']))
text(x=as.numeric(quant_stds[,'Standard_concentration']), y=as.numeric(quant_stds[,'Fucose']) , labels=as.character(quant_stds[,'nr']), 
     cex=0.65, pos=2,col="red")
##### some code to test different models

model_stds<-lm(as.numeric(Fucose) ~ Standard_concentration+  Fucose:nr_squared -1, quant_stds)

par(mfrow=c(2,1))
plot(x=as.numeric(quant_stds[,'Standard_concentration']),y=as.numeric(quant_stds[,'Fucose']))
text(x=as.numeric(quant_stds[,'Standard_concentration']), y=as.numeric(quant_stds[,'Fucose']) , labels=as.character(quant_stds[,'nr']), 
     cex=0.65, pos=2,col="red")

plot(x=as.numeric(quant_stds[,'Standard_concentration']),
     y=as.numeric(quant_stds[,'Fucose'])*quant_stds[,'nr']*quant_stds[,'Standard_concentration']/model_stds$coefficients[3]+model_stds$coefficients[1])
plot(x=quant_stds[,'Standard_concentration'], 
     y=quant_stds[,'Fucose']-(quant_stds[,'Fucose']*quant_stds[,'nr_squared']*model_stds$coefficients[2]))

text(x=as.numeric(quant_stds[,'Standard_concentration']), y=as.numeric(quant_stds[,'Fucose']) , labels=as.character(quant_stds[,'nr']), 
     cex=0.65, pos=2,col="red")



# Correct for decrease over time ------------------------------------------

# a. calculate factors from quant standards
quant_stds$nr_squared<-(quant_stds$nr^2)
correction<-data.frame(monomer=character(),
                       concentration_factor=numeric(),
                       time_concentration_factor=numeric(),
                       pvalue_conc=numeric(),
                       pvalue_time_conc=numeric(),
                       R_squared=numeric()
                       )

for (i in 4:18){
  monomer_temp<-colnames(quant_stds)[i]
  model_temp<-lm(quant_stds[,i]~quant_stds$Standard_concentration+quant_stds[,i]:quant_stds$nr_squared -1)
  conc_factor_temp<-model_temp$coefficients[1]
  time_conc_factor_temp<-model_temp$coefficients[2]
  pvalue_conc_temp<-summary(model_temp)$coefficients[1,4]
  pvalue_time_conc_temp<-summary(model_temp)$coefficients[2,4]
  R_squared_temp<-summary(model_temp)$adj.r.squared
  correction<-rbind(correction, c(monomer=monomer_temp,
                      concentration_factor=conc_factor_temp,
                      time_concentration_factor=time_conc_factor_temp,
                      pvalue_conc=pvalue_conc_temp,
                      pvalue_time_conc=pvalue_time_conc_temp,
                      R_squared=R_squared_temp))
  colnames(correction)<-c("monomer", "concentration_factor","time_concentration_factor",
                          "pvalue_conc", "pvalue_time_conc", "R_squared")
}


##### some code to test correction on samples

triplicate_data$Fucose_corr<-triplicate_data$Fucose-(triplicate_data$Fucose*triplicate_data$nr*triplicate_data$nr*model_stds$coefficients[2])

quant_stds<-triplicate_data[complete.cases(triplicate_data$Standard_concentration),]
plot(x=quant_stds$Standard_concentration, y=quant_stds$Fucose_corr)

triplicate_data_noA<-triplicate_data[!grepl(pattern = "^A", triplicate_data$Sample),]
triplicate_data_noA_nostds<-triplicate_data_noA[triplicate_data_noA$Sample!="Standard",]


# b. correction of signals in all runs

triplicate_data_corr<-data.frame(nr=numeric(),
                                 Standard_concentration=numeric(),
                                 Sample=character(),
                                 Monomer=character(),
                                 Intensity=numeric(),
                                 Intensity_corr=numeric())
for (i in 1:length(triplicate_data[,1])){
    nr_temp<-triplicate_data$nr[i]
    nr_squared_temp<-nr_temp*nr_temp
    Standard_concentration_temp<-triplicate_data$Standard_concentration[i]
    Sample_temp<-triplicate_data$Sample[i]
    
    for (j in 4:18){
      Monomer_temp<-colnames(triplicate_data)[j]
      Intensity_temp<-triplicate_data[i,j]
      Intensity_corr_temp<-Intensity_temp-(Intensity_temp*nr_squared_temp*as.numeric(correction[correction$monomer==Monomer_temp,3]))
      triplicate_data_corr<-rbind(triplicate_data_corr, c(nr_temp, nr_squared_temp, 
                                                          Standard_concentration_temp, Sample_temp,
                                                          Monomer_temp, Intensity_temp, Intensity_corr_temp))
      colnames(triplicate_data_corr)<-c("nr", "nr_squared", "Standard_concentration", "Sample",
                                        "Monomer", "Intensity", "Intensity_corr")
    }

}



# calculate means of triplicates ------------------------------------------


triplicate_data_corr_means<-triplicate_data_corr %>%
  group_by(Standard_concentration, Sample, Monomer)%>%
  summarize(mean_int_corr=mean(as.numeric(Intensity_corr)), std_dev_corr=sd(as.numeric(Intensity_corr)), n=n())
  


# Turn sample names into metadata -----------------------------------------
meta_temp=data.frame(plant=character(),
                     timepoint=character(),
                     treatment=character(), 
                     time_treat=character())

plant_temp=character()
timepoint_temp=character()
treatmeant_temp=character()
time_treat_temp=character()

for (i in 1:length(triplicate_data_corr_means$Standard_concentration)){
  if (grepl(triplicate_data_corr_means$Sample[i], pattern = "-")){
  name_parts_temp<-unlist(str_split(triplicate_data_corr_means$Sample[i], "-"))
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
  else {
    plant_temp<-"none"
    timepoint_temp<-"none"
    treatment_temp<-"none"
    time_treat_temp<-"none"
  }
  meta_temp<-rbind(meta_temp, c(as.character(plant_temp), as.character(timepoint_temp), 
                                as.character(treatment_temp), as.character(time_treat_temp)))
  colnames(meta_temp)<-c("plant", "timepoint", "treatment", "time_treat_temp")
}

triplicate_data_corr_means$Plant<-meta_temp$plant
triplicate_data_corr_means$Timepoint<-meta_temp$timepoint  
triplicate_data_corr_means$Treatment<-meta_temp$treatment 
triplicate_data_corr_means$Time_treat<-meta_temp$time_treat_temp

# plotting corrected means by plant ---------------------------------------
time_cat_colors<-c(scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.8),
                   scico(1, palette = "davos",begin=0.1,end=0.2, alpha=0.8),
                   scico(1, palette = "lajolla",begin=0.2,end=0.3, alpha=0.8))

triplicate_data_corr_means$Monomer<-factor(triplicate_data_corr_means$Monomer,
              levels=c("Fucose", "Rhamnose","Galactosamine","Arabinose","Glucosamine",
                       "Galactose","Glucose","Mannose","Xylose",
                       "Gluconic_acid","Muramic_acid","Galacturonic_acid",
                       "Glucuronic_acid","Mannuronic_acid","Iduronic_acid"))

ggplot(data=triplicate_data_corr_means[triplicate_data_corr_means$Plant!="none",])+
  geom_bar(aes(x=Monomer, y=mean_int_corr, fill=Time_treat),
           position = position_dodge(), stat = "summary", colour="black")+
  geom_point(aes(x=Monomer, y=mean_int_corr, fill=Time_treat), 
             position = position_dodge(.9), pch=21, size=2)+
  scale_x_discrete(name="Monosaccharide", 
                   labels=c("fucose", "rhamnose","galactosamine","arabinose","glucosamine",
                            "galactose","glucose","mannose\n(bad separation)","xylose",
                            "gluconic\nacid","muramic\nacid","galacturonic\nacid",
                            "glucuronic\nacid","mannuronic\nacid","iduronic\nacid"))+
  scale_y_continuous(name=expression(paste(Concentration~(microg~L^-1))))+
  scale_fill_manual(name=NULL,  breaks=c("initial", "end_dark", "end_light"),
                    labels=c("before incubation",
                                         "incubation in dark",
                                         "incubation in light"),
                    values = time_cat_colors)+
  facet_grid(rows="Plant")+
  theme_bw()


#Seagrass WC
Samples <- rbind(SG_WC<-triplicate_data[grepl(pattern = "^S-", triplicate_data$Sample),],
                 triplicate_data[grepl(pattern = "^A", triplicate_data$Sample),])


#Algae WC
A_WC<-triplicate_data[grepl(pattern = "^A", triplicate_data$Sample),]
A_WC_<-


plot(x=triplicate_data$nr, y=triplicate_data$Fucose)
plot(x=triplicate_data_noA_nostds$nr, y=triplicate_data_noA_nostds$Glucosamine)
text(x=triplicate_data_noA_nostds$nr, y=triplicate_data_noA_nostds$Glucosamine, labels = triplicate_data_noA_nostds$Sample, pos=3, cex = 0.7)

