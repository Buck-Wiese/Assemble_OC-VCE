### method dev. anionic extraction

#packages ####
require(ggplot2)
require(tidyr)
require(reshape2)
require(wesanderson)
require(dplyr)
require(tidyverse)
require(scico)
require(Rmisc)


# I. method developent ####
#hydrolysate from 4 different column types on HPAEC-PAD Nov 2020

raw_data<-read.table("C:/Users/admin/ownCloud/Assemble/Data/AEX/20201201_Finland_ANXextracts-meth-dev_NaClelution_only-ug-per-L_sampleIDs_dbatches.txt", 
                     h=T, stringsAsFactors = F, colClasses = character())

raw_data[raw_data=="n.a."]<-0

temp_raw<-sapply(raw_data[,7:21], as.numeric)
raw_data$Amount_ugperL_sum<-rowSums(temp_raw)


# Data treatment: substract signal in dialysis blanks per batch from samples

dblank_batch1<-colMeans(sapply(raw_data[raw_data$Sample=="Dialysis_blank" & raw_data$Dialysis_batch==1,7:21], as.numeric))
dblank_batch2<-colMeans(sapply(raw_data[raw_data$Sample=="Dialysis_blank" & raw_data$Dialysis_batch==2,7:21], as.numeric))
dblank_batch3<-colMeans(sapply(raw_data[raw_data$Sample=="Dialysis_blank" & raw_data$Dialysis_batch==3,7:21], as.numeric))
dblank_batch4<-colMeans(sapply(raw_data[raw_data$Sample=="Dialysis_blank" & raw_data$Dialysis_batch==4,7:21], as.numeric))

dialysis_corrected_data<-raw_data[raw_data$Environment!="Lab",]
for (line in 1:length(dialysis_corrected_data[,1])){
  if (dialysis_corrected_data$Dialysis_batch[line]==1){
    for (column in 1:14){
      dialysis_corrected_data[line,column+6]<-as.numeric(dialysis_corrected_data[line,column+6])-dblank_batch1[column]
    }
  }
  if (dialysis_corrected_data$Dialysis_batch[line]==2){
    for (column in 1:14){
      dialysis_corrected_data[line,column+6]<-as.numeric(dialysis_corrected_data[line,column+6])-dblank_batch2[column]
    }
  }
  if (dialysis_corrected_data$Dialysis_batch[line]==3){
    for (column in 1:14){
      dialysis_corrected_data[line,column+6]<-as.numeric(dialysis_corrected_data[line,column+6])-dblank_batch3[column]
    }
  }
  if (dialysis_corrected_data$Dialysis_batch[line]==4){
    for (column in 1:14){
      dialysis_corrected_data[line,column+6]<-as.numeric(dialysis_corrected_data[line,column+6])-dblank_batch4[column]
    }
  }
  
}

#take out algae anxff that failed

dialysis_corrected_data<-dialysis_corrected_data[dialysis_corrected_data$Sample!="ANXff_Algae_A_NaCl",]

data_dc_melted<-melt(dialysis_corrected_data, id.vars = c("Number", "Injection", "Dialysis_batch","Sample","Environment","Column"))

# in dialysis blanks, there is 
column_cols<-c(wes_palette("Chevalier1", n=3), wes_palette("Zissou1", n=1), wes_palette("Darjeeling1", n=1))
column_cols<-c(column_cols[1],column_cols[3],column_cols[2],column_cols[4], column_cols[5])

data_dc_melted$value[data_dc_melted$value<=0]<-0

tiff("C:/Users/hbuck/ownCloud/Assemble/Data/AEX/Fig_methdev_ANXFF-dialysisblanks.tiff",
     width = 36, height = 12, units = "cm", res = 150)
ggplot(data=data_dc_melted[data_dc_melted$Column=="ANXff",])+
  geom_boxplot(aes(x=variable, y=as.numeric(value), fill=Environment), 
               position = position_dodge2(width=1), pch=21, alpha=0.6)+
  geom_point(aes(x=variable, y=as.numeric(value), fill=Environment), 
             position = position_dodge2(width=0.8), pch=21, size=2)+
  scale_x_discrete(name=expression(paste("Monosaccharide")), 
                   labels=c("Fucose", "Rhamnose","Galactosamine","Arabinose","Glucosamine","Galactose","Glucose",
                            "Mannose","Xylose","Gluconic\nacid"," Muramic\nacid","Galacturonic\nacid","Glucuronic\nacid",
                            "Mannuronic\nacid","Iduronic\nacid"))+
  scale_fill_manual(values=column_cols)+
  scale_y_continuous(name=expression(paste("Concentration (",mu,"g L"^"-1",")")))+
  theme_bw()
dev.off()


raw_data$Amount_ugperL_sum<-
  column_cols<-c(wes_palette("GrandBudapest1", n=3),wes_palette("GrandBudapest2", n=1))
tiff("C:/Users/hbuck/ownCloud/Assemble/Data/AEX/Fig_prelim_sum_extracted_wspiked.tiff", 
     width = 17.5, height = 12, units = "cm", res = 150)
ggplot(data=raw_data[which(raw_data$Column %in% c("ANXff") ),])+
  geom_boxplot(aes(x=Environment, y=as.numeric(Amount_ugperL_sum), fill=Environment), 
               position = position_dodge2(width=1), pch=21, alpha=0.8)+
  scale_fill_manual(values=column_cols)+
  scale_y_continuous(name=expression(paste("Concentration (",mu,"g L"^"-1",")")))+
  theme_bw()
dev.off()

column_cols<-c(wes_palette("GrandBudapest1", n=3),wes_palette("GrandBudapest2", n=1))
tiff("C:/Users/hbuck/ownCloud/Assemble/Data/AEX/Fig_prelim_sum_extracted_WOspike.tiff", 
     width = 17.5, height = 12, units = "cm", res = 150)
ggplot(data=raw_data[which(raw_data$Column %in% c("ANXseph", "ANXff","Source15Q","HP") &raw_data$Environment!="spikedSand"),])+
  geom_boxplot(aes(x=Environment, y=as.numeric(Amount_ugperL_sum), fill=Column), 
               position = position_dodge2(width=1), pch=21, alpha=0.8)+
  scale_fill_manual(values=column_cols)+
  theme_bw()
dev.off()