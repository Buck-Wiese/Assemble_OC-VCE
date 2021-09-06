##### Benthic chamber raw samples algae, seawater dialyzed, acid hydrolyzed, HPAEC-PAD


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
require(ggpmisc)
require(grid)
require(FSA)
require(lme4)
require(lmerTest)
require(emmeans)
require(cowplot)
require(vegan)

fonts <- list(
  sans = "Helvetica",
  mono = "Consolas",
  `Times New Roman` = "DejaVu Serif"
)
font_family<-"sans"
font_size<-10
stats_text_size<-3

time_cat_colors<-c(scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.7),
                   scico(2, palette = "vikO",begin=0.25,end=0.75, alpha=0.7))
time_cat_4colors<-c("grey",scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.7),
                    scico(2, palette = "vikO",begin=0.25,end=0.75, alpha=0.7))
time_cat_5colors<-c("white", "grey",scico(1, palette = "cork",begin=0.6,end=0.7, alpha=0.7),
                    scico(2, palette = "vikO",begin=0.25,end=0.75, alpha=0.7))
plant_colors<-scico(2, palette="tofino", begin = 0.6, end=0.9, alpha=0.7)
'%!in%' <- function(x,y)!('%in%'(x,y))

# V: RAW SAMPLES --------------------------------------

raw_samples<-read_xls("C:/Users/admin/ownCloud/Assemble/Data/raw samples AH/20210802_monosaccharides_raw_samples_AH_results_for_R.xls")

ggplot(raw_samples[raw_samples$Sample_group=="raw_sample",])+
  geom_point(aes(x=Concentration, y=Fucose))

raw_samples$Group_for_plot<-factor(raw_samples$Group_for_plot,
                                   levels=c("dialysis_blank", "procedure_blank",
                                            "initial", "dark_end", "light_end"),
                                   ordered = T)

Fig_rawsamples_fucose<-ggplot(data=raw_samples[which(raw_samples$Sample_group %in% 
                                                c("raw_sample", "Procedure_blank", "Dialysis_blank")),])+
  geom_bar(aes(x=Group_for_plot, y=Fucose, fill=Group_for_plot), stat="summary", color="#262626")+
  geom_point(aes(x=Group_for_plot, y=Fucose, fill=Group_for_plot, group=Group_for_plot),
             pch=21, position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size=2)+
  scale_fill_manual(name="", labels=c("dialysis\nblank","procedure\nblank", 
                                      "initial\nseawater", "dark\nincubation", "light\nincubation"),
                    values=time_cat_5colors)+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Fucose signal")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        axis.text.x = element_blank(),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.ticks.x = element_blank(),
        strip.background.x = element_rect(fill = "#262626"))

Fig_flowthrough_fucose<-ggplot(data=raw_samples[which(raw_samples$Sample_group %in% 
                                                c("ANX_flowthrough", "raw_sample")),])+
  geom_bar(aes(x=Group_for_plot, y=Fucose, group=Sample_group, fill=Group_for_plot), 
           stat="summary", color="#262626", position="dodge")+
  geom_point(aes(x=Group_for_plot, y=Fucose, group=Sample_group, fill=Group_for_plot),
             pch=21, position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size=2)+
  scale_fill_manual(name="", labels=c("initial\nseawater", "dark\nincubation", "light\nincubation"),
                    values=time_cat_colors)+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Fucose signal")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        axis.text.x = element_blank(),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.ticks.x = element_blank(),
        strip.background.x = element_rect(fill = "#262626"))


# melt monosaccharides

melted_raw_samples<-melt(raw_samples, measure.vars = colnames(raw_samples)[7:21])

Fig_rawsamples_allmonos<-ggplot(data=melted_raw_samples[which(melted_raw_samples$Sample_group %in% 
                                                                c("raw_sample", "Procedure_blank", "Dialysis_blank")),])+
  geom_bar(aes(x=variable, y=value, fill=Group_for_plot), 
           stat="summary", color="#262626", position=position_dodge2(width = 0.5))+
  geom_point(aes(x=variable, y=value, fill=Group_for_plot, group=Group_for_plot),
             pch=21, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), size=2)+
  scale_fill_manual(name="", labels=c("dialysis\nblank","procedure\nblank", 
                                      "initial\nseawater", "dark\nincubation", "light\nincubation"),
                    values=time_cat_5colors)+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Fucose signal")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.x = element_rect(fill = "#262626"))

selected_monos_raw_samples<-melted_raw_samples[which(melted_raw_samples$variable %in%
                                                       c("Fucose", "Galactose", "Glucose",
                                                         "Rhamnose", "Glucuronic_acid")),]

# calibration using dialyzed, hydrolyzed fucoidan #####
# fucoidan standard

formula<-y~x
Fig_std_fucoidan<-ggplot(data=selected_monos_raw_samples[which(selected_monos_raw_samples$Sample_group %in%
                                                                 "Fucoidan_std"),])+
  geom_smooth(method = "lm",aes(x=as.numeric(Concentration), y=value/20, fill=variable), 
           color="#262626", position=position_dodge2(width = 0.5),  formula = formula)+
  geom_point(aes(x=as.numeric(Concentration), y=value/20, fill=variable, group=variable),
             pch=21, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), size=2)+
  stat_fit_glance(method="lm", method.args = list(formula=formula), geom="text", 
                  aes(y=value/20, x=as.numeric(Concentration),  
                      label = sprintf('R²=%.3f p=%.2g',stat(r.squared), stat(p.value))), 
                  color=scico(1, alpha=1, end=0.8, palette = "vikO"),
                  label.x=0, label.y=seq(200, 110, -30), size=stats_text_size)+
  scale_x_continuous(name=expression(paste(Fucoidan~(microg~L^-1))))+
  scale_y_continuous(name=expression(paste("Concentration (", mu, "g",L^-1,")")))+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.x = element_rect(fill = "#262626"))

tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210802_Fig_raw samples monosaccharides.tiff",
     width=44, height=12, units = "cm", res=300)
plot_grid(Fig_std_fucoidan, Fig_rawsamples_selectedmonos, Fig_flowthrough_selectedmonos, nrow = 1)
dev.off()

#data set for fucoidan calibration
calibration_fucoidan_data<-selected_monos_raw_samples[which(selected_monos_raw_samples$Sample_group %in%
                                   "Fucoidan_std"),]
fucoidan_std_cal_all<-lm(value/20~as.numeric(Concentration)*variable, 
                     data=calibration_fucoidan_data)
summary(fucoidan_std_cal_all)
plot(fucoidan_std_cal_all)
calibration_fucoidan_data_fucose<-calibration_fucoidan_data[calibration_fucoidan_data$variable=="Fucose",]
fucoidan_std_cal_fucose<-lm(value/20~as.numeric(Concentration)-1, 
                            data=calibration_fucoidan_data_fucose)
summary(fucoidan_std_cal_fucose)
slope_fucoidan_cal<-fucoidan_std_cal_fucose$coefficients

raw_samples$Fucoidan_calc<-raw_samples$Fucose/20/slope_fucoidan_cal # 20 is dialysis concentration factor, slope is from fucoidan standard series

Fig_rawsamples_fucoidan<-ggplot(data=raw_samples[which(raw_samples$Sample_group %in% 
                                                       c("raw_sample", "Procedure_blank", "Dialysis_blank")),])+
  geom_bar(aes(x=Group_for_plot, y=Fucoidan_calc, fill=Group_for_plot), stat="summary", color="#262626")+
  geom_point(aes(x=Group_for_plot, y=Fucoidan_calc, fill=Group_for_plot, group=Group_for_plot),
             pch=21, position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size=2)+
  scale_fill_manual(name="", labels=c("dialysis\nblank","procedure\nblank", 
                                      "seawater", "dark\ninc.", "light\ninc."),
                    values=time_cat_5colors)+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Fucoidan signal")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        axis.text.x = element_blank(),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.ticks.x = element_blank(),
        strip.background.x = element_rect(fill = "#262626"))

raw_samples$FucoidanC_calc<-raw_samples$Fucoidan_calc*0.245 # 0.245 is C content of Sigma pure F. vesiculosus fucoidan

Fig_rawsamples_fucoidan_carbon<-ggplot(data=raw_samples[which(raw_samples$Sample_group %in% 
                                                         c("raw_sample", "Procedure_blank", "Dialysis_blank")),])+
  geom_bar(aes(x=Group_for_plot, y=FucoidanC_calc/1000, fill=Group_for_plot), stat="summary", color="#262626")+
  geom_point(aes(x=Group_for_plot, y=FucoidanC_calc/1000, fill=Group_for_plot, group=Group_for_plot),
             pch=21, position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size=2)+
  scale_fill_manual(name="", labels=c("dialysis\nblank","procedure\nblank", 
                                      "seawater", "dark\ninc.", "light\ninc."),
                    values=time_cat_5colors)+
  scale_x_discrete(name="")+
  scale_y_continuous(name=expression(paste("Fucoidan carbon (g ", m^-3,")")))+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        axis.text.x = element_blank(),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.ticks.x = element_blank(),
        strip.background.x = element_rect(fill = "#262626"))

# substract flowthrough #####

raw_samples_subsetted<-raw_samples[which(raw_samples$Sample_group %in% c("raw_sample", "ANX_flowthrough")),]
raw_samples_blanks<-raw_samples[which(raw_samples$Sample_group %in% c("Procedure_blank", "Dialysis_blank")),]
raw_samples_subsetted$Sample_name[13:24]<-raw_samples_subsetted$Sample_name[1:12]
samples_substracted<-data.frame(list(raw_samples_subsetted$Sample_name[1:12], 
                                     raw_samples_subsetted$Group_for_plot[1:12]))
colnames(samples_substracted)<-c("Sample_name", "Group_for_plot")
samples_substracted$raw_sample_FucoidanC_calc<-raw_samples_subsetted$FucoidanC_calc[1:12]
samples_substracted$flowthrough_FucoidanC_calc<-raw_samples_subsetted$FucoidanC_calc[13:24]
samples_substracted$net_FucoidanC_calc<-samples_substracted$raw_sample_FucoidanC_calc-
  samples_substracted$flowthrough_FucoidanC_calc
blanks_for_substracted<-data.frame(list(raw_samples_blanks$Sample_name, 
                                        raw_samples_blanks$Group_for_plot))
colnames(blanks_for_substracted)<-c("Sample_name", "Group_for_plot")
blanks_for_substracted$raw_sample_FucoidanC_calc<-raw_samples_blanks$FucoidanC_calc

blanks_for_substracted$flowthrough_FucoidanC_calc<-NA
blanks_for_substracted$net_FucoidanC_calc<-blanks_for_substracted$raw_sample_FucoidanC_calc

net_fucoidan<-rbind(samples_substracted, blanks_for_substracted)


# fucoidan net figure #####

Fig_substracted_fucoidan_carbon<-ggplot(data=net_fucoidan)+
  geom_bar(aes(x=Group_for_plot, y=net_FucoidanC_calc/1000, fill=Group_for_plot), stat="summary", color="#262626")+
  geom_point(aes(x=Group_for_plot, y=net_FucoidanC_calc/1000, fill=Group_for_plot, group=Group_for_plot),
             pch=21, position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size=2)+
  scale_fill_manual(name="", labels=c("dialysis\nblank","procedure\nblank", 
                                      "seawater", "dark\ninc.", "light\ninc."),
                    values=time_cat_5colors)+
  scale_x_discrete(name="")+
  scale_y_continuous(name=expression(paste("Fucoidan carbon (g ", m^-3,")")))+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        axis.text.x = element_blank(),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.ticks.x = element_blank(),
        strip.background.x = element_rect(fill = "#262626"))

pdf("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210819_Fig_netfucoidanC.pdf", 
    height=3, width = 5)
Fig_substracted_fucoidan_carbon
dev.off()

# raw samples selected monos

selected_monos_raw_samples$Group_for_plot<-factor(selected_monos_raw_samples$Group_for_plot, 
                                                  levels=c("dialysis_blank","procedure_blank",
                                                           "initial", "dark_end", "light_end"),
                                                  ordered = T)
selected_monos_raw_samples$variable<-factor(selected_monos_raw_samples$variable, 
                                            levels=c("Fucose", "Galactose","Glucuronic_acid",
                                                     "Rhamnose", "Glucose"), 
                                            labels=c("Fucose", "Galactose","Glucuronic acid",
                                                     "Rhamnose", "Glucose"),
                                            ordered = T)

# try models on raw samples -----

raw_samples_fucose<-selected_monos_raw_samples[selected_monos_raw_samples$variable=="Fucose",]
raw_samples_fucose <- raw_samples_fucose[which(raw_samples_fucose$Sample_group %in% 
                                                 c("raw_sample")),]

model_raw_fucose<-lm(log(value)~Group_for_plot, raw_samples_fucose)
anova(model_raw_fucose)
plot(model_raw_fucose)
summary(model_raw_fucose)
tukey_hsd(model_raw_fucose)

# assumptions for lm okay when using log and no blanks
# result: end signf different from initial, dark vs. light pval at 0.0535


# non-parametric tests
kruskal.test(x=raw_samples_fucose$value, g=raw_samples_fucose$Group_for_plot)
dunnTest(raw_samples_fucose$value~raw_samples_fucose$Group_for_plot, method = )
# non-parametric test results are ridiculous

# try some multivariate on sugars
require(vegan)

multi_data<-raw_samples[which(raw_samples$Sample_group %in% c("raw_sample", "Procedure_blank", "Dialysis_blank")),]
sugars<-as.data.frame(multi_data[,7:21])
rownames(sugars)<-multi_data$Sample_name

sugars_groups<-as.data.frame(multi_data[,4:6])
rownames(sugars_groups)<-multi_data$Sample_name

sugars_rel<-decostand(sugars, method = "total")
sugars_distmatrix<-vegdist(sugars_rel, method="bray")

cluster_sugars<-hclust(sugars_distmatrix)
plot(cluster_sugars, hang = -1)

svg("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210805_raw_dendrogramm.svg",
     height = 10, width = 14)
plot(cluster_sugars, hang = -1)
dev.off()

pdf("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210819_raw_dendrogramm.pdf",
    height = 10, width = 14)
plot(cluster_sugars, hang = -1)
dev.off()

# remove blanks as too many zeros
multi_data<-raw_samples[which(raw_samples$Sample_group %in% c("raw_sample")),]
sugars<-as.data.frame(multi_data[,7:21])
rownames(sugars)<-multi_data$Sample_name

sugars_groups<-as.data.frame(multi_data[,4:6])
rownames(sugars_groups)<-multi_data$Sample_name

sugars_rel<-decostand(sugars, method = "total")
sugars_distmatrix<-vegdist(sugars_rel, method="bray")
sugars_distmatrix_mat<-as.matrix(sugars_distmatrix, labels=T)

sugars_NMS <-
  metaMDS(sugars_distmatrix_mat,
          distance = "bray",
          k = 2)


plot(sugars_NMS, "sites")
orditorp(sugars_NMS, "sites")



colvec <- c("gray0", "gray24", "gray49", "gray60")   # Identifies colors for group assignments
pchvec <- c(21, 22, 23, 24)   # Identifies character symbols for group assignments

plot(sugars_NMS)
with(sugars_groups,
     points(sugars_NMS,
            display = "sites",
            col = "black",
            pch = pchvec[Group_for_plot],
            bg = colvec[Group_for_plot]))

#Create convex hulls that highlight point clusters based on grouping dataframe
ordihull(
  sugars_NMS,
  sugars_groups$Group_for_plot,
  display = "sites",
  draw = c("polygon"),
  col = NULL,
  border = c("gray0", "gray20","gray30", "gray48"),
  lty = c(1, 2, 1, 2),
  lwd = 2.5
)
scrs <-
  scores(sugars_NMS, display = "sites", "species")
cent <-
  aggregate(scrs ~ Group_for_plot, data = sugars_groups, FUN = "mean")
names(cent) [-1] <- colnames(scrs)
points(cent [,-1],
       pch = c( 8 , 8 , 8),
       col = c("gray0", "gray0", "gray48"),
       bg = c("black"),
       lwd = 3.0,
       cex = 2.0 # Plots centroids as points on ordination
)


# RDA
sugars_rda<-rda(sugars)
#droplevels(sugars_groups$Group_for_plot)
#sugars_groups$Group_for_plot<-factor(sugars_groups,
#                                        levels=c("light_end", "initial", "dark_end"))
sugars_groups$Group_for_plot<-as.character(sugars_groups$Group_for_plot)
for (i in 1:length(sugars_groups$Concentration)){
  if (sugars_groups$Concentration[i]=="initial"){
    sugars_groups$Group_for_plot[i]<-"Seawater"
  }
  if (sugars_groups$Concentration[i]=="dark_end"){
    sugars_groups$Group_for_plot[i]<-"Dark inc."
  }
  if (sugars_groups$Concentration[i]=="light_end"){
    sugars_groups$Group_for_plot[i]<-"Light inc."
  }
  
}
sugars_groups$Group_for_plot
tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210818_eigenvalues_RDAinsert.tiff", 
     height=12, width = 12, units = "cm", res=300)
biplot(sugars_rda,
       display = c("sites", 
                   "species"),
       type = c("points"),
       xlim=c(-110,30)
)

ordihull(sugars_rda,
  sugars_groups$Group_for_plot,
  display = "sites",
  draw = c("polygon"),
  col = time_cat_colors[c(2,3,1)],
  border = time_cat_colors[c(2,3,1)],
  lty = 1,
  lwd = 2.5,
  label = T
)
dev.off()

# eigenval figure insert #####
pdf("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210819_eigenvalues_RDAinsert.pdf", 
     height=4, width = 4)
biplot(sugars_rda,
       display = c("sites", 
                   "species"),
       type = c("points"),
       xlim=c(-110,30)
)

ordihull(sugars_rda,
         sugars_groups$Group_for_plot,
         display = "sites",
         draw = c("polygon"),
         col = time_cat_colors[c(2,3,1)],
         border = time_cat_colors[c(2,3,1)],
         lty = 1,
         lwd = 2.5,
         label = T
)

dev.off()


# Eigenval figure ---------------------------------------------------------
sug <- summary(sugars_rda)$species %>% as.data.frame()
sug$sug <- rownames(sug)

eigenvals<-sug[,c(1, 7)]
eigenvals_other<-eigenvals[which(eigenvals$sug %in% c("Fucose", "Glucuronic_acid", "Xylose", 
                                                "Galactose", "Glucose")),]
other_sum<-sum(abs(eigenvals$PC1[which(eigenvals$sug %!in% c("Fucose", "Glucuronic_acid", "Xylose", 
                                                             "Galactose", "Glucose"))]))

eigenvals_other<-rbind(eigenvals_other, c(other_sum, "others"))
eigenvals_other$sug<-factor(eigenvals_other$sug,
                            levels=c("Fucose", "Glucuronic_acid", "Xylose", 
                                     "Galactose", "Glucose", "others"),
                            labels=c("Fucose", "Glucuronic\nacid", "Xylose", 
                                     "Galactose", "Glucose", "others\nsummed"))
Fig_eigenvalues<-ggplot(dat=eigenvals_other,aes(x=sug,y=abs(as.numeric(PC1)), 
                                                fill=sug))+
  geom_bar(stat = "identity")+
  scale_x_discrete(name="")+
  scale_y_continuous(name="absolute eigenvalue")+
  scale_fill_manual(name="", values = scico(6, alpha=0.7,end=0.8, palette = "vikO"))+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background.x = element_rect(fill = "#262626"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "#262626"))

tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210818_eigenvalues.tiff", 
     height=10, width = 12, units = "cm", res=300)
Fig_eigenvalues
dev.off()
pdf("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210819_eigenvalues.pdf", 
     height=4, width = 4.8)
Fig_eigenvalues
dev.off()


axesimp<-summary(sugars_rda)$cont$importance

Fig_rawsamples_selectedmonos<-ggplot(data=selected_monos_raw_samples[which(selected_monos_raw_samples$Sample_group %in% 
                                                                             c("raw_sample", "Procedure_blank", "Dialysis_blank")),])+
  geom_bar(aes(x=variable, y=value/20, fill=Group_for_plot), 
           stat="summary", color="#262626", position=position_dodge2(width = 0.5))+
  geom_point(aes(x=variable, y=value/20, fill=Group_for_plot, group=Group_for_plot),
             pch=21, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), size=2)+
  scale_fill_manual(name="", labels=c("dialysis\nblank","procedure\nblank", 
                                      "initial\nseawater", "dark\nincubation", "light\nincubation"),
                    values=time_cat_5colors)+
  scale_x_discrete(name="")+
  scale_y_continuous(name=expression(paste(Concentration~(microg~L^-1))))+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.x = element_rect(fill = "#262626"))

fucose_raw_samples<-selected_monos_raw_samples[selected_monos_raw_samples$variable=="Fucose",] 
                                            

Fig_rawsamples_fucose<-ggplot(data=fucose_raw_samples[which(fucose_raw_samples$Sample_group %in% 
                                                                c("raw_sample", "Procedure_blank", "Dialysis_blank")),])+
  geom_bar(aes(x=variable, y=value/20/1000, fill=Group_for_plot), 
           stat="summary", color="#262626", position=position_dodge2(width = 0.5))+
  geom_point(aes(x=variable, y=value/20/1000, fill=Group_for_plot, group=Group_for_plot),
             pch=21, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), size=2)+
  scale_fill_manual(name="", labels=c("dialysis\nblank","procedure\nblank", 
                                      "initial\nseawater", "dark\nincubation", "light\nincubation"),
                    values=time_cat_5colors)+
  scale_x_discrete(name="")+
  scale_y_continuous(name=expression(paste(Fucoidan~concentration~(mg*C~L^-1))))+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.x = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.background.x = element_rect(fill = "#262626"))

tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210805_Fig_fucoidan.tiff", 
     height=6, width = 14, units = "cm", res=300)
Fig_rawsamples_fucose
dev.off()



# quatsch -----------------------------------------------------------------

dat <- fucose_raw_samples[which(fucose_raw_samples$Sample_group %in% 
                           c("raw_sample")),]

require(nlme)

ctrl <- lmeControl(opt='optim',msMaxIter = 100)


m0 <- lm(log(value+1)~Group_for_plot,data=dat)
summary(m0)
plot(m0)
fixvar<-varFixed(~value)

dat1<-dat[-c(2,4, 6), ] %>% droplevels()
dat1$Group_for_plot <- as.character(dat1$Group_for_plot )
m1 <- gls(value~Group_for_plot,data=dat1,weights=fixvar,control = list(singular.ok=T), verbose = T)
summary(m1)
plot(m1)
anova(m1)



# looking at flowthrough: ######
# older stuff ######
Fig_flowthrough_allmonos<-ggplot(data=
                                        melted_raw_samples[which(melted_raw_samples$Sample_group %in% 
                                                                           c("raw_sample", "ANX_flowthrough")),])+
  geom_bar(aes(x=variable, y=value, fill=Sample_group), 
           stat="summary", color="#262626", position=position_dodge2(width = 0.5))+
  geom_point(aes(x=variable, y=value, fill=Sample_group, group=Sample_group),
             pch=21, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), size=2)+
  scale_fill_manual(name="", labels=c("ANX flowthrough", "Raw_sample"),
                    values=time_cat_5colors[c(2:3)])+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Fucose signal")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.y = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.y = element_rect(fill = "#262626"))+
  facet_wrap(~Group_for_plot, nrow=3, strip.position = "right")
Fig_flowthrough_selectedmonos<-ggplot(data=
                                        selected_monos_raw_samples[which(selected_monos_raw_samples$Sample_group %in% 
                                                                           c("raw_sample", "ANX_flowthrough")),])+
  geom_bar(aes(x=variable, y=value, fill=Sample_group), 
           stat="summary", color="#262626", position=position_dodge2(width = 0.5))+
  geom_point(aes(x=variable, y=value, fill=Sample_group, group=Sample_group),
             pch=21, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), size=2)+
  scale_fill_manual(name="", labels=c("Raw_sample", "ANX flowthrough"),
                    values=scico(2,alpha = 0.7, palette = "vikO", begin=0, end=0.5))+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Fucose signal")+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.y = element_text(size=font_size, family = font_family, color = "white"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.y = element_rect(fill = "#262626"))+
  facet_wrap(~Group_for_plot, nrow=3, strip.position = "right", scales = "free")


# latest #####
selected_monos_raw_samples<-melted_raw_samples[which(melted_raw_samples$variable %in%
                                                       c("Fucose", "Galactose", "Glucose",
                                                         "Xylose", "Glucuronic_acid")),]


flowthough_fig<-selected_monos_raw_samples[which(selected_monos_raw_samples$Sample_group %in% 
                                                   c("raw_sample", "ANX_flowthrough")),]
flowthough_fig$Sample_group<-factor(flowthough_fig$Sample_group, 
                                    levels=c("raw_sample", "ANX_flowthrough"),
                                    ordered = T)
flowthough_fig$variable<-factor(flowthough_fig$variable, 
                                levels=c("Fucose", "Glucuronic_acid", "Xylose","Galactose",
                                          "Glucose"), 
                                labels=c("Fucose", "Glucuronic\nacid","Xylose","Galactose",
                                          "Glucose"),
                                ordered = T)
flowthough_fig$Group_for_plot<-factor(flowthough_fig$Group_for_plot, 
                                      levels=c("initial", "dark_end", "light_end"),
                                      labels=c("Seawater", "Dark inc.", "Light inc."),
                                      ordered=T)


Fig_flowthrough<-ggplot(data=flowthough_fig)+
  geom_bar(aes(x=variable, y=value/20, fill=Sample_group), 
           stat="summary", color="#262626", position=position_dodge2(width = 0.3))+
  geom_point(aes(x=variable, y=value/20, fill=Sample_group, group=Sample_group),
             pch=21, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), size=2)+
  scale_fill_manual(name="", labels=c("Sample prior to AEX   ", "AEX flowthrough"),
                    values=scico(2,alpha = 0.8, palette = "oleron", begin=0.68, end=0.98))+
  scale_x_discrete(name="")+
  scale_y_continuous(name=expression(paste("Concentration (",mu,"g",~L^-1,")")))+
  theme_bw()+  
  theme(text=element_text(size=font_size, family = font_family, colour = "#262626"), 
        strip.text.y = element_text(size=font_size, family = font_family, color = "#262626"),
        legend.text = element_text(size=font_size, family = font_family),
        axis.text = element_text(size=font_size, family=font_family, colour = "#262626"),
        panel.grid = element_blank(),
        legend.position = "top",
        strip.background.y = element_rect(fill = c(time_cat_colors[1], 
                                                   time_cat_colors[2], time_cat_colors[3])),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "#262626"))+
  facet_wrap(~Group_for_plot, nrow=3, strip.position = "right", scales = "free")

g <- ggplot_gtable(ggplot_build(Fig_flowthrough))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- time_cat_colors
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210824_extraction.tiff", 
     height=10, width = 12, units = "cm", res=300)
grid.draw(g)
dev.off()
pdf("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210824_extraction.pdf", 
     height=4, width = 4.8)
grid.draw(g)
dev.off()


tiff("C:/Users/admin/ownCloud/Assemble_OC-VCE/Figures/20210802_Fig_raw samples monosaccharides and flowthrough.tiff",
     width=30, height=12, units = "cm", res=300)
plot_grid(Fig_rawsamples_selectedmonos, g, nrow = 1)
dev.off()
