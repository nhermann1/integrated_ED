
#This script is created to input data files both of raw records and summarized averages
#Additional quality assurances are included for checking data completeness, see Metadata document included with publication
#Finally, a figure comparing distributions and medians by phyla is generated

#Code operated with R version 4.3.2 "Eye Holes"
#Created by NTH
#Last updated 12/19/24 (m/d/y)



#insert file path to the downloaded zip folder 
setwd("")

#Packages needed for this script to function
library(tidyverse)  #Version 1.3.1
library(ggplot2)    #Version 3.5.1
library(readr)      #Version 2.0.1
library(gridExtra)  #Version 2.3

#Personalization
theme_set(theme_bw(base_size=25))
options(dplyr.summarise.inform = FALSE,scipen=999)

#Read the raw database file into the environment
preyEDs<-read_csv("preyEDs_integrated_12_24.csv",col_types = list(Indigestible_perc=col_number()))

#Verification against the database supplied
colSums(preyEDs[,sapply(preyEDs,is.numeric)],na.rm=T)


#Viewing summaries for each phylum
preyEDs%>%group_by(Phylum)%>%mutate(Wet_KJg_min=ifelse(is.na(Wet_KJg_min),Wet_KJg_mean,Wet_KJg_min),
                                    Wet_KJg_max=ifelse(is.na(Wet_KJg_max),Wet_KJg_mean,Wet_KJg_max))%>%
  reframe(n=n_distinct(Taxa),
          m=mean(Wet_KJg_mean,na.rm=T),
          s=sd(Wet_KJg_mean,na.rm=T),
          l=min(Wet_KJg_min,na.rm=T),
          u=max(Wet_KJg_max,na.rm=T))%>%arrange(-n)

#### Plotting the Phyla energy density distributions ####

#Creating a summary dataframe for each phylum
phylaEDs<-preyEDs%>%
  filter(Phylum!="N/A")%>%
  group_by(Phylum)%>%
  reframe(n_spec=n_distinct(Taxa),
          n_samp=sum(N_replicates,na.rm=T),
          Wet_KJg=mean(Wet_KJg_mean,na.rm=T),
          Wet_SD=sd(Wet_KJg_mean,na.rm=T),
          max=max(Wet_KJg_mean,na.rm=T),
          median=median(Wet_KJg_mean,na.rm=T))%>%
  arrange(Wet_KJg)%>%
  mutate(height=seq(87,100,length=17),
         Phylum=fct_reorder(factor(Phylum),n_spec))%>%
  filter(Wet_KJg>0)


#Creating a violin plot with each phylum distribution and datapoints
phylaViolins<-left_join(preyEDs,phylaEDs)%>%mutate(Phylum=fct_reorder(factor(Phylum),n_spec,.desc=T))%>%
  filter(Phylum%in%phylaEDs$Phylum)%>%
  mutate(tissue_shape=case_when(grepl("whole",Sample_Type,ignore.case=T)~"Whole",
                                grepl("liver",Sample_Type,ignore.case=T)~"Liver",
                                grepl("muscle",Sample_Type,ignore.case=T)~"Muscle",
                                Sample_Type=="Roe"~"Egg",
                                Sample_Type=="Pooled"~"Whole",
                                Sample_Type=="Tail"~"Muscle",
                                T~"Unknown"),
         tissue_shape=factor(tissue_shape,levels=c("Whole","Muscle","Liver","Egg","Unknown")))%>%
  ggplot()+
  geom_jitter(aes(x=Phylum, y=Wet_KJg_mean,fill=tissue_shape,shape=tissue_shape),
              position=position_jitterdodge(0.2),color="black",alpha=0.85,size=2.5)+
  geom_violin(aes(x=Phylum, y=Wet_KJg_mean),fill="transparent",color="red",
              show.legend=F,draw_quantiles = c(0.5),scale="width",linewidth=2)+
  geom_violin(aes(x=Phylum, y=Wet_KJg_mean),fill="transparent",color="black",
              show.legend=F,scale="width",linewidth=2)+
  geom_text(data=filter(phylaEDs,max>0),aes(Phylum,y=max+0.5,label=paste0("N = ",n_spec)),
            size=5,fontface="bold",hjust=0)+
  scale_y_continuous(name="Energy Density (kJ/g WW)",limits=c(0,20),expand=expansion(add=c(0.1,0.5)))+
  scale_x_discrete(expand=expansion(add=c(0.5,1)))+
  scale_shape_manual(values=c(21:25),name="Tissue Type\nfor Energy Density\nEstimation")+
  scale_fill_manual(values=c(viridis::viridis_pal()(4),"black"),
                    name="Tissue Type\nfor Energy Density\nEstimation")+
  guides(shape=guide_legend(override.aes=list(alpha=1,size=10)))+
  theme(axis.text.x=element_text(size=25),panel.border=element_blank(),
        axis.line.x.bottom=element_line(linewidth=2),axis.line.x.top=element_line(linewidth=2),
        axis.line.y.right = element_blank(),axis.line.y.left=element_line(linewidth=2),
        plot.margin=unit(c(-1,1,1,1),"cm"),
        legend.position=c(0.7,0.8),legend.background=element_rect(color="black"))+
  coord_flip(clip="off")


totalViolin<-preyEDs%>%
  mutate(tissue_shape=case_when(grepl("whole",Sample_Type,ignore.case=T)~"Whole",
                                grepl("liver",Sample_Type,ignore.case=T)~"Liver",
                                grepl("muscle",Sample_Type,ignore.case=T)~"Muscle",
                                Sample_Type=="Roe"~"Egg",
                                Sample_Type=="Pooled"~"Whole",
                                Sample_Type=="Tail"~"Muscle",
                                T~"Unknown"),
         tissue_shape=factor(tissue_shape,levels=c("Whole","Muscle","Liver","Egg","Unknown")))%>%
  ggplot()+
  geom_jitter(aes(1, y=Wet_KJg_mean,fill=tissue_shape,shape=tissue_shape), 
              position=position_jitterdodge(0.25),show.legend=F,alpha=0.85,size=2.5)+
  geom_violin(aes(1, y=Wet_KJg_mean),color="red", fill="transparent",show.legend=F,
              draw_quantiles = c(0.5),scale="width",linewidth=2)+
  geom_violin(aes(1, y=Wet_KJg_mean),color="black", fill="transparent",show.legend=F,
              scale="width",linewidth=2)+
  geom_text(aes(1,y=18.5,label=paste0("N = ",n_distinct(preyEDs$Taxa))),
            size=5,fontface="bold",hjust=0)+
  scale_shape_manual(values=c(21:25),name="Tissue Type\nfor Energy Density\nEstimation")+
  scale_fill_manual(values=c(viridis::viridis_pal()(4),"black"),
                    name="Tissue Type\nfor Energy Density\nEstimation")+
  scale_y_continuous(limits=c(0,20),expand=expansion(add=c(0.1,0.5)))+
  scale_x_continuous(name=" ",breaks=c(1),labels=c("            Total"))+
  theme(panel.border = element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),
        axis.title=element_blank(),axis.text.y=element_text(size=30,face="bold"),
        axis.line.y.left=element_line(linewidth=2),axis.line.x.top=element_line(linewidth=2),
        axis.line.x.bottom = element_line(size=1,linetype="dashed"),axis.line.y.right=element_line(linewidth=2),
        plot.margin=unit(c(1,1,-0.5,1.04),"cm"))+
  coord_flip(clip="off")

grid.arrange(phylaViolins,totalViolin,ncol=1,layout_matrix = t(t(c(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))))




# Summarized database

#Read in the summarized database
tableSummaryEDs<-read_csv("preyEDs_taxaSummarized.7.24.csv")

#Verification against the database supplied
colSums(tableSummaryEDs[,sapply(tableSummaryEDs,is.numeric)],na.rm=T)



