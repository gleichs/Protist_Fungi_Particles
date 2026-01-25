### Protist + Fungi MetaT ###
### Fungi gene expression on particles ###

# Load libraries
library(edgeR)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(stringr)
library(KEGGREST)
library(ggpubr)

# Load in taxonomic IDs (dfEuk), functional IDs (dfEgg), and transcript abundances (dfSalmon)
dfEukProfun <- read.delim("./Data/eukulele_profun_eukprot.out",header=TRUE,sep="\t")
dfEgg <- read.delim("./Data/eggnog.annotations",header=TRUE,skip=3)
dfSalmon <- read.csv("./Data/salmon.csv",header=TRUE,row.names=1)

# Get taxonomic ID with the max pid.
dfEukProfun <- dfEukProfun[c("transcript_name","full_classification","max_pid")]

dfEukProfun <- dfEukProfun %>% arrange(desc(max_pid)) %>% distinct(transcript_name,.keep_all = TRUE) %>% as.data.frame() 
dfEukProfun$max_pid <- NULL

# Get the eggnog-mapper columns we'll be using for functional annotations: KEGG ko and CAZy
dfEgg <- dfEgg[c(1,9,16)]
dfEgg <- dfEgg[c("X.query_name","KEGG_ko","CAZy")]
colnames(dfEgg) <- c("transcript_name","KEGG","CAZy")

# Join EUKulele and eggnog-mapper tables and modify names of proteinIDs
dfEukProfunEgg <- left_join(dfEukProfun,dfEgg)

colz <- colsplit(dfEukProfunEgg$transcript_name,"\\.p",c("name","p"))
dfEukProfunEgg$transcript_name <- colz$name

colnames(dfSalmon)[1] <- "transcript_name"

# Join EUKulele + eggnog-mapper table and salmon (transcript counts) table.
dfProfun <- left_join(dfEukProfunEgg,dfSalmon)

# Grab all eukaryotic proteinIDs - we don't care about anything that is not eukaryotic.
dfProfun <- subset(dfProfun,grepl("Fungi",dfProfun$full_classification))

# Parse ProFun KEGG ko data
dfProfun$KEGG <- ifelse(is.na(dfProfun$KEGG),"",dfProfun$KEGG)
dfProfun$CAZy <- ifelse(is.na(dfProfun$CAZy),"",dfProfun$CAZy)
dfProfun$KEGG <- str_remove_all(dfProfun$KEGG,"ko:")
dfProfun$KEGG <- ifelse(is.na(dfProfun$KEGG),"",dfProfun$KEGG)
dfProfun$KEGG <- str_split(dfProfun$KEGG,",")
dfProfun$CAZy <- str_split(dfProfun$CAZy,",")
dfWideProfun <- unnest(dfProfun,KEGG)
dfWideProfun <- unnest(dfWideProfun,CAZy)
colnames(dfWideProfun)[1:4] <- c("Name","Taxonomy","KEGG","CAZy")
dfWideProfun  <-  dfWideProfun %>% distinct(.,.keep_all = TRUE) %>% as.data.frame()

# Separate Exp #1 and Exp #2 from Exp #3 (Experiment #3 was sequenced separately)
dfExp12P <- dfWideProfun[,c("Name","Taxonomy","KEGG","CAZy","RotT0_Exp1","RotT3_Exp1","RotT6_Exp1","RotT0_Exp2","RotT3_Exp2","RotT6_Exp2","WaterColumn1_Exp1","WaterColumn2_Exp1","WaterColumn1_Exp2","WaterColumn2_Exp2")]
dfExp3P <- dfWideProfun[,c("Name","Taxonomy","KEGG","CAZy","RotT0_Exp3","RotT3_Exp3","RotT6_Exp3","WaterColumn1_Exp3","WaterColumn2_Exp3")] 

# Normalize Exp 1&2 Libraries and Calculate CPM
dgeExp12_listP <- DGEList(counts=dfExp12P[5:14],genes=dfExp12P[1:4],group=c(rep("T0_Exp1",1),rep("T3_Exp1",1),rep("T6_Exp1",1),rep("T0_Exp2",1),rep("T3_Exp2",1),rep("T6_Exp2",1),rep("Water_Exp1",2),rep("Water_Exp2",2)))
dgeExp12_listP$samples # Visualize groups

Exp12NormP <- calcNormFactors(dgeExp12_listP,method = "TMM")
Exp12CpmP <- cpm(Exp12NormP, normalized.lib.sizes=TRUE, log=FALSE) 
Exp12CpmP <-as.data.frame(Exp12CpmP)  
Exp12CpmP <- data.frame(Exp12NormP$genes,Exp12CpmP)

# Normalize Exp 3 Libraries and Calculate CPM 
dgeExp3_listP <- DGEList(counts=dfExp3P[5:9],genes=dfExp3P[1:4],group=c(rep("T0_Exp3",1),rep("T3_Exp3",1),rep("T6_Exp3",1),rep("Water_Exp3",2)))
dgeExp3_listP$samples # Visualize groups

Exp3NormP <- calcNormFactors(dgeExp3_listP,method = "TMM")
Exp3CpmP <- cpm(Exp3NormP, normalized.lib.sizes=TRUE, log=FALSE) 
Exp3CpmP <-as.data.frame(Exp3CpmP)  
Exp3CpmP <-data.frame(Exp3NormP$genes,Exp3CpmP)

# Combine Normalized Exp 1&2 data with Normalized Exp 3 data
dfFinP <- left_join(Exp12CpmP,Exp3CpmP)

# Remove rows without KO terms
dfKO <- subset(dfFinP,KEGG !="")
dfKO$Name <- NULL
dfKO$Taxonomy <- NULL
dfKO$CAZy <- NULL
dfKO <- dfKO %>% group_by(KEGG) %>% summarize_all(sum)
dfKO <- dfKO %>% filter(if_all(everything(), ~ . != 0))
dfK0Melt <- melt(dfKO, id.vars="KEGG")
dfK0Melt <- subset(dfK0Melt,!grepl("Water",dfK0Melt$variable))

zscore_fxn <- function(df,exp){
  df <- subset(df,grepl(exp,df$variable))
  meanz <- df %>% group_by(KEGG) %>% summarize(m=mean(value),s=sd(value))
  df <- left_join(df,meanz)
  df$z <- (df$value-df$m)/df$s
  return(df)
}

z_exp1 <- zscore_fxn(dfK0Melt,"Exp1")
z_exp2 <- zscore_fxn(dfK0Melt,"Exp2")
z_exp3 <- zscore_fxn(dfK0Melt,"Exp3")
all <- rbind(z_exp1,z_exp2,z_exp3)

kegg <- unique(all$KEGG)
out <- NULL
for(i in 1:length(kegg)){
  s <- subset(all,KEGG==kegg[i])
  cols <- colsplit(s$variable,"_",c("Sample","Exp"))
  s <- cbind(s,cols)
  s$Sample <- ifelse(grepl("Water",s$Sample),"Water Column",s$Sample)
  s <- s %>% group_by(Sample,Exp) %>% summarize(m=mean(z))
  s1 <- subset(s,Exp=="Exp1")
  s2 <- subset(s,Exp=="Exp2")
  s3 <- subset(s,Exp=="Exp3")
  s1max <- subset(s1,m==max(s1$m))
  s2max <- subset(s2,m==max(s2$m))
  s3max <- subset(s3,m==max(s3$m))
  s1min <- subset(s1,m==min(s1$m))
  s2min <- subset(s2,m==min(s2$m))
  s3min <- subset(s3,m==min(s3$m))
  
  if(s1max$Sample==s2max$Sample & s2max$Sample==s3max$Sample & s1min$Sample==s2min$Sample & s2min$Sample==s3min$Sample){
    vec <- data.frame(v1=kegg[i],v2=s1max$Sample)
    out <- rbind(out,vec)  
  }
}


all <- subset(all,KEGG %in% out$v1)
cols <- colsplit(all$variable,"_",c("Sample","Exp"))
all <- cbind(all,cols)
all$Sample <- ifelse(grepl("Water",all$Sample),"Water Column",all$Sample)
all <- all %>% group_by(Sample,Exp,KEGG) %>% summarize(m=mean(z))

all$Sample <- ifelse(all$Sample=="RotT0","Net Trap T0",all$Sample)
all$Sample <- ifelse(all$Sample=="RotT3","Net Trap T3",all$Sample)
all$Sample <- ifelse(all$Sample=="RotT6","Net Trap T6",all$Sample)
all$Sample <- factor(all$Sample,levels=c("Water Column","Net Trap T0","Net Trap T3","Net Trap T6"))

all$Exp <- ifelse(all$Exp=="Exp1","Experiment #1 - 2021",all$Exp)
all$Exp <- ifelse(all$Exp=="Exp2","Experiment #2 - 2021",all$Exp)
all$Exp <- ifelse(all$Exp=="Exp3","Experiment #3 - 2022",all$Exp)

names <- data.frame(KEGG=c(out$v1),names=c("small subunit ribosomal protein S6e","heat shock 70kDa protein 1/6/8","ubiquitin B","clathrin heavy chain"),var=c(out$v2))
all <- left_join(all,names)
all$var <- factor(all$var,levels=c("Water Column","RotT0","RotT3","RotT6"))
all$var <- as.numeric(as.factor(all$var))

exp1 <-  all %>% filter(Exp=="Experiment #1 - 2021")%>%ggplot(aes(x=Sample,y=reorder(names,var),fill=m))+
  geom_tile(color="grey")+
  scale_fill_gradient2(name="Z-score Standardized Expression",low="#aa3900",mid="white",high="#6b1461",midpoint=0)+
  theme_bw(base_size=14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("")+
  ggtitle("Experiment #1 - 2021")

exp2 <-  all %>% filter(Exp=="Experiment #2 - 2021")%>%ggplot(aes(x=Sample,y=reorder(names,var),fill=m))+
  geom_tile(color="grey")+
  scale_fill_gradient2(name="Z-score Standardized Expression",low="#aa3900",mid="white",high="#6b1461",midpoint=0)+
  theme_bw(base_size=14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("")+
  ggtitle("Experiment #2 - 2021")

exp3 <-  all %>% filter(Exp=="Experiment #3 - 2022")%>%ggplot(aes(x=Sample,y=reorder(names,var),fill=m))+
  geom_tile(color="grey")+
  scale_fill_gradient2(name="Z-score Standardized Expression",low="#aa3900",mid="white",high="#6b1461",midpoint=0)+
  theme_bw(base_size=14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("")+
  ggtitle("Experiment #3 - 2022")



ggarrange(exp1,exp2,exp3,common.legend=TRUE,nrow=1,ncol=3)
ggsave("../../Figure7.png",width=18,height=6)

