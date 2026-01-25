### Protist + Fungi MetaT ###
### Comparison of MMTESP vs. ProFun ###

# Load libraries
library(edgeR)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(stringr)
library(patchwork)

# Load in taxonomic IDs (dfEuk), functional IDs (dfEgg), and transcript abundances (dfSalmon)
dfEuk <- read.delim("./Data/eukulele_mmetsp.out",header=TRUE,sep="\t")
dfEukProfun <- read.delim("./Data/eukulele_profun_eukprot.out",header=TRUE,sep="\t")
dfEgg <- read.delim("./Data/eggnog.annotations",header=TRUE,skip=3)
dfSalmon <- read.csv("./Data/salmon.csv",header=TRUE,row.names=1)

# Get taxonomic ID with the max pid.
dfEuk <- dfEuk[c("transcript_name","full_classification","max_pid")]
dfEukProfun <- dfEukProfun[c("transcript_name","full_classification","max_pid")]

dfEuk <- dfEuk %>% arrange(desc(max_pid)) %>% distinct(transcript_name,.keep_all = TRUE) %>% as.data.frame() 
dfEuk$max_pid <- NULL
dfEukProfun <- dfEukProfun %>% arrange(desc(max_pid)) %>% distinct(transcript_name,.keep_all = TRUE) %>% as.data.frame() 
dfEukProfun$max_pid <- NULL

# Get the eggnog-mapper columns we'll be using for functional annotations: KEGG ko and CAZy
dfEgg <- dfEgg[c(1,9,16)]
dfEgg <- dfEgg[c("X.query_name","KEGG_ko","CAZy")]
colnames(dfEgg) <- c("transcript_name","KEGG","CAZy")

# Join EUKulele and eggnog-mapper tables and modify names of proteinIDs
dfEukEgg <- left_join(dfEuk,dfEgg)
dfEukProfunEgg <- left_join(dfEukProfun,dfEgg)

colz <- colsplit(dfEukEgg$transcript_name,"\\.p",c("name","p"))
dfEukEgg$transcript_name <- colz$name

colz <- colsplit(dfEukProfunEgg$transcript_name,"\\.p",c("name","p"))
dfEukProfunEgg$transcript_name <- colz$name

colnames(dfSalmon)[1] <- "transcript_name"

# Join EUKulele + eggnog-mapper table and salmon (transcript counts) table.
dfMmetsp <- left_join(dfEukEgg,dfSalmon)
dfProfun <- left_join(dfEukProfunEgg,dfSalmon)

# Grab all eukaryotic proteinIDs - we don't care about anything that is not eukaryotic.
dfMmetsp <- subset(dfMmetsp,grepl("Eukaryot",dfMmetsp$full_classification))
dfProfun <- subset(dfProfun,grepl("Eukaryot",dfProfun$full_classification))

# Only look at seqs that are annotated below "Eukaryota"
dfMmetsp <- subset(dfMmetsp,full_classification!="Eukaryota")
dfProfun <- subset(dfProfun,full_classification!="Eukaryota")

# Parse MMETSP KEGG ko data
dfMmetsp$KEGG <- ifelse(is.na(dfMmetsp$KEGG),"",dfMmetsp$KEGG)
dfMmetsp$CAZy <- ifelse(is.na(dfMmetsp$CAZy),"",dfMmetsp$CAZy)
dfMmetsp$KEGG <- str_remove_all(dfMmetsp$KEGG,"ko:")
dfMmetsp$KEGG <- ifelse(is.na(dfMmetsp$KEGG),"",dfMmetsp$KEGG)
dfMmetsp$KEGG <- str_split(dfMmetsp$KEGG,",")
dfMmetsp$CAZy <- str_split(dfMmetsp$CAZy,",")
dfWideMmetsp <- unnest(dfMmetsp,KEGG)
dfWideMmetsp <- unnest(dfWideMmetsp,CAZy)
colnames(dfWideMmetsp)[1:4] <- c("Name","Taxonomy","KEGG","CAZy")
dfWideMmetsp  <- dfWideMmetsp %>% distinct(.,.keep_all = TRUE) %>% as.data.frame()

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
dfExp12M <- dfWideMmetsp[,c("Name","Taxonomy","KEGG","CAZy","RotT0_Exp1","RotT3_Exp1","RotT6_Exp1","RotT0_Exp2","RotT3_Exp2","RotT6_Exp2","WaterColumn1_Exp1","WaterColumn2_Exp1","WaterColumn1_Exp2","WaterColumn2_Exp2")]
dfExp3M <- dfWideMmetsp[,c("Name","Taxonomy","KEGG","CAZy","RotT0_Exp3","RotT3_Exp3","RotT6_Exp3","WaterColumn1_Exp3","WaterColumn2_Exp3")] 

dfExp12P <- dfWideProfun[,c("Name","Taxonomy","KEGG","CAZy","RotT0_Exp1","RotT3_Exp1","RotT6_Exp1","RotT0_Exp2","RotT3_Exp2","RotT6_Exp2","WaterColumn1_Exp1","WaterColumn2_Exp1","WaterColumn1_Exp2","WaterColumn2_Exp2")]
dfExp3P <- dfWideProfun[,c("Name","Taxonomy","KEGG","CAZy","RotT0_Exp3","RotT3_Exp3","RotT6_Exp3","WaterColumn1_Exp3","WaterColumn2_Exp3")] 

# Normalize Exp 1&2 Libraries and Calculate CPM
dgeExp12_listM <- DGEList(counts=dfExp12M[5:14],genes=dfExp12M[1:4],group=c(rep("T0_Exp1",1),rep("T3_Exp1",1),rep("T6_Exp1",1),rep("T0_Exp2",1),rep("T3_Exp2",1),rep("T6_Exp2",1),rep("Water_Exp1",2),rep("Water_Exp2",2)))
dgeExp12_listM$samples # Visualize groups

dgeExp12_listP <- DGEList(counts=dfExp12P[5:14],genes=dfExp12P[1:4],group=c(rep("T0_Exp1",1),rep("T3_Exp1",1),rep("T6_Exp1",1),rep("T0_Exp2",1),rep("T3_Exp2",1),rep("T6_Exp2",1),rep("Water_Exp1",2),rep("Water_Exp2",2)))
dgeExp12_listP$samples # Visualize groups

Exp12NormM <- calcNormFactors(dgeExp12_listM,method = "TMM")
Exp12CpmM <- cpm(Exp12NormM, normalized.lib.sizes=TRUE, log=FALSE) 
Exp12CpmM <-as.data.frame(Exp12CpmM)  
Exp12CpmM <- data.frame(Exp12NormM$genes,Exp12CpmM)

Exp12NormP <- calcNormFactors(dgeExp12_listP,method = "TMM")
Exp12CpmP <- cpm(Exp12NormP, normalized.lib.sizes=TRUE, log=FALSE) 
Exp12CpmP <-as.data.frame(Exp12CpmP)  
Exp12CpmP <- data.frame(Exp12NormP$genes,Exp12CpmP)

# Normalize Exp 3 Libraries and Calculate CPM 
dgeExp3_listM <- DGEList(counts=dfExp3M[5:9],genes=dfExp3M[1:4],group=c(rep("T0_Exp3",1),rep("T3_Exp3",1),rep("T6_Exp3",1),rep("Water_Exp3",2)))
dgeExp3_listM$samples # Visualize groups

dgeExp3_listP <- DGEList(counts=dfExp3P[5:9],genes=dfExp3P[1:4],group=c(rep("T0_Exp3",1),rep("T3_Exp3",1),rep("T6_Exp3",1),rep("Water_Exp3",2)))
dgeExp3_listP$samples # Visualize groups

Exp3NormM <- calcNormFactors(dgeExp3_listM,method = "TMM")
Exp3CpmM <- cpm(Exp3NormM, normalized.lib.sizes=TRUE, log=FALSE) 
Exp3CpmM <-as.data.frame(Exp3CpmM)  
Exp3CpmM <-data.frame(Exp3NormM$genes,Exp3CpmM)

Exp3NormP <- calcNormFactors(dgeExp3_listP,method = "TMM")
Exp3CpmP <- cpm(Exp3NormP, normalized.lib.sizes=TRUE, log=FALSE) 
Exp3CpmP <-as.data.frame(Exp3CpmP)  
Exp3CpmP <-data.frame(Exp3NormP$genes,Exp3CpmP)

# Combine Normalized Exp 1&2 data with Normalized Exp 3 data
dfFinM <- left_join(Exp12CpmM,Exp3CpmM)
dfFinP <- left_join(Exp12CpmP,Exp3CpmP)

# Label datasets and combine
dfFinM$type <- "MMETSP"
dfFinP$type <- "ProFun"
dfAll <- rbind(dfFinM,dfFinP)
dfAll$Name <- NULL
dfAll$KEGG <- NULL
dfAll$CAZy <- NULL

# Determine taxonomic groups
dfAll$Tax <- ifelse(grepl("Fungi",dfAll$Taxonomy),"Fungi","Other Eukaryote")
dfAll$Tax <- ifelse(grepl("Metazoa",dfAll$Taxonomy),"Metazoa",dfAll$Tax)
dfAll$Taxonomy <- NULL
dfMelt <- melt(dfAll,id.vars=c("Tax","type"))

# Specify sample names
dfMelt$Sample <- ifelse(grepl("RotT0",dfMelt$variable),"Net Trap T0",NA)
dfMelt$Sample <- ifelse(grepl("RotT3",dfMelt$variable),"Net Trap T3",dfMelt$Sample)
dfMelt$Sample <- ifelse(grepl("RotT6",dfMelt$variable),"Net Trap T6",dfMelt$Sample)
dfMelt$Sample <- ifelse(grepl("Water",dfMelt$variable),"Water Column",dfMelt$Sample)
dfMelt$Exp <- ifelse(grepl("Exp1",dfMelt$variable),"Experiment #1",NA)
dfMelt$Exp <- ifelse(grepl("Exp2",dfMelt$variable),"Experiment #2",dfMelt$Exp)
dfMelt$Exp <- ifelse(grepl("Exp3",dfMelt$variable),"Experiment #3",dfMelt$Exp)

# Add reads across all sample types
dfMeltSum <- dfMelt %>% group_by(Tax,Sample,,Exp,type) %>% summarize(s=sum(value))
dfMeltSum$Tax <- factor(dfMeltSum$Tax,levels=c("Fungi","Metazoa","Other Eukaryote")) 

# Step 1: Calculate percentages for each Sample-type-Exp combo
temp <- dfMeltSum %>% 
  group_by(Sample, type, Exp) %>%
  mutate(percent = s / sum(s) * 100)

# Step 2: Average across Exp column
new <- temp %>%
  group_by(Sample, type, Tax) %>%
  summarise(
    avg_s = mean(s),
    avg_percent = mean(percent),
    .groups = "drop"
  ) %>%
  group_by(Sample, type) %>%  # Changed: removed Tax here
  mutate(
    label = paste0(round(avg_percent), "%"),
    ypos = cumsum(avg_s) - 0.5 * avg_s - 1  # Now cumsum works across all Tax
  )
new$ypos <- new$ypos-1

# Step 3: Plot fxn
plot_fxn <- function(samp,db,main){
  new %>% filter(Sample==samp & type==db) %>%
  ggplot(aes_string(x="1",y="avg_percent",fill="Tax"))+
  geom_bar(width = 1, stat = "identity",color="black") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = main)+
  scale_fill_manual(name="Taxonomy",values=c("grey85","grey30","white"),breaks=c("Fungi","Metazoa","Other Eukaryote"))+
  geom_text(data=subset(new,Tax=="Other Eukaryote" & Sample==samp & type==db),aes(y = 100, label = label), color = "black",vjust=5.2)
}

# Generate pie charts
m_wc <- plot_fxn("Water Column","MMETSP","Water Column\nMMETSP")
m_t0 <- plot_fxn("Net Trap T0","MMETSP","Net Trap T0\nMMETSP")
m_t3 <- plot_fxn("Net Trap T3","MMETSP","Net Trap T3\nMMETSP")
m_t6 <- plot_fxn("Net Trap T6","MMETSP","Net Trap T6\nMMETSP")
p_wc <- plot_fxn("Water Column","ProFun","Water Column\nProFun")
p_t0 <- plot_fxn("Net Trap T0","ProFun","Net Trap T0\nProFun")
p_t3 <- plot_fxn("Net Trap T3","ProFun","Net Trap T3\nProFun")
p_t6 <- plot_fxn("Net Trap T6","ProFun","Net Trap T6\nProFun")

# Combine plots
m_wc + m_t0 +m_t3+m_t6+p_wc+p_t0+p_t3+p_t6+plot_layout(ncol=4,nrow=2,guides="collect")+plot_annotation(tag_levels="a")
ggsave("../../Figure2.png",width=12,height=4)
