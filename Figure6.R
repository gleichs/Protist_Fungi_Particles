### Protist + Fungi MetaT ###
### Protist rel ab + gene expression ###

# Load libraries
library(edgeR)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(stringr)
library(KEGGREST)
library(patchwork)

norm_transcripts <- function(tax){  
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
  dfProfun <- subset(dfProfun,grepl(tax,dfProfun$full_classification))
  
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
  return(dfFinP)
}

getKOs <- function(df,tax){
  # Subset dataframe for taxon of interest
  dfTax <- subset(df,KEGG!="")
  dfTax$Name <- NULL
  dfTax$Taxonomy <- NULL
  dfTax$CAZy <- NULL
  
  # Sum CPM valuNULL# Sum CPM values for each KO
  dfTax <- dfTax %>% group_by(KEGG) %>% summarize_all(sum) %>% as.data.frame()
  
  # Melt dataframe
  dfMelt <- melt(dfTax,id.vars="KEGG")
  
  # Separate sample name and experiment number into different columns; Rename water column samples 
  dfNames <- colsplit(dfMelt$variable,"_",c("Sample","Experiment"))
  dfNames$Sample <- ifelse(grepl("Water",dfNames$Sample),"Water Column",dfNames$Sample)
  dfMelt <- cbind(dfMelt,dfNames)
  dfMelt$variable <- NULL
  
  # Find the KO terms that are not found in any samples from an entire experiment
  zeroKEGG <- dfMelt %>% group_by(Experiment,KEGG) %>% summarize(s=sum(value)) %>% filter(s==0) %>% as.data.frame()
  zeroKEGG <- unique(zeroKEGG$KEGG) 
  
  # Subset dataframe so that all KO terms are found in at least one sample per experiment
  `%ni%` <- Negate(`%in%`)
  dfMelt <- subset(dfMelt, KEGG %ni% zeroKEGG)
  
  # Calculate mean and sd for each KO for each experiment; will be used in the z-score standardization
  dfMean <- dfMelt %>% group_by(KEGG,Experiment) %>% summarize(mz=mean(value),sdz=sd(value)) %>% as.data.frame()
  dfMelt <- left_join(dfMelt,dfMean)
  
  # Calculate z-score
  dfMelt$z <- (dfMelt$value - dfMelt$mz)/dfMelt$sdz
  # range(dfMelt$z)
  
  # Calculate mean z-score for each KEGG, Sample, Exp (i.e., average water column duplicates)
  dfMelt <- dfMelt[c(1,3:4,7)]
  dfMeltAvg <- dfMelt %>% group_by(Experiment,Sample,KEGG) %>% summarize(mz=mean(z)) %>% as.data.frame()
  
  # List of all KO terms in dataframe to loop through
  keggTerms <- unique(dfMeltAvg$KEGG)
  
  # Establish empty dataframe
  out <- NULL
  
  # Loop through dataframe. For each KO term and experiment, find the Sample with the highest z-score.
  for (kegg in 1:length(keggTerms)){
    s <- subset(dfMeltAvg,KEGG==keggTerms[kegg])
    maxVals <- s %>% group_by(Experiment) %>% filter(mz==max(mz))
    if(length(unique(maxVals$Sample))==1){
      out <- rbind(out,maxVals)
    }
  }
  
  # Find all KO terms with max z-scores at specific time points (dependent on the temporal profile of taxa)
  if(tax=="Discob"){
    out <- subset(out,Sample=="RotT3") # Max discobid signal at T3
  }
  
  if(tax=="Cilioph"){
    out <- subset(out,Sample=="RotT6") # Max ciliate signal at T6
  }
  
  if(tax=="Laby"){
    out <- subset(out,Sample=="RotT3") # Max laby signal at T3
  }
  
  if(tax=="Foram"){
    out <- subset(out,Sample=="RotT0") # Max laby signal at T3
  }
  
  if(tax=="Amoeboz"){
    out <- subset(out,Sample=="RotT6") 
  }
  
  # Subset full dataset for only KO terms in out df
  subOut <- subset(dfMeltAvg,KEGG %in% out$KEGG)
  return(subOut)
}


calc_relab <- function(tax){
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
  dfAll <- subset(dfProfun,grepl("Eukaryot",dfProfun$full_classification))
  dfAll$tax <- ifelse(grepl(tax,dfAll$full_classification),"Tax","Other")
  dfAll$transcript_name <- NULL
  dfAll$full_classification <- NULL
  dfAll$KEGG <- NULL
  dfAll$CAZy <- NULL
  dfSum <- dfAll %>% group_by(tax) %>% summarize_all(sum) %>% as.data.frame()
  rownames(dfSum) <- dfSum$tax
  dfSum$tax <- NULL
  dfSum<- as.data.frame(t(dfSum))
  dfSum$frac <- dfSum$Tax/(dfSum$Tax+dfSum$Other)
  cols <- colsplit(rownames(dfSum),"_",c("Sample","Experiment"))
  dfFin <- cbind(dfSum,cols)
  dfFinMean <- dfFin %>% group_by(Experiment) %>% summarize(m=mean(frac),s=sd(frac))
  dfFin <- left_join(dfFin,dfFinMean)
  dfFin$z <- (dfFin$frac-dfFin$m)/dfFin$s
  dfFin$Sample <- ifelse(grepl("Water",dfFin$Sample),"Water Column",dfFin$Sample)
  dfFin <- dfFin %>% group_by(Experiment,Sample)%>% summarize(mz=mean(z))
  dfFin$day <- ifelse(dfFin$Sample=="Water Column",0,NA)
  dfFin$day <- ifelse(dfFin$Sample=="RotT0",1,dfFin$day)
  dfFin$day <- ifelse(dfFin$Sample=="RotT3",2,dfFin$day)
  dfFin$day <- ifelse(dfFin$Sample=="RotT6",3,dfFin$day)
  return(dfFin)
}

plot_fxn <- function(ko,relab,kegg,title){
  ko$day <- ifelse(ko$Sample=="Water Column",0,NA)
  ko$day <- ifelse(ko$Sample=="RotT0",1,ko$day)
  ko$day <- ifelse(ko$Sample=="RotT3",2,ko$day)
  ko$day <- ifelse(ko$Sample=="RotT6",3,ko$day)
  
  ko %>% filter(day!=0 & KEGG==kegg) %>% ggplot(aes(x=day,y=mz,color=Experiment,fill=Experiment))+
  # geom_smooth(method = "glm", formula = y ~ exp(x), se = FALSE)+
  theme_bw(base_size=14)+
  xlab("Sample")+
  ylab("Z-Score Standardized\nExpression")+
  scale_color_manual(name="Experiment",labels=c("Experiment #1 - 2021","Experiment #2 - 2021","Experiment #3 - 2022"),values=c("dodgerblue","darkgoldenrod2","indianred"))+
  scale_fill_manual(name="Experiment",labels=c("Experiment #1 - 2021","Experiment #2 - 2021","Experiment #3 - 2022"),values=c("dodgerblue","darkgoldenrod2","indianred"))+
  scale_x_continuous(labels=c("Water Column","Net Trap T0","Net Trap T3","Net Trap T6"),breaks=c(0,1,2,3))+
  geom_line(data=relab %>% filter(day!=0),aes(x=day,y=mz),size=1)+
  geom_jitter(size=3,color="black",shape=21,height = 0,width=0.2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
  ggtitle(title)
}


# Apply fxns

# Ciliates
normCil <- norm_transcripts("Cilio")
koOutCiliate <- getKOs(normCil,"Cilioph")
relCil <- calc_relab("Cilioph")
cil1 <- plot_fxn(koOutCiliate,relCil,'K00128',"Ciliate\nK00128 - aldehyde dehydrogenase (NAD+)")+theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10))
cil2 <- plot_fxn(koOutCiliate,relCil,'K01265',"Ciliate\nK01265 - methionyl aminopeptidase")+theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10))
cil3 <- plot_fxn(koOutCiliate,relCil,'K01369',"Ciliate\nK01369 - legumain")+theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10))

# Labyrinthulomycetes
normLaby <- norm_transcripts("Labyrinthu")
koOutLaby <- getKOs(normLaby,"Laby")
relLaby <- calc_relab("Labyrinthu")
laby1 <- plot_fxn(koOutLaby,relLaby,'K01692',"Labyrinthulomycetes\nK01692 - enoyl-CoA hydratase")+theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10))
laby2 <- plot_fxn(koOutLaby,relLaby,'K01647',"Labyrinthulomycetes\nK01647 - citrate synthase")+theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10))
laby3 <- plot_fxn(koOutLaby,relLaby,'K02133',"Labyrinthulomycetes\nK02133 - F-type H+-transporting ATPase")+theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10))

# Amoebozoans
normAmoe <- norm_transcripts("Amoeboz")
koOutAmoe <- getKOs(normAmoe,"Amoeboz")
relAmoe <- calc_relab("Amoeboz")
amoe1 <- plot_fxn(koOutAmoe,relAmoe,'K01279',"Amoebozoa\nK01279 - tripeptidyl-peptidase I")+theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10))
amoe2 <- plot_fxn(koOutAmoe,relAmoe,'K04646',"Amoebozoa\nK04646 - clathrin heavy chain")+theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10))
amoe3 <- plot_fxn(koOutAmoe,relAmoe,'K04409',"Amoebozoa\nK04409 - p21-activated kinase 1")+theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10))

# Foraminifera
normForam <- norm_transcripts("Foramin")
koOutForam <- getKOs(normForam,"Foramin")
relForam <- calc_relab("Foramin")
foram1 <- plot_fxn(koOutForam,relForam,'K07374',"Foraminifera\nK07374 - tubulin alpha")+theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10))
foram2 <- plot_fxn(koOutForam,relForam,'K01304',"Foraminifera\nK01304 - pyroglutamyl-peptidase")+theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10))
foram3 <- plot_fxn(koOutForam,relForam,'K19942',"Foraminifera\nK19942 - dynein regulatory complex subunit 4")+theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10))

foram1+foram2+foram3+dis1+dis2+dis3+laby1+laby2+laby3+cil1+cil2+cil3+amoe1+amoe2+amoe3+plot_layout(ncol=3,nrow=5,guides="collect",axes="collect")+plot_annotation(tag_levels = "a")
ggsave("../../Figure6.png",width=20,height=14)
