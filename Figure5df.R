### Protist + Fungi MetaT ###
### Taxa barplot fungi: Net T0-T6 ###

# Load libraries
library(edgeR)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(stringr)
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

# Remove anything that can't be classified below "Eukaryota"
dfTax <- subset(dfFinP,Taxonomy!="Eukaryota")

# Group by taxonomic units
taxz <- colsplit(dfTax$Taxonomy,";",c("d","k","p","c","o","f","g","s"))
taxz$c <- str_remove_all(taxz$c,"^ ")
taxz$c <- ifelse(taxz$c=="","Unknown Fungi",taxz$c)

# Melt data and summarize
dfTax$tax <- taxz$c
dfTax$Taxonomy <- NULL
dfTax <- dfTax[c(4:19)]
dfTaxMelt <- melt(dfTax,id.vars="tax")
cols <- colsplit(dfTaxMelt$variable,"_",c("Sample","Exp"))
dfTaxMelt <- cbind(dfTaxMelt,cols)
dfTaxMelt$variable <- NULL
dfTaxMelt <- dfTaxMelt %>% group_by(Sample,Exp,tax)%>% summarize(s=sum(value))

# Average WC together
dfTaxMelt$Sample <- ifelse(grepl("Water",dfTaxMelt$Sample),"Water Column",dfTaxMelt$Sample)
dfTaxMelt <- dfTaxMelt %>% group_by(Sample,Exp,tax)%>% summarize(s=mean(s))


# Colors for taxonomy
colrs <- randomcoloR::distinctColorPalette(length(unique(taxz$c)))

colrs <- setNames(colrs, unique(taxz$c))

# Plot function
taxPltFxnNets <- function(df,exp,title){
  subs <- subset(df, Exp==paste(exp))
  subsWater <- subset(subs,grepl("Water",subs$Sample))
  subsWater <- subsWater %>% group_by(tax) %>% summarize(s=mean(s))
  subsWater$Sample <- "Water Column"
  subs <- subset(subs,!grepl("Water",subs$Sample))
  subs$variable <- NULL
  subs$Exp <- NULL
  subs$Sample <- ifelse(subs$Sample=="RotT0",0,subs$Sample)
  subs$Sample <- ifelse(subs$Sample=="RotT3",3,subs$Sample)
  subs$Sample <- ifelse(subs$Sample=="RotT6",6,subs$Sample)
  
  rotPlot <- ggplot(subs,aes(x=as.numeric(Sample),y=s,fill=tax))+geom_area(alpha = 0.6, position = 'fill')+geom_col(width = 1.5, color = 'gray20', position = 'fill')+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+theme_classic(base_size = 16)+xlab("")+ylab("Relative Transcript Abundance")+theme(axis.text.x = element_text(angle = 45, hjust =1))+scale_x_continuous(breaks=c(0,3,6),labels=c("Net Trap T0","Net Trap T3","Net Trap T6"))
  
  out <- rotPlot+plot_layout(guides = "collect",nrow=1)+plot_annotation(title=paste(title))
  return(out)}

# Execute plot function
pExp1 <- taxPltFxnNets(dfTaxMelt,"Exp1"," Experiment #1\n 2021")
pExp2 <- taxPltFxnNets(dfTaxMelt,"Exp2"," Experiment #2\n 2021")
pExp3 <- taxPltFxnNets(dfTaxMelt,"Exp3"," Experiment #3\n 2022")

# Make combined ggarrange plot
ggarrange(pExp1,pExp2,pExp3,common.legend = TRUE,nrow=1,ncol=3,labels=c("d","e","f"),font.label = list(size=12,color="black",face="plain"))
ggsave("../../Figure5df.png",width=13,height=7)
