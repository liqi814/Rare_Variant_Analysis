rm(list = ls())
library(data.table)
library(openxlsx)
library(tidyverse)

## setting 
setwd("D:/0000_LIST-GENO-VAR/2022-05_All_WGS_WES")
pop <- c("afr", "amr", "asj", "eas", "sas", "fin", "nfe")
pop2 <- pop[!pop == "asj"]
sequence_type <- c("WGS", "WES")

# write dataset
wb <- createWorkbook()
# define style
Color_code <- c("QVs meet at least one criteria (REVEL, PrimateAI and Polyphen)",
                "Dominant none variants that encoded fs, stop, or splice site variants",
                "QVs meet rareEnsemble model")
yellow_style <- createStyle(fgFill="#FFFF00")  # at leaste one criteria (REVEL, PrimateAI and Polyphen) 
# the dominant none variants that encoded fs, stop, or splice site variants (these I coded in orange as they are probably more deleterious)
orange_style <- createStyle(fgFill = "#FFA500") 
green_style <- createStyle(fgFill = "#90ee90") # RE2 variants

##################################################
geneList <- fread("FPF_candidate_genelist.txt", col.names = "List of gene checked in this study")
addWorksheet(wb, sheetName="Gene List")
writeData(wb, sheet="Gene List", x=geneList)

##################################################
addWorksheet(wb, sheetName="Color Code")
writeData(wb, sheet="Color Code", x=as.data.frame(Color_code))
setColWidths(wb, sheet="Color Code", cols = 1, widths = "auto")
addStyle(wb, sheet="Color Code", style=yellow_style, rows=2, cols=1, gridExpand=TRUE) # +1 for header line
addStyle(wb, sheet="Color Code", style=orange_style, rows=3, cols=1, gridExpand=TRUE) # +1 for header line
addStyle(wb, sheet="Color Code", style=green_style, rows=4, cols=1, gridExpand=TRUE) # +1 for header line

##################################################
## adding whole genome/exome sequencing samples
##################################################
for (sampletype in sequence_type) {
  caseData <- fread(list.files(pattern = paste0(sampletype, "_CasesAll.txt")))
  
  QV_RareEnsemble2 <- read.csv(list.files(pattern = paste0(sampletype, "-RareEnsemble2_genotypes.csv$"), recursive = TRUE))
  QV_dominantNone <- read.csv(list.files(pattern = paste0(sampletype, "-dominantNone_genotypes.csv$"), recursive = TRUE))
  
  ## remove redunant columns 
  filtdataDomNone <- QV_dominantNone %>% select(colnames(QV_dominantNone)[1:17], 
                                                gnomAD.Genome.global_AF, paste0("gnomAD.Exome.", pop, "_AF"),
                                                paste0("ExAC.", pop2, ".af"), Polyphen.Humdiv.Prediction, REVEL, PrimateAI, Sample.Name)
  filtdataRE2 <- QV_RareEnsemble2 %>% select(colnames(QV_RareEnsemble2)[1:17], 
                                             gnomAD.Genome.global_AF, paste0("gnomAD.Exome.", pop, "_AF"),
                                             paste0("ExAC.", pop2, ".af"), Polyphen.Humdiv.Prediction, REVEL, PrimateAI, Sample.Name)

  ##################################################
  addWorksheet(wb, sheetName=paste0(sampletype, " RareEnsembl QVs"))
  writeData(wb, sheet=paste0(sampletype, " RareEnsembl QVs"), x=filtdataRE2)
  
  ##################################################
  addWorksheet(wb, sheetName=paste0(sampletype, " DominantNone QVs"))
  writeData(wb, sheet=paste0(sampletype, " DominantNone QVs"), x=filtdataDomNone)
  
  # Polyphen.Humdiv.Prediction
  y <- which(colnames(filtdataDomNone)=="Polyphen.Humdiv.Prediction")
  x <- which(filtdataDomNone$Polyphen.Humdiv.Prediction == "probably")
  addStyle(wb, sheet=paste0(sampletype, " DominantNone QVs"), style=yellow_style, rows=x+1, cols=y, gridExpand=TRUE) # +1 for header line
  
  # REVEL
  y <- which(colnames(filtdataDomNone)=="REVEL")
  x <- which(filtdataDomNone$REVEL>=0.5)
  addStyle(wb, sheet=paste0(sampletype, " DominantNone QVs"), style=yellow_style, rows=x+1, cols=y, gridExpand=TRUE) # +1 for header line
  
  # PrimateAI
  y <- which(colnames(filtdataDomNone)=="PrimateAI")
  x <- which(filtdataDomNone$PrimateAI>=0.8)
  addStyle(wb, sheet=paste0(sampletype, " DominantNone QVs"), style=yellow_style, rows=x+1, cols=y, gridExpand=TRUE) # +1 for header line
  
  # the dominant none variants that encoded fs, stop (these I coded in orange as they are probably more deleterious)
  y <- which(colnames(filtdataDomNone)=="HGVS_p")
  x <- which(str_detect(filtdataDomNone$Effect, "frameshift|stop_gained"))
  addStyle(wb, sheet=paste0(sampletype, " DominantNone QVs"), style=orange_style, rows=x+1, cols=y, gridExpand=TRUE) # +1 for header line
  
  # the dominant none variants that encoded splice site variants (these I coded in orange as they are probably more deleterious)
  y <- which(colnames(filtdataDomNone)=="HGVS_c")
  x <- which(str_detect(filtdataDomNone$Effect, "splice_donor|splice_acceptor"))
  addStyle(wb, sheet=paste0(sampletype, " DominantNone QVs"), style=orange_style, rows=x+1, cols=y, gridExpand=TRUE) # +1 for header line
  
  # variant ID meet RareEnsemble criteria
  y <- which(colnames(filtdataDomNone)=="Variant.ID")
  x <- which(filtdataDomNone$Polyphen.Humdiv.Prediction == "probably" | filtdataDomNone$REVEL>=0.5 | filtdataDomNone$PrimateAI>=0.8)
  addStyle(wb, sheet=paste0(sampletype, " DominantNone QVs"), style=yellow_style, rows=x+1, cols=y, gridExpand=TRUE) # +1 for header line
  x <- which(str_detect(filtdataDomNone$Effect, "frameshift|stop_gained|splice_donor|splice_acceptor"))
  addStyle(wb, sheet=paste0(sampletype, " DominantNone QVs"), style=orange_style, rows=x+1, cols=y, gridExpand=TRUE) # +1 for header line
  x <- which(filtdataDomNone$Variant.ID %in% filtdataRE2$Variant.ID)
  addStyle(wb, sheet=paste0(sampletype, " DominantNone QVs"), style=green_style, rows=x+1, cols=y, gridExpand=TRUE) # +1 for header line
  
  ##################################################
  addWorksheet(wb, sheetName=paste0(sampletype, " No QVs"))
  NoQV_case <- caseData$V1[!caseData$V1 %in% filtdataDomNone$Sample.Name]
  writeData(wb, sheet=paste0(sampletype, " No QVs"), x=as.data.frame(NoQV_case))
}

##################################################
# write final result into a excel file
##################################################
saveWorkbook(wb, paste0(Sys.Date(), "_WGS_WES_genotypes.xlsx"), overwrite=TRUE)
