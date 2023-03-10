library(vroom)
library(minfi)
library(illuminaio)
library(devtools)
library(dplyr)

file <- "EPIC-8v2-0_A1.csv"

e1 <- vroom(file)

control.line <- grep("Controls",e1$Illumina)+1
assay.line <- grep("\\[Assay",e1$Illumina)+1

rm(e1)
gc()

manifest <- vroom(file,skip=assay.line,n_max =control.line-assay.line-2) 

# remove mgmt dupliated probes
manifest <- manifest[-which(manifest$IlmnID %in% c("cg12434587_BO11","cg12981137_BO11","cg12981137_TC21","cg12981137_BC21")),]

# check remaining duplicates and select those in manifest$manifest_probe_match == TRUE when possible
duplicated_probes <- manifest$Name[duplicated(manifest$Name)]

remove_duplicated <- function(probe){
  sel <- which(manifest$Name == probe)
  if(any(manifest[sel,"Manifest_probe_match"])){ 
    sel_ID <- manifest$IlmnID[sel[which(manifest[sel,"Manifest_probe_match"]%>%pull())]][1]
  }else{
    sel_ID <- manifest$IlmnID[sel][1]
  }
  rmv_ID <- setdiff(manifest$IlmnID[sel],sel_ID)
  return(rmv_ID)
}

rmv_IlmnID <- unlist(sapply(duplicated_probes,remove_duplicated))

manifest <- manifest[-which(manifest$IlmnID %in% rmv_IlmnID),]

any(duplicated(manifest$Name))

manifest$AddressA_ID <- gsub("^0*", "", manifest$AddressA_ID)
manifest$AddressB_ID <- gsub("^0*", "", manifest$AddressB_ID)

manifest$AddressA_ID[is.na(manifest$AddressA_ID)] <- ""
manifest$AddressB_ID[is.na(manifest$AddressB_ID)] <- ""

TypeI <- manifest[manifest$Infinium_Design_Type == "I",
                  c("Name", "AddressA_ID", "AddressB_ID", "Color_Channel", "Next_Base",
                    "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")]

names(TypeI)[c(2, 3, 4, 5, 6, 7)] <- c("AddressA", "AddressB", "Color", "NextBase", "ProbeSeqA", "ProbeSeqB")

TypeI <- as(TypeI, "DataFrame")
TypeI$ProbeSeqA <- DNAStringSet(TypeI$ProbeSeqA)
TypeI$ProbeSeqB <- DNAStringSet(TypeI$ProbeSeqB)
TypeI$NextBase <- DNAStringSet(TypeI$NextBase)
TypeI$nCpG <- as.integer(
  oligonucleotideFrequency(TypeI$ProbeSeqB, width = 2)[, "CG"] - 1L)
TypeI$nCpG[TypeI$nCpG < 0] <- 0L
TypeSnpI <- TypeI[grep("^rs", TypeI$Name), ]
TypeI <- TypeI[-grep("^rs", TypeI$Name), ]
TypeI <- TypeI[-grep("^nv", TypeI$Name), ]

TypeII <- manifest[
  manifest$Infinium_Design_Type == "II",
  c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
names(TypeII)[c(2,3)] <- c("AddressA", "ProbeSeqA")
TypeII <- as(TypeII, "DataFrame")
TypeII$ProbeSeqA <- DNAStringSet(TypeII$ProbeSeqA)
TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeqA, letters = "R"))
TypeII$nCpG[TypeII$nCpG < 0] <- 0L
TypeSnpII <- TypeII[grep("^rs", TypeII$Name), ]
TypeII <- TypeII[-grep("^rs", TypeII$Name), ]
TypeII <- TypeII[-grep("^nv", TypeII$Name), ]

controls <- read.table(
  file = file,
  skip = control.line,
  sep = ",",
  comment.char = "",
  quote = "",
  colClasses = c(rep("character", 5)))[, 1:5]

TypeControl <- controls[, 1:4]
names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
TypeControl <- as(TypeControl, "DataFrame")


maniTmp <- list(
  manifestList = list(
    TypeI = TypeI,
    TypeII = TypeII,
    TypeControl = TypeControl,
    TypeSnpI = TypeSnpI,
    TypeSnpII = TypeSnpII),
  manifest = manifest,
  controls = controls)


# checks
manifest <- maniTmp$manifest

epic <- readIDAT("./206891110005/206891110005_R02C01_Grn.idat")

address.epic <- as.character(epic$MidBlock)

dropCpGs <- manifest$Name[manifest$AddressB_ID != "" & !manifest$AddressB_ID %in% address.epic]
table(substr(dropCpGs, 1,2))

dropCpGs <- manifest$Name[manifest$AddressA_ID != "" & !manifest$AddressA_ID %in% address.epic]
table(substr(dropCpGs, 1,2))

## Controls ok
maniTmp$controls[!maniTmp$manifestList$TypeControl$Address %in% address.epic,]

## Manifest package
maniList <- maniTmp$manifestList

IlluminaHumanMethylationEPICv2manifest <- IlluminaMethylationManifest(TypeI = maniList$TypeI,
                                                              TypeII = maniList$TypeII,
                                                              TypeControl = maniList$TypeControl,
                                                              TypeSnpI = maniList$TypeSnpI,
                                                              TypeSnpII = maniList$TypeSnpII,
                                                              annotation = "IlluminaHumanMethylationEPICv2")

use_data(IlluminaHumanMethylationEPICv2manifest, internal=TRUE, overwrite = T)


