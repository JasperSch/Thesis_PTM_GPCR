#script to further filter dataset of good files
#maps PTMs to the files
#Provides tables with information regarding these PTMs

library(readr)
library(rpart)
library(rpart.plot)
library(hash)
library(stringr)
library(vcd)
library(doBy)
library(readxl)
library(xlsx)

DS <- read_csv("~/School/2e_master_bioinformatica/Thesis/tables/Overview_GPROT_raw", 
                col_types = cols(`ENGINEERED_MUTATION` = col_logical(), 
                                 RESOLUTION = col_double(),
                                 GAP_AVERAGE = col_double(),
                                 FUSIONCHAINS = col_logical(), 
                                 MONOMERIC = col_logical(),
                                 TRANSMEMBRANE = col_logical(), 
                                 RESOLUTION_FINE= col_logical(),
                                 HAS_PROTEIN_STRUCTURE = col_logical(),
                                 COVERAGE_AVERAGE = col_double(),
                                 COVERAGE_GOOD = col_logical()))

PhosphorylationDS <- read_delim("~/School/2e_master_bioinformatica/Thesis/tables/Phosphorylation_dbPTM_20181122.txt", 
                                "\t", escape_double = FALSE, col_names = FALSE, 
                                trim_ws = TRUE)
colnames(PhosphorylationDS) = c('name',  "Uniprot Accession","Position", "Type","PMIDs","Sequence")
Homo_sapiensDS<- read_delim("~/School/2e_master_bioinformatica/Thesis/tables/Homo_sapiens_PLMD_20181122.elm", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
AcetylationDS <- subset(Homo_sapiensDS, Homo_sapiensDS$Type == 'Acetylation')


#make dataframe of database identifiers + UNP ID's
  DB_ID = str_split_fixed(DS$DB_IDENTIFIERS, " ", length(DS$ID) )
  DB_ID = as.data.frame(DB_ID)

  #Check from which collumn empty collumns begin and remove all empty ones
  for (i in length(DB_ID[,1]):1){
   colNumber = i
    if (paste(DB_ID[,i], collapse = '') != "")
     break
  }                             
  DB_ID = DB_ID[1:colNumber]
  DB_ID = cbind(DB_ID, ID = DS$ID)

  #Specialize dataframe to only include uniprot identifiers and count the amount of chains
  DB_ID_UNP <- as.data.frame(lapply(DB_ID, function(x) gsub("[0-9a-zA-Z_]+\\((?!UNP)[0-9a-zA-Z_]+\\)", "", x, perl=TRUE))) #removes al ID's except UNP ID's
  DB_ID_UNP <- as.data.frame(lapply(DB_ID_UNP, function(x) gsub("\\(UNP\\)", "", x, perl=TRUE)))
  DB_ID_UNP <- data.frame(lapply(DB_ID_UNP, as.character), stringsAsFactors=FALSE)
  DB_ID_UNP$Entry_count <- apply(DB_ID_UNP, 1, function(x) sum(x!= "")-1) #Minus 1 to not include the ID-cell. *****************BEWARE! DO NOT RUN THIS TWICE***********


#Get counts for PTMs
  i=1
  ac_count = vector()
  ph_count = vector()
  hasPTM = vector()
  for (i in 1:length(DB_ID_UNP[,1])){
    id_unique =apply(DB_ID_UNP[i,],1, function(x) unique(x))
    ac =lapply(id_unique[,1], function(x) sum(AcetylationDS$`Uniprot Accession` == x))
    ph =lapply(id_unique[,1], function(x) sum(PhosphorylationDS$`Uniprot Accession` == x))
    ac_count[i] = sum(unlist(ac))
    ph_count[i] = sum(unlist(ph))
  }

  hasPTM= ac_count+ ph_count > 0
  DB_ID_UNP = cbind(DB_ID_UNP, ac_count)
  DB_ID_UNP = cbind(DB_ID_UNP, ph_count)
  DB_ID_UNP = cbind(DB_ID_UNP, hasPTM)
  DS= cbind(DS, ac_count)
  DS = cbind(DS, ph_count)
  DS = cbind(DS, hasPTM)
  
#Get subset of pdb files with PTMs
  DS_PTM = subset(DS, DS$hasPTM == TRUE)
  DB_ID_UNP_PTM = subset(DB_ID_UNP, DB_ID_UNP$hasPTM == TRUE)

#Get the overlaplist for final selection
  #Making the list with overlaps
  Overlaplist = data.frame(matrix(ncol = 3),stringsAsFactors = FALSE)
  Overlaplist= Overlaplist[!is.na(Overlaplist)]
  
  #Searching every unique UNP identifier
  chains = DB_ID_UNP_PTM[,1]
  for (i in 2:colNumber){
    chains = c(as.character(chains), as.character(DB_ID_UNP_PTM[,i]))
  }
  chains = unique(chains[chains !=""])

  #For every unique UNP > search lines with UNP > extract structures > Compare structures and put in overlaplist
  total <- length(chains)
  
  for (chain in chains){
    #Find all unqiue structures in DB_ID_UNP_PTM that have the chain
    matches = DB_ID_UNP_PTM[grep(chain, DB_ID_UNP_PTM[,1]),]
    for (i in 1:colNumber){
      matches = rbind(matches, DB_ID_UNP_PTM[grep(chain, DB_ID_UNP_PTM[,i+1]),])
    }
    matches = unique(matches)
    
    #If only 1 structure has the UNP, no need to start comparing structures
    if (length(matches$V1)==1){
    }
    else{
      #compare structures pairwaise and calculate overlap
      for (i in 1:(length(matches$V1)-1)){
        x = matches[i,]
        IDx = toString(x$ID)
        UNPx = x[1,1:colNumber]
        UNPx = UNPx[UNPx!=""]
        for (j in (i+1):length(matches$V1)){
          y = matches[j,]
          IDy = toString(y$ID)
          UNPy = y[1,1:colNumber]
          UNPy = UNPy[UNPy!=""]
          
          s1 = sum(is.element(UNPx,UNPy))
          s2 = sum(is.element(UNPy,UNPx))
          count = min(s1, s2) #This way 1 chain in 1 structure gets mapped to only 1 identical chain in the outher structure
          
          Totalchains = x$Entry_count + y$Entry_count
          overlap = count/Totalchains*2
          newrow = c(IDx, IDy, overlap)
          Overlaplist = rbind(Overlaplist, newrow)
        }
      }
    }
  }
  #sort collumns alphabetically and remove redundant information
  x = paste(Overlaplist[,1], Overlaplist[,2], sep= " ")
  Overlaplist[,1] = as.vector(unlist(lapply(x, function(x){ str_sort(str_split(x, " ")[[1]])[1]})))
  Overlaplist[,2] = as.vector(unlist(lapply(x, function(x){ str_sort(str_split(x, " ")[[1]])[2]})))
  Overlaplist = unique(Overlaplist)

#Manual futher selection (see powerpoint slideset GPCR_strcutre_selection)
  sel = c("1DS6", "1J1B","1NR4","1OW3","2ACX","2GTP","2J7Z","3BBC","4IP8","4XO6","5L9E",
          "1FL7","1KV3","1O7Z","1Z92","2GDZ","2IK8","2ODE","3VZD","4L23","5HK1","6F2U")
  DS_sel = subset(DS, DS$ID %in% sel, select= c(1,3:9,15,16,19,20))
  
#meta data about selection
  PDBs = length(sel)
  AA_count = 0
  AA_range = 0
  for(i in 1:length(DS_sel$ID)){
  AA_count = AA_count + sum(as.integer(unlist(str_split(DS_sel$PRESENT_AA[i], " "))))
  AA_range = AA_range + sum(as.integer(unlist(str_split(DS_sel$PRESENT_AA_RANGE[i], " "))))
  }
  DB_ID_UNP_sel = subset(DB_ID_UNP, ID %in% sel)
  x = which(colnames(DB_ID_UNP_sel)=="ID" )
  x= unique(unlist(as.list(DB_ID_UNP_sel[,1:x-1])))
  UNP_IDs = x[x != ""]
  UniqueProteins = length(UNP_IDs)
  PhosphorylationDS_sel = subset(PhosphorylationDS, select= 2:6, subset = PhosphorylationDS$`Uniprot Accession` %in% UNP_IDs)
  AcetylationDS_sel = subset(AcetylationDS, select= c(2:5,7), subset = AcetylationDS$`Uniprot Accession` %in% UNP_IDs)
  Ac= length(AcetylationDS_sel$`Uniprot Accession`)
  Ph= length(PhosphorylationDS_sel$`Uniprot Accession`)
  PTMs = Ac + Ph
  PMIDs= vector()
  for(i in UNP_IDs){
    x = subset(AcetylationDS, select = PMIDs, `Uniprot Accession` == i)
    y = subset(PhosphorylationDS, select = 5, `Uniprot Accession` == i)
    PMIDs = c(PMIDs, unlist(str_split(paste(unlist(x), collapse= ";"), ";")))
    PMIDs = c(PMIDs, unlist(str_split(paste(unlist(y), collapse= ";"), ";")))
  }
  PMIDs = unique(PMIDs[PMIDs !=""])
  PTM_PMIDs = length(PMIDs)
  DS_sel_metadata = data.frame(PDBs, UniqueProteins, PTMs, Ac, Ph, AA_count, AA_range, PTM_PMIDs,
                               UNP_IDs =paste(UNP_IDs,collapse=" "),
                               PMIDs =paste(PMIDs,collapse=" "))
  
  PMIDs
 
#Making PTM data frame after PMID search of PTMS and updating metadata
  PMIDs_stressNormal <- read_excel("~/School/2e_master_bioinformatica/Thesis/tables/PMIDs_stressNormal.xlsx")
  PMIDs_normal = subset(PMIDs_stressNormal, PMIDs_stressNormal$Condition %in% c("Normal", "Both"))
  PMIDs_treated = subset(PMIDs_stressNormal, PMIDs_stressNormal$Condition %in% c("Treated", "Both"))
  DS_sel_PTMdata = merge(AcetylationDS_sel, PhosphorylationDS_sel, all= TRUE)
  

  isNormal <- function(PMIDs, PMIDs_normal){
    x = unlist(str_split(PMIDs, ";"))
    if (length(intersect(x, PMIDs_normal))>0){
      return (TRUE)
    }
    else{
      return(FALSE)
    }
  }
  

  isTreated <- function(PMIDs, PMIDs_treated){
    x = unlist(str_split(PMIDs, ";"))
    if (length(intersect(x, PMIDs_treated))>0){
      return (TRUE)
    }
    else{
      return(FALSE)
    }
  }
  
  condition <- function(treated, normal){
    if (treated & normal){
      return ("both")
    }
    if(treated){
      return ("treated")
    }
    if(normal){
      return ("normal")
    }
    else{
      return("unknown")
    }
  }
  
  DS_sel_PTMdata$normal = sapply(DS_sel_PTMdata$PMIDs, isNormal,PMIDs_normal$PMID)
  DS_sel_PTMdata$treated = sapply(DS_sel_PTMdata$PMIDs, isTreated,PMIDs_treated$PMID)
  DS_sel_PTMdata$condition = mapply(condition, DS_sel_PTMdata$treated, DS_sel_PTMdata$normal)
  
  #Manully checked PTMs on phosphositePlust, who could not be assigned to a condition by using PMIDs. Otherwise I had 92 unkown PTMs
  PTMs_unknown_annotation <- read_excel("~/School/2e_master_bioinformatica/Thesis/tables/PTMs_unknown_annotation.xlsx")
  DS_sel_PTMdata = subset(DS_sel_PTMdata, condition != "unknown")
  DS_sel_PTMdata= rbind(DS_sel_PTMdata, PTMs_unknown_annotation)
  DS_sel_PTMdata = DS_sel_PTMdata[order(DS_sel_PTMdata$`Uniprot Accession`, DS_sel_PTMdata$Position),]
  
  #adding some metadata
  x = as.data.frame.matrix(t(table(DS_sel_PTMdata$condition)))
  colnames(x) = paste("PTMs_", colnames(x), sep="")                   
  DS_sel_metadata = cbind(DS_sel_metadata, x)
  
  #Seeing whitch UNP identifiers also PDB files have treated and normal PTMs
  UNP_normal = unique(subset(DS_sel_PTMdata, select='Uniprot Accession', subset = condition %in% c("normal", "both")))
  UNP_treated= unique(subset(DS_sel_PTMdata, select='Uniprot Accession', subset = condition %in% c("treated", "both")))
  
  PDB_normal = vector()
  PDB_treated = vector()
  PDB_onlyUNK = vector()
  x = which(colnames(DB_ID_UNP_sel)=="ID" )
  for(i in 1:length(DB_ID_UNP_sel$ID)){
  y = paste(DB_ID_UNP_sel[i,1:x-1])
  PDB_normal[i] =   length(intersect(y, UNP_normal$`Uniprot Accession`))> 0
  PDB_treated[i] =  length(intersect(y, UNP_treated$`Uniprot Accession`)) >0
  PDB_onlyUNK[i] =   PDB_normal[i]==0 &  PDB_treated[i]==0
  }
  
  DS_sel$hasNormalPTMs = PDB_normal
  DS_sel$hasTreatedPTMs = PDB_treated
  DS_sel$hasOnlyUNK = PDB_onlyUNK
  DS_sel$hasBoth = PDB_normal + PDB_treated == 2
  DS_sel_metadata$PDBs_treated = sum(PDB_treated)
  DS_sel_metadata$PDBs_normal = sum(PDB_normal)
  DS_sel_metadata$PDBs_both = sum(DS_sel$hasBoth)
  DS_sel_metadata$PDBs_onlyUNK = sum(DS_sel$hasOnlyUNK)
  
  check = subset(DS_sel_PTMdata, DS_sel_PTMdata$condition== "unknown")
  
  
#Code to search where ptms map within the structures.
#Meant to be interactively used with the pdb4amber file, uniprot_seq_DB and other extra information if needed. 
#
#chainEnd = the residue number of the last residue in the previous chain
#start = #of amino acids according to uniprot database before the amino acids starting in the pdb file
#gapb = beginning of a gap
#gape = end of a gap
x = subset(DS_sel_PTMdata, 
       select= Position, 
       subset= DS_sel_PTMdata$`Uniprot Accession` =="P42330" & condition %in% c("normal","both")  & Type == "Phosphorylation")
chainEnd = 314
start = 5
gap = 0 
x = x -start+chainEnd - gap
gapb = 228 -start+chainEnd - gap
gape = 232 -start+chainEnd - gap
gapb
gape
gape - gapb -1
as.character(x)



#code to get dataframe with PTM counting information
  Condition = c("total","normal", "treated")
  DS_sel_PTMCounts_UNP =  data.frame(matrix(ncol = 4, nrow = 0))
  colnames(DS_sel_PTMCounts_UNP) <- c('UNP', 'Condition', "AcUni", "PhUni")
  for(unp in unlist(strsplit(as.character(DS_sel_metadata$UNP_IDs), " "))){
  x = subset(DS_sel_PTMdata, subset= DS_sel_PTMdata$`Uniprot Accession` == unp)
  a =length(which(x$Type == "Acetylation"))
  b =length(which(x$Type == "Phosphorylation"))
  c =length(which(x$condition %in% c("normal", "both") & x$Type == "Acetylation"))
  d =length(which(x$condition %in% c("normal", "both") & x$Type == "Phosphorylation"))
  e =length(which(x$condition %in% c("normal", "both")& x$Type == "Acetylation"))
  f =length(which(x$condition %in% c("treated", "both") & x$Type == "Phosphorylation"))
  y =cbind(Condition, rbind(c(a,b),  c(c,d),   c(e,f)))
  y =cbind(rep(unp, 3),y)
  colnames(y) = colnames(DS_sel_PTMCounts_UNP)
  DS_sel_PTMCounts_UNP= rbind(DS_sel_PTMCounts_UNP, y)
  }
 
  Condition = c("total","normal", "treated")
  DS_sel_PTMCounts_PDB =  data.frame(matrix(ncol = 4, nrow = 0))
  colnames(DS_sel_PTMCounts_PDB) <- c('PDB', 'Condition', "AcUni", "PhUni")
  for(pdb in DS_sel$ID){
    x_unp = subset(DB_ID_UNP_sel, ID == pdb)[1,1:12]
    x = subset(DS_sel_PTMdata, subset= DS_sel_PTMdata$`Uniprot Accession` %in% x_unp)
    a =length(which(x$Type == "Acetylation"))
    b =length(which(x$Type == "Phosphorylation"))
    c =length(which(x$condition %in% c("normal", "both") & x$Type == "Acetylation"))
    d =length(which(x$condition %in% c("normal", "both") & x$Type == "Phosphorylation"))
    e =length(which(x$condition %in% c("normal", "both")& x$Type == "Acetylation"))
    f =length(which(x$condition %in% c("treated", "both") & x$Type == "Phosphorylation"))
    y =cbind(Condition, rbind(c(a,b),  c(c,d),   c(e,f)))
    y =cbind(rep(pdb, 3),y)
    colnames(y) = colnames(DS_sel_PTMCounts_PDB)
    DS_sel_PTMCounts_PDB= rbind(DS_sel_PTMCounts_PDB, y)
  }
  
  # write.xlsx(DS_sel_PTMCounts_UNP, "/Users/Jasper/Documents/School/2e_master_bioinformatica/Thesis/tables/PTM_counts.xlsx", sheetName = 'UNP_1', row.names=FALSE, append=TRUE)
  # write.xlsx(DS_sel_PTMCounts_PDB, "/Users/Jasper/Documents/School/2e_master_bioinformatica/Thesis/tables/PTM_counts.xlsx", sheetName = "PDB_1", row.names=FALSE, append=TRUE)
  # write.xlsx(DB_ID_UNP_sel, "/Users/Jasper/Documents/School/2e_master_bioinformatica/Thesis/tables/PTM_counts.xlsx", sheetName = "PDB_build", row.names=FALSE, append=TRUE)
  # write.xlsx(DS_sel_PTMdata, "/Users/Jasper/Documents/School/2e_master_bioinformatica/Thesis/tables/UNP_PTM_data.xlsx", row.names=FALSE, append=TRUE)
  

  c=c("19,46,50,76,77,105,113,124,157,161,191,314,333,360,364,390,391,419,427,438,471,475,505,628")
  c = unlist(strsplit(c, ","))
  as.character(length(c))
c

