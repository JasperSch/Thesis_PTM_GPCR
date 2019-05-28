#script to explore the GO annotation of pdb structures
#Used on excell table downloaded from uniprot after querying for good pdb structures
#Extra information can be added with script PDB_subset_overview
library(readxl)
df <- read_excel("~/School/2e_master_bioinformatica/Thesis/Uniprot/GoodStructGO_annotation.xlsx")
names(df)
names(df)[13] <- "PDBs"


plot(as.factor(df$Organism))
#apperently there are a couple of non-human proteins

#Get all GO terms
GO_str  = gsub(";","",paste(df$`Gene ontology IDs`, collapse="; "))
GO_vec = strsplit(GO_str, " ")[[1]]
GO_ids = unique(GO_vec)

GO_counts =as.data.frame(table(GO_vec))

#take subset of GO terms that occur more then 30 times to reduce computationial time (Because this will be at least 30 PDB structures)
GO_freqGood = as.vector(GO_counts$GO_vec[GO_counts$Freq>30])


checkGO <- function(x, GO_freqGood){
  for (term in GO_freqGood){
    if (grepl(term, x)){
      return(TRUE)
    }
  }
  return(FALSE)
}


df$Sel = sapply(X = df$`Gene ontology IDs`, checkGO, GO_freqGood)

#make extra collumns for each good GO term saying wether or not a protein has this term
  for (term in GO_freqGood){
    df[term] = sapply(X = df$`Gene ontology IDs`, checkGO, as.vector(term))
}

#Make overview of amount of proteins and PDBs for each GO term and make df for each GO term
GO_counts = data.frame(GO = character(), Protein_count = integer(), PDB_count = integer(), ratio= double(), PDBs=character(), UNPs = character(), stringsAsFactors=FALSE)
i=1

for (term in GO_freqGood){
  x = subset(df, select = 1:14, subset= df[term] == TRUE)
  Protein_count = length(x[[1]])
  PDBs = unique(strsplit(paste(x$PDBs, collapse=","), ",")[[1]])
  UNPs = unique(strsplit(paste(x$Entry, collapse=","), ",")[[1]])
  PDB_count  =length(PDBs)
  if (PDB_count>40 & PDB_count <100){
    assign(gsub(":", "_",term), x)
    GO_counts[i,]=c(term, Protein_count, PDB_count, Protein_count/PDB_count, paste(PDBs, collapse = " "), paste(UNPs, collapse = " "))
    i=i+1
  }
}
names(GO_counts) = c("GO", 'Protein_Count', 'PDB_Count', 'ratioProtein_PDB', 'PDBs', "Proteins")
GO_counts$Protein_Count = as.integer(GO_counts$Protein_Count)
GO_counts$PDB_Count = as.integer(GO_counts$PDB_Count)
GO_counts$ratioProtein_PDB = as.double(GO_counts$ratioProtein_PDB)

GO_counts[,1]

