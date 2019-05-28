#OVerlaplist revisited
#13/02/2018
library(readr)
library(rpart)
library(rpart.plot)
library(hash)
library(stringr)



#Enter FILE to translate it to dataset
DS <- read_csv("~/School/2e_master_bioinformatica/Thesis/tables/Overview_GPROT_27", 
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
DS$METHOD = as.factor(DS$METHOD)




#make dataframe of database identifiers + structure ID's
DB_Identifiers = str_split_fixed(DS$DB_IDENTIFIERS, " ", length(DS$ID) )
DB_Identifiers = as.data.frame(DB_Identifiers)

#Check from which collumn empty collumns begin and remove all empty ones
for (i in length(DB_Identifiers[,1]):1){
  colNumber = i
  if (paste(DB_Identifiers[,i], collapse = '') != "")
    break
}                             

DB_Identifiers = DB_Identifiers[1:colNumber]
DB_Identifiers = cbind(DB_Identifiers, ID = DS$ID)

#Specialize dataframe to only include uniprot identifiers and count the amount of chains
DB_Identifiers_UNP <- as.data.frame(lapply(DB_Identifiers, function(x) gsub("[0-9a-zA-Z_]+\\((?!UNP)[0-9a-zA-Z_]+\\)", "", x, perl=TRUE))) #removes al ID's except UNP ID's
DB_Identifiers_UNP <- as.data.frame(lapply(DB_Identifiers_UNP, function(x) gsub("\\(UNP\\)", "", x, perl=TRUE)))
DB_Identifiers_UNP <- data.frame(lapply(DB_Identifiers_UNP, as.character), stringsAsFactors=FALSE)
DB_Identifiers_UNP$Entry_count <- apply(DB_Identifiers_UNP, 1, function(x) sum(x!= "")-1) #Minus 1 to not include the ID-cell. *****************BEWARE! DO NOT RUN THIS TWICE***********

#Making the list with overlaps
Overlaplist = data.frame(matrix(ncol = 3),stringsAsFactors = FALSE)
Overlaplist= Overlaplist[!is.na(Overlaplist)]


#Searching every unique UNP identifier
chains = DB_Identifiers_UNP[,1]
for (i in 2:colNumber){
  chains = c(as.character(chains), as.character(DB_Identifiers_UNP[,i]))
}
print(length(chains))
chains = unique(chains[chains !=""])
print(length(chains))

#For every unique UNP > search lines with UNP > extract structures > Compare structures and put in overlaplist

total <- length(chains)

for (chain in chains){
  #Find all unqiue structures in DB_Identifiers_UNP that have the chain
  matches = DB_Identifiers_UNP[grep(chain, DB_Identifiers_UNP[,1]),]
  for (i in 1:colNumber){
    matches = rbind(matches, DB_Identifiers_UNP[grep(chain, DB_Identifiers_UNP[,i+1]),])
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
check = unique(Overlaplist)
x = paste(Overlaplist[,1], Overlaplist[,2], sep= " ")
Overlaplist[,1] = as.vector(unlist(lapply(x, function(x){ str_sort(str_split(x, " ")[[1]])[1]})))
Overlaplist[,2] = as.vector(unlist(lapply(x, function(x){ str_sort(str_split(x, " ")[[1]])[2]})))
Overlaplist = unique(Overlaplist)

plot(Overlaplist[,3])
#At the time when I didn't include the coverage,  one outlier spotted: 1YDI has overlap of zero. This is because of '-2' standing behind the UNP identifier in the PDB file...

