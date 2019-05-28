#PTM conservation
#Jasper Schelfhout
#Run together with Conservation_analysis_shiny to get interactive plots
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(readr)
library(readxl)
library(stringr)
library("RColorBrewer")
library(evobiR)

setwd("~/School/2e_master_bioinformatica/Thesis/tables")


#Loading data
Table_Complex_Build <- read_excel("~/School/2e_master_bioinformatica/Thesis/tables/PTM_counts.xlsx", 
                                  sheet = "PDB_build")
Table_UNP_EMBL <- read_excel("~/School/2e_master_bioinformatica/Thesis/tables/PTM_counts.xlsx", 
                                  sheet = "UNP_EMBL")
Table_chain_range <- read.csv("Table_chain_range", sep="")
UNP_PTM_data <- read_excel("~/School/2e_master_bioinformatica/Thesis/tables/UNP_PTM_data.xlsx")
colnames(UNP_PTM_data)[1:2]=c("UNP", "RESN")
UNP_length <- read_csv("~/School/2e_master_bioinformatica/Thesis/Uniprot/uniprot_proteome_ID_Length")
Table_DSSP<-read_delim("~/School/2e_master_bioinformatica/Thesis/tables/DSSP_ENS/Table_DSSP.csv", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE)

#path of folder with alignments
Al_path = "~/School/2e_master_bioinformatica/Thesis/tables/Conservation_scores/"

#Amino acids
AA = c("M", "C", "R", "T", "L", "A", "F", "P", "E", "K", "G", "I", "H", "S", "D", "W", "N", "V", "Q", "Y")


#Function that returns a matrix with conservation scores for PTM residues
#data is the matrix to do the calculations upon. Each row is a residue position in the alignment, each column is a different protein
#Seq is the column of the data matrix to calculate the conservation scores for
ConsRates <- function(Seq, data){
  ConsLys = vector(mode="double", length=nrow(data))
  ConsSer = vector(mode="double", length=nrow(data))
  ConsTyr = vector(mode="double", length=nrow(data))
  ConsThr = vector(mode="double", length=nrow(data))
  ConsGen = vector(mode="double", length=nrow(data))
  for(res in which(Seq == 'K')){
    a=data[res,] =='K'
    b=data[res-1,]=='K'
    c=data[res+1,]=='K'
    d = a+b+c
    ConsLys[res]= sum(d!=0)/(ncol(data))
  }
  for(res in which(Seq == 'Y')){
    a=sum(data[res,] =='Y')
    ConsTyr[res]= a/(ncol(data))
  }
  for(res in which(Seq == 'S')){
    a=sum(data[res,] =='S' )
    b=sum(data[res,] =='T' )
    ConsSer[res]= (a+b)/(ncol(data))
  }
  for(res in which(Seq == 'T')){
    a=sum(data[res,] =='S')
    b=sum(data[res,]=='T')
    ConsThr[res]= (a+b)/(ncol(data))
  }
  for(res in 1:length(Seq)){
    a=sum(data[res, ]== Seq[res])
    ConsGen[res]= a/ncol(data)
  }
  ConsRates = data.frame(ConsLys, ConsTyr, ConsThr, ConsSer, ConsGen)
  return(ConsRates)
}


#function to calculate Shannon entropy
#Shannon entropy = - sum over all AA p*log(p)
ShaEnt = function(Seq, data){
  sh = vector()
  for (i in 1:length(Seq)){
    s=0
    for (res in unique(data[i,])){
      p = sum(data[i,]== res) / ncol(data)
      s = s + p*log(p)
    }
    sh[i] = -s
  }
  return(sh)
}

#Making conservation dataframe Table_conservation  
  cons= vector()
  for(i in 1:nrow(Table_UNP_EMBL)){
    unp = Table_UNP_EMBL$UNP[i]
    ens = Table_UNP_EMBL$ENSEMBL[i]
    path = paste(Al_path, unp, "_al_t", sep="")
    x <- read_delim(path, " ", escape_double = FALSE, trim_ws = TRUE)
    d= t(str_split_fixed(x$Seq,"",nchar(x$Seq[1])))
    colnames(d)= x$Id
    coln = which(colnames(d)== paste("9606.",ens, sep=""))
    RES = d[,coln]
    sh = ShaEnt(RES,d)
    d = as.data.frame(cbind(RES,ConsRates(RES,d)))
    d$ConsSha = sh
    d= subset(d, RES != "-")
    d$RESN = 1:nrow(d)
    d$UNP = rep(unp,nrow(d))
    cons = rbind(cons,d)
  }
  Table_conservation = merge(cons, UNP_PTM_data, all.x=TRUE)
  attach(Table_conservation)
  Table_conservation = Table_conservation[order(UNP, RESN),]
  rownames(Table_conservation) = 1:nrow(Table_conservation)
  detach(Table_conservation)
  
 #add a collumn which combines the ConsGen but adjusts it for residues K,Y,S,T
  mapAdjCons = function(Lys, Ser, Thr, Tyr, Gen, Res){
    if (Res == "K"){
      return(Lys)
    }
    if (Res == "S"){
      return(Ser)
    }
    if (Res == "T"){
      return(Thr)
    }
    if (Res == "Y"){
      return(Tyr)
    }
    else{
      return(Gen)
    }
  }
  
  d= Table_conservation
  Table_conservation$ConsAdj = mapply(FUN = mapAdjCons, Lys = d$ConsLys,
                                                        Ser = d$ConsSer, 
                                                        Thr = d$ConsThr, 
                                                        Tyr = d$ConsTyr, 
                                                        Gen = d$ConsGen, 
                                                        Res = as.character(d$RES))
  
  #Add a collumn with PTM type
  mapPTMtype = function(Res){
    if (Res %in% c("S","T","Y")){
      return("Phosphorylation")
    }
    if (Res == "K"){
      return("Acetylation")
    }
    else{
      return(NA)
    }
  }
  Table_conservation$PTMtype = unlist(lapply(Table_conservation$RES, function(x) as.character(mapPTMtype(x))))
 
  #Add a collumn with background conservation (mean of windowsize)
  #ConsBck is an adjusted shannon entropy score, scaling it to [0,1] for each protein
  attach(Table_conservation)
  Table_conservation$ConsBck = 1-ConsSha/max(ConsSha) #scales shannon entropy
  ConsBck_loess = vector()
  for (unp in unique(UNP)){
    x= subset(Table_conservation, UNP==unp)
    lo = loess(ConsBck ~ RESN, data = x, span=20/nrow(x)) #always use 20 neighbourings residues
    p = predict(lo)
    ConsBck_loess = c(ConsBck_loess, p)
  }
  Table_conservation$ConsBck_loess = ConsBck_loess
  
  #correct for mean conservation score of each protein
  BckMean = vector()
  BckLowBound=vector()
  i=1
  for (unp in Table_conservation$UNP){
    x=subset(Table_conservation, UNP== unp)
    BckMean[i] = mean(x$ConsBck)
    BckLowBound[i] = quantile(x$ConsBck_loess, 0.4)
    i=i+1
  }
  Table_conservation$BckMean = BckMean
  Table_conservation$ConsRel = ConsAdj / BckMean
  detach(Table_conservation)
  attach(Table_conservation)
  Table_conservation$ConsRel = (ConsRel - min(ConsRel)) / max(ConsRel - min(ConsRel))
  detach(Table_conservation)
  


  #Get structural information out of DSSP dataframe to Table_conservation
  mapUNP <- function (ens){
    return(Table_UNP_EMBL$UNP[Table_UNP_EMBL$ENSEMBL == ens])
  }
  mapStructure <- function(s){
    if(s %in% c("H", "G", "I")){
      return("helix")
    }
    if(s %in% c("B", "E")){
      return("sheet")
    }
    else{
      return("loop")
    }
  }
  Table_DSSP$UNP = unlist(lapply(Table_DSSP$ENS, function(x) mapUNP(x)))
  Table_conservation= merge(Table_conservation, Table_DSSP[c(1,3,4,5,6)], sort=FALSE)
  Table_conservation$Structure = unlist(lapply(Table_conservation$S, function(x) mapStructure(x)))
  Table_conservation$Conserved_DSSP = Table_conservation$Structure != "loop"
  Table_conservation$Conserved_REG = Table_conservation$ConsBck_loess > BckLowBound

  Table_conservation$inTemplate = rep(FALSE, nrow(Table_conservation))
  attach(Table_UNP_EMBL)
  for (unp in UNP){
    x =unlist(strsplit(Table_UNP_EMBL$Range[UNP == unp], ","))
    for(r in x){
      b = as.integer(unlist(strsplit(r,"-")))[1] + min(which(Table_conservation$UNP==unp)) -1
      e = as.integer(unlist(strsplit(r,"-")))[2] + min(which(Table_conservation$UNP==unp)) -1
      Table_conservation$inTemplate[b:e]=TRUE
    }
  }
  detach(Table_UNP_EMBL)
  
  # Bad working LOESS
  # c= vector()
  # i=1
  # for (unp in unique(Table_conservation$UNP)){
  #   x = subset(Table_conservation, inTemplate == TRUE & UNP==unp)
  #   c[i] =  table(x$Conserved_DSSP)[1] / nrow(x)
  #   print(unp)
  #   print(c[i])
  #   i=i+1
  # 
  # }
  # ggplot() + aes(c)+ geom_histogram() + GGtheme
  # 
  # x = subset(Table_conservation, inTemplate == TRUE)
  # table(x$Conserved_DSSP)/ nrow(x)
  # table(x$Conserved_DSSP)
  # table(x$Conserved_DSSP, x$Conserved_REG)
  

  
#################
#Plots          #
#################
  colPal = brewer.pal(n= length(unique(Table_conservation$Type)), name = "Set2")
  GGtheme = theme_bw() + 
    theme(axis.title = element_text(size=20), 
          legend.title = element_text(size=20), 
          legend.text = element_text(size=15))
  
  #plot of background scores
  smoothPlot_Bck_Cons = function(df){
    df_PTM = subset(df, RES %in% c("Y", "K", "T", "S"))
    ggplot(df, aes(x=RESN, y= ConsBck)) + 
      geom_point(size=2) + 
      #geom_point(data = df_PTM, aes(x=RESN, y= ConsBck, color=PTMtype, shape=is.na(Type)), alpha =0.5, size=5) +
      #stat_smooth(aes(group=1),se = FALSE,span = 20/nrow(df), color="red") +
      geom_line(data=df, aes(x=RESN, y=ConsBck_loess), color="red", size=1.2) +
      geom_area(data=df, aes(x=RESN, y=as.integer(Conserved_DSSP)), color="black", alpha=0.1)+
      scale_color_manual(values = c("Phosphorylation"="#FC8D62","Acetylation"="#66C2A5"))  +
      GGtheme
  }
  
  #plot of background scores for all combinded
  smoothPlot_Bck_Cons_all = function(df){
    ggplot(df, aes(x=1:nrow(df), y= ConsBck)) + 
      geom_point(size=1) + 
      geom_line(data = df, aes(x=1:nrow(df), y= ConsBck_loess), color="red") +
      GGtheme
  }
  
  #general boxplots of conservation scores factored by residue
  # df must be a subset of Table_conservation
  df=Table_conservation
  boxPlot_Res_Cons = function(df){
    ggplot(df, aes(x= RES, y=ConsGen, fill=RES)) +
      scale_fill_manual(values = c(c("Y"="#FC8D62", "S"="#FC8D62", "T"="#FC8D62", "K"="#66C2A5"), rep("grey",16)), guide=FALSE)  +
      geom_boxplot() +
      theme_bw() +
      labs(x="residue", y= "conservation")
  }
  boxPlot_Res_Cons(df) 

  #split boxplots of conservation scores factored by PTM residue and modification
  # df must be a subset of Table_conservation
  df=Table_conservation
  boxPlot_PTM_ConsAdj = function(df){
    df = subset(df, RES %in% c("Y","K","S","T"))
    ggplot(df, aes(x= Type, y=ConsAdj, fill =PTMtype)) +
      geom_boxplot() +
      facet_grid(.~RES, labeller = label_both, drop=TRUE) +
      theme_bw() +
      theme(axis.title= element_text(size=20), legend.title = element_text(size=20),
            axis.text=element_text(size=15), legend.text=element_text(size=15),strip.text = element_text(size=15)) +
      scale_x_discrete(limits = unique(Table_conservation$Type), labels=c("No", "Ph", "Ac"))+
      scale_fill_manual(values = c("#66C2A5", "#FC8D62")) + 
      labs(x="", y= "conservation adjusted")
  }
  boxPlot_PTM_ConsAdj(df)

  
  #split boxplots of conservation scores factored by PTM residue and modification
  # df must be a subset of Table_conservation
  # This uses the general conservation score instead of the modified one for specific residues
  df=Table_conservation
  boxPlot_PTM_ConsGen = function(df){
    df = subset(df, RES %in% c("Y","K","S","T"))
    ggplot(df, aes(x= Type, y=ConsGen, fill =PTMtype)) +
      geom_boxplot() +
      facet_grid(.~RES, labeller = label_both, drop=TRUE) +
      theme_bw() +
      theme(axis.title= element_text(size=20), legend.title = element_text(size=20),
            axis.text=element_text(size=15), legend.text=element_text(size=15), strip.text = element_text(size=15)) +
      scale_x_discrete(limits = unique(Table_conservation$Type), labels=c("No", "Ph", "Ac"))+
      scale_fill_manual(values = c("#66C2A5", "#FC8D62")) + 
      labs(x="", y= "conservation general")
  }
  boxPlot_PTM_ConsGen(df)
  
  #split boxplots of conservation scores factored by PTM residue and modification
  # df must be a subset of Table_conservation
  # This uses the relative conservation score instead of the modified one for specific residues
  df=Table_conservation
  boxPlot_PTM_ConsRel = function(df){
    df = subset(df, RES %in% c("Y","K","S","T"))
    ggplot(df, aes(x= Type, y=ConsRel, fill =PTMtype)) +
      geom_boxplot() +
      facet_grid(.~RES, labeller = label_both, drop=TRUE) +
      theme_bw() +
      theme(axis.title= element_text(size=20), legend.title = element_text(size=20),
            axis.text=element_text(size=15), legend.text=element_text(size=15), strip.text = element_text(size=15)) +
      scale_x_discrete(limits = unique(Table_conservation$Type), labels=c("No", "Ph", "Ac"))+
      scale_fill_manual(values = c("#66C2A5", "#FC8D62")) + 
      labs(x="", y= "conservation relative (Adj)")
  }
  boxPlot_PTM_ConsRel(df)
  
  df=Table_conservation
  Table_conservation$ConsRel_Bck = df$ConsBck - df$ConsBck_loess - BckMean
  df=Table_conservation
  boxPlot_PTM_ConsRel_Bck = function(df){
    df = subset(df, RES %in% c("Y","K","S","T"))
    ggplot(df, aes(x= Type, y=df$ConsRel_Bck, fill =PTMtype)) +
      geom_boxplot() +
      facet_grid(.~RES, labeller = label_both, drop=TRUE) +
      theme_bw() +
      theme(axis.title= element_text(size=20), legend.title = element_text(size=20),
            axis.text=element_text(size=15), legend.text=element_text(size=15), strip.text = element_text(size=15)) +
      scale_x_discrete(limits = unique(Table_conservation$Type), labels=c("No", "Ph", "Ac"))+
      scale_fill_manual(values = c("#66C2A5", "#FC8D62")) + 
      labs(x="", y= "conservation relative (Sha)")
  }
  boxPlot_PTM_ConsRel_Bck(df)

  #facetted boxplot ConsAdj based on DSSP.
  df = subset(Table_conservation, inTemplate==TRUE)
  df$Structure_grouped = df$Structure == "loop"
  Boxplot_PTM_ConsDSSP <- function(df){
    df = subset(df, RES %in% c("K","Y","T","S")) 
    ggplot(df, aes(x=Structure_grouped,y=ConsAdj, fill=Type)) + 
      geom_boxplot() +
      GGtheme + 
      stat_compare_means(method = "t.test", label="p.format") +
      expand_limits(x = unique(df$Structure_grouped, fill=unique(Table_conservation$Type)))+
      facet_grid(rows= vars(RES), scales="free")  +
      labs(x="", y="CA", fill="") +
      scale_fill_manual(labels= c("Acetylated", "Phosphorylated", "Unmodified"), values =c("#0072B2", "#D55E00")) +
      scale_x_discrete(labels=c("Structured", "Unstructured")) +
      theme(axis.text.x   = element_text(size = 15),
            axis.text.y   = element_text(size =10),
            strip.text=element_text(size=15),
            strip.background=element_rect(fill="white"),
            legend.position = "top")
  } 
  Boxplot_PTM_ConsDSSP(df)
  
  df = subset(Table_conservation, inTemplate==TRUE)
  Boxplot_PTM_ACC <- function(df){
    df = subset(df, RES %in% c("K","Y","T","S")) 
    ggplot(df, aes(x=1,y=ACC, fill=Type)) + 
      geom_boxplot() +
      GGtheme + 
      stat_compare_means(method = "t.test", label='p.format') +
      facet_grid(cols= vars(RES), scales="free") +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.position="bottom",
            legend.title = element_blank()) +
      scale_fill_manual(labels= c("Acetylated", "Phosphorylated", "Unmodified"), values =c("#0072B2", "#D55E00")) +
      theme(strip.text=element_text(size=15),
          strip.background=element_rect(fill="white"))
    }


#play ground analysis
  
  table(df$Structure, df$Type)
  table(df$Type)
  table(UNP_PTM_data$Type, UNP_PTM_data$condition) 
  table(Table_conservation$Type, Table_conservation$UNP)
  table(df$Structure) / nrow (df)
  
  
 x=subset(Table_conservation, UNP=="P57771" & RESN >63 & RESN <196)
 smoothPlot_Bck_Cons(x) + 
   geom_line(y= 0.35, size=1) + 
   theme(legend.text = element_blank(), 
         legend.title = element_blank(),
         axis.text = element_text(size=15)) +
   labs(y= "CS", x="Residue")
 print(as.integer(x$ConsBck_loess>0.35))

 for (i in 1:nrow(Table_UNP_EMBL)){
   y = subset(Table_conservation, select=c("RESN", "RES", "ConsBck", "ConsBck_loess", "Conserved_DSSP", "Conserved_REG"), 
              UNP== Table_UNP_EMBL$UNP[i]) #& RESN >63 & RESN <195)
  # write.table(as.integer(y$Conserved_DSSP),  file=paste("C:/Users/Jasper/Documents/School/2e_master_bioinformatica/Thesis/tables/Pymol_cons_colors/",Table_UNP_EMBL$ENSEMBL[i],".csv", sep=""), col.names=FALSE, row.names = FALSE, quote=FALSE)
 }
 

 
 

 
 
 
 
 # for(r in c("K","Y","T","S")) {
 #   df = subset(Table_conservation, RES ==r & inTemplate == TRUE)
 #   assign(paste("p",r, sep=""),
 #          ggplot(df, aes(x=Structure,y=ConsAdj, fill=Type)) + geom_boxplot() +GGtheme + labs(title=paste(r, "in template"))
 #   )
 # }
 # pK
 # pS
 # pT
 # pY