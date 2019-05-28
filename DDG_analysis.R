#Jasper Schelfhout
#Script to analyze Delta G values
#Run together with DDG_analyis_shiny to interactively play with plots

library(ggplot2)
library(gridExtra)
library(readr)
library(readxl)
library(stringr)
library("RColorBrewer")


##############
#Loading data#
##############
setwd("~/School/2e_master_bioinformatica/Thesis/tables")

Table_DG_complex <- read.csv("Table_DG_complex", sep="")
Table_DG_residues <- read.csv("Table_DG_residues", sep="")
Table_chain_range <- read.csv("Table_chain_range", sep="")
Table_PTM_position <- read_delim("~/School/2e_master_bioinformatica/Thesis/tables/PTM_pos_normal.csv", " ", escape_double = FALSE, trim_ws = TRUE)
Table_PTM_position <- Table_PTM_position[1:22,]
Table_PTM_position$PDB <- tolower(Table_PTM_position$PDB)
Table_Complex_Build <- read_excel("~/School/2e_master_bioinformatica/Thesis/tables/PTM_counts.xlsx", 
                                  sheet = "PDB_build")
colnames(Table_Complex_Build)= c("A", "B", "C", "D", "E", "F", "PDB", "Count_raw", "Count_clean")
Table_Complex_Build$PDB <- tolower(Table_Complex_Build$PDB)

#path of folder with alignments
Al_path = "~/School/2e_master_bioinformatica/Thesis/tables/Conservation_scores/"



###############################################
# Combining dataframes to master table        #
###############################################

#Making dataframes with delta delta G
  x = subset(Table_DG_complex,select = DG, COND== "no" & !(PDB %in% c("1nr4", "2gtp")))$DG #exclude those whithout NC PTMs
  y = Table_DG_complex$DG[Table_DG_complex$COND == "NC"]
  Table_DDG_complex = cbind(subset(Table_DG_complex, select=c(1,3,5), COND=="NC"), DDG =y-x)
  
  x = subset(Table_DG_residues,select = DG, COND== "no" & !(PDB %in% c("1nr4", "2gtp")))$DG #exclude those whithout NC PTMs
  y = Table_DG_residues$DG[Table_DG_residues$COND == "NC"]
  Table_DDG_residues= cbind(subset(Table_DG_residues, select=c(1,3,4,5,7), COND=="NC"), DDG =y-x)
  rownames(Table_DDG_residues) = c(1:length(Table_DDG_residues$PDB))


#adding collumn with Chain specification to residue delta delta G mastertable
  mapChain <- function(pdb, res, resn){
    x = subset(Table_chain_range, PDB == pdb)
    for(i in (1:length(x[,1]))){
      b = x$BEG[i]
      e = x$END[i]
      if(resn %in% c(b:e)){
      return(x$CHAIN[i])
      }
    }
    return(res)
  }
  
  Table_DDG_residues$CHAIN = mapply(FUN= mapChain, Table_DDG_residues$PDB,
                                                   Table_DDG_residues$RES,
                                                   Table_DDG_residues$RESN)

#Adding collumn with PTM information to DDG_residue master table
  # Table_DDG_residues =cbind(Table_DDG_residues, PTM = rep(NA,length(Table_DDG_residues$PDB)))
  # Table_DDG_residues$PTM = as.character(Table_DDG_residues$PTM)
  # for(s in unique(Table_DDG_residues$PDB)){
  #   ph = unlist(strsplit(Table_PTM_position$PH[Table_PTM_position$PDB== s], ","))
  #   ac = unlist(strsplit(Table_PTM_position$AC[Table_PTM_position$PDB==s], ","))
  #   for (pos in ph){
  #     x = subset(Table_DDG_residues, PDB== s & RESN == pos)
  #     for (r in as.integer(rownames(x))){
  #       Table_DDG_residues$PTM[r] = "Ph"
  #     }
  #   }
  #   for (pos in ac){
  #     x = subset(Table_DDG_residues, PDB== s & RESN == pos)
  #     for (r in as.integer(rownames(x))){
  #       Table_DDG_residues$PTM[r] = "Ac"
  #     }
  #   }
  # }
  # Table_DDG_residues$PTM = as.factor(Table_DDG_residues$PTM)
  # 
  # x= subset(Table_DDG_residues, PTM=="Ac" & RES=="LYS")
 
  mapPTMsites <- function (res){
    if (res == "ALY"){
      return("Ac")
    }
    if (res %in% c("PTR", "TPO", "SEP")){
      return("Ph")
    }
    else{
      return('unmodified')
    }
  }
  Table_DDG_residues$PTM = unlist(lapply(Table_DDG_residues$RES, function(x) mapPTMsites(x)))
  
#exclude gap residues
  Table_DDG_residues = subset(Table_DDG_residues, RES != "ACE" & RES != "NME")
  
  #fix numeration
  for (pdb in unique(Table_DDG_residues$PDB)){
    for(comp in unique(Table_DDG_residues$COMP[Table_DDG_residues$PDB == pdb])){
      x = subset(Table_DDG_residues, COMP == comp & PDB == pdb)
      Table_DDG_residues$RESN[Table_DDG_residues$PDB ==pdb & Table_DDG_residues$COMP == comp]= c(1: nrow(x))
    }
  }
  
#Conservation scores mapping to dataframes
  # cons= vector()
  # for(pdb in unique(Table_DDG_residues$PDB)){ 
  #   z=vector()
  #   print(pdb)
  #   for(chain in unique(Table_chain_range$CHAIN[Table_chain_range$PDB == pdb])){
  #     path = paste("~/School/2e_master_bioinformatica/Thesis/tables/Conservation_scores/", pdb, "_", chain, "_cons.txt", sep="")
  #     x <- read_delim(path, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 4 )
  #     y = apply(x, 1, function(x) str_sub(x,-1,-1))
  #     x$X3 = y[3,]
  #     x = subset(x, X3 != "-", c(2,3))
  #     colnames(x) = c("Shannon_Adj", "RES_short")
  #     x$CHAIN = rep(chain, length(x[,1]))
  #     x$PDB = rep(pdb, length(x[,1]))
  #     z = rbind(z, x)
  #     print(paste("Check", chain))
  #   }
  #   z$RESN = c(1:length(z$PDB))
  #   cons = rbind(cons,z)
  # }
  # Table_DDG_residues$rown <- 1:nrow(Table_DDG_residues) 
  # x = merge(Table_DDG_residues,cons, all.x=TRUE)
  # Table_DDG_residues <- x[order(x$rown), ]
 

  #Function that returns a matrix with conservation scores for PTM residues
  #data is the matrix to do the calculations upon. Each row is a residue position in the alignment, each column is a different protein
  #Seq is the column of the data matrix to calculate the conservation scores for
  ConsRates <- function(Seq, data){
    ConsLys = vector(mode="double", length=nrow(data))
    ConsSer = vector(mode="double", length=nrow(data))
    ConsTyr = vector(mode="double", length=nrow(data))
    ConsThr = vector(mode="double", length=nrow(data))
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
    ConsRates = cbind(ConsLys, ConsTyr, ConsThr, ConsSer)
    colnames(ConsRates)= c("ConsLys", "ConsTyr", "ConsThr", "ConsSer")
    return(ConsRates)
  }
  
  #Determining conservation scores
  cons= vector()
  for(pdb in unique(Table_chain_range$PDB)){ 
    z=vector()
    print(pdb)
    i=1
    for(chain in unique(Table_chain_range$CHAIN[Table_chain_range$PDB == pdb])){
      path = paste(Al_path, pdb, "_", chain, "_cons.txt", sep="")
      x <- read_delim(path, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 4 )
      d= str_split_fixed(x$X3,"",nchar(x$X3[1]))
      y = apply(x, 1, function(x) str_sub(x,-1,-1))
      x$X3 = y[3,]
      x= cbind(x, ConsRates(d[,ncol(d)],d))
      x = subset(x, X3 != "-", c(2,3:ncol(x)))
      colnames(x)[1:2] = c("Shannon_Adj", "RES_short")
      x$CHAIN = rep(chain, length(x[,1]))
      x$PDB = rep(pdb, length(x[,1]))
      z = rbind(z, x)
    }
    z$RESN = c(1:nrow(z))
    cons = rbind(cons,z)
  }
  Table_DDG_residues$rown <- 1:nrow(Table_DDG_residues) 
  x = merge(Table_DDG_residues,cons, all.x=TRUE)
  Table_DDG_residues <- x[order(x$rown), ]
  attach(Table_DDG_residues)
  Table_DDG_residues$consPTM = ConsLys + ConsTyr + ConsThr + ConsSer
  detach(Table_DDG_residues)

  #mapping UNP ids of Table_Complex_Build to cons dataframe
  mapUNPcons <- function(pdb, chain){
  x = subset(Table_Complex_Build, PDB == pdb)
  col = which(colnames(x)==chain)
  return(x[col])
  }
  cons$UNP = mapply(FUN= mapUNPcons, cons$PDB, cons$CHAIN)
  
#map amino acid properties
  mapCharge <- function(Res){
    if (Res %in% c("K", "R", "H")){
      return("positive")
    }
    if (Res %in% c("D", "E")){
      return("negative")
    }
    else{
      return ("no")
    }
  }
  Table_DDG_residues$Charge = unlist(lapply(Table_DDG_residues$RES_short, function(x) as.character(mapCharge(x))))
  
#Get dataframes with ranges of structures and decomposed structures (for plotting purposes)
  i=1
  e = vector()
  b = vector()
  for (s in unique(Table_DDG_residues$PDB)){
    x = subset(Table_DDG_residues, subset = PDB == s & GDP==FALSE)
    b[i] = min(as.numeric(rownames(x)))
    e[i] = max(as.numeric(rownames(x)))
    i=i+1
  }
  ranges_struct <- data.frame(begin = b, end = e, col = unique(Table_DDG_residues$PDB))
  
  i=1
  e = vector()
  b = vector()
  for (s in unique(paste(Table_DDG_residues$PDB, Table_DDG_residues$COMP))){
      x = subset(Table_DDG_residues, PDB == unlist(strsplit(s, " "))[1] & COMP == unlist(strsplit(s, " "))[2])
      b[i] = min(as.numeric(rownames(x)))
      e[i] = max(as.numeric(rownames(x)))
      i=i+1
    }
  ranges_comp<- data.frame(begin = b, end = e, col = unique(paste(Table_DDG_residues$PDB, Table_DDG_residues$COMP)))
  ranges_comp$PDB = do.call(rbind, strsplit(as.character(ranges_comp$col), " "))[,1]
  ranges_comp$CHAIN = do.call(rbind, strsplit(as.character(ranges_comp$col), " "))[,2]
  

  
#######################
# PLOTS AND ANALYSIS  #     #Written as functions for easy use in shiny
#######################

colPal  = brewer.pal(n= length(unique(Table_DDG_residues$CHAIN)), name = "Set2")
GGtheme = theme_bw() + 
          theme(axis.title = element_text(size=20), 
          legend.title = element_text(size=20), 
          legend.text = element_text(size=15))

#One big plot of all residue delta G's in all structures, but only for PTMs 
  ResidueDDGplot_all= function(s, PTM){
    if (PTM== "all"){
      d = subset(Table_DDG_residues, GDP == FALSE)
      sd = subset(d, PDB == s)
      dotsize = 0.5
    }
    if (PTM== "PTM"){
      d= subset(Table_DDG_residues, GDP==FALSE & !(is.na(PTM)))
      sd= subset(d, !(is.na(PTM)) & PDB == s)
      dotsize = 1.5
    }
    p = ggplot() +
      geom_point(data = d, aes(x=as.numeric(rownames(d)), y=DDG), size = dotsize, shape=20) +
      geom_point(data = sd, aes(x=as.numeric(rownames(sd)), y=DDG), size = 2, color="red", shape =20) +
      GGtheme +
      theme(axis.text.x= element_blank())+
      geom_rect(data=ranges_struct, aes(xmin= begin, xmax=end, ymin = -Inf, ymax= +Inf, fill = col), alpha=0.1) +
      geom_rect(data=subset(ranges_struct, col == s), aes(xmin= begin, xmax=end, ymin = -Inf, ymax= +Inf, fill = col), alpha=0.4) +
      geom_vline(xintercept=c(1,ranges_comp$end), linetype="dotted", size=0.2) +
      geom_vline(xintercept=c(1,ranges_struct$end), linetype=2, size=0.3)+
      scale_x_discrete(limits=c(1,length(d$PDB))) +
      labs(x="", fill="PDB", y=expression(Delta*Delta*"G(kcal/mol"))
    return(p)
  }

#seperate plots
  ResidueDDGplot = function(s, gdp){
    if (gdp ==FALSE){
      d= subset(Table_DDG_residues, PDB == s & GDP == gdp)
    }
    if (gdp ==TRUE){
      d= subset(Table_DDG_residues, PDB == s & GDP == gdp & COMP %in% c("GDP", "GDPA", "GDPC"))
    }
    conservation= d$Shannon_Adj*max(abs(d$DDG))
    p= ggplot(data = d, aes(x=1:length(PDB), y=DDG, fill=CHAIN, color=CHAIN)) +
             geom_col(width = 0.1) +
             geom_point(size = 2) +
             geom_point(data = d, aes(fill=d$CHAIN, shape=d$PTM, color=d$CHAIN), size = 7, alpha = 0.5) +
             geom_area(aes(y=conservation), alpha = 0.2, colour = "grey", fill = "grey") +
             scale_color_manual(values = colPal) +
             scale_fill_manual(values = colPal) +
             scale_shape_discrete(limits=c(levels(Table_DDG_residues$PTM))) +
             scale_y_continuous(sec.axis = sec_axis(~./max(abs(d$DDG)), name = "conservation rate")) +
             GGtheme + 
             expand_limits(y=c(-1.5,1)) +
             labs(x="Residue", shape="PTM", y=expression(Delta*Delta*"G (kcal/mol)")) 
    
    return(p)
  }
  
  
  # saving residue decompositions
 # for (pdb in unique (subset(Table_DDG_complex, GDP==TRUE)$PDB)){
 #   png(paste(pdb,"_resDecompGDP.png", sep=""), width = 1500, height = 350)
 #   x = pdb
 #   p = ResidueDDGplot(x, TRUE) + 
 #     theme(axis.text= element_text(size=15), legend.position = "right",
 #           title=element_text(size=13)) +
 #     labs(color="Chain", fill="Chain", shape="PTM", title=paste("Residue decomposition of",toupper(x)))
 #   print(p)
 #   dev.off()
 #   }

 p = ResidueDDGplot("2ik8", TRUE) + 
    theme(axis.text= element_text(size=15), legend.position = "top",
          title=element_text(size=13)) +
    labs(color="", fill="", shape="")
 p
 
#DDG plot for all
 x    = subset(Table_DDG_complex, GDP==FALSE)
 x_ac = subset(Table_DDG_residues, select= c("DDG", "PDB", "COMP"), GDP==FALSE &  PTM == "Ac")
 x_ph = subset(Table_DDG_residues,  select= c("DDG", "PDB", "COMP"),GDP==FALSE & PTM == "Ph")
 
 compTrans <- function(comp){
   if (comp == 'Dimer'){
     return('A')
   }
   else{
     return(substr(comp, 6, 6))
   }
 }
 
 x$chain    = unlist(lapply(x$COMP, function(x) compTrans(x)))
 x_ac$chain = unlist(lapply(x_ac$COMP, function(x) compTrans(x)))
 x_ph$chain = unlist(lapply(x_ph$COMP, function(x) compTrans(x)))
 
 x$component    = paste(x$PDB, x$chain)
 x_ac$component = paste(x_ac$PDB, x_ac$chain)
 x_ph$component = paste(x_ph$PDB, x_ph$chain)
 
 x_ac_aggr = aggregate(x_ac$DDG, FUN = sum, by=list(x_ac$component))
 x_ph_aggr = aggregate(x_ph$DDG, FUN = sum, by=list(x_ph$component))
 x_ac_aggr$cnt_ac = table(x_ac$component)
 x_ph_aggr$cnt_ph = table(x_ph$component)
 
 colnames(x_ac_aggr)=c("component", "DDG_ac", "cnt_ac")
 colnames(x_ph_aggr)=c("component", "DDG_ph", "cnt_ph")
 
 x=merge(x, x_ac_aggr, all.x=TRUE)
 x=merge(x, x_ph_aggr, all.x=TRUE)

 x$component = with(x, reorder(toupper(component), DDG))
 x$cnt_ac    = with(x, reorder(cnt_ac, DDG))
 x$cnt_ph    = with(x, reorder(cnt_ph, DDG))
 x$cnt_ph[is.na(x$cnt_ph)] <- ""


 ggplot(x, aes(x=component, y=DDG)) + 
   GGtheme +
   geom_col(alpha=0.5) +
   geom_col(data= x, aes(x=component, y=DDG_ac, fill="Acetylation "), alpha=0.8, width=0.4, position="dodge") +
   geom_col(data= x, aes(x=component, y=DDG_ph, fill= "Phosphorylation"), alpha=0.8, width=0.4, position="dodge") +
   geom_text(data= x, aes(label=cnt_ac), y=-33, col="#0072B2") +
   geom_text(data= x, aes(label=cnt_ph), y=-40, col="#D55E00") +
   theme(axis.text.x = element_text(angle=90, size=10),
         axis.text.y = element_text( size=10),
         legend.direction = "vertical",
         legend.box = "horizontal",
         legend.position = c(0.005,0.995),
         legend.justification = c(0, 1),
         legend.title = element_blank()) + 
   labs(x="",  y=expression(Delta*Delta*"G (kcal/mol)"))+
   # scale_y_continuous(
   #   sec.axis = sec_axis(~ . * 1, name = expression(Delta*Delta*"G contribution (kcal/mol)"))
   # ) +
   scale_fill_manual(labels=c(expression("Acetylation        "*Delta*Delta*"G contribution (kcal/mol)"),
                              expression("    Phosphorylation "*Delta*Delta*"G contribution (kcal/mol)")),
                     values=c("#0072B2", "#D55E00"))
 
#local binding influence of P and A
  x=droplevels(subset(Table_DDG_residues,PTM %in% c("Ph", "Ac") & GDP ==FALSE))
  cnt=table(abs(x$DDG)<0.5, x$RES)
  tot=table(x$RES)
  ggplot(subset(x, abs(DDG)>0.5), aes(x=RES_short, y=DDG, fill=PTM)) + 
     geom_violin(trim=TRUE, scale="width") + #use adjust for bandwith
     geom_boxplot(width=0.1) +
     GGtheme +
     labs(x="Residue",  y=expression(Delta*Delta*"G (kcal/mol)")) +
     stat_summary(geom = 'text', label = cnt[1,], fun.y = max, vjust = -0.75, hjust = 1.5)+
     stat_summary(geom = 'text', label = tot, fun.y = max, vjust = -0.75, hjust = -0.4) +
     stat_summary(geom = 'text', label = rep("/",4), fun.y = max, vjust = -0.75, hjust = -0.1)+
     scale_fill_manual(values=c("#0072B2","#D55E00")) +
     theme(axis.text = element_text(size=10), 
           legend.position="none")
     
#conservation vs DDG
  x=subset(Table_DDG_residues, GDP == FALSE & PTM %in% c("Ac", "Ph"))
  x$group <- cut(x$DDG, 
                 #breaks = c(-Inf, -5, -3, -1, 1, 3, 5, Inf), 
                 #labels = c("<-5", "[-5,-3]", "[-3,-1]", "[-1,1]", "[1,3]", "[3,5]",">5"), 
                 breaks = c(-Inf, -5, -1, 1, 5, Inf), 
                 labels = c("<-5", "[-5,-1]", "[-1,1]", "[1,5]",">5"), 
                 right = FALSE)
  ggplot(x, aes(x=group, y=consPTM, fill=PTM))+
    #geom_violin(draw_quantiles = TRUE) +
    geom_dotplot() +
    geom_boxplot()+
    facet_grid(rows=x$RES) +
    GGtheme +
    scale_fill_manual(labels= c("Acetylated", "Phosphorylated"), values =c("#0072B2", "#D55E00")) +
    theme(strip.text=element_text(size=15),
          strip.background=element_rect(fill="white")) +
    labs(x= expression(Delta*Delta*"G contribution (kcal/mol)") , y= "CA")
  

  

#################
#Playground code#
#################

  
#trying to find correlations
  d= subset(Table_DDG_residues, DDG < 10 & abs(DDG)>2)
  plot(d$Shannon_Adj,d$DDG)

#some conservation against PTM plotting
  d=subset(Table_DDG_residues, RES_short=="S" & COMP %in% c("Dimer", "chainA"))
  ggplot(d, aes(x=1:nrow(d), y=ConsSer, col=PTM)) + geom_point()

#PTM DDG vs no PTM DDG
  x = Table_DDG_residues$DDG[is.na(Table_DDG_residues$PTM) & abs(Table_DDG_residues$DDG) >2 ]
  y = Table_DDG_residues$DDG[!is.na(Table_DDG_residues$PTM) & abs(Table_DDG_residues$DDG) >2]
  t.test(x,y)

  d = subset(Table_DDG_residues, abs(Table_DDG_residues$DDG) >1 & abs(Table_DDG_residues$DDG)<20 & RES_short %in% c("K", "S", "T", "Y"))
  ggplot(d, aes(x=RES_short, y=d$DDG, fill=!(is.na(d$PTM)) )) + 
    geom_boxplot() +
    theme_bw() +
    labs(x= "residue", y="DDG (kcal/mol)", fill= "PTM")

#PTM ShannonAdj vs noPTM Shannon Adj
  d = subset(Table_DDG_residues, RES_short %in% c("K", "S", "T", "Y"))
  ggplot(d, aes(x=RES_short, y=d$Shannon_Adj, fill=!(is.na(d$PTM)) )) + 
    geom_boxplot() +
    theme_bw() +
    labs(x= "residue", y="Conservation", fill= "PTM")
  
#plot of DDG influence factored by residue
  d = subset(Table_DDG_residues, abs(DDG)>2 & DDG<100)
  ggplot(d, aes(x=RES_short, y=DDG, fill=Charge)) + geom_boxplot()
  
#plot of resdiue conservation factored by residue
  d = Table_DDG_residues
  ggplot(d, aes(x=RES_short, y=Shannon_Adj, fill=Charge)) + geom_boxplot()

  x= subset(Table_DDG_residues, PDB=="2ode" & COMP == "GDP")
  sum(x$DDG)
  

  

   # x=subset(Table_DDG_residues, select = c("PDB", "COMP","CHAIN", "RESN", "RES_short","PTM", "DDG"), GDP==TRUE & abs(DDG) > 5 & COMP == "GDP")
   # write.table(x,  file=paste("C:/Users/Jasper/Documents/School/2e_master_bioinformatica/Thesis/tables/","DDG_large.csv", sep=""), col.names=FALSE, row.names = FALSE, quote=FALSE)
    
  
#Old conservation loading code, not deled in case it becomes usefull someday
  # Cons = data.frame(PDB = character(),
  #                   CHAIN = character() ,
  #                   RESN = integer(), 
  #                   RES = character(), 
  #                   SHANNON = double())
  # pdb = '1fl7'
  # chain = 'A'
  # path = paste("~/School/2e_master_bioinformatica/Thesis/tables/Conservation_scores/", pdb, "_", chain,"_cons",sep="" )
  # x <- read_csv(path, "\t", trim_ws = TRUE, skip = 3)
  # y = t(apply(x, 1 , function(x) unlist(strsplit(as.character(x), " "))[c(1,12,20)] ))
  # y = cbind(rep(chain, length(y[,1])), y)
  # y = cbind(rep(pdb, length(y[,1])), y)
  # colnames(y)= c("PDB", "CHAIN","RESN", "RES", "SHANNON")
  # Cons = rbind(Cons, y)
  # Cons = subset(Cons, Cons$RES != "-")
  