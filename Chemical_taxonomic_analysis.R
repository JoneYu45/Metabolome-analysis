#Library activated
library(tcltk)
library(rvest)

#Data input
workbook <- "D:/Metabolome/Data_Process/HMDB/SignificantMeatboltes.csv" #Input the significant different metabolites data
SelectedMetabolites <- read.csv(workbook,header = T)
remove(workbook)
database <- " D:/Metabolome/Data_Process/HMDB/hmdb_metabolites.xml" #Input the HMDB database file

#Information collection
########################
divlinnum = 50000 #The database file is divided into fractions, and the divlinnum represents how many lines each fractions should contain
close(HMDB)
HMDB <- file(database,"r")
line0 <- NA
Chemicaltaxa <- array(NA, dim = c(length(SelectedMetabolites[,1]),5)) #The chemical taxonomy information and the KEGGID will be stored in this array.
colnames(Chemicaltaxa) <- c("Metabolites","direct_parent","super_class","class","KEGGID")
kegg_info <- NA
x = 1 #Judge whether program should stop
y = 1 #No. of the fragment

#Progress report
plot.new()
pb <- tkProgressBar("HMDB Mapping","Started",  0, 100)
info <- sprintf("HMDB Mapping Started")
begin_time <- Sys.time()
setTkProgressBar(pb, 0, sprintf("HMDB Mapping Started at (%s)", begin_time), info)

while (x != 0) {
  #Database divided 
  line = readLines(HMDB,n=divlinnum)
  #Collect the information of interest from the database fragment
  Lname1 <- grep("  <name>", line, ignore.case = T)
  Lname2 <- grep("   <name>", line, ignore.case = T)
  Lname1 <- Lname1[which(!(Lname1 %in% Lname2) == T)]
  
  Lsynonym1 <- grep("    <synonym>", line, ignore.case = T)
  Lsynonym2 <- grep("     <synonym>", line, ignore.case = T)
  Lsynonym1 <- Lsynonym1[which(!(Lsynonym1 %in% Lsynonym2) == T)]
  
  Ldirect_parent <- grep("    <direct_parent>", line, ignore.case = T)
  Lsuper_class <- grep("    <super_class>", line, ignore.case = T)
  Lclass <- grep("    <class>", line, ignore.case = T)
  Lkeggmap <- grep("        <kegg_map_id>", line, ignore.case = T)
  Lkeggid <- grep("  <kegg_id>", line, ignore.case = T)
  
  Wanteddata <- c(line0,line[sort(c(Lname1,Lsynonym1,Ldirect_parent,Lsuper_class,Lclass,
                                    Lkeggmap,Lkeggid))])
  line0 <- line[sort(c(Lname1,Lsynonym1,Ldirect_parent,Lsuper_class,Lclass,
                       Lkeggmap,Lkeggid))]
  
  #Matching the metabolites by name and output the corresponding taxonomy information and KEGGID
  for (i in which(is.na(Chemicaltaxa[,4]))) {
    #Progress report
    info <- sprintf("Finished by %f",round((1 - length(which(is.na(Chemicaltaxa[,1])))/length(SelectedMetabolites[,1]))*100,6))
    info2 <- sprintf("Database fragment No %d",y)
    setTkProgressBar(pb, 
                     round((1 - length(which(is.na(Chemicaltaxa[,1])))/length(SelectedMetabolites[,1]))*100,6), 
                     sprintf("HMDB Mapping Started at (%s)", begin_time), 
                     paste(info, info2, sep = "\t"))
    
    try({
      n <- grep(paste("<name>",SelectedMetabolites[i,2],"</name>",sep = ""),Wanteddata,ignore.case = T)
      if (length(n) == 0){
        n <- grep(paste("<synonym>",SelectedMetabolites[i,2],sep = ""),Wanteddata,ignore.case = T)
      }
      m <- grep("<name>",Wanteddata)
      narrowdown <- Wanteddata[n:(m[which(m > n)][1]-1)]
      Chemicaltaxa[i,1] <- as.character(SelectedMetabolites[i,2])
      Chemicaltaxa[i,2] <- narrowdown[grep("    <direct_parent>", narrowdown)]
      Chemicaltaxa[i,3] <- narrowdown[grep("    <super_class>", narrowdown)]
      Chemicaltaxa[i,4] <- narrowdown[grep("    <class>", narrowdown)]
      Chemicaltaxa[i,5] <- narrowdown[grep("  <kegg_id>", narrowdown)]
      kegg_info1 <- list(narrowdown[grep("        <kegg_map_id>", narrowdown)])
      names(kegg_info1) <- as.character(SelectedMetabolites[i,2])
      kegg_info <- c(kegg_info,kegg_info1)
    }
    ,silent = T)
  }
  
  #Next
  x = length(line[length(line)])
  y = y+1
}


