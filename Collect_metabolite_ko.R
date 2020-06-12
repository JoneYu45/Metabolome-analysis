#Library activated
library(tcltk)
library(rvest)

# Collect pathyway ID
mapIDcollection <- c(NA,NA,NA)
n = 0
for (i in which(!is.na(Chemicaltaxa[,1]))) {
  #Progress repsort
  info <- sprintf("KEGG mapID Mapping finished by  %f", round(n/length(which(!is.na(Chemicaltaxa[,1])))*100,6))
  setTkProgressBar(pb, 
                   round(n/length(which(!is.na(Chemicaltaxa[,1])))*100,6), 
                   sprintf("KEGG Mapping Started at (%s)", begin_time), 
                   info)
  
  try({
    targetID <- gsub("  <kegg_id>","",Chemicaltaxa[i,5])
    targetID <- gsub("</kegg_id>","",targetID)
    targetID
    Chemicaltaxa[i,5]
    url <- paste("https://www.genome.jp/dbget-bin/www_bget?cpd:",
                 targetID,
                 sep = "")
    webpage <- read_html(url)
    mapID_html <- html_nodes(webpage,"nobr a")
    mapID <- html_text(mapID_html)
    mapIDcollection0 <- mapID[grep("map",mapID)]
    mapIDcollection0 <- cbind(array(Chemicaltaxa[i,1],dim = c(length(mapIDcollection0),1)),
                              mapIDcollection0,
                              array(NA,dim = c(length(mapIDcollection0),1)))
    mapIDcollection <- rbind(mapIDcollection, mapIDcollection0)
  },silent = T)
  
  n = n + 1
  Sys.sleep(5) #According to the Crawl delay
}
mapIDcollection <- mapIDcollection[-1,]

# Collect KO pathyway ID
for (i in 1:length(mapIDcollection)) {
  #Progress repsort
  info <- sprintf("KEGG koID Mapping finished by %f",round(i/length(mapIDcollection)*100,6))
  setTkProgressBar(pb, 
                   round(i/length(mapIDcollection)*100,6), 
                   sprintf("KEGG Mapping Started at (%s)", begin_time), 
                   info)
  
  try({
    url <- paste("https://www.kegg.jp/dbget-bin/www_bget?",
                 mapIDcollection[i,2],
                 sep = "")
    webpage <- read_html(url)
    koID_html <- html_nodes(webpage,"nobr a")
    koID <- html_text(koID_html)
    mapIDcollection[i,3] <- koID[grep("ko",koID)]
  },silent = T)
  
  Sys.sleep(5) #According to the Crawl delay
}

