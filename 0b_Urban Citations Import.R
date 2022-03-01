# 0b_Urban Citations Import ####

# Taken from https://github.com/viralemergence/virionette/blob/master/04_predictors/Citations.R

# Takes an hour or two ####

library(easyPubMed)
library(tidyverse)

CiteCounter <- function(Name) {
  
  Query <- Name %>% str_trim %>% #str_replace_all(" ", "-") %>% 
    str_replace_all("_", " ") %>% 
    paste0("\"", ., "\"") %>% 
    get_pubmed_ids
  
  Original <- Query$OriginalQuery
  
  Translation <- Query$QueryTranslation
  
  if(str_count(Original) == (str_count(Translation) - 12)){
    
    Query %>% 
      extract2("Count") %>% 
      as.character %>% 
      as.numeric %>% 
      return
    
  }else{
    
    print(paste0(Name, ": Maybe failure?"))
    
    return(0)
    
  }
  
}

UrbanDF$Citations <- NA

i <- 1

for(i in i:nrow(UrbanDF)){
  
  Sp <- UrbanDF$Sp[i]
  
  print(Sp)
  
  CiteTemp <- CiteCounter(Sp)
  
  UrbanDF$Citations[i] <- CiteTemp
  
  print(paste0("Citations: ", UrbanDF$Citations[i]))
  
}

UrbanDF %<>% mutate(LogCites = log10(Citations + 1))

UrbanDF %>% select(Sp, Citations) -> CitationDF

saveRDS(CitationDF, file = "Intermediate/CitationDF.rds")

# Comparing new and old ####

New <- readRDS("Intermediate/CitationDF.rds")

Old <- readRDS("Intermediate/OldCitationDF.rds")


