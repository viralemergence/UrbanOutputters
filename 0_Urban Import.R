
# 0_Urban Import #####

{
  
  library(tidyverse); library(data.table); library(magrittr); library(ggregplot); library(cowplot)
  library(geiger);library(ape);library(picante)
  
  theme_set(theme_cowplot())
  
}

"Urban/Urban.csv" %>% 
  fread %>% 
  as.data.frame() %>% 
  rename(UrbanStatus = Urban) %>%
  mutate_at("Species", ~.x %>% str_trim %>% str_replace_all(" ", "_")) %>%
  mutate_at("UrbanStatus", ~.x %>% CamelConvert) ->
  
  UrbanDFInitial

Panth1 <- read.delim("Data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename_all(~str_replace(.x, "MSW05_", "h")) %>%
  rename(Sp = hBinomial) %>%
  mutate_at("Sp", ~.x %>% str_trim %>% str_replace_all(" ", "_"))

Panth1 %>% 
  filter(hOrder%in%c("Cetacea", "Sirenia")|
           hFamily%in%c("Phocidae", "Odobenidae", "Otariidae")) %>% 
  pull(Sp) ->
  MarineSp

# Parasite stuff ####

# GMPD ####

GMPD_main <- read.csv("GMPD_datafiles/GMPD_main.csv")

GMPD_main %<>% 
  #filter(ParType == "Virus") %>% 
  select(HostCorrectedName, ParasiteReportedName, ParasiteCorrectedName, 
         ParasiteType = ParType) %>%
  mutate_all(~str_replace_all(.x, " ", "_"))

Remove1String <- function(a, Sep = "_"){
  
  String <- a %>% str_split(Sep) %>% unlist()
  
  String <- String[-1]
  
  String %>% paste(collapse = Sep)
  
}

GMPD_main %<>% 
  mutate_at(c("ParasiteReportedName", "ParasiteCorrectedName"), 
            ~ifelse(ParasiteType == "Virus", sapply(.x, Remove1String), .x))

# HP3 ####

HP3 <- 
  read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/associations.csv") %>% 
  data.frame()

HP3 %<>% 
  rename(Parasite = vVirusNameCorrected,
         Host = hHostNameFinal) %>% 
  mutate_at(c("Host", "Parasite"), ~str_replace_all(.x, " ", "_"))

# EID2 ####

EIDSpecies <- read.csv("Data/EID/SpeciesInteractions_EID2.csv", header = T)

EIDSpecies %<>% 
  mutate_at(c("Carrier", "Cargo"),  ~str_replace(.x, " ", "_") %>% CamelConvert)# %>%
  # filter(#Cargo.classification == "Virus", 
  #   Carrier.classification %in%c("Mammal", 
  #                                "Primate", "Rodent",
  #                                "Human"))

EIDLocations <- read.csv("Data/EID/LocationInteractions_EID2.csv", header = T)

# Shaw ####

Shaw <- read.csv("Data/Shaw.csv", header = T)

Shaw %<>% 
  mutate_at(c("HostSpecies", "Species"),  
            ~.x %>% str_trim %>% str_replace_all(" ", "_") %>% 
              tolower %>% CamelConvert) %>% 
  mutate(Dataset = "Shaw") %>% 
  rename(Parasite = Species, Host = HostSpecies,
         ParasiteType = Type)

Shaw %>% 
  filter(HostGroup %in% c("Mammalia", 
                          "Primates", "Rodentia",
                          "Human", "Ungulates", "Carnivora", "Cetacea", "Chiroptera"))

# Combining ####

HP3 %>% select(Host, Parasite) %>% 
  mutate(Dataset = "HP3", ParasiteType = "Virus") %>% 
  bind_rows(GMPD_main %>% select(Host = HostCorrectedName, 
                                 #Virus = ParasiteReportedName) %>% 
                                 Parasite = ParasiteCorrectedName,
                                 ParasiteType) %>% 
              mutate(Dataset = "GMPD2")) %>% 
  bind_rows(EIDSpecies %>% select(Host = Carrier, 
                                  Parasite = Cargo,
                                  ParasiteType = Cargo.classification) %>% 
              mutate(Dataset = "EID2")) %>% 
  mutate_at("Parasite", ~.x %>% tolower %>% CamelConvert %>% str_replace_all(" ", "_")) %>% 
  bind_rows(Shaw %>% select(Host, Parasite, Dataset, ParasiteType)) %>% 
  
  arrange(Host, Parasite) -> 
  
  FullAssociations

FullAssociations %<>% 
  mutate_at("ParasiteType", ~.x %>% str_replace_all("Fungi", "Fungus"))

# Cleaning Names Initially ####

NameCorrects <- c("corona_virus" = "coronavirus",
                  "Andes_hantavirus" = "Andes_virus",
                  "Bovine_fever_ephemerovirus" = "Bovine_ephemeral_fever_virus",
                  "Rift_valley_phlebovirus" = "Rift_valley_fever_virus",
                  "Rift_valley_fever_phlebovirus" = "Rift_valley_fever_virus",
                  "Carnivore_protoparvovirus_1" = "Feline_panleukopenia_virus",
                  "Junin_mammarenavirus" = "Junin_virus",
                  "Laguna_negra_hantavirus" = "Laguna_negra_virus",
                  "Akabane_orthobunyavirus" = "Akabane_virus",
                  "Crimean-congo_hemorrhagic_fever_nairovirus" = "Crimean-congo_hemorrhagic_fever_virus")

FullAssociations %>% filter(Host %in% Panth1$Sp) %>% filter(Host == unique(Host)[23])

FullAssociations$Parasite %<>% str_replace_all(NameCorrects)

Omit <- c("Not_identified_to_genus",
          "Abolished",
          "_sp[.]",
          "Splitted_in_ictv",
          "Cjd_agent",
          "Bse_agent",
          "No_binomial_name")

Omit %<>% paste0("$") %>% paste0(collapse = "|")

FullAssociations %<>% filter(!str_detect(Parasite, Omit))

FullAssociations %<>% 
  mutate_at("Parasite", as.character) %>% 
  mutate_at("Parasite", ~.x %>% 
              str_replace_all(c("[()]" = '',
                                "\\]" = '',
                                "\\[" = '',
                                "\\'" = '')) %>% 
              str_trim %>% str_replace_all(" ", "_") %>% 
              CamelConvert
  ) 

if(!file.exists("Intermediate/Translations.csv")){
  
  # Host Names ####
  
  library(taxize)
  library(rentrez)
  
  set_entrez_key("acde4a6a16c824df046d750a42aa0b563208")
  Sys.getenv("ENTREZ_KEY")
  
  # Parasite Names ####
  
  FullAssociations %>% filter(Host %in% Panth1$Sp) %>% 
    pull(Parasite) %>% unique %>% sort -> ParasiteNames
  
  ParasiteNames %<>% 
    str_replace_all(c("_" = " ",
                      "[()]" = '',
                      "\\]" = '',
                      "\\[" = '',
                      "\\'" = '')) %>% 
    CamelConvert %>% sort
  
  ParasiteNames %>% lapply(function(x){
    
    print(x)
    
    ## classify
    cset <- classification(x, db = "ncbi", rows = 1)
    
    cset2 <- cset[[1]]
    
    if(!is.na(cset2)){
      
      cset2 %>% slice(n()) %>% pull(name)
      
    }else{
      
      return(x)
      
    }
    
  }) -> CleanedParasiteNames
  
  NameTranslate <- data.frame(
    
    ParasiteNames,
    CleanedParasiteNames = CleanedParasiteNames %>% unlist
    
  )
  
  NameTranslate %>% write.csv("Intermediate/Translations.csv", row.names = F)
  
}else{
  
  NameTranslate <- read.csv("Intermediate/Translations.csv", header = T)
  
}

NameTranslate %<>% 
  mutate_all(
    ~.x %>% str_replace_all(c("[()]" = '',
                              "\\]" = '',
                              "\\[" = '',
                              "\\'" = '')) %>% 
      str_trim %>% str_replace_all(" ", "_") %>% 
      CamelConvert
  )

Replace <- NameTranslate$ParasiteNames %>% as.character
names(Replace) <- NameTranslate$CleanedParasiteNames %>% as.character

Replace <- Replace[str_length(Replace)>0]

FullAssociations %>% 
  mutate_at("Parasite", ~str_replace_all(.x, Replace)) %>% 
  arrange(Host, Parasite) %>% 
  select(Host, Parasite, ParasiteType, Dataset) -> FullAssociations1

Replace2 <- NameTranslate$CleanedParasiteNames %>% as.character
names(Replace2) <- NameTranslate$ParasiteNames %>% as.character

Replace2 <- Replace2[str_length(Replace2)>0]

FullAssociations %>% 
  mutate_at("Parasite", ~str_replace_all(.x, Replace2)) %>% 
  arrange(Host, Parasite) %>% 
  select(Host, Parasite, ParasiteType, Dataset) -> 
  FullAssociations2

FullAssociations <- FullAssociations2

# Adding Clover ####

"Data/clover" %>% list.files(full.names = T, pattern = "FlatFile") %>%
  map_dfr(read.csv) -> Clover

Clover %<>%
  mutate(ParasiteType = PathogenType, Parasite = Pathogen) %>%
  mutate_at("ParasiteType", ~.x %>% str_split("/") %>% map_chr(1) %>%
              CamelConvert %>% str_replace("Fungi", "Fungus")) %>%
  dplyr::select(Host,
                Parasite,
                ParasiteType,
                Dataset = Database) %>%
  mutate_at(c("Host", "Parasite"), ~.x %>% str_trim %>%
              str_replace_all(" ", "_") %>%
              CamelConvert)

# FullAssociations %<>% filter(ParasiteType == "Arthropod") %>%
#   bind_rows(Clover)

"Data/Virion.csv.gz" %>% vroom::vroom() -> Virion

Virion %>% mutate(ParasiteType = "Virus",
                  Parasite = Virus,
                  Database = "Virion") %>%
  filter(ICTVRatified == T,
         HostNCBIResolved == T,
         VirusNCBIResolved == T) %>%
  dplyr::select(Host,
                Parasite,
                ParasiteType,
                Dataset = Database) %>%
  mutate_at(c("Host", "Parasite"), ~.x %>% str_trim %>%
              str_replace_all(" ", "_") %>%
              CamelConvert) -> Virion

FullAssociations %<>% filter(ParasiteType == "Arthropod") %>%
  # bind_rows(Clover) %>%
  bind_rows(Clover %>% filter(!ParasiteType == "Virus")) %>%
  bind_rows(Virion)

FullAssociations$Parasite %<>% 
  str_remove("\\[") %>% 
  str_remove("\\]") %>% 
  # head(20) %>% 
  CamelConvert

# Putting together ####

FullAssociations %<>% filter(Host %in% Panth1$Sp)

FullAssociations %<>% select(Host, Parasite, ParasiteType, Dataset) %>% unique

FullAssociations %>% 
  filter(Host == "Homo_sapiens") %>% 
  pull(Parasite) %>% unique ->
  
  ZoonoticParasites

FullAssociationsNoEID2 <- 
  FullAssociations %>% filter(!Dataset == "EID2")

FullAssociations %>% 
  dplyr::select(Host, Parasite) %>% 
  unique %>%
  mutate(Zoonotic = as.numeric(Parasite%in%ZoonoticParasites)) %>%
  group_by(Host) %>% 
  summarise(Richness = nunique(Parasite),
            ZoonoticRichness = sum(Zoonotic)) %>%
  mutate(NonZoonoticRichness = Richness - ZoonoticRichness) %>%
  mutate(PropZoonotic = ZoonoticRichness/Richness) -> 
  
  AssociationRichness

ParasiteTypes <- c(
  
  "Helminth", "Arthropod", "Bacteria", "Protozoa", "Virus", "Fungus"
  
)

ParasiteTypes %>% 
  map(~FullAssociations %>% 
        filter(ParasiteType == .x) %>% 
        dplyr::select(Host, Parasite) %>% 
        unique %>%
        mutate(Zoonotic = as.numeric(Parasite%in%ZoonoticParasites)) %>%
        group_by(Host) %>% 
        summarise(Richness = nunique(Parasite),
                  ZoonoticRichness = sum(Zoonotic)) %>%
        mutate(NonZoonoticRichness = Richness - ZoonoticRichness) %>%
        mutate(PropZoonotic = ZoonoticRichness/Richness) %>% 
        mutate_at(vars(matches("Richness|PropZoonotic")), ~ifelse(is.na(.x), 0, .x)) %>% 
        rename_at(vars(matches("Richness|PropZoonotic")), 
                  function(a) paste0(a, ".", .x))) %>% 
  reduce(~full_join(.x, .y, by = "Host")) %>%
  full_join(AssociationRichness, ., by = "Host") ->
  
  AssociationRichness

Resps <- c("ZoonoticRichness", "Richness", "PropZoonotic", "NonZoonoticRichness")

AssociationRichness %<>% rename_at(Resps, ~paste0(.x, ".Overall"))

# Uniting urban data and parasite data ####

UrbanDFInitial %>% 
  select(Sp = Species, UrbanStatus) %>% 
  full_join(Panth1) %>% 
  mutate(Marine = as.numeric(Sp %in% MarineSp)) %>% 
  filter(!Marine) %>% 
  mutate_at("UrbanStatus", ~ifelse(is.na(.x), "Non-urban", .x)) %>% 
  mutate_at("UrbanStatus", ~factor(as.character(.x), 
                                   levels = c("Non-urban", "Visitor", "Visitor/dweller",
                                              "Dweller"))) %>%
  full_join(AssociationRichness, by = c("Sp" = "Host")) %>%
  mutate(Pseudoabsence = as.numeric(is.na(Richness.Overall))) %>%
  mutate_at(vars(contains("Richness")), 
            ~ifelse(is.na(.x), 0, .x)) -> 
  
  UrbanDF

UrbanDF %<>% arrange(Sp)

UrbanDF %<>% filter(!Sp == "Homo_sapiens")

UrbanDF %<>% mutate(UrbanBinary = as.numeric(!UrbanStatus == "Non-urban") %>% as.factor)

Domestic <- read.csv("Data/Domestic.csv")

Domestic$Species.and.subspecies %>% 
  str_split("[(]") -> Bracket1 # Detecting species names

Bracket1 %>% map(~.x[str_detect(.x, "[)]")]) -> Bracket2 # detecting incomplete species names

Bracket2 %>% map(~.x %>% str_split("[)]") %>% map(1)) -> Bracket3

1:length(Bracket3) %>% map(function(b){
  
  # print(b)
  
  a <- Bracket3[[b]]
  
  Genus <- a[[1]] %>% str_split(" ") %>% extract2(1) %>% extract2(1)
  
  if((a[[1]] %>% str_split(" ") %>% extract2(1) %>% length)>1){
    
    a %>% map(~.x %>% str_split(" ") %>% map_chr(2) %>% paste0(Genus, "_", .))
    
  }else NA
  
}) %>% unlist %>% unique %>% sort -> DomesticSpecies

AddDomestic <- c("Capra_hircus", 
                 "Neovison_vison")

DomesticSpecies %<>% c(AddDomestic) %>% unique %>% sort

UrbanDF %>% 
  mutate(DomesticBinary = as.numeric(Sp %in% DomesticSpecies) %>% as.factor) -> 
  UrbanDF

CitationDF <- readRDS("Intermediate/CitationDF.rds")

# CitationDF <- readRDS("Intermediate/OldCitationDF.rds")

UrbanDF %<>% left_join(CitationDF, by = "Sp") %>% 
  mutate(LogCites = log10(Citations + 1)) %>% 
  mutate(LogRichness = log(Richness.Overall + 1))# %>% 
#mutate_at(vars(matches("Richness")), ~log10(.x + 1))

# Phylogenetic data ####

STFull <- read.nexus("Data/ele_1307_sm_sa1.tre")[[1]]
FullSTMatrix <- as.data.frame(cophenetic(STFull)) %>% as.matrix

FullSTMatrix["Homo_sapiens",] %>% as.data.frame() %>% 
  rownames_to_column %>% rename("HumanDistance" = ".") %>% 
  left_join(UrbanDF, ., by = c("Sp" = "rowname")) -> 
  UrbanDF

# Fixing phenotypic stuff ####

UrbanDF %>% 
  rename(Mass = X5.1_AdultBodyMass_g, 
         HomeRange = X22.1_HomeRange_km2) %>% 
  mutate(LogMass = log(Mass)) ->
  
  UrbanDF

DegreeGet <- function(a){
  
  a[a>0] <- 1
  rowSums(a) %>% return
  
}

# Adding Spatial Data ####

load("Data/FullRangeOverlapMercator.Rdata") # From Albersnet ####

FullRangeAdj %>% colSums -> SympatryStrengthVector
FullRangeAdj %>% DegreeGet -> SympatryDegreeVector

SympatryStrengthVector[UrbanDF$Sp] -> UrbanDF$SympatryStrength
SympatryDegreeVector[UrbanDF$Sp] -> UrbanDF$SympatryDegree

UrbanDF %>% 
  select(contains("Sympatry")) %>% 
  mutate_all(log10) ->
  UrbanDF[,paste0("Log", c("SympatryStrength", "SympatryDegree"))]

# Comparing with HP3 urban variables ####

HP3Hosts <- read.csv("Data/Hosts.csv")

logp = function(x){   # Fn to take log but make zeros less 10x less than min
  # x[is.na(x)] <- 0
  m = min(x[ x > 0], na.rm=T)
  x = log( x + m )
  return(x)
}

HP3Hosts %>% mutate(
  UrbRurPopRatioLn = logp((urbc_2005AD)/(rurc_2005AD)),
  UrbRurPopRatioChg = logp((urbc_2005AD)/(rurc_2005AD)) - logp((urbc_1970AD)/(rurc_1970AD))) %>% 
  select(hHostNameFinal, UrbRurPopRatioLn, UrbRurPopRatioChg) %>% 
  left_join(UrbanDF, ., by = c("Sp" = "hHostNameFinal")) -> 
  
  UrbanDF

# Importing Eskew Data ####

EskewPCA <- read.csv("EskewFiles/reservoir_data_for_greg.csv")

EskewPCA %<>% mutate_at("PC1", ~ -.x)

EskewPCA %>% 
  mutate(Sp = binom %>% str_trim %>%  str_replace_all(" ", "_")) %>% 
  select(PC1:Sp) %>% 
  left_join(UrbanDF, ., by = "Sp") ->
  UrbanDF

# Eliminating non-Eutherians ####

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia", 
                   "Notoryctemorphia",
                   "Monotremata")

NonEutherianSp <- Panth1[Panth1$hOrder%in%NonEutherians,"Sp"]

UrbanDF %<>% filter(!Sp %in% NonEutherianSp)

UrbanDF %>% dim

# Adding centroids etc ####

CentroidList <- readRDS("Data/CentroidList.rds")
ContinentsInhabitedList <- readRDS("Data/ContinentsInhabitedList.rds")

ContinentsInhabitedList %>% unlist %>% unique -> ContinentVars

CentroidList %<>% 
  map(~.x %>% 
        matrix(nrow = 2) %>% t %>% as.data.frame %>% rename(X = V1, Y = V2) %>% 
        summarise_all(mean)
  )

ContinentDF <- 
  data.frame(Sp = names(ContinentsInhabitedList))

ContinentDF$Continents <- ContinentsInhabitedList

ContinentDF %<>% 
  unnest(Continents) %>% 
  mutate(Presence = 1) %>% 
  tidyr::pivot_wider(#id_cols = "Sp", 
    names_from = "Continents", values_from = "Presence") %>% 
  mutate_at(-1, ~ifelse(is.na(.x), 0, .x))

CentroidList %>% bind_rows(.id = "Sp") %>% 
  left_join(ContinentDF, .) %>% 
  left_join(UrbanDF, .) -> UrbanDF

UrbanDF %>% 
  select(Eurasia:`NAm`) %>% names -> ContinentVar

UrbanDF %>% 
  dplyr::select(all_of(ContinentVar)) %>% rowSums ->
  UrbanDF[,"NContinents"]

UrbanDF %<>% mutate_at("NContinents", ~ifelse(.x == 0, 1, .x))

AreaList <- read.csv("Data/RangeAreas.csv")

UrbanDF %<>% 
  
  left_join(AreaList, by = c("Sp" = "Species"))

(UrbanDF$Area+1) %>% 
  log -> UrbanDF$LogArea

# Elton ####

Elton <- read.delim("Data/Defunct/MamFuncDat.txt")

Elton %<>% 
  mutate(Sp = Scientific %>% str_trim %>% str_replace_all(" ", "_"))

Elton %>% 
  dplyr::select(contains("Diet.")) %>% 
  mutate_all(AsBinary) %>% 
  rowSums(na.rm = T) -> Elton$DietNo

Elton %>% dplyr::select(Sp, DietNo) %>% left_join(UrbanDF, .) -> UrbanDF

library(vegan)

Elton %>% 
  dplyr::select(contains("Diet.")) %>% 
  dplyr::select(-matches("Source|Certainty")) -> Proportions

Proportions %>% diversity -> Diversities

Elton$DietDiversity <- Diversities

# Elton %>% ggplot(aes(DietNo, DietDiversity)) + geom_point() + geom_smooth(method = lm)

Elton %>% dplyr::select(Sp, DietDiversity) %>% 
  left_join(UrbanDF, .) -> 
  UrbanDF

# Wildlife Trade ####

if(0){
  
  Trade <- read.csv("Data/WildlifeTrade.csv")
  
  Traded <- 
    Trade %>% filter(Category == c("Traded")) %>% 
    mutate_at("Species.Name", 
              ~.x %>% 
                str_trim %>% 
                str_replace_all(" ", "_")) %>% 
    pull(Species.Name)
  
  UrbanDF %<>% 
    mutate(Traded = as.factor(as.numeric(Sp %in% Traded)))
  
}

# Endangered ####

IUCNClov <- read.csv("Data/IUCNCLOV.csv")

IUCNClov %<>% mutate_at("Host",
                        ~.x %>% 
                          str_trim %>% 
                          str_replace_all(" ", "_") %>% 
                          CamelConvert)

UrbanDF <- IUCNClov %>% dplyr::select(1:2) %>% left_join(UrbanDF, ., by = c("Sp" = "Host"))

UrbanDF %<>% mutate(RedList = redlistCategory%>% 
                      str_trim %>% 
                      str_replace_all(" ", "_") %>% 
                      str_remove_all("_in_the_Wild") %>% 
                      CamelConvert)

UrbanDF %<>% 
  mutate(Endangered = as.numeric(str_detect(RedList, "Endangered"))) %>% 
  mutate(Data_Deficient = as.numeric(is.na(RedList) | RedList == "Data_Deficient"))

UrbanDF %<>% 
  mutate_at("RedList", 
            ~factor(.x,
                    levels = c("Least_Concern", "Near_Threatened", 
                               "Vulnerable", "Endangered", "Critically_Endangered",
                               "Extinct",
                               "Data_Deficient"
                    )
            ))
