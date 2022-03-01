
# 1_Urban Models ####

library(tidyverse); library(INLA); library(ggregplot); library(glue); library(fs)

# source("0_Urban Import.R")

dir_create("Output")

# UrbanDF %<>% filter(!Pseudoabsence)

Resps <- c("Richness", "ZoonoticRichness")

FamilyList <- rep("nbinomial", 2)

names(FamilyList) <- Resps

FullCovar <- c("hOrder", 
               "LogCites", 
               "DomesticBinary", "HumanDistance",
               # "Traded",
               # "UrbRurPopRatioLn",
               "LogArea", "DietDiversity",
               "LogMass", "PC1", "PC2", "UrbanBinary")

PhyloList <- 
  IMList <- 
  TestDFList <- 
  list()

ParasiteTypes <- c(
  
  "Overall", "Virus", "Helminth", "Arthropod", "Bacteria", "Protozoa", "Fungus"
  
)

q <- 1

for(q in q:length(ParasiteTypes)){
  
  print(ParasiteTypes[q])
  
  r <- 1
  
  for(r in r:length(Resps[1:2])){
    
    FocalResp <- paste0(Resps[r], ".", ParasiteTypes[q])
    
    print(FocalResp)
    
    UrbanDF %>% 
      select(FocalResp,
             X, Y, Eurasia:`NA`,
             FullCovar, Sp) %>% 
      na.omit %>% droplevels ->
      TestDF
    
    if(ParasiteTypes[q] == "Fungus"){
      
      TestDF <- 
        TestDF[TestDF[,paste0(Resps[r], ".", ParasiteTypes[q])] <= 10,]
      
    }
    
    TestDF[,"Response"] <- TestDF[,FocalResp]
    
    TestDF %>% 
      group_by(hOrder) %>% 
      summarise(N = n(), 
                Prevalence = Prev(Response), 
                NInf = sum(Response>0)) %>% 
      filter(N>20, Prevalence > 0.001, NInf > 1) %>%
      pull(hOrder) ->
      
      OrderInclude
    
    TestDF %>% 
      filter(hOrder %in% OrderInclude) %>% 
      droplevels ->
      
      TestDF
    
    print(nrow(TestDF))
    
    TestDF %>% 
      select(FullCovar) %>% 
      mutate_if(is.numeric, ~.x %>% scale %>% c) ->
      TestDF[,FullCovar]
    
    PhyloMatrix <- 1 - FullSTMatrix[TestDF$Sp, TestDF$Sp] %>% 
      scales::rescale(to = c(10^-6, 1-(10^-6)))
    
    PhyloList[[FocalResp]] <- PhyloMatrix
    
    TestDF$Phylogeny <- match(TestDF$Sp, colnames(PhyloMatrix))
    
    IM1 <- INLAModelAdd(
      Response = FocalResp,
      Explanatory = FullCovar %>% setdiff("LogCites"), # %>% setdiff("UrbanBinary"),
      # Add = c(#"UrbanBinary:hOrder", 
      #         # "LogCites",
      #         # "Traded:hOrder",
      #         #"UrbanBinary:LogCites"
      #   # "UrbanBinary", "UrbRurPopRatioLn"
      #   ),
      Add = c("LogCites"),
      Family = FamilyList[[Resps[r]]],
      Data = TestDF,
      Base = T,
      AllModels = T,
      AddSpatial = T,
      Groups = T,
      GroupVar = "UrbanBinary"
      
    )
      
    # IM1$FinalModel %>% Efxplot
    
    IMList[[FocalResp]] <- IM1
    
    IM1$FinalModel %>% Efxplot %>% plot
    
    if(0){ 
      
      inla.setOption(num.threads = 8)
      
      glue("{FocalResp} ~ {as.character(IM1$FinalFormula %>% str_split(' [+] ') %>% extract2(3) %>% setdiff(c('hOrder')) %>% paste0(collapse = ' + '))} + {last(FullCovar)}") %>% 
        as.formula -> F1
      
      IM2 <- inla(
        data = TestDF,
        F1,
        family = FamilyList[[Resps[r]]],
        control.compute = list(dic = TRUE)
        
      )
      
      IMList[[FocalResp]]$Phylogenetic <- IM2
      
    }
    
    IMList[[FocalResp]] %>%
      saveRDS(file = paste0("Output/Models_", FocalResp, ".rds"))
    
  }
}

# Continent slopes ####

UrbanDF %>% 
  select(Eurasia:`NAm`) %>% names -> ContinentVar

q <- 1

for(q in q:length(ParasiteTypes)){
  
  print(ParasiteTypes[q])
  
  r <- 1
  
  for(r in r:length(Resps[1:2])){
    
    FocalResp <- paste0(Resps[r], ".", ParasiteTypes[q])
    
    print(FocalResp)
    
    UrbanDF %>% 
      select(FocalResp,
             X, Y, Eurasia:`NA`, NContinents,
             FullCovar, Sp) %>% 
      na.omit %>% droplevels ->
      TestDF
    
    TestDF %<>% mutate(AbsLat = abs(Y) %>% kader:::cuberoot())
    
    if(ParasiteTypes[q] == "Fungus"){
      
      TestDF <- 
        TestDF[TestDF[,paste0(Resps[r], ".", ParasiteTypes[q])]<20,]
      
    }
    
    TestDF[,"Response"] <- TestDF[,FocalResp]
    
    TestDF %>% 
      group_by(hOrder) %>% 
      summarise(N = n(), 
                Prevalence = Prev(Response), 
                NInf = sum(Response>0)) %>% 
      filter(N>20, Prevalence > 0.001, NInf > 1) %>% 
      pull(hOrder) ->
      
      OrderInclude
    
    TestDF %>% 
      filter(hOrder %in% OrderInclude) %>% 
      droplevels ->
      
      TestDF
    
    TestDF %>% 
      select(FullCovar) %>% 
      mutate_if(is.numeric, ~.x %>% scale %>% c) ->
      TestDF[,FullCovar]
    
    print(nrow(TestDF))
    
    PhyloMatrix <- 1 - FullSTMatrix[TestDF$Sp, TestDF$Sp] %>% 
      scales::rescale(to = c(10^-6, 1-(10^-6)))
    
    PhyloList[[FocalResp]] <- PhyloMatrix
    
    TestDF$Phylogeny <- match(TestDF$Sp, 
                              colnames(PhyloMatrix))
    
    IM1 <- INLAModelAdd(
      Response = FocalResp,
      Explanatory = FullCovar %>% c(ContinentVar),
      Add = list(#ContinentVar, 
        "AbsLat",
        paste0(ContinentVar, ":UrbanBinary"),
        paste0("AbsLat", ":UrbanBinary"),
        # paste0("NContinents", ":UrbanBinary"),
        paste0("LogArea", ":UrbanBinary")),
      Family = FamilyList[[Resps[r]]],
      Groups = T, GroupVar = "UrbanBinary", 
      Data = TestDF,
      Delta = 5,
      AddSpatial = T
      
    )
    
    IMList[[FocalResp]] <- IM1
    
    IM1$FinalModel %>% Efxplot %>% plot
    
    if(0){ 
      
      inla.setOption(num.threads = 8)
      
      glue("{FocalResp} ~ {as.character(IM1$FinalFormula %>% str_split(' [+] ') %>% extract2(3) %>% setdiff(c('hOrder')) %>% paste0(collapse = ' + '))} + {last(FullCovar)}") %>% 
        as.formula -> F1
      
      IM2 <- inla(
        data = TestDF,
        F1,
        family = FamilyList[[Resps[r]]],
        control.compute = list(dic = TRUE)
        
      )
      
      IMList[[FocalResp]]$Phylogenetic <- IM2
      
    }
    
    IMList[[FocalResp]] %>%
      saveRDS(file = paste0("Output/ContinentModels_", FocalResp, ".rds"))
    
  }
}

# Multivariate ####

MultivResps <- ParasiteTypes[2:7] %>% paste0("Richness.", .)

UrbanDF %>% filter(Sp %in% IMList[[1]]$Data$Sp) %>% 
  select(FocalResp,
         all_of(MultivResps),
         FullCovar, Sp) %>% 
  na.omit %>% droplevels ->
  TestDF

N <- length(MultivResps)

Prior <- list(R = list(V = diag(N), nu = 6.002))

Formula <- 
  glue::glue("cbind({paste(MultivResps, collapse = ', ')}) ~ trait - 1 + trait:({paste(FullCovar, collapse = ' + ')})") %>% 
  as.formula

M1 <- 
  MCMCglmm(
    
    Formula,
    data = TestDF,
    prior = Prior,
    rcov = ~us(trait):units,
    family = rep("poisson", N)
    
  )

M1$Sol %>% as.data.frame() %>% dplyr::select(contains("UrbanBinary")) %>% names

M1 %>% MCMCFactorComp2(MultivResps, "UrbanBinary") %>% dim

M1 %>% saveRDS("Output/Multivariate.rds")

# Zoonotic Multivariate ####

MultivResps <- ParasiteTypes[2:7] %>% paste0("ZoonoticRichness.", .)

UrbanDF %>% filter(Sp %in% IMList[[1]]$Data$Sp) %>% 
  select(FocalResp,
         all_of(MultivResps),
         FullCovar, Sp) %>% 
  na.omit %>% droplevels ->
  TestDF

N <- length(MultivResps)

Prior <- list(R = list(V = diag(N), nu = 6.002))

Formula <- 
  glue::glue("cbind({paste(MultivResps, collapse = ', ')}) ~ trait - 1 + trait:({paste(FullCovar, collapse = ' + ')})") %>% 
  as.formula

M1 <- 
  MCMCglmm(
    
    Formula,
    data = TestDF,
    prior = Prior,
    rcov = ~us(trait):units,
    family = rep("poisson", N)
    
  )

M1$Sol %>% as.data.frame() %>% dplyr::select(contains("UrbanBinary")) %>% names

M1 %>% MCMCFactorComp2(MultivResps, "UrbanBinary") %>% dim

M1 %>% saveRDS("Output/ZoonoticMultivariate.rds")

# Path Analysis ####

# UrbanDF %<>% filter(!Pseudoabsence)

Resps <- c("Richness", "ZoonoticRichness")

FamilyList <- rep("nbinomial", 2)

names(FamilyList) <- Resps

FullCovar <- c("hOrder", 
               "LogCites", 
               "DomesticBinary", "HumanDistance",
               "LogArea", "DietDiversity",
               "LogMass", "PC1", "PC2", "UrbanBinary")

AddCovar <- c("HumanDistance", "LogMass", "PC1")

Covar <- FullCovar %>% setdiff(AddCovar)

PhyloList <- 
  IMList <- 
  TestDFList <- 
  list()

ParasiteTypes <- c(
  
  "Overall", "Virus", "Helminth", "Arthropod", "Bacteria", "Protozoa", "Fungus"
  
)

q <- 1

# for(q in q:length(ParasiteTypes[1])){

print(ParasiteTypes[q])

r <- 2

# for(r in r:length(Resps[1])){

FocalResp <- paste0(Resps[r], ".", ParasiteTypes[q])

print(FocalResp)

UrbanDF %>% 
  mutate(LogRichness = log(Richness.Overall + 1)) %>% 
  select(FocalResp,
         X, Y, Eurasia:`NA`,
         FullCovar, Sp, LogRichness) %>% 
  na.omit %>% droplevels ->
  TestDF

TestDF[,"Response"] <- TestDF[,FocalResp]

TestDF %>% 
  group_by(hOrder) %>% 
  summarise(N = n(), 
            Prevalence = Prev(Response), 
            NInf = sum(Response>0)) %>% 
  filter(N>20, Prevalence > 0.001, NInf > 1) %>% 
  pull(hOrder) ->
  
  OrderInclude

TestDF %>% 
  filter(hOrder %in% OrderInclude) %>% 
  droplevels ->
  
  TestDF

CouldBeNumeric <- function(X){
  
  if(is.numeric(X)|is.integer(X)|is.numeric.difftime(X)){
    
    return(TRUE)
    
  }else return(FALSE)
  
}

TestDF %>% 
  select(FullCovar, LogRichness) %>% 
  mutate_if(CouldBeNumeric, ~.x %>% as.numeric %>% scale %>% c) ->
  TestDF[,FullCovar %>% c("LogRichness")]

print(nrow(TestDF))

PhyloMatrix <- 1 - FullSTMatrix[TestDF$Sp, TestDF$Sp] %>% 
  scales::rescale(to = c(10^-6, 1-(10^-6)))

PhyloList[[FocalResp]] <- PhyloMatrix

TestDF$Phylogeny <- match(TestDF$Sp, colnames(PhyloMatrix))

IM1 <- INLAModelAdd(
  Response = FocalResp,
  Explanatory = FullCovar %>% c("LogRichness") %>% setdiff("LogCites"), #[1:8],
  # Add = c("UrbanBinary:hOrder"),# %>% c(FullCovar[9:10]),
  Add = "LogCites",
  # Clashes = list(FullCovar[9:10]),
  Family = FamilyList[[Resps[r]]],
  Data = TestDF,
  Base = T,
  AllModels = T,
  Groups = T, GroupVar = "UrbanBinary", 
  AddSpatial = T
  
)

IM1$FinalModel %>% Efxplot %>% plot

saveRDS(IM1, file = paste0("Output/PathModels_", FocalResp, ".rds"))

list(IM1$FinalModel, IM1$Spatial$Model, IM1$Spatial$SpatiotemporalModel) %>% 
  INLADICFig()

IM1$Spatial$Model %>% ggField(Mesh = IM1$Spatial$Mesh) + 
  scale_fill_discrete_sequential(AlberPalettes[[1]]) +
  geom_sf(data = world,inherit.aes = F, fill = NA, colour = "black", size = 1)

IM1$Spatial$SpatiotemporalModel %>% 
  ggField(Mesh = IM1$Spatial$Mesh,
          GroupLabels = c("Non-urban", "Urban"),
          Groups = 2, GroupVar = "UrbanBinary") + 
  scale_fill_discrete_sequential(AlberPalettes[[1]]) +
  geom_sf(data = world,inherit.aes = F, fill = NA, colour = "black", size = 1)

# Zero-inflated ####

Models_Richness.Overall <- readRDS("C:/Users/gfalb/Documents/Script/UrbanOutputters/Output/Models_Richness.Overall.rds")

ziPrior <- list(R =list(V = diag(2), nu = 0, fix = 2))

F1 <- 
  glue::glue("Richness.Overall ~ trait - 1 + trait:({paste(FullCovar, collapse = ' + ')})") %>% 
  as.formula

mf <- 10

ziModel <- MCMCglmm(F1,
                    rcov =~ idh(trait):units, 
                    family = "zipoisson",
                    prior = ziPrior, 
                    data = Models_Richness.Overall$Data, 
                    # pr = TRUE, 
                    pl = TRUE, 
                    nitt = 13000*mf,
                    thin = 10*mf,burnin = 3000*mf)

ziModel %>% summary
ziModel %>% Efxplot

ziModel %>% saveRDS("Output/ziModel.rds")

HUModel <- MCMCglmm(F1,
                    rcov =~ idh(trait):units, 
                    family = "hupoisson",
                    prior = ziPrior, 
                    data = Models_Richness.Overall$Data, 
                    # pr = TRUE, 
                    pl = TRUE, 
                    nitt = 13000*mf,
                    thin = 10*mf,burnin = 3000*mf)

HUModel %>% summary

HUModel %>% saveRDS("Output/HUModel.rds")

# Zero-inflated Zoonoses ####

ziPrior <- list(R =list(V = diag(2), nu = 0, fix = 2))

F1 <- 
  glue::glue("ZoonoticRichness.Overall ~ trait - 1 + trait:({paste(FullCovar, collapse = ' + ')})") %>% 
  as.formula

mf <- 10

ziModel <- MCMCglmm(F1,
                    rcov =~ idh(trait):units, 
                    family = "zipoisson",
                    prior = ziPrior, 
                    data = Models_ZoonoticRichness.Overall$Data, 
                    # pr = TRUE, 
                    pl = TRUE, 
                    nitt = 13000*mf,
                    thin = 10*mf,burnin = 3000*mf)

ziModel %>% summary

ziModel %>% saveRDS("Output/ziModelZoonotic.rds")

HUModel <- MCMCglmm(F1,
                    rcov =~ idh(trait):units, 
                    family = "hupoisson",
                    prior = ziPrior, 
                    data = Models_ZoonoticRichness.Overall$Data, 
                    # pr = TRUE, 
                    pl = TRUE, 
                    nitt = 13000*mf,
                    thin = 10*mf,burnin = 3000*mf)

HUModel %>% summary

HUModel %>% saveRDS("Output/HUModelZoonotic.rds")
