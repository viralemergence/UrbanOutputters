
# Urban Figures ####

{
  
  library(tidyverse); library(ggregplot); library(cowplot); library(patchwork); library(colorspace)
  library(magrittr); library("rnaturalearth"); library("rnaturalearthdata"); library(INLA)
  library(fs); library(MCMCglmm)
  
  dir_create("Figures")
  
  UrbanColours <- c("#2e74b5", "#28bcae")
  
  theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))
  
  Resps <- c("ZoonoticRichness", "Richness", "PropZoonotic", "NonZoonoticRichness")
  
  ParasiteTypes <- c(
    
    "Overall", "Virus", "Helminth", "Arthropod", "Bacteria", "Protozoa", "Fungus"
    
  )
  
  M1 <- readRDS("Output/Multivariate.rds")
  M2 <- readRDS("Output/ZoonoticMultivariate.rds")
  
  Models_Richness.Overall <- 
    readRDS("Output/Models_Richness.Overall.rds")
  
  Models_ZoonoticRichness.Overall <- 
    readRDS("Output/Models_ZoonoticRichness.Overall.rds")
  
  PathModels_ZoonoticRichness.Overall <- 
    readRDS("Output/PathModels_ZoonoticRichness.Overall.rds")
  
  ContinentModels_ZoonoticRichness.Overall <- readRDS("Output/ContinentModels_ZoonoticRichness.Overall.rds")
  ContinentModels_Richness.Overall <- readRDS("Output/ContinentModels_Richness.Overall.rds")
  
}

# Figure 1 ####

colnames(M2$Sol) <- colnames(M2$Sol) %>% str_remove("Zoonotic")

Models_Richness.Overall$Spatial$Model %>% 
  GetEstimates("UrbanBinary1") %>% 
  paste0(
    "; P = ",
    Models_Richness.Overall$Spatial$Model %>% 
      INLAPValue("UrbanBinary1") %>% unlist %>% round(3)
  ) -> Label

YMax <- Models_Richness.Overall$Data[,paste0("Richness.", "Overall")] %>% max

Models_ZoonoticRichness.Overall$Spatial$Model %>% 
  GetEstimates("UrbanBinary1") %>% 
  paste0(
    "; P = ",
    Models_ZoonoticRichness.Overall$Spatial$Model %>% 
      INLAPValue("UrbanBinary1") %>% unlist %>% round(4)
  ) -> Label2

YMax2 <- Models_ZoonoticRichness.Overall$Data[,paste0("ZoonoticRichness.", "Overall")] %>% max

(TopRow <- Models_Richness.Overall$Data %>%
    mutate_at(paste0("Richness.", "Overall"), ~.x + 1) %>% 
    mutate_at("UrbanBinary", ~as.factor(c("Non-Urban", "Urban")[as.numeric(as.character(.x))+1])) %>%
    SinaGraph("UrbanBinary", paste0("Richness.", "Overall"),
              Scale = "width", Alpha = 0.6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Overall parasite richness") +
    # geom_text(data = NULL, aes(x = 1.5, y = YMax + YMax*0.1, label = Label),
    #           colour = UrbanColours[[1]]) +
    scale_y_log10(breaks = c(1, 11, 101, 1001), 
                  labels = c("0", "10", "100", "1000"),
                  limits = c(1, 600)) +
    scale_colour_manual(values = UrbanColours) +
    ggtitle("Overall parasite richness") +
    
    Models_ZoonoticRichness.Overall$Data %>%
    mutate_at(paste0("ZoonoticRichness.", "Overall"), ~.x + 1) %>% 
    mutate_at("UrbanBinary", ~as.factor(c("Non-Urban", "Urban")[as.numeric(as.character(.x))+1])) %>%
    SinaGraph("UrbanBinary", paste0("ZoonoticRichness.", "Overall"),
              Scale = "width", Alpha = 0.6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Zoonotic richness") +
    # geom_text(data = NULL, aes(x = 1.5, y = YMax2 + YMax2*0.1, label = Label2),
    #           colour = UrbanColours[[1]]) +
    scale_y_log10(breaks = c(1, 11, 101, 1001), 
                  labels = c("0", "10", "100", "1000"),
                  limits = c(1, 600)) +
    scale_colour_manual(values = UrbanColours) +
    ggtitle("Zoonotic parasite richness")) 

EstimateList <- 
  
  list(
    
    Models_Richness.Overall$AllModels[[1]],
    Models_Richness.Overall$FinalModel,
    Models_Richness.Overall$Spatial$Model,
    
    Models_ZoonoticRichness.Overall$AllModels[[1]],
    Models_ZoonoticRichness.Overall$FinalModel,
    Models_ZoonoticRichness.Overall$Spatial$Model
    
    
  )

names(EstimateList) <- c("Base", "+ Citations", "+ Spatial") %>% 
  rep(2)

EstimateDF <- 
  EstimateList %>% 
  map(~INLAPValue(.x, "UrbanBinary1", Method = "MCMC") %>% unlist) %>% 
  as.data.frame %>% 
  gather("Model", "PValue") %>% 
  mutate_at("Model", ~str_replace(.x, "X..", "+ ")) %>% 
  mutate_at("Model", ~str_remove(.x, ".1")) %>% 
  mutate(Response = rep(c("Overall", "Zoonotic"), each = 3))

EstimateDF %<>% 
  mutate_at("PValue", ~round(.x, 3)) %>% 
  mutate(Asterisks = ifelse(PValue < 0.001, "***", 
                            ifelse(PValue < 0.01, "**", 
                                   ifelse(PValue < 0.05, "*", "")))) %>% 
  # mutate_at("PValue", ~str_replace(.x, "^0$", "P < 0.000001")) %>% 
  # mutate_at("PValue", ~str_replace(.x, "^0.0", "P = 0.0")) %>% 
  mutate_at("PValue", ~str_replace(.x, "^0$", "< 0.000001")) %>% 
  mutate_at("PValue", ~paste0(.x, Asterisks))

EstimateDF <- 
  EstimateList %>% 
  map(c("summary.fixed", as.data.frame, rownames_to_column)) %>% 
  bind_rows(.id = "Model") %>% 
  filter(rowname == "UrbanBinary1") %>% 
  rename(Estimate = 3, Lower = 5, Upper = 7) %>% 
  mutate(Response = rep(c("Overall", "Zoonotic"), each = 3)) %>% 
  left_join(EstimateDF, .)

EstimateDF %<>% 
  mutate_at("Model", ~factor(.x, levels = c("Base", "+ Citations", "+ Spatial")))

Estimates <- 
  EstimateDF$Response %>% unique %>% 
  map(function(a){
    
    EstimateDF %>% filter(Response == a) %>% 
      ggplot(aes(Model, Estimate, colour = Model)) +#
      geom_hline(aes(yintercept = 0), lty = 2, alpha = 0.3) +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                    width = 0.2) +
      geom_point(size = 5) +
      geom_text(aes(y = Upper + 0.1, label = PValue), 
                # hjust = 0, 
                show.legend = F) +
      scale_colour_manual(values = c(UrbanColours[[2]], 
                                     "black",
                                     UrbanColours[[1]])) +
      # scale_x_discrete(limits  = rev(c("Base", "+ Citations", "+ Spatial"))) +
      # lims(y = c(NA, 2)) +
      labs(y = "Urban effect estimate", 
           x = NULL) +
      theme(legend.position = "none") +
      # coord_flip() +
      # facet_grid(~Response) +
      NULL
    
  }) %>% ArrangeCowplot()

(TopRow  / 
    Estimates) +
  plot_annotation(tag_levels = "A")

# (M1 %>% list(M2) %>% 
#    Efxplot(ModelNames = c("Overall", "Zoonotic"),
#            VarOrder = M1$Sol %>% as.data.frame() %>% dplyr::select(contains("UrbanBinary")) %>% names %>% rev,
#            VarNames = M1$Sol %>% as.data.frame() %>% dplyr::select(contains("UrbanBinary")) %>% names %>% 
#              str_remove("traitRichness.") %>% 
#              str_remove(":UrbanBinary1") %>% rev,
#            PointOutline = T) +
#    lims(y = c(-0.75, 1.75)) +
#    scale_y_continuous(breaks = c(-2):2/2) +
#    scale_colour_manual(values = c(UrbanColours[[1]], UrbanColours[[2]])) +
#    labs(x = NULL, y = "Estimate for urban effect")) +
# plot_layout(heights = c(2, 1)) +
# plot_annotation(tag_levels = "A")

ggsave("Figures/Figure1.jpeg", units = "mm",
       width = 220, height = 180)

ggsave("Figures/Figure1.pdf", units = "mm",
       width = 220, height = 180)

# Figure 2 ####

(UrbanDF %>% 
   mutate_at("UrbanBinary", ~as.numeric(.x) - 1) %>% 
   ggplot(aes(LogCites, UrbanBinary)) + 
   geom_point(alpha = 0.1) + 
   geom_smooth(method = glm, 
               method.args = list(family = "binomial"), 
               colour = UrbanColours[[2]],
               fill = UrbanColours[[2]], alpha = 0.2) +
   # lims(y = c(NA, 1.2)) +
   labs(y = "Urban      ", 
        x = "log(Citations + 1)") +
   scale_y_continuous(breaks = c(0:5/5), limits = c(NA, 1.2)) +
   ggpubr::stat_cor(method = "spearman") +
   
   UrbanDF %>% 
   ggplot(aes(LogCites, log(Richness.Overall + 1))) + 
   geom_point(alpha = 0.1, colour = UrbanColours[[1]]) + 
   # geom_smooth(method = lm, 
   #             # method.args = list(family = "poisson"), 
   #             colour = UrbanColours[[2]],
   #             fill = UrbanColours[[2]], alpha = 0.2) +
   labs(y = "log(Parasite richness + 1)", 
        x = "log(Citations + 1)") +
   lims(y = c(NA, 7.5)) +
   ggpubr::stat_cor(method = "spearman") +
   
   UrbanDF %>% 
   ggplot(aes(LogCites,
              #log(Richness.Overall + 1), 
              log(ZoonoticRichness.Overall + 1))) + 
   geom_point(alpha = 0.1, colour = UrbanColours[[1]]) + 
   # geom_smooth(method = lm, 
   #             # method.args = list(family = "poisson"), 
   #             colour = UrbanColours[[2]],
   #             fill = UrbanColours[[2]], alpha = 0.2) +
   labs(y = "log(Zoonotic richness + 1)", 
        x = "log(Citations + 1)")  +
   lims(y = c(NA, 7.5)) +
   ggpubr::stat_cor(method = "spearman")) +
  plot_annotation(tag_levels = "A")

ggsave("Figures/Figure2.jpeg", units = "mm", 
       height = 100, width = 240, dpi = 300)

ggsave("Figures/Figure2.pdf", units = "mm", 
       height = 100, width = 240, dpi = 300)

# Figure 2b ####

(
  UrbanDF %>%
    mutate_at(paste0("Citations"), ~.x + 1) %>% 
    mutate_at("UrbanBinary", ~as.factor(c("Non-Urban", "Urban")[as.numeric(as.character(.x))+1])) %>%
    SinaGraph("UrbanBinary", "Citations",
              Scale = "width", Alpha = 0.6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Citations") +
    # geom_text(data = NULL, aes(x = 1.5, y = YMax + YMax*0.1, label = Label),
    #           colour = UrbanColours[[1]]) +
    scale_y_log10(breaks = c(1, 11, 101, 1001, 10001), 
                  labels = c("0", "10", "100", "1000", "10000"),
                  limits = c(1, 50000)) +
    scale_colour_manual(values = UrbanColours) +
    
    UrbanDF %>% 
    ggplot(aes(LogCites, log(Richness.Overall + 1))) + 
    geom_point(alpha = 0.1, colour = UrbanColours[[1]]) + 
    # geom_smooth(method = lm, 
    #             # method.args = list(family = "poisson"), 
    #             colour = UrbanColours[[2]],
    #             fill = UrbanColours[[2]], alpha = 0.2) +
    labs(y = "log(Parasite richness + 1)", 
         x = "log(Citations + 1)") +
    lims(y = c(NA, 7.5)) +
    ggpubr::stat_cor(method = "spearman") +
    
    UrbanDF %>% 
    ggplot(aes(LogCites,
               #log(Richness.Overall + 1), 
               log(ZoonoticRichness.Overall + 1))) + 
    geom_point(alpha = 0.1, colour = UrbanColours[[1]]) + 
    # geom_smooth(method = lm, 
    #             # method.args = list(family = "poisson"), 
    #             colour = UrbanColours[[2]],
    #             fill = UrbanColours[[2]], alpha = 0.2) +
    labs(y = "log(Zoonotic richness + 1)", 
         x = "log(Citations + 1)")  +
    lims(y = c(NA, 7.5)) +
    ggpubr::stat_cor(method = "spearman")) +
  plot_annotation(tag_levels = "A")

ggsave("Figures/Figure2b.jpeg", units = "mm", 
       height = 100, width = 240, dpi = 300)

ggsave("Figures/Figure2b.pdf", units = "mm", 
       height = 100, width = 240, dpi = 300)

# Figure 3 ####

NDraws <- 1000

Path1 <- Models_Richness.Overall$Spatial$Model %>% GetEstimates("UrbanBinary1", Draw = T, NDraws = NDraws)
Path2 <- PathModels_ZoonoticRichness.Overall$Spatial$Model %>% GetEstimates("LogRichness", Draw = T, NDraws = NDraws)
Path3 <- PathModels_ZoonoticRichness.Overall$Spatial$Model %>% GetEstimates("UrbanBinary1", Draw = T, NDraws = NDraws)

PathOutputs <- data.frame(
  
  Urban_Richness = Path1,
  Richness_ZoonoticRichness = Path2,
  Urban_ZoonoticRichness = Path3
  
) %>% mutate(Indirect = Path1*Path2)

MCMCSummary <- function(X){
  
  data.frame(
    Mean = mean(X),
    Lower = HPDinterval(as.mcmc(X))[1],
    Upper = HPDinterval(as.mcmc(X))[2],
    P = PCalc(X)
    
  ) %>% list
}

PathOutputs %>% 
  summarise_all(~.x %>% MCMCSummary %>% Unlist1 %>% Unlist1 %>% round(3) %>% list) %>% 
  # map(bind_rows) %>% 
  map(c(bind_cols, rownames_to_column)) %>% 
  bind_rows(.id = "Model") %>% 
  rename(Parameter = 2, Estimate = 3) %>% 
  mutate_at("Parameter", ~c("Mean", "Lower", "Upper", "P")[as.numeric(.x)]) %>% 
  pivot_wider(names_from = "Parameter", values_from = "Estimate") -> 
  PathEstimates

PathEstimates %>% write.csv("PathCoefficients.csv")

# Figure 4 ####

Models_Richness.Overall$FinalModel$summary.fixed %>% 
  rownames %>% setdiff("(Intercept)") ->
  Vars

Vars <- 
  Vars[!str_detect(Vars, "hOrder")]

Vars %>% str_remove("hOrder|Binary1|Log|Strength") ->
  VarLabels

EffectNames <- c("Overall", #"+ Spatial", 
                 "Zoonotic", "+ Spatial ")

(ModelEffects <- 
    list(#Models_Richness.Overall$FinalModel,
      Models_Richness.Overall$Spatial$Model,
      # Models_ZoonoticRichness.Overall$FinalModel, 
      Models_ZoonoticRichness.Overall$Spatial$Model) %>% 
    Efxplot(PointOutline = T, Intercept = F, 
            ModelNames = EffectNames, VarOrder = rev(Vars), 
            VarNames = rev(VarLabels)) +
    theme(legend.position = c(0.15, 0.15), legend.justification = c(0, 0)) +
    # scale_colour_discrete_sequential(AlberPalettes[[1]], nmax = 6, order = c(2,5)) +
    scale_colour_manual(values = c(UrbanColours[[1]], UrbanColours[[2]])))

# Spatial

world <- ne_countries(scale = "medium", returnclass = "sf")

Map <- 
  ggField(Models_Richness.Overall$Spatial$Model,
          Models_Richness.Overall$Spatial$Mesh) +  
  scale_fill_discrete_sequential(AlberPalettes[[3]]) +
  geom_sf(data = world,inherit.aes = F, fill = NA, colour = "black", size = 1) +
  #labs(x = "Longitude", y = "Latitude", fill = "Richness") +
  labs(x = NULL, y = NULL, fill = "Richness") +
  coord_sf(xlim = c(min(UrbanDF$X, na.rm = T), max(UrbanDF$X, na.rm = T)), 
           ylim = c(min(UrbanDF$Y, na.rm = T), max(UrbanDF$Y, na.rm = T) + 20))

# Continents 

UrbanDF %>% 
  select(Eurasia:`NAm`) %>% names -> ContinentVar

ContinentModels_Richness.Overall$FinalModel$summary.fixed %>% 
  rownames %>% setdiff("(Intercept)") ->
  Vars2

Vars2 %>% 
  str_replace("NAm", "N.Am") %>% 
  str_replace("SAm", "S.Am") %>% 
  str_replace("AbsLat", "Latitude") %>% 
  str_remove("hOrder|Binary1|Log|Strength") ->
  VarLabels2

EffectNames <- c("Overall", "+ Spatial", "Zoonotic", "+ Spatial ")

ModelEffects2 <- 
  list(ContinentModels_Richness.Overall$FinalModel,
       # ContinentModels_Richness.Overall$Spatial$Model,
       ContinentModels_ZoonoticRichness.Overall$FinalModel) %>%  
  # ContinentModels_ZoonoticRichness.Overall$Spatial$Model) %>%
  Efxplot(PointOutline = T, Intercept = F, ModelNames = EffectNames[c(1,3)],
          VarOrder = ContinentVar %>% c(., paste0("UrbanBinary1:", .)) %>% 
            # Vars2 %>% 
            rev %>% c(., "AbsLat"),
          VarNames = ContinentVar %>% c(., paste0("Urban:", .)) %>% 
            # VarLabels2 %>% 
            rev %>% c(., "AbsLat") %>% 
            str_replace("NAm", "N.Am") %>% 
            str_replace("SAm", "S.Am") %>% 
            str_replace("AbsLat", "Latitude")) +
  theme(legend.position = "none") +
  # theme(legend.position = c(0.15, 0.15), legend.justification = c(0, 0)) +
  # scale_colour_discrete_sequential(AlberPalettes[[1]], nmax = 6, order = c(2, 5)) +
  scale_colour_manual(values = c(UrbanColours[[1]], UrbanColours[[2]])) +
  # scale_colour_discrete_sequential(AlberPalettes[[1]]) +
  scale_y_continuous(breaks = c((-1):1))

(((ModelEffects | 
     ModelEffects2) +
    plot_layout(guides = "collect")) / Map) +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "A") 

ggsave("Figures/Figure4.jpeg", units = "mm", 
       height = 200, width = 200, dpi = 300)

ggsave("Figures/Figure4.pdf", units = "mm", 
       height = 200, width = 200, dpi = 300)

# ~~~~~ Supplement ####

# Figure SI1-2 ####

# LabelList <- YMaxList <- list()

# ParasiteTypes %>% map(function(a){
#   
#   print(a)
#   
#   Models_Richness.Overall <-
#     readRDS(paste0("Output/Models_Richness.", a, ".rds"))
# 
#   UrbanEstimates <- Models_Richness.Overall$FinalModel %>%
#     GetEstimates("UrbanBinary1", Round = 2)
#   
#   if(!UrbanEstimates == "NA (NA, NA)"){
#     
#     UrbanEstimates %>% 
#       paste0(
#         
#         "; P = ",
#         Models_Richness.Overall$FinalModel %>% 
#           INLAPValue("UrbanBinary1") %>% unlist %>% round(4)
#         
#       ) ->> LabelList[[a]]
#     
#   }else{
#     
#     LabelList[[a]] <<- "NS"
#     
#   }
#   
#   YMaxList[[a]] <<- 
#     Models_Richness.Overall$Data[,paste0("Richness.", a)] %>% max
#   
# })

ParasiteTypes[2:7] %>% map(function(a){
  
  
  Models_Richness.Overall <-
    readRDS(paste0("Output/Models_Richness.", a, ".rds"))
  
  Models_Richness.Overall$Data %>%
    mutate_at("UrbanBinary", as.factor) %>%
    SinaGraph("UrbanBinary", paste0("Richness.", a),
              Scale = "width", Alpha = 0.3) +
    theme(legend.position = "none") +
    labs(x = "Urban", y = paste0(a)) +
    # geom_text(data = NULL, 
    #           aes(x = 1.5, y = YMaxList[[a]] + YMaxList[[a]]*0.1, 
    #               label = LabelList[[a]]),
    #           colour = UrbanColours[[1]]) +
    scale_colour_manual(values = UrbanColours)
  
}) %>% 
  ArrangeCowplot()

ggsave(filename = "Figures/RichnessSina.jpeg",
       units = "mm", height = 150, width = 250)

# LabelList <- YMaxList <- list()

# ParasiteTypes[2:7-1] %>% map(function(a){
#   
#   print(a)
#   
#   Models_Richness.Overall <-
#     readRDS(paste0("C:/Users/gfalb/Documents/Script/UrbanOutputters/Output/Models_ZoonoticRichness.", a, ".rds"))
#   
#   UrbanEstimates <- Models_Richness.Overall$FinalModel %>% 
#     GetEstimates("UrbanBinary1", Round = 2)
#   
#   if(!UrbanEstimates == "NA (NA, NA)"){
#     
#     UrbanEstimates %>% 
#       paste0(
#         
#         "; P = ",
#         Models_Richness.Overall$FinalModel %>% 
#           INLAPValue("UrbanBinary1") %>% unlist %>% round(4)
#         
#       ) ->> LabelList[[a]]
#     
#   }else{
#     
#     LabelList[[a]] <<- "NS"
#     
#   }
#   
#   YMaxList[[a]] <<- 
#     Models_Richness.Overall$Data[,paste0("ZoonoticRichness.", a)] %>% max
#   
# })

ParasiteTypes[2:7] %>% map(function(a){
  
  # Models_Richness.Overall <-
  #   readRDS(paste0("Output/Models_ZoonoticRichness.", a, ".rds"))
  
  # Models_Richness.Overall$Data %>%
  UrbanDF %>% 
    mutate_at("UrbanBinary", as.factor) %>%
    SinaGraph("UrbanBinary", paste0("ZoonoticRichness.", a),
              Scale = "width", Alpha = 0.3) +
    theme(legend.position = "none") +
    labs(x = "Urban", y = paste0(a)) +
    # geom_text(data = NULL, 
    #           aes(x = 1.5, y = YMaxList[[a]] + YMaxList[[a]]*0.1, 
    #               label = LabelList[[a]]),
    #           colour = UrbanColours[[1]]) +
    scale_colour_manual(values = UrbanColours)
  
}) %>%
  ArrangeCowplot()

ggsave(filename = "Figures/ZoonoticRichnessSina.jpeg",
       units = "mm", height = 150, width = 250)

# Figure SI3 ####

PathModels_ZoonoticRichness.Overall$FinalModel %>% 
  list(., PathModels_ZoonoticRichness.Overall$Spatial$Model) %>% 
  Efxplot(ModelNames = c("Base", "SPDE"))

ggsave("Figures/SIFigure3.jpeg", 
       units = "mm", height = 150, width = 200)

# Figure SI4 ####

"Output" %>% 
  list.files(full.names = T, pattern = "^Models_Richness") %>% 
  map(readRDS) %>% 
  map("FinalModel") -> ModelList1

"Output" %>% 
  list.files(full.names = T, pattern = "^Models_Richness") %>% 
  map(readRDS) %>% 
  map(c("Spatial", "Model")) -> ModelList1

names(ModelList1) <- 
  "Output" %>% list.files(pattern = "^Models_Richness") %>% 
  str_remove_all("Models_Richness.|.rds")

ModelList1 <-
  ModelList1[c("Overall", ParasiteTypes)]

ModelList1[c("Overall", ParasiteTypes[2:6])] %>% 
  Efxplot(ModelNames = c("Overall", ParasiteTypes[2:6]))

ggsave("Figures/SIFigure4.jpeg", 
       units = "mm", height = 150, width = 200)

# Figure SI5 ####

# "Output" %>% list.files(full.names = T, pattern = "^Models_ZoonoticRichness") %>% 
#   map(readRDS) %>% map("FinalModel") -> ModelList

"Output" %>% list.files(full.names = T, pattern = "^Models_ZoonoticRichness") %>% 
  map(readRDS) %>% map(c("Spatial", "Model")) -> ModelList

names(ModelList) <- 
  "Output" %>% list.files(pattern = "^Models_ZoonoticRichness") %>% 
  str_remove_all("Models_ZoonoticRichness.|.rds")

ModelList <-
  ModelList[c("Overall", ParasiteTypes)]

ModelList[c("Overall", ParasiteTypes[1:6])] %>% 
  Efxplot(ModelNames = names(ModelList)) 

ggsave("Figures/SIFigure5.jpeg", 
       units = "mm", height = 150, width = 200)

# Figure SI6 ####

"Output" %>% list.files(pattern = "ziModel", full.names = T) %>% map(readRDS) ->
  XIModels

XIModels %>% lapply(function(a){
  
  a %>% summary %>% extract2("solutions") %>% 
    as.data.frame %>% rownames_to_column("Var") %>% 
    separate(Var, ":", into = c("Resp", "Explanatory")) %>% 
    mutate_at("Resp", ~c("Count", "Inflation")[as.numeric(str_detect(.x, "_")) + 1])
  
}) %>% bind_rows(.id = "Model") %>% 
  mutate_at("Explanatory", ~ifelse(is.na(.x), "Intercept", .x)) %>% 
  mutate_at("Explanatory", ~factor(.x, levels = rev(unique(.x)))) %>% 
  mutate_at("Model", ~c("Overall", "Zoonotic")[as.numeric(.x)]) %>% 
  rename(Mean = 4, Lower = 5, Upper = 6) %>% 
  filter(!str_detect(Explanatory, "hOrder")) %>% 
  ggplot(aes(Explanatory, Mean, shape = Resp, group = Resp, colour = Resp)) + 
  geom_hline(yintercept = 0, lty = 2, alpha = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3, position = position_dodge(w = 0.5)) + 
  geom_point(width = 0.3, size = 2.5, colour = "black", position = position_dodge(w = 0.5)) + 
  geom_point(width = 0.3, position = position_dodge(w = 0.5)) + 
  scale_colour_manual(values = c(UrbanColours[[1]], UrbanColours[[2]])) +
  labs(y = "Effect estimate") +
  facet_wrap(Resp~Model, scales = "free_x") + coord_flip()

ggsave("Figures/SIFigure6.jpeg", units = "mm",
       height = 200, width = 200)

# Multivariate

Multivariate <- readRDS("Output/Multivariate.rds")

MultivResps <- ParasiteTypes[2:7] %>% paste0("Richness.", .)

MultivOutput <- Multivariate %>% 
  MCMCFactorComp2(MultivResps, "UrbanBinary1")

MultivOutput %>% rownames %>% 
  str_remove("traitRichness.") %>% 
  str_remove(":UrbanBinary1") ->
  rownames(MultivOutput) -> colnames(MultivOutput)

colnames(MultivOutput)[1] <- "Intercept"

DimNames <- dimnames(MultivOutput)

data.frame(X = 0, Y = 4:20) %>%
  rowwise %>% do(x = rep(.$X, .$Y) %>% paste0(collapse = "") %>% paste0("1")) %>% 
  unnest(x) %>% pull(x) %>% paste0(collapse = "|") %>% 
  str_remove(MultivOutput, .) %>% 
  matrix(ncol = ncol(MultivOutput), nrow = nrow(MultivOutput)) ->
  MultivOutput

dimnames(MultivOutput) <- DimNames

# BES Figures ####


IMList$Richness.Overall$AllModels[[1]] %>% list %>% append(list(IMList$Richness.Overall$FinalModel)) %>% 
  Efxplot

IMList$Richness.Overall$FinalModel %>% INLAPValue("UrbanBinary1")
IMList$ZoonoticRichness.Overall$FinalModel %>% INLAPValue("UrbanBinary1")

IMList$Richness.Overall$Spatial$Model %>% INLAPValue("UrbanBinary1")
IMList$ZoonoticRichness.Overall$Spatial$Model %>% INLAPValue("UrbanBinary1")

IMList$Richness.Overall$AllModels[[1]] %>% list %>% append(list(IMList$Richness.Overall$Spatial$Model)) %>% 
  Efxplot(ModelNames = c("Base", "With Citations")) + 
  scale_x_discrete(limits = c("UrbanBinary1", "LogCites")) +
  scale_colour_manual(values = c(UrbanColours[1], UrbanColours[2])) +
  scale_y_continuous(limits = c(-0.25, 2.5))

ggsave("CitationsFigure.jpeg", units = "mm", height = 120, width = 150)


world <- ne_countries(scale = "medium", returnclass = "sf")

# Map <- 

IMList$Richness.Overall$Spatial$Model %>% 
  ggField(IMList$Richness.Overall$Spatial$Mesh) +
  # coord_sf()  +
  scale_fill_discrete_sequential(AlberPalettes[[3]]) +
  geom_sf(data = world,inherit.aes = F, fill = NA, colour = "black", size = 1) +
  #labs(x = "Longitude", y = "Latitude", fill = "Richness") +
  labs(x = NULL, y = NULL, fill = "Richness") +
  coord_sf(xlim = c(min(UrbanDF$X, na.rm = T), max(UrbanDF$X, na.rm = T)), 
           ylim = c(min(UrbanDF$Y, na.rm = T), max(UrbanDF$Y, na.rm = T) + 20))

ggsave("Map.jpeg", units = "mm", height = 120, width = 180)

IMList$Richness.Overall$Spatial$Model %>% Efxplot

TestDF %>% 
  ggplot(aes(X, Y)) +
  geom_sf(data = world,inherit.aes = F, fill = NA, colour = "black", size = 1) +
  #labs(x = "Longitude", y = "Latitude", fill = "Richness") +
  labs(x = NULL, y = NULL, fill = "Richness") +
  coord_sf(xlim = c(min(UrbanDF$X, na.rm = T), max(UrbanDF$X, na.rm = T)), 
           ylim = c(min(UrbanDF$Y, na.rm = T), max(UrbanDF$Y, na.rm = T) + 20)) +
  geom_point(colour = UrbanColours[[1]], alpha = 0.3)

ggsave("SampleMap.jpeg", units = "mm", height = 120, width = 180)
