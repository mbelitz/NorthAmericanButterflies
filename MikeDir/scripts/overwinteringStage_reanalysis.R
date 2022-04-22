#load libraries
library(lme4)
library(lmerTest)
library(car)
library(tidyverse)
library(MuMIn)
library(ape)
library(stringr)
library(sjPlot)

bfly_trends = read.table("data/butterfly_traits_envars_trends_50km_NoMergeJul_m5_trim_6traits.txt",sep='\t',as.is=T,check.names=F,header=T)

# build first model

m <- lmer(Abundance.trend ~ Cropland.trend + Built.2005_2015 + Precip.1993_2018 +
          Temp.1993_2018 + LarvalColor + LarvalHair + AdultSize + AdultColor +  
            Diet.breadth.families + (1 | grid_id) + (1 | Species),
          data = bfly_trends, REML = F, lmerControl(optimizer = "bobyqa"))

ms <- step(m)
ms

## so our results are similar to crossley et all, which is assuring  
cross_top <- lmer(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                (1 | grid_id) + (1 | Species),
              data = bfly_trends, REML = F, lmerControl(optimizer = "bobyqa"))

summary(cross_top)
r.squaredGLMM(cross_top) ## 0.00657 R2m -- 0.588 R2c

# first we shall look at butterflies across north america results with traits we think are important
#read in traits
st <- read.csv("MikeDir/data/TraitsButterflyCrossley.csv")
# join to butterfly trend data
bfly_trends_st <- left_join(bfly_trends, st)

# do any species not have OWS?
dplyr::filter(bfly_trends_st, is.na(overwintering.stage))
# voltinism?
dplyr::filter(bfly_trends_st, is.na(voltinism))
# we have traits for every species
# what sort of unique traits do we have?
unique(bfly_trends_st$voltinism) #F for facultative, "M" for multivoltine, and "U" for univoltine "" needs to be removed
unique(bfly_trends_st$overwintering.stage) # so lots of stuff, going to make "A (migrate)" into "M" and then 

# great, we have traits for all species so let's make into mdf
mdf <- bfly_trends_st %>% 
  mutate(overwintering.stage = if_else(condition = overwintering.stage == "A (migrate)",
                                       true = "M",
                                       false = overwintering.stage)) %>% 
  filter(overwintering.stage %in% c("L", "A", "P", "M", "E")) %>% 
  filter(voltinism != "")
# so by removing species with more ambiguous or complicated overwintering stages, 
# we now have 9360 rows from the original dataset that had 9535 rows

# run linear model to get top model
m <- lmer(Abundance.trend ~ Cropland.trend + Built.2005_2015 + Precip.1993_2018 +
            Temp.1993_2018 + LarvalColor + LarvalHair + AdultSize + AdultColor +  
            Diet.breadth.families + overwintering.stage + voltinism + 
            Temp.1993_2018:overwintering.stage + 
            Precip.1993_2018:overwintering.stage +
            (1 | grid_id) + (1 | Species),
          data = mdf, REML = F, lmerControl(optimizer = "bobyqa"))

ms <- step(m)
ms

m_top <- lmer(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                overwintering.stage + 
                Temp.1993_2018:overwintering.stage +
                Precip.1993_2018:overwintering.stage +
                (1 | grid_id) + (1 | Species),
              data = mdf, REML = F, lmerControl(optimizer = "bobyqa"))


summary(m_top)                
vif(m_top)                
r.squaredGLMM(m_top) ## 0.0379 # 0.547 --> greatly improving the variation explained in model

# plot results
plot_model(m_top, type = "pred", terms = "overwintering.stage")
plot_model(m_top, type = "pred", terms = c("Precip.1993_2018", "overwintering.stage"))
plot_model(m_top, type = "pred", terms = c("Temp.1993_2018", "overwintering.stage"))

## now take that top model and make it a phylogenetic linear mixed model
# read in phylogeny
t <- read.tree(file = "MikeDir/data/iScience/SupDryad_treepl.tre")
length(t$tip.label)
head(t$tip.label)
# read in tiplabels that will match Spp of mdf
tip_labels_df <- read.csv("MikeDir/data/iScience/tipLabel.csv") # all tip labels
head(tip_labels_df$validName)

t$tip.label <- tip_labels_df$validName
head(tip_labels_df$validName)
head(t$tip.label)

# drop species not in phylogeny
mdf_phylo <- mdf %>% 
  filter(Species %in% t$tip.label)

# now we drop our data from 9535 rows of Crossley to 9092 
# total of 4.6% loss

## now add in a phylogenetic effect
library(phyr)

mdf_phylo <- mdf_phylo %>% 
  mutate(Precip.1993_2018 = scale(Precip.1993_2018),
         Temp.1993_2018 = scale(Temp.1993_2018))

tree_sp <- t$tip.label
sppNotInAnalysis <- data.frame(Species = tree_sp) %>% 
  filter(!Species %in% mdf_phylo$Species)

t2 <- ape::drop.tip(t, tip = sppNotInAnalysis$Species)

pglmm <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                 overwintering.stage + 
                 Temp.1993_2018:overwintering.stage +
                 Precip.1993_2018:overwintering.stage +
                 (1 | grid_id) + (1 | Species__) ,
               data = mdf_phylo,
               cov_ranef = list(Species = t2), 
               bayes = T, verbose = T)

# look into model summary
summary(pglmm)
# look at posterior distributions
p <- phyr::plot_bayes(x = pglmm)
ggsave(plot = p, filename = "MikeDir/figures/pglmm_effects_entireExtent.png",
       dpi = 400, width = 12, height = 8)
# now plot interactions
pred_vals <- -3:3

inla_mdf <- mdf_phylo %>% 
  dplyr::select(Abundance.trend, Temp.1993_2018, 
                Precip.1993_2018, overwintering.stage,
                Species, grid_id)

term1 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Temp.1993_2018 = pred_vals,
  Precip.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("E", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))

term2 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Temp.1993_2018 = pred_vals,
  Precip.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("L", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))

term3 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Temp.1993_2018 = pred_vals,
  Precip.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("P", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))

term4 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Temp.1993_2018 = pred_vals,
  Precip.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("A", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))

term5 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Temp.1993_2018 = pred_vals,
  Precip.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("M", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))


pred.df.1 <- rbind(inla_mdf, term1)
pred.df.2 <- rbind(inla_mdf, term2)
pred.df.3 <- rbind(inla_mdf, term3)
pred.df.4 <- rbind(inla_mdf, term4)
pred.df.5 <- rbind(inla_mdf, term5)

pglmm1 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.1,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

pglmm2 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.2,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

pglmm3 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.3,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

pglmm4 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.4,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

pglmm5 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.5,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

## make df of results
rdf1 <- pglmm1$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Egg")

rdf2 <- pglmm2$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Larvae")

rdf3 <- pglmm3$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Pupae")

rdf4 <- pglmm4$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Adult")

rdf5 <- pglmm5$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Migratory")

rdf_total <- rbind(rdf1, rdf2, rdf3, rdf4, rdf5)

temp_ows <- ggplot(rdf_total, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2), size = 1.05) +
  labs(x = "Temperature", y = "Abundance trend", 
       fill = "Overwintering stage", color = "Overwintering stage") + 
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_classic()

temp_ows
ggsave(plot = temp_ows, filename = "MikeDir/figures/temp_OWS_interaction_fullExtent.png",
       dpi = 450, width = 6, height = 3.5)



## same thing but do this for the precipitation x ows interaction
# now plot interactions
term1 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Precip.1993_2018 = pred_vals,
  Temp.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("E", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))

term2 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Precip.1993_2018 = pred_vals,
  Temp.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("L", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))

term3 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Precip.1993_2018 = pred_vals,
  Temp.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("P", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))

term4 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Precip.1993_2018 = pred_vals,
  Temp.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("A", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))

term5 <- data.frame(
  Abundance.trend = rep(NA, length(pred_vals)),
  Precip.1993_2018 = pred_vals,
  Temp.1993_2018 = rep(NA, length(pred_vals)),
  overwintering.stage = rep("M", length(pred_vals)),
  grid_id = rep(NA, length(pred_vals)),
  Species = rep(NA, length(pred_vals)))


pred.df.1 <- rbind(inla_mdf, term1)
pred.df.2 <- rbind(inla_mdf, term2)
pred.df.3 <- rbind(inla_mdf, term3)
pred.df.4 <- rbind(inla_mdf, term4)
pred.df.5 <- rbind(inla_mdf, term5)

pglmm1 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.1,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

pglmm2 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.2,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

pglmm3 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.3,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

pglmm4 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.4,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

pglmm5 <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                  overwintering.stage + 
                  Temp.1993_2018:overwintering.stage +
                  Precip.1993_2018:overwintering.stage +
                  (1 | grid_id) + (1 | Species__) ,
                data = pred.df.5,
                cov_ranef = list(Species = t2), 
                bayes = T, verbose = T)

## make df of results
rdf1 <- pglmm1$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Egg")

rdf2 <- pglmm2$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Larvae")

rdf3 <- pglmm3$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Pupae")

rdf4 <- pglmm4$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Adult")

rdf5 <- pglmm5$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Migratory")

rdf_total <- rbind(rdf1, rdf2, rdf3, rdf4, rdf5)

prec_ows <- ggplot(rdf_total, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2), size = 1.05) +
  labs(x = "Precipitation", y = "Abundance trend", 
       fill = "Overwintering stage", color = "Overwintering stage") + 
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_classic()

prec_ows


ggsave(plot = prec_ows, filename = "MikeDir/figures/prec_OWS_interaction_fullExtent.png",
       dpi = 450, width = 6, height = 3.5)


cp <- cowplot::plot_grid(temp_ows, prec_ows, labels = c("A","B"))
cp

ggsave(plot = cp, filename = "MikeDir/figures/OWS_interaction_fullExtent.png",
       dpi = 450, width = 12, height = 4)