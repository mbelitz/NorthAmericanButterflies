#load libraries
library(lme4)
library(lmerTest)
library(car)
library(tidyverse)
library(MuMIn)
library(ape)
library(stringr)

bfly_trends = read.table("data/butterfly_traits_envars_trends_50km_NoMergeJul_m5_trim_6traits.txt",sep='\t',as.is=T,check.names=F,header=T)

# build first model

m <- lmer(Abundance.trend ~ Cropland.trend + Built.2005_2015 + Precip.1993_2018 +
          Temp.1993_2018 + LarvalColor + LarvalHair + AdultSize + AdultColor +  
            Diet.breadth.families + (1 | grid_id) + (1 | Species),
          data = bfly_trends, REML = F, lmerControl(optimizer = "bobyqa"))

ms <- step(m)

## so our results are similar to crossley et all, which is cool 
cross_top <- lmer(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                (1 | grid_id) + (1 | Species),
              data = bfly_trends, REML = F, lmerControl(optimizer = "bobyqa"))

summary(cross_top)
r.squaredGLMM(cross_top) ## 0.00657 R2m -- 0.588 R2c

# let's filter butterflies to only our study extent
# (32° N, 50° N) and (95° W, 60° W), 
bfly_trends_crop <- bfly_trends %>% 
  filter(Lon > -95 & Lon < -60 ) %>% 
  filter(Lat > 32 & Lat < 50)


## let's get a list of species
spp_list <- bfly_trends$Species %>% unique()

write.csv(spp_list, "MikeDir/data/spp_list.csv", row.names = F)

#simpleton traits
st <- read.csv("MikeDir/data/TraitsButterflyCrossley.csv")

bfly_trends_st <- left_join(bfly_trends_crop, st)

# do any species not have OWS?
dplyr::filter(bfly_trends_st, is.na(overwintering.stage))
# voltinism?
dplyr::filter(bfly_trends_st, is.na(voltinism))

# great, we have traits for all species so let's make into mdf
mdf <- bfly_trends_st %>% 
  filter(overwintering.stage %in% c("L", "A", "P", "M", "E")) %>% 
  filter(voltinism != "")

# read in phylogeny
t <- read.tree(file = "MikeDir/data/iScience/SupDryad_treepl.tre")

# read in tiplabels that will match Spp of mdf
tip_labels_df <- read.csv("MikeDir/data/iScience/tipLabel.csv")

t$tip.label <- tip_labels_df$validName

# run linear model to get top model
m <- lmer(Abundance.trend ~ Cropland.trend + Built.2005_2015 + Precip.1993_2018 +
            Temp.1993_2018 + LarvalColor + LarvalHair + AdultSize + AdultColor +  
            Diet.breadth.families + overwintering.stage + voltinism + 
            Temp.1993_2018:overwintering.stage + 
            overwintering.stage:voltinism + 
            Precip.1993_2018:overwintering.stage +
            (1 | grid_id) + (1 | Species),
          data = mdf, REML = F, lmerControl(optimizer = "bobyqa"))

ms <- step(m)
ms

m_top <- lmer(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                AdultColor + 
                overwintering.stage + 
                Precip.1993_2018:overwintering.stage +
                (1 | grid_id) + (1 | Species),
              data = mdf, REML = F, lmerControl(optimizer = "bobyqa"))


summary(m_top)                
vif(m_top)                
r.squaredGLMM(m_top) ## 0.0379 # 0.547

library(sjPlot)
plot_model(m_top, type = "pred", terms = "overwintering.stage")
plot_model(m_top, type = "pred", terms = "AdultColor")
plot_model(m_top, type = "pred", terms = c("Precip.1993_2018", "overwintering.stage"))

## Remove OWS that aren't in our study
mdf2 <- filter(mdf, overwintering.stage %in% c("E", "L", "P"))

m <- lmer(Abundance.trend ~ Cropland.trend + Built.2005_2015 + Precip.1993_2018 +
            Temp.1993_2018 + LarvalColor + LarvalHair + AdultSize + AdultColor +  
            Diet.breadth.families + overwintering.stage + voltinism + 
            Temp.1993_2018:overwintering.stage + 
            overwintering.stage:voltinism + 
            Precip.1993_2018:overwintering.stage +
            (1 | grid_id) + (1 | Species),
          data = mdf2, REML = F, lmerControl(optimizer = "bobyqa"))

ms <- step(m)
ms

m_top <- lmer(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                overwintering.stage + 
                Precip.1993_2018:overwintering.stage +
                (1 | grid_id) + (1 | Species),
              data = mdf2, REML = F, lmerControl(optimizer = "bobyqa"))


summary(m_top)                
vif(m_top)                
r.squaredGLMM(m_top) ## 0.0379 # 0.547

plot_model(m_top, type = "pred", terms = c("Precip.1993_2018", "overwintering.stage"))

# drop species not in phylogeny
mdf3 <- mdf2 %>% 
  filter(!Species %in% c("Brephidium isophthalma", "Copaeodes minimus", 
                         "Nymphalis vaualbum", "Papilio palamedes", "Satyrium edwardsii",
                         "Speyeria aphrodite", "Pieris napi"))

missing <- mdf3 %>% 
  filter(!mdf3$Species %in% t$tip.label)

## now add in a phylogenetic effect
library(phyr)

mdf3 <- mdf3 %>% 
  mutate(Precip.1993_2018 = scale(Precip.1993_2018),
         Temp.1993_2018 = scale(Temp.1993_2018))

tree_sp <- t$tip.label
sppNotInAnalysis <- data.frame(Species = tree_sp) %>% 
  filter(!Species %in% mdf3$Species)

t2 <- ape::drop.tip(t, tip = sppNotInAnalysis$Species)


pglmm <- pglmm(Abundance.trend ~ Precip.1993_2018 + Temp.1993_2018 + 
                 overwintering.stage + 
                 Precip.1993_2018:overwintering.stage +
                 (1 | grid_id) + (1 | Species__),
               data = mdf3,
               cov_ranef = list(Species = t2), 
               bayes = T, verbose = T)


summary(pglmm)

phyr::plot_bayes(x = pglmm)
